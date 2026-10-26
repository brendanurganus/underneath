"""Underneath data refresh.

Rebuilds the page's earthquake data and volcano alert levels from live public sources, adds real
seismograms for new big quakes, and drops seismograms for quakes that have left the past-year window.
Plate boundaries, the Earth model and the relief texture never change, so they are not touched.

Usage: python3 refresh.py CURRENT_DIR OUT_DIR
  CURRENT_DIR  folder holding the currently published files at their published paths:
               data/quakes.json, data/volcanoes.json, data/plates.json, data/waves/index.json
               (data/waves/skipped.json is optional)
  OUT_DIR      receives the new files at their published paths, plus OUT_DIR/manifest.json:
               {"files": {"<published path>": "<absolute local path>", ...},
                "remove": ["<published path>", ...], "summary": {...}}
Exit 0: refreshed; manifest written.
Exit 2: the USGS catalog could not be downloaded, or failed a sanity check; nothing written.
Needs: Python 3 with numpy. ObsPy is pip-installed automatically, only when a new big quake needs seismograms.
Network: earthquake.usgs.gov, volcanoes.usgs.gov, service.earthscope.org (and pypi.org for ObsPy).
"""
import base64
import datetime as dt
import json
import math
import os
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

import numpy as np

USGS = 'https://earthquake.usgs.gov/fdsnws/event/1'
HANS = 'https://volcanoes.usgs.gov/hans-public/api/volcano/getElevatedVolcanoes'
SVC = 'https://service.earthscope.org/fdsnws'
NETS = 'IU,II,IC,CU,GT,G,GE,MN,CN,AU,IM'
UA = {'User-Agent': 'Underneath-refresh/1.0'}
MU = 255.0
DIGITAL_ERA_MS = 1072915200000  # 2004-01-01 UTC
LATENCY_H = 6       # wait this long after a quake's recording window ends before fetching seismograms
GIVE_UP_DAYS = 3    # stop retrying seismograms for a quake this old that still lacks usable stations


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def http(url, data=None, tries=4, timeout=300, as_json=False):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read()
            return json.loads(body) if as_json else body
        except urllib.error.HTTPError as e:
            if e.code == 204:
                return [] if as_json else b''
            log(f'  HTTP {e.code} on attempt {k + 1}: {url[:140]}')
        except Exception as e:
            log(f'  attempt {k + 1} failed ({e}): {url[:140]}')
        time.sleep(4 * (k + 1))
    return None


def iso(t):
    return t.strftime('%Y-%m-%dT%H:%M:%S')


# ---------------------------------------------------------------- quakes
def fetch_year(end):
    start = end - dt.timedelta(days=365)
    feats, a = {}, start
    while a < end:
        b = min(a + dt.timedelta(days=30), end)
        d = http(f'{USGS}/query?format=geojson&starttime={iso(a)}&endtime={iso(b)}&minmagnitude=4.5'
                 '&eventtype=earthquake&orderby=time-asc', as_json=True)
        if d is None:
            return None, start
        for f in d['features']:
            feats[f['id']] = f
        a = b
    return sorted(feats.values(), key=lambda f: f['properties']['time']), start


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def plate_mids(plates):
    r = np.array([row[:4] for row in plates['rows']], dtype=np.float64)
    m = unit(r[:, 1], r[:, 0]) + unit(r[:, 3], r[:, 2])
    return m / np.linalg.norm(m, axis=1, keepdims=True)


def snap(lat, lon, mids):
    """Index of the nearest plate-boundary step, or -1 when it is more than 650 km away."""
    q = unit(np.asarray(lat, float), np.asarray(lon, float))
    out = np.empty(len(q), dtype=np.int64)
    for i in range(0, len(q), 2000):
        dots = q[i:i + 2000] @ mids.T
        k = dots.argmax(1)
        km = np.arccos(np.clip(dots[np.arange(len(k)), k], -1, 1)) * 6371.0
        out[i:i + 2000] = np.where(km < 650, k, -1)
    return out


def is_preliminary_duplicate(f):
    """Tsunami-warning-center entries with only an initial magnitude (Mi/Mwp) that USGS never merged into its own
    solution. In this catalog they duplicate real events (sometimes with the wrong date), so they are left out."""
    p = f['properties']
    return (p.get('net') or '').lower() in ('pt', 'at') and (p.get('magType') or '').lower() in ('mi', 'mwp')


DROPPED = []


def rows_from(feats, mids):
    for f in feats:
        if is_preliminary_duplicate(f):
            DROPPED.append(f"{f['id']} M{f['properties'].get('mag')} {f['properties'].get('place')}")
    keep = [f for f in feats if f['properties'].get('mag') is not None and not is_preliminary_duplicate(f)]
    if not keep:
        return []
    segs = snap([f['geometry']['coordinates'][1] for f in keep], [f['geometry']['coordinates'][0] for f in keep], mids)
    rows = []
    for f, s in zip(keep, segs):
        p = f['properties']
        lon, lat, dep = f['geometry']['coordinates']
        mmi = p.get('mmi')
        rows.append([f['id'], p['time'], round(lat, 3), round(lon, 3), round(dep if dep is not None else 10.0, 1),
                     round(p['mag'], 1), p.get('magType') or '', p.get('place') or '', p.get('alert') or '',
                     p.get('felt') or 0, None if mmi is None else round(float(mmi), 1), int(s)])
    return rows


def us_context(start, end):
    base = f'{USGS}/count?starttime={iso(start)}&endtime={iso(end)}&eventtype=earthquake'
    allq = http(base, timeout=180)
    l48 = http(base + '&minlatitude=24&maxlatitude=50&minlongitude=-125&maxlongitude=-66', timeout=180)
    ak = http(base + '&minlatitude=50&maxlatitude=72&minlongitude=-180&maxlongitude=-129', timeout=180)
    try:
        return {'all': int(allq), 'us': int(l48) + int(ak)}
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------- seismograms
def ensure_obspy():
    try:
        import obspy  # noqa: F401
        return True
    except ImportError:
        log('installing ObsPy…')
        r = subprocess.run([sys.executable, '-m', 'pip', 'install', '--break-system-packages', '-q', 'obspy'],
                           capture_output=True, text=True, timeout=900)
        if r.returncode != 0:
            log(r.stderr[-800:])
            return False
        try:
            import obspy  # noqa: F401,F811
            return True
        except ImportError:
            return False


def compand(x):
    y = np.sign(x) * np.log1p(MU * np.abs(x)) / math.log1p(MU)
    return np.clip(np.round(y * 127), -127, 127).astype(np.int8)


def process_waves(ev, raw_dir, out_dir):
    """Download and process one quake's long-period vertical seismograms. Returns {id, n, bytes} or None."""
    from obspy import UTCDateTime, read
    from obspy.geodetics import gps2dist_azimuth, locations2degrees
    from obspy.taup import TauPyModel
    model = TauPyModel('iasp91')
    eid, lat, lon, dep, mag = ev['id'], ev['lat'], ev['lon'], max(ev['depth'], 0.0), ev['mag']
    t0 = UTCDateTime(ev['time'] / 1000.0)
    giant = mag >= 8.5
    win = 240 * 60 if giant else 100 * 60
    sps = 0.5 if giant else 1.0
    t1 = t0 + win
    raw = http(f'{SVC}/station/1/query?net={NETS}&cha=LHZ&level=channel&format=text'
               f'&starttime={t0.strftime("%Y-%m-%dT%H:%M:%S")}&endtime={t1.strftime("%Y-%m-%dT%H:%M:%S")}')
    if not raw:
        return None
    sts = {}
    for line in raw.decode().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        p = line.split('|')
        rank = {'00': 0, '10': 1, '': 2, '20': 3}.get(p[2], 4)
        key = (p[0], p[1])
        if key not in sts or rank < sts[key]['rank']:
            sts[key] = {'net': p[0], 'sta': p[1], 'loc': p[2], 'lat': float(p[4]), 'lon': float(p[5]), 'rank': rank}
    body = '\n'.join(f'{s["net"]} {s["sta"]} {s["loc"] or "--"} LHZ {(t0 - 120).strftime("%Y-%m-%dT%H:%M:%S")} '
                     f'{(t1 + 120).strftime("%Y-%m-%dT%H:%M:%S")}' for s in sts.values())
    mseed = http(f'{SVC}/dataselect/1/query', data=body.encode(), timeout=900)
    if not mseed:
        return None
    path = os.path.join(raw_dir, eid + '.mseed')
    open(path, 'wb').write(mseed)
    st = read(path)
    picks = []
    for (net, sta), s in sts.items():
        sel = st.select(network=net, station=sta, channel='LHZ', location=s['loc'] or '')
        if not len(sel):
            continue
        sel = sel.copy()
        try:
            sel.merge(method=1, fill_value=None)
        except Exception:
            continue
        tr = sel[0]
        if abs(tr.stats.sampling_rate - 1.0) > 1e-3:
            continue
        tr.trim(t0 - 60, t1 + 60, pad=True, fill_value=None)
        data = tr.data
        if np.ma.isMaskedArray(data):
            mask = np.ma.getmaskarray(data)
            if mask.mean() > 0.03:
                continue
            good = np.flatnonzero(~mask)
            tr.data = np.interp(np.arange(len(data)), good, np.asarray(data)[good]).astype(np.float64)
        tr.data = np.asarray(tr.data, dtype=np.float64)
        if tr.data.std() == 0:
            continue
        tr.detrend('demean'); tr.detrend('linear'); tr.taper(0.03)
        tr.filter('bandpass', freqmin=0.008, freqmax=0.1, corners=4, zerophase=True)
        if giant:
            tr.decimate(2, no_filter=True)
        tr.trim(t0, t1, pad=True, fill_value=0)
        x = tr.data[: int(win * sps)]
        if len(x) < int(win * sps) * 0.98:
            continue
        dist = locations2degrees(lat, lon, s['lat'], s['lon'])
        _, az, _ = gps2dist_azimuth(lat, lon, s['lat'], s['lon'])
        arrs = model.get_travel_times(source_depth_in_km=dep, distance_in_degree=dist,
                                      phase_list=['p', 'P', 'Pdiff', 'PKP', 'PKIKP'])
        if not arrs:
            continue
        ip = int(min(a.time for a in arrs) * sps)
        if ip > 150 * sps:
            noise = x[: max(int(ip - 40 * sps), 10)]
            snr = np.abs(x[ip:]).max() / (np.sqrt(np.mean(noise ** 2)) + 1e-30)
        else:
            snr = 50.0
        if int(np.argmax(np.abs(x))) < ip - 20 * sps or snr < 8:
            continue
        picks.append({'net': net, 'sta': sta, 'lat': round(s['lat'], 3), 'lon': round(s['lon'], 3),
                      'dist': round(dist, 2), 'az': round(az, 1), 'snr': float(snr), 'x': x / np.abs(x).max()})
    best = {}
    for p in picks:  # best signal in each distance bin, so the record section spans 0-180 degrees
        b = min(int(p['dist'] / 180 * 56), 55)
        if b not in best or p['snr'] > best[b]['snr']:
            best[b] = p
    chosen = sorted(best.values(), key=lambda p: p['dist'])
    log(f'  {eid}: {len(sts)} stations listed, {len(picks)} usable, {len(chosen)} chosen')
    if len(chosen) < 12:
        return None
    rec = {'id': eid, 'time': ev['time'], 'lat': lat, 'lon': lon, 'depth': dep, 'mag': mag, 'place': ev['place'],
           'sps': sps, 'seconds': win, 'filter': '10-125 s band-pass, vertical (LHZ)', 'mu': MU,
           'stations': [{k: p[k] for k in ('net', 'sta', 'lat', 'lon', 'dist', 'az')} for p in chosen],
           'traces': [base64.b64encode(compand(p['x']).tobytes()).decode() for p in chosen],
           'source': 'EarthScope (IRIS) FDSN dataselect'}
    out = os.path.join(out_dir, eid + '.json')
    json.dump(rec, open(out, 'w'), separators=(',', ':'))
    return {'id': eid, 'n': len(chosen), 'bytes': os.path.getsize(out)}


# ---------------------------------------------------------------- main
def load(path, default=None):
    try:
        return json.load(open(path, encoding='utf-8'))
    except (OSError, ValueError):
        return default


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    cur, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    if os.path.isfile(cur):  # a path to one of the files was passed: use the folder that holds data/
        cur = os.path.dirname(cur)
    while not os.path.isfile(os.path.join(cur, 'data/quakes.json')) and os.path.basename(cur) in ('data', 'waves', 'tools'):
        cur = os.path.dirname(cur)
    old_q = load(os.path.join(cur, 'data/quakes.json'))
    old_v = load(os.path.join(cur, 'data/volcanoes.json'))
    plates = load(os.path.join(cur, 'data/plates.json'))
    if not (old_q and old_v and plates):
        log('missing current data files in ' + cur)
        sys.exit(2)
    index = load(os.path.join(cur, 'data/waves/index.json'), []) or []
    skipped = load(os.path.join(cur, 'data/waves/skipped.json'), {}) or {}
    now = dt.datetime.now(dt.timezone.utc).replace(tzinfo=None, microsecond=0)
    mids = plate_mids(plates)
    summary = {}

    # 1) the past year of M4.5+ quakes
    feats, start = fetch_year(now)
    if feats is None:
        log('could not download the USGS catalog; keeping the published data')
        sys.exit(2)
    year = rows_from(feats, mids)
    if len(year) < max(1000, 0.85 * len(old_q['year'])):
        log(f'only {len(year)} quakes arrived (published: {len(old_q["year"])}); keeping the published data')
        sys.exit(2)
    old_ids = {r[0] for r in old_q['year']}
    fresh = sorted((r for r in year if r[0] not in old_ids and r[1] > Date_ms(old_q['end']) - 86400000 * 2),
                   key=lambda r: -r[5])
    summary['quakes'] = {'count': len(year), 'previous': len(old_q['year']), 'new_since_last': len(fresh),
                         'biggest_new': [f'M{r[5]} {r[7]}' for r in fresh[:3]]}

    # 2) giants since 1900 (rarely changes; keep the published list if USGS does not answer)
    g = http(f'{USGS}/query?format=geojson&starttime=1900-01-01&minmagnitude=8.5&eventtype=earthquake&orderby=time-asc',
             as_json=True)
    giants = rows_from(g['features'], mids) if g and g.get('features') else old_q['giants']
    summary['giants'] = len(giants)

    # 3) how US-heavy the full catalog is (the About panel quotes this)
    ctx = us_context(start, now) or old_q.get('context')
    summary['context'] = ctx

    quakes = {'source': old_q.get('source', 'USGS ANSS Comprehensive Catalog (ComCat), earthquake.usgs.gov'),
              'fetched': now.strftime('%Y-%m-%dT%H:%M:%SZ'), 'start': start.strftime('%Y-%m-%dT%H:%M:%SZ'),
              'end': now.strftime('%Y-%m-%dT%H:%M:%SZ'), 'fields': old_q['fields'], 'year': year, 'giants': giants}
    if ctx:
        quakes['context'] = ctx

    # 4) volcano alert levels (US volcanoes only)
    vol = old_v
    h = http(HANS, as_json=True, timeout=60)
    if isinstance(h, list):
        F = {f: i for i, f in enumerate(vol['fields'])}
        alerts = {str(a['vnum']): [a['color_code'], a['alert_level'], a['sent_utc'], a['obs_fullname']] for a in h}
        for row in vol['rows']:
            row[F['alert']] = alerts.get(str(row[F['num']]))
        vol['alertsFetched'] = quakes['fetched']
        names = {str(row[F['num']]): row[F['name']] for row in vol['rows']}
        summary['volcano_alerts'] = [f'{names.get(k, k)} {v[0].lower()}' for k, v in alerts.items()]
        unknown = [k for k in alerts if k not in names]
        if unknown:
            summary['volcano_alerts_not_on_map'] = unknown
    else:
        summary['volcano_alerts'] = 'USGS volcano alerts did not answer; kept the previous levels'

    # 5) seismograms
    os.makedirs(os.path.join(out, 'data/waves'), exist_ok=True)
    targets = {}
    for r in year:
        if r[5] >= 7.0:
            targets[r[0]] = r
    for r in giants:
        if r[1] >= DIGITAL_ERA_MS:
            targets[r[0]] = r
    have = {e['id'] for e in index}
    removed = sorted(have - set(targets))
    index = [e for e in index if e['id'] in targets]
    skipped = {k: v for k, v in skipped.items() if k in targets}
    todo, waiting = [], []
    for eid, r in targets.items():
        if eid in have or eid in skipped:
            continue
        win_h = (240 if r[5] >= 8.5 else 100) / 60
        age_h = ((now - dt.datetime(1970, 1, 1)).total_seconds() - r[1] / 1000) / 3600
        (todo if age_h >= win_h + LATENCY_H else waiting).append((eid, r, age_h))
    added, deferred = [], [f'{r[7]} M{r[5]} (too recent, {age:.0f} h old)' for _, r, age in waiting]
    if todo:
        if not ensure_obspy():
            deferred += [f'{r[7]} M{r[5]} (ObsPy could not be installed)' for _, r, _ in todo]
        else:
            with tempfile.TemporaryDirectory() as raw_dir:
                for eid, r, age_h in todo:
                    ev = {'id': eid, 'time': r[1], 'lat': r[2], 'lon': r[3], 'depth': r[4], 'mag': r[5], 'place': r[7]}
                    try:
                        res = process_waves(ev, raw_dir, os.path.join(out, 'data/waves'))
                    except Exception as e:
                        log(f'  {eid} failed: {e}')
                        res = None
                    if res:
                        index.append(res)
                        added.append(eid)
                    elif age_h / 24 >= GIVE_UP_DAYS:
                        skipped[eid] = f'fewer than 12 usable stations as of {quakes["fetched"]}'
                    else:
                        deferred.append(f'{r[7]} M{r[5]} (not enough station data yet)')
    summary['seismograms'] = {'bundled': len(index), 'added': added, 'removed': removed, 'deferred': deferred,
                              'gave_up': sorted(skipped)}
    summary['left_out_preliminary_duplicates'] = DROPPED

    # 6) write everything and the manifest
    files = {}

    def put(rel, obj):
        p = os.path.join(out, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        json.dump(obj, open(p, 'w', encoding='utf-8'), separators=(',', ':'), ensure_ascii=False)
        files[rel] = p

    put('data/quakes.json', quakes)
    put('data/volcanoes.json', vol)
    put('data/waves/index.json', index)
    put('data/waves/skipped.json', skipped)
    for eid in added:
        files[f'data/waves/{eid}.json'] = os.path.join(out, 'data/waves', eid + '.json')
    manifest = {'files': files, 'remove': [f'data/waves/{eid}.json' for eid in removed], 'summary': summary}
    json.dump(manifest, open(os.path.join(out, 'manifest.json'), 'w'), indent=1)
    print(json.dumps(summary, ensure_ascii=False))


def Date_ms(s):
    return int(dt.datetime.fromisoformat(s.rstrip('Z')).replace(tzinfo=dt.timezone.utc).timestamp() * 1000)


if __name__ == '__main__':
    main()
