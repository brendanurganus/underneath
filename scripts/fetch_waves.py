"""Download and process real long-period seismograms (record sections) for Underneath's Wave Race.

Usage: python3 fetch_waves.py YEAR_JSON GIANTS_GEOJSON RAW_DIR OUT_DIR
  YEAR_JSON       merged USGS features (fetch_quakes.py output); its M7+ events are used
  GIANTS_GEOJSON  USGS M8.5+ since 1900; events from 2004 onward are used (digital era)
  RAW_DIR         where raw miniSEED is cached
  OUT_DIR         writes OUT_DIR/<usgs id>.json per event and OUT_DIR/index.json

Data: EarthScope (IRIS) FDSN station + dataselect services, LHZ channels (1 sample/s vertical) from the
GSN and partner networks. Traces are demeaned, detrended, tapered and band-passed 10-125 s, normalised
per trace and stored as 8-bit mu-law (base64).
"""
import base64
import io
import json
import math
import os
import sys
import time
import urllib.request

import numpy as np
from obspy import UTCDateTime, read
from obspy.geodetics import gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel

SVC = 'https://service.earthscope.org/fdsnws'
NETS = 'IU,II,IC,CU,GT,G,GE,MN,CN,AU,IM'
MU = 255.0
MODEL = TauPyModel('iasp91')


def http(url, data=None, tries=4, timeout=600):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'User-Agent': 'Underneath-build/1.0'})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 204:
                return b''
            print(f'  HTTP {e.code} attempt {k + 1}', file=sys.stderr)
        except Exception as e:
            print(f'  attempt {k + 1} failed: {e}', file=sys.stderr)
        time.sleep(4 * (k + 1))
    return None


def stations(t0, t1):
    url = (f'{SVC}/station/1/query?net={NETS}&cha=LHZ&level=channel&format=text'
           f'&starttime={t0.strftime("%Y-%m-%dT%H:%M:%S")}&endtime={t1.strftime("%Y-%m-%dT%H:%M:%S")}')
    raw = http(url)
    out = {}
    if not raw:
        return out
    for line in raw.decode().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        p = line.split('|')
        net, sta, loc, cha = p[0], p[1], p[2], p[3]
        key = (net, sta)
        rank = {'00': 0, '10': 1, '': 2, '20': 3}.get(loc, 4)
        if key not in out or rank < out[key]['rank']:
            out[key] = {'net': net, 'sta': sta, 'loc': loc, 'lat': float(p[4]), 'lon': float(p[5]), 'rank': rank}
    return out


def first(arrs, names):
    t = [a.time for a in arrs if a.name in names]
    return min(t) if t else None


def compand(x):
    y = np.sign(x) * np.log1p(MU * np.abs(x)) / math.log1p(MU)
    return np.clip(np.round(y * 127), -127, 127).astype(np.int8)


def process(ev, raw_dir, out_dir):
    eid, lat, lon, dep, mag = ev['id'], ev['lat'], ev['lon'], ev['depth'], ev['mag']
    t0 = UTCDateTime(ev['time'] / 1000.0)
    giant = mag >= 8.5
    win = 240 * 60 if giant else 100 * 60
    sps = 0.5 if giant else 1.0
    t1 = t0 + win
    cache = os.path.join(raw_dir, eid + '.mseed')
    sts = stations(t0, t1)
    print(f'{eid} M{mag} {ev["place"]}: {len(sts)} stations', flush=True)
    if not os.path.exists(cache):
        loc = lambda s: s['loc'] if s['loc'] else '--'
        body = '\n'.join(f'{s["net"]} {s["sta"]} {loc(s)} LHZ {(t0 - 120).strftime("%Y-%m-%dT%H:%M:%S")} '
                         f'{(t1 + 120).strftime("%Y-%m-%dT%H:%M:%S")}' for s in sts.values())
        raw = http(f'{SVC}/dataselect/1/query', data=body.encode())
        if not raw:
            print('  no waveform data', flush=True)
            return None
        open(cache, 'wb').write(raw)
    st = read(cache)
    picks = []
    for (net, sta), s in sts.items():
        sel = st.select(network=net, station=sta, channel='LHZ',
                        location=s['loc'] if s['loc'] else '')
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
            miss = np.ma.getmaskarray(data).mean()
            if miss > 0.03:
                continue
            tr.data = data.filled(0) if miss == 0 else np.interp(
                np.arange(len(data)), np.flatnonzero(~np.ma.getmaskarray(data)),
                np.asarray(data)[~np.ma.getmaskarray(data)]).astype(np.float64)
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
        arrs = MODEL.get_travel_times(source_depth_in_km=max(dep, 0.0), distance_in_degree=dist,
                                      phase_list=['p', 'P', 'Pdiff', 'PKP', 'PKIKP'])
        tp = first(arrs, {'p', 'P', 'Pdiff', 'PKP', 'PKIKP'})
        if tp is None:
            continue
        ip = int(tp * sps)
        if ip > 150 * sps:
            noise = x[: max(int(ip - 40 * sps), 10)]
            sig = x[ip:]
            snr = np.abs(sig).max() / (np.sqrt(np.mean(noise ** 2)) + 1e-30)
        else:
            snr = 50.0
        peak_i = int(np.argmax(np.abs(x)))
        if peak_i < ip - 20 * sps or snr < 8:
            continue
        x = x / np.abs(x).max()
        picks.append({'net': net, 'sta': sta, 'lat': round(s['lat'], 3), 'lon': round(s['lon'], 3),
                      'dist': round(dist, 2), 'az': round(az, 1), 'snr': float(snr), 'x': x})
    # spread across distance: best signal in each bin
    nb = 56
    best = {}
    for p in picks:
        b = min(int(p['dist'] / 180 * nb), nb - 1)
        if b not in best or p['snr'] > best[b]['snr']:
            best[b] = p
    chosen = sorted(best.values(), key=lambda p: p['dist'])
    print(f'  usable {len(picks)}, chosen {len(chosen)}', flush=True)
    if len(chosen) < 12:
        return None
    rec = {'id': eid, 'time': ev['time'], 'lat': lat, 'lon': lon, 'depth': dep, 'mag': mag, 'place': ev['place'],
           'sps': sps, 'seconds': win, 'filter': '10-125 s band-pass, vertical (LHZ)', 'mu': MU,
           'stations': [{k: p[k] for k in ('net', 'sta', 'lat', 'lon', 'dist', 'az')} for p in chosen],
           'traces': [base64.b64encode(compand(p['x']).tobytes()).decode() for p in chosen],
           'source': 'EarthScope (IRIS) FDSN dataselect'}
    path = os.path.join(out_dir, eid + '.json')
    json.dump(rec, open(path, 'w'), separators=(',', ':'))
    print(f'  wrote {path} {os.path.getsize(path) / 1e6:.2f} MB', flush=True)
    return {'id': eid, 'n': len(chosen), 'bytes': os.path.getsize(path)}


def main():
    year_json, giants_json, raw_dir, out_dir = sys.argv[1:5]
    evs = []
    for f in json.load(open(year_json))['features']:
        p = f['properties']
        if p['mag'] is not None and p['mag'] >= 7.0:
            evs.append(f)
    for f in json.load(open(giants_json))['features']:
        if f['properties']['time'] >= UTCDateTime(2004, 1, 1).timestamp * 1000:
            evs.append(f)
    todo = []
    for f in evs:
        lon, lat, dep = f['geometry']['coordinates']
        p = f['properties']
        todo.append({'id': f['id'], 'time': p['time'], 'lat': lat, 'lon': lon, 'depth': dep or 10.0,
                     'mag': p['mag'], 'place': p['place']})
    print(len(todo), 'events', flush=True)
    index = []
    for ev in todo:
        path = os.path.join(out_dir, ev['id'] + '.json')
        if os.path.exists(path):
            d = json.load(open(path))
            index.append({'id': ev['id'], 'n': len(d['stations']), 'bytes': os.path.getsize(path)})
            continue
        try:
            r = process(ev, raw_dir, out_dir)
        except Exception as e:
            print('  failed', ev['id'], e, flush=True)
            r = None
        if r:
            index.append(r)
    json.dump(index, open(os.path.join(out_dir, 'index.json'), 'w'))
    print('done', len(index), 'events,', sum(i['bytes'] for i in index) / 1e6, 'MB')


if __name__ == '__main__':
    main()
