"""Fetch USGS M4.5+ earthquakes for the past 365 days in monthly chunks.

Usage: python3 fetch_quakes.py OUT_DIR END_ISO
Writes OUT_DIR/year45.json (merged GeoJSON features, de-duplicated by id).
"""
import datetime as dt
import json
import os
import sys
import time
import urllib.request

API = 'https://earthquake.usgs.gov/fdsnws/event/1/query'


def get(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Underneath-build/1.0'})
            with urllib.request.urlopen(req, timeout=240) as r:
                return json.loads(r.read())
        except Exception as e:
            print(f'  attempt {k + 1} failed: {e}', file=sys.stderr)
            time.sleep(3 * (k + 1))
    raise SystemExit('USGS fetch failed: ' + url)


def main():
    out, end = sys.argv[1], dt.datetime.fromisoformat(sys.argv[2])
    start = end - dt.timedelta(days=365)
    feats = {}
    a = start
    while a < end:
        b = min(a + dt.timedelta(days=30), end)
        url = (f'{API}?format=geojson&starttime={a:%Y-%m-%dT%H:%M:%S}&endtime={b:%Y-%m-%dT%H:%M:%S}'
               '&minmagnitude=4.5&eventtype=earthquake&orderby=time-asc')
        d = get(url)
        for f in d['features']:
            feats[f['id']] = f
        print(f'{a:%Y-%m-%d} .. {b:%Y-%m-%d}: {len(d["features"])}')
        a = b
    arr = sorted(feats.values(), key=lambda f: f['properties']['time'])
    json.dump({'start': start.strftime('%Y-%m-%dT%H:%M:%SZ'), 'end': end.strftime('%Y-%m-%dT%H:%M:%SZ'),
               'features': arr}, open(os.path.join(out, 'year45.json'), 'w'))
    print('total', len(arr))


if __name__ == '__main__':
    main()
