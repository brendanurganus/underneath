"""Build Underneath's data files from the raw downloads.

Usage: python3 make_data.py RAW_DIR SITE_DIR
Writes SITE_DIR/data/quakes.json, plates.json, volcanoes.json, earth.json.
"""
import json
import math
import os
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

CLASSES = ['SUB', 'OCB', 'CCB', 'OSR', 'CRB', 'OTF', 'CTF']  # convergent x3, divergent x2, transform x2


def unit(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def r(v, n):
    return None if v is None else round(float(v), n)


def main():
    raw, site = sys.argv[1], sys.argv[2]
    out = os.path.join(site, 'data')
    os.makedirs(out, exist_ok=True)

    # ---------- plates ----------
    names = {f['properties']['Code']: f['properties']['PlateName']
             for f in json.load(open(os.path.join(raw, 'plates/PB2002_plates.json')))['features']}
    steps = json.load(open(os.path.join(raw, 'plates/PB2002_steps.json')))['features']
    steps.sort(key=lambda f: (f['properties']['PLATEBOUND'], f['properties']['SEQNUM']))
    bounds = sorted({f['properties']['PLATEBOUND'] for f in steps})
    bidx = {b: i for i, b in enumerate(bounds)}
    rows = []
    for f in steps:
        p = f['properties']
        c = f['geometry']['coordinates']
        rows.append([round(c[0][0], 3), round(c[0][1], 3), round(c[-1][0], 3), round(c[-1][1], 3),
                     CLASSES.index(p['STEPCLASS']), r(p['VELOCITYLE'], 1), r(p['VELOCITYAZ'], 0),
                     r(p['VELOCITYDI'], 1), r(p['VELOCITYRI'], 1), bidx[p['PLATEBOUND']]])
    bnames = []
    for b in bounds:
        a, c = re.split(r'[-/\\]', b)[:2] if re.search(r'[-/\\]', b) else (b, '')
        bnames.append([b, names.get(a, a), names.get(c, c)])
    plates = {'source': 'Bird (2003) PB2002 plate model, via github.com/fraxen/tectonicplates',
              'classes': ['Subduction zone', 'Oceanic convergent boundary', 'Continental collision',
                          'Spreading ridge', 'Continental rift', 'Oceanic transform fault', 'Continental transform fault'],
              'fields': ['lon1', 'lat1', 'lon2', 'lat2', 'cls', 'vel', 'velAz', 'div', 'right', 'bound'],
              'bounds': bnames, 'rows': rows}
    json.dump(plates, open(os.path.join(out, 'plates.json'), 'w'), separators=(',', ':'))

    mid = unit(np.array([(s[1] + s[3]) / 2 for s in rows]), np.array([(s[0] + s[2]) / 2 for s in rows]))
    tree = cKDTree(mid)

    # ---------- quakes ----------
    def quake_rows(feats):
        out_rows = []
        for f in feats:
            p = f['properties']
            lon, lat, dep = f['geometry']['coordinates']
            if p['mag'] is None:
                continue
            d, k = tree.query(unit(np.array(lat), np.array(lon)))
            km = 2 * math.asin(min(1.0, d / 2)) * 6371.0
            out_rows.append([f['id'], p['time'], round(lat, 3), round(lon, 3), round(dep if dep is not None else 10.0, 1),
                             round(p['mag'], 1), p.get('magType') or '', p.get('place') or '', p.get('alert') or '',
                             p.get('felt') or 0, r(p.get('mmi'), 1), int(k) if km < 650 else -1])
        return out_rows

    year = json.load(open(os.path.join(raw, 'usgs/year45.json')))
    giants = json.load(open(os.path.join(raw, 'usgs/giants.geojson')))
    q = {'source': 'USGS ANSS Comprehensive Catalog (ComCat), earthquake.usgs.gov',
         'fetched': year['end'], 'start': year['start'], 'end': year['end'],
         'fields': ['id', 'time', 'lat', 'lon', 'depth', 'mag', 'magType', 'place', 'alert', 'felt', 'mmi', 'seg'],
         'year': quake_rows(year['features']), 'giants': quake_rows(giants['features'])}
    json.dump(q, open(os.path.join(out, 'quakes.json'), 'w'), separators=(',', ':'), ensure_ascii=False)

    # ---------- volcanoes ----------
    alerts = {str(v['vnum']): v for v in json.load(open(os.path.join(raw, 'volc/hans_elevated.json')))}
    vrows = []
    for f in json.load(open(os.path.join(raw, 'volc/gvp.json')))['features']:
        p = f['properties']
        summ = (p.get('Geological_Summary') or '').strip()
        first = re.split(r'(?<=[a-z0-9\)])\.\s', summ, maxsplit=1)[0]
        if len(first) > 220:
            first = first[:217].rsplit(' ', 1)[0] + '…'
        elif first and not first.endswith('.'):
            first += '.'
        a = alerts.get(str(p['Volcano_Number']))
        vrows.append([p['Volcano_Number'], p['Volcano_Name'], r(p['Latitude'], 3), r(p['Longitude'], 3),
                      p.get('Elevation'), p.get('Primary_Volcano_Type') or '', p.get('Last_Eruption_Year'),
                      p.get('Country') or '', first,
                      [a['color_code'], a['alert_level'], a['sent_utc'], a['obs_fullname']] if a else None])
    vol = {'source': 'Smithsonian Institution Global Volcanism Program, Volcanoes of the World (Holocene list); '
                     'alert levels from the USGS Volcano Hazards Program (HANS)',
           'fetched': year['end'],
           'fields': ['num', 'name', 'lat', 'lon', 'elev', 'type', 'lastYear', 'country', 'summary', 'alert'],
           'rows': vrows}
    json.dump(vol, open(os.path.join(out, 'volcanoes.json'), 'w'), separators=(',', ':'), ensure_ascii=False)

    for f in ('quakes.json', 'plates.json', 'volcanoes.json'):
        print(f, round(os.path.getsize(os.path.join(out, f)) / 1e3), 'KB')
    print('quakes', len(q['year']), 'giants', len(q['giants']), 'on boundaries',
          sum(1 for x in q['year'] if x[-1] >= 0), 'volcanoes', len(vrows), 'alerts',
          [v[1] for v in vrows if v[-1]])


if __name__ == '__main__':
    main()
