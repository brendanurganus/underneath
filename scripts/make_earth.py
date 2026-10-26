"""Earth model and travel-time tables for Underneath's Wave Race.

Usage: python3 make_earth.py SITE_DIR
Writes SITE_DIR/data/earth.json:
  model: iasp91 velocity profile (depth km, vp, vs km/s) as shipped with ObsPy's TauP
  tt: first-arrival times (seconds) by phase, for source depths x distances 0..180 degrees, from TauP (iasp91)
"""
import json
import os
import sys

import obspy
from obspy.taup import TauPyModel

DEPTHS = [0, 10, 20, 35, 50, 70, 100, 150, 200, 300, 400, 500, 600, 700]
PHASES = {
    'P': ['p', 'P'], 'Pdiff': ['Pdiff'], 'PKP': ['PKIKP', 'PKP', 'PKiKP'], 'PP': ['PP'],
    'S': ['s', 'S'], 'Sdiff': ['Sdiff'], 'SKS': ['SKS', 'SKIKS'], 'SS': ['SS'],
}


def main():
    site = sys.argv[1]
    tvel = os.path.join(os.path.dirname(obspy.__file__), 'taup', 'data', 'iasp91.tvel')
    prof = []
    for line in open(tvel).read().splitlines()[2:]:
        p = line.split()
        if len(p) >= 3:
            prof.append([float(p[0]), float(p[1]), float(p[2])])
    m = TauPyModel('iasp91')
    allp = sorted({x for v in PHASES.values() for x in v})
    tt = {k: [] for k in PHASES}
    for d in DEPTHS:
        rows = {k: [] for k in PHASES}
        for dist in range(0, 181):
            arrs = m.get_travel_times(source_depth_in_km=d, distance_in_degree=max(dist, 0.01), phase_list=allp)
            for k, names in PHASES.items():
                t = [a.time for a in arrs if a.name in names]
                rows[k].append(round(min(t), 1) if t else None)
        for k in PHASES:
            tt[k].append(rows[k])
        print('depth', d, 'done', flush=True)
    out = {'model': 'iasp91 (Kennett & Engdahl 1991), via ObsPy TauP', 'profile': prof,
           'layers': [
               {'name': 'Crust', 'top': 0, 'bottom': 35},
               {'name': 'Upper mantle', 'top': 35, 'bottom': 410},
               {'name': 'Transition zone', 'top': 410, 'bottom': 660},
               {'name': 'Lower mantle', 'top': 660, 'bottom': 2889},
               {'name': 'Outer core (liquid)', 'top': 2889, 'bottom': 5153.9},
               {'name': 'Inner core (solid)', 'top': 5153.9, 'bottom': 6371}],
           'depths': DEPTHS, 'phases': list(PHASES), 'tt': tt}
    os.makedirs(os.path.join(site, 'data'), exist_ok=True)
    path = os.path.join(site, 'data', 'earth.json')
    json.dump(out, open(path, 'w'), separators=(',', ':'))
    print('wrote', path, os.path.getsize(path) // 1000, 'KB')


if __name__ == '__main__':
    main()
