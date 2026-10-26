"""Download NOAA ETOPO1 (ice surface) from ERDDAP at every 2nd grid point (2 arc-minutes), in bands.

Usage: python3 fetch_etopo.py OUT_DIR
Writes OUT_DIR/band_XX.nc and OUT_DIR/etopo_2min.npy (int16, rows south->north, 5401 x 10801).
"""
import os
import sys
import time
import urllib.request

import numpy as np

BASE = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?altitude'
NLAT = 10801


def get(url, path, tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Underneath-build/1.0'})
            with urllib.request.urlopen(req, timeout=600) as r, open(path, 'wb') as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            return
        except Exception as e:
            print(f'  attempt {k + 1} failed: {e}', file=sys.stderr)
            time.sleep(5 * (k + 1))
    raise SystemExit('ETOPO download failed: ' + url)


def main():
    out = sys.argv[1]
    from scipy.io import netcdf_file
    rows = []
    step = 900
    for i0 in range(0, NLAT, step):
        i1 = min(i0 + step - 2, NLAT - 1)
        path = os.path.join(out, f'band_{i0 // step:02d}.nc')
        if not os.path.exists(path) or os.path.getsize(path) < 1000:
            url = f'{BASE}%5B{i0}:2:{i1}%5D%5B0:2:21600%5D'
            print('fetching rows', i0, i1, flush=True)
            get(url, path)
        with netcdf_file(path, 'r', mmap=False) as nc:
            a = np.array(nc.variables['altitude'].data, dtype=np.int16)
            lat = np.array(nc.variables['latitude'].data)
        print('  band', i0, a.shape, lat[0], lat[-1], flush=True)
        rows.append(a)
    z = np.concatenate(rows, 0)
    np.save(os.path.join(out, 'etopo_2min.npy'), z)
    print('saved', z.shape, z.min(), z.max())


if __name__ == '__main__':
    main()
