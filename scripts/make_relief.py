"""Bake Underneath's globe texture from NOAA ETOPO1 elevations.

Usage: python3 make_relief.py ETOPO_NPY OUT_DIR
Writes OUT_DIR/relief_8k.jpg (8192x4096) and OUT_DIR/relief_4k.jpg (4096x2048).
Slate-blue seafloor shaded by depth, warm umber land shaded by height, both hill-shaded from the north-west.
"""
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

R = 6371000.0


def ramp(v, stops):
    xs = np.array([s[0] for s in stops], dtype=np.float32)
    cs = np.array([[int(s[1][i:i + 2], 16) for i in (1, 3, 5)] for s in stops], dtype=np.float32)
    out = np.empty(v.shape + (3,), np.float32)
    for k in range(3):
        out[..., k] = np.interp(v, xs, cs[:, k])
    return out


def hillshade(z, lat, px_w, px_h, exag, az=315.0, alt=38.0):
    dzdy, dzdx = np.gradient(z * exag)
    dx = (2 * np.pi * R * np.cos(np.radians(lat)))[:, None] / px_w
    dy = np.pi * R / px_h
    gx = dzdx / np.maximum(dx, 2000.0)
    gy = -dzdy / dy  # rows run north -> south
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    a, h = np.radians(360 - az + 90), np.radians(alt)
    s = np.sin(h) * np.cos(slope) + np.cos(h) * np.sin(slope) * np.cos(a - aspect)
    return np.clip(s, 0, 1)


def main():
    src, out = sys.argv[1], sys.argv[2]
    z0 = np.load(src).astype(np.float32)  # south -> north, 5401 x 10801, -180..180
    W, H = 8192, 4096
    lat = 90 - (np.arange(H) + 0.5) * 180.0 / H
    lon = -180 + (np.arange(W) + 0.5) * 360.0 / W
    rows = (lat + 90) / 180.0 * (z0.shape[0] - 1)
    cols = (lon + 180) / 360.0 * (z0.shape[1] - 1)
    rr, cc = np.meshgrid(rows, cols, indexing='ij')
    z = ndimage.map_coordinates(z0, [rr, cc], order=1, mode='nearest').astype(np.float32)
    del z0, rr, cc
    print('resampled', z.shape, z.min(), z.max(), flush=True)

    land = z > 0
    # shading: fine relief plus a broad pass so ridges and trenches read from far away
    zs = ndimage.gaussian_filter(z, 1.0)
    fine = hillshade(zs, lat, W, H, exag=np.where(land, 14.0, 22.0))
    broad = hillshade(ndimage.gaussian_filter(z, 6.0), lat, W, H, exag=np.where(land, 26.0, 70.0))
    flat = np.sin(np.radians(38.0))
    shade = 0.55 * fine + 0.45 * broad
    shade = shade / flat  # 1.0 on flat ground

    sea = ramp(-z, [(0, '#3c5462'), (150, '#334857'), (600, '#2a3c4b'), (2500, '#1f2f3e'),
                    (4200, '#182635'), (5800, '#121d2b'), (8000, '#0b1320'), (11000, '#060b14')])
    gnd = ramp(z, [(0, '#3f2e22'), (300, '#4a3628'), (900, '#5a4231'), (1800, '#6c503c'),
                   (3000, '#82634b'), (4500, '#98795f'), (6500, '#b09478'), (9000, '#c8ae93')])
    # ice sheets read as pale ash rather than umber
    latg = lat[:, None]
    lng = lon[None, :]
    ice = land & ((latg < -60) | ((latg > 60) & (lng > -75) & (lng < -10) & (z > 900)))
    gnd = np.where(ice[..., None], gnd * 0.35 + ramp(z, [(0, '#5e5852'), (4000, '#8f877e')]) * 0.65, gnd)

    s_sea = np.clip(1.0 + 0.75 * (shade - 1.0), 0.40, 1.55)
    s_gnd = np.clip(0.50 + 0.62 * (shade - 1.0) + 0.50, 0.35, 1.45)
    img = np.where(land[..., None], gnd * s_gnd[..., None], sea * s_sea[..., None])

    # faint warm coastline
    edge = land ^ ndimage.binary_erosion(land, iterations=1)
    img[edge] = img[edge] * 0.6 + np.array([150, 112, 80], np.float32) * 0.4
    img = np.clip(img, 0, 255).astype(np.uint8)
    os.makedirs(out, exist_ok=True)
    im = Image.fromarray(img)
    im.save(os.path.join(out, 'relief_8k.jpg'), quality=84, optimize=True, progressive=True)
    im.resize((4096, 2048), Image.LANCZOS).save(os.path.join(out, 'relief_4k.jpg'), quality=86, optimize=True,
                                                progressive=True)
    im.resize((1600, 800), Image.LANCZOS).save(os.path.join(out, 'relief_preview.jpg'), quality=88)
    for f in ('relief_8k.jpg', 'relief_4k.jpg'):
        print(f, os.path.getsize(os.path.join(out, f)) / 1e6, 'MB')


if __name__ == '__main__':
    main()
