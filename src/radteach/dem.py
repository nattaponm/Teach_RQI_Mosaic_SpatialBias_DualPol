"""DEM: mosaic tiles, crop/resample, hillshade and a fast terrain interpolator."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.merge import merge
from scipy.interpolate import RegularGridInterpolator


def mosaic_tiles(tiles, bounds=None, res=None, out_path=None, resampling="average", nodata=-32768):
    """Merge DEM tiles into one array.

    tiles  : list of GeoTIFF paths (same CRS)
    bounds : (lon_min, lat_min, lon_max, lat_max) to crop, or None for the union of the tiles
    res    : output pixel size in degrees (None = keep native)
    Returns (elevation array [m] with NaN for no-data, affine transform, profile)
    """
    srcs = [rasterio.open(t) for t in tiles]
    arr, tr = merge(srcs, bounds=bounds, res=res, nodata=nodata,
                    resampling=getattr(Resampling, resampling))
    prof = srcs[0].profile.copy()
    for s in srcs:
        s.close()
    z = arr[0].astype("float32")
    z[z == nodata] = np.nan
    z[z < -500] = np.nan
    prof.update(height=z.shape[0], width=z.shape[1], transform=tr, count=1, dtype="float32", nodata=np.nan,
                compress="deflate", tiled=True, blockxsize=256, blockysize=256)
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(out_path, "w", **prof) as dst:
            dst.write(z, 1)
    return z, tr, prof


def read(path):
    with rasterio.open(path) as s:
        z = s.read(1).astype("float32")
        if s.nodata is not None and np.isfinite(s.nodata):
            z[z == s.nodata] = np.nan
        return z, s.transform, s.profile


def extent(z, tr):
    """[left, right, bottom, top] for imshow."""
    return [tr.c, tr.c + tr.a * z.shape[1], tr.f + tr.e * z.shape[0], tr.f]


def coords(z, tr):
    lon = tr.c + tr.a * (np.arange(z.shape[1]) + 0.5)
    lat = tr.f + tr.e * (np.arange(z.shape[0]) + 0.5)
    return lon, lat


def interpolator(z, tr, sea_level_fill=0.0):
    """Function f(lat, lon) -> terrain height (m); bilinear, NaN (sea/no data) replaced by 0 m."""
    lon, lat = coords(z, tr)
    zz = np.where(np.isfinite(z), z, sea_level_fill)
    if lat[0] > lat[-1]:
        lat, zz = lat[::-1], zz[::-1]
    f = RegularGridInterpolator((lat, lon), zz, bounds_error=False, fill_value=sea_level_fill)

    def terrain(lat_, lon_):
        lat_ = np.asarray(lat_, float); lon_ = np.asarray(lon_, float)
        return f(np.column_stack([lat_.ravel(), lon_.ravel()])).reshape(lat_.shape)
    return terrain


def slope_deg(z, tr):
    """Terrain slope (degrees) from a geographic DEM (pixel size converted to metres)."""
    lon, lat = coords(z, tr)
    dy = abs(tr.e) * 111_320.0
    dx = tr.a * 111_320.0 * np.cos(np.radians(lat))[:, None]
    gy, gx = np.gradient(np.nan_to_num(z))
    return np.degrees(np.arctan(np.hypot(gx / dx, gy / dy)))
