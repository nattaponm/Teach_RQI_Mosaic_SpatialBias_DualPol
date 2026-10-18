"""Standard maps: basin boundaries, hillshade, radars, graticule, scale bar, north arrow, inset, source note.

Every map in the course is built with `basemap(ax, extent)` so that all figures share the same look.
Coordinates are geographic (EPSG:4326); the y/x aspect is set to 1/cos(latitude) so that distances
look correct (equivalent to a local equirectangular projection).
"""
from __future__ import annotations

from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from pyproj import Geod

from . import plotting as P

GEOD = Geod(ellps="WGS84")
RADAR_COLORS = {"CRI": "#e67e22", "LMP": "#c0392b", "MHS": "#8e44ad", "PHS": "#27ae60", "TAK": "#2e86c1"}


@lru_cache(maxsize=4)
def boundaries(root=None):
    """(main basins, sub-basins, Thailand outline) as GeoDataFrames in EPSG:4326."""
    import geopandas as gpd
    from .config import paths
    p = paths(root).boundaries
    mb = gpd.read_file(p / "onwr_main_basins.gpkg")
    sb = gpd.read_file(p / "onwr_sub_basins.gpkg")
    th = gpd.read_file(p / "thailand_outline.gpkg")
    return mb, sb, th


def radar_sites(root=None):
    import pandas as pd
    from .config import paths
    return pd.read_csv(paths(root).static / "radar_sites.csv").set_index("radar")


def set_extent(ax, extent):
    x0, x1, y0, y1 = extent
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect(1 / np.cos(np.radians((y0 + y1) / 2)))


def graticule(ax, step=1.0, labels=True, grid=True):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    ax.set_xticks(np.arange(np.ceil(x0 / step) * step, x1 + 1e-9, step))
    ax.set_yticks(np.arange(np.ceil(y0 / step) * step, y1 + 1e-9, step))
    fmt = "{:.0f}" if step >= 1 else "{:.1f}"
    if labels:
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: fmt.format(v) + "°E"))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: fmt.format(v) + "°N"))
    else:
        ax.set_xticklabels([]); ax.set_yticklabels([])
    if grid:
        ax.grid(ls=":", lw=0.4, color="0.45", zorder=0.5)
    ax.set_xlabel(""); ax.set_ylabel("")


def scalebar(ax, km=50, loc=(0.05, 0.05), lw=3):
    """Scale bar in kilometres at axes-fraction position loc (lower-left corner)."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    lon = x0 + loc[0] * (x1 - x0); lat = y0 + loc[1] * (y1 - y0)
    lon2, _, _ = GEOD.fwd(lon, lat, 90, km * 1e3)
    half = (lon + lon2) / 2
    ax.plot([lon, half], [lat, lat], color="k", lw=lw, solid_capstyle="butt", zorder=20)
    ax.plot([half, lon2], [lat, lat], color="w", lw=lw, solid_capstyle="butt", zorder=20,
            path_effects=None)
    ax.plot([lon, lon2], [lat, lat], color="k", lw=lw + 1.2, solid_capstyle="butt", zorder=19)
    ax.text((lon + lon2) / 2, lat + 0.012 * (y1 - y0), f"{km} km", ha="center", va="bottom",
            fontsize=plt.rcParams["font.size"] * 0.85, zorder=21,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))


def north_arrow(ax, loc=(0.94, 0.90), size=0.07):
    ax.annotate("N", xy=(loc[0], loc[1]), xytext=(loc[0], loc[1] - size), xycoords="axes fraction",
                ha="center", va="center", fontweight="bold", fontsize=plt.rcParams["font.size"] * 1.1,
                arrowprops=dict(facecolor="k", edgecolor="k", width=2.5, headwidth=8, headlength=7), zorder=21)


def source_note(ax, text, loc=(0.99, 0.01)):
    ax.text(loc[0], loc[1], text, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=plt.rcParams["font.size"] * 0.7, color="0.3", zorder=21,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))


def inset_thailand(fig, rect, extent):
    """Small locator map of Thailand with a red box for `extent`. rect = [x, y, w, h] in figure fraction."""
    mb, _, th = boundaries()
    ins = fig.add_axes(rect)
    th.plot(ax=ins, color="0.82", lw=0)
    th.boundary.plot(ax=ins, color="0.4", lw=0.3)
    x0, x1, y0, y1 = extent
    ins.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="red", lw=1))
    ins.set_xlim(97.2, 106); ins.set_ylim(5.5, 20.6)
    ins.set_aspect(1 / np.cos(np.radians(13)))
    ins.set_xticks([]); ins.set_yticks([]); ins.set_xlabel(""); ins.set_ylabel("")
    ins.set_facecolor("white")
    for s in ins.spines.values():
        s.set_lw(0.6)
    return ins


def range_ring(lon, lat, km, n=361):
    az = np.linspace(0, 360, n)
    lo, la, _ = GEOD.fwd(np.full(n, lon), np.full(n, lat), az, np.full(n, km * 1e3))
    return lo, la


def plot_radars(ax, radars=None, rings_km=None, labels=True, ms=7):
    sites = radar_sites()
    for name, r in sites.iterrows():
        if radars is not None and name not in radars:
            continue
        c = RADAR_COLORS.get(name, "k")
        ax.plot(r.lon, r.lat, "^", color=c, mec="k", mew=0.6, ms=ms, zorder=15)
        if rings_km:
            for km in np.atleast_1d(rings_km):
                x, y = range_ring(r.lon, r.lat, km)
                ax.plot(x, y, color=c, lw=0.7, ls="--", zorder=14)
        if labels:
            ax.annotate(name, (r.lon, r.lat), xytext=(4, 4), textcoords="offset points", fontweight="bold",
                        fontsize=plt.rcParams["font.size"] * 0.85, zorder=16,
                        bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75))


def hillshade(ax, dem, extent_dem, alpha=0.35, vert_exag=0.05, dx_m=275.0):
    from matplotlib.colors import LightSource
    ls = LightSource(azdeg=315, altdeg=40)
    hs = ls.hillshade(np.nan_to_num(dem), vert_exag=vert_exag, dx=dx_m, dy=dx_m)
    ax.imshow(hs, extent=extent_dem, cmap="gray", alpha=alpha, zorder=0, origin="upper")


def basemap(ax, extent, basins=("01", "02", "06", "07", "08", "09"), subbasins=None, all_basins=True,
            radars=True, rings_km=None, grid_step=1.0, scale_km=50, north=True, title=None, lw=0.9):
    """Draw the standard background: basins, optional sub-basins, radars, graticule, scale bar, north arrow."""
    mb, sb, _ = boundaries()
    if all_basins:
        mb.boundary.plot(ax=ax, color="0.6", lw=0.35, zorder=5)
    if basins:
        mb[mb.MB_CODE.isin(basins)].boundary.plot(ax=ax, color="k", lw=lw, zorder=6)
    if subbasins is not None:
        sb[sb.SB_CODE.isin(subbasins)].boundary.plot(ax=ax, color="k", lw=0.7, ls="--", zorder=6)
    if radars:
        plot_radars(ax, None if radars is True else radars, rings_km)
    set_extent(ax, extent)
    graticule(ax, grid_step)
    if scale_km:
        scalebar(ax, scale_km)
    if north:
        north_arrow(ax)
    if title:
        ax.set_title(title, loc="left", fontweight="bold")
    return ax


def grid_extent(glon, glat):
    r = abs(glon[1] - glon[0])
    return [glon[0] - r / 2, glon[-1] + r / 2, min(glat) - r / 2, max(glat) + r / 2]


def show_grid(ax, arr, glon, glat, cmap, norm, zorder=2, alpha=1.0):
    """imshow a lon/lat grid whose first row is the NORTHERN edge (glat decreasing)."""
    origin = "upper" if glat[0] > glat[-1] else "lower"
    return ax.imshow(arr, extent=grid_extent(glon, glat), origin=origin, cmap=cmap, norm=norm,
                     zorder=zorder, alpha=alpha, interpolation="nearest")


def sweep_lonlat(sweep):
    """lon, lat of every gate of a one-sweep Py-ART radar object."""
    from .geometry import polar_lonlat
    return polar_lonlat(float(sweep.longitude["data"][0]), float(sweep.latitude["data"][0]),
                        sweep.azimuth["data"].astype(float), sweep.range["data"].astype(float),
                        float(sweep.fixed_angle["data"][0]))


def ppi(ax, sweep, values, cmap, norm, zorder=2, lonlat=None):
    """Plot one sweep (azimuth × range array) on a geographic map."""
    lon, lat = lonlat if lonlat is not None else sweep_lonlat(sweep)
    o = np.argsort(sweep.azimuth["data"])
    return ax.pcolormesh(lon[o], lat[o], np.ma.masked_invalid(np.asarray(values, float)[o]), cmap=cmap, norm=norm,
                         shading="nearest", zorder=zorder, rasterized=True)
