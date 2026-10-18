"""Polar → grid, hybrid scans and multi-radar RQI mosaics.

Products (all on the 0.01° grid, 15-min depth and hourly accumulation):
  LS_MP   lowest sweep, Z–R Marshall–Palmer, beams above the melting-layer bottom discarded (baseline)
  LS_DP   lowest sweep, dual-pol blended estimate
  HYZ_MP  reflectivity hybrid scan (Fulton et al. 1998): for every polar bin the lowest tilt with
          CBB < 0.5 and beam bottom ≥ 150 m above terrain; partial-blockage correction
          ΔZ = −10 log10(1 − CBB); VPR correction above the melting layer
  PHY_DP  polarimetric hybrid scan: R(KDP)/R(A) from the lowest tilt with CBB < 0.7 that is below the
          melting layer and has a phase estimate (phase is immune to partial blockage); HYZ elsewhere
Mosaic:  R = Σ RQIᵢ Rᵢ / Σ RQIᵢ  over radars i with RQIᵢ > 0 (weighted) — or the radar with max RQI.
"""
from __future__ import annotations

import time as _time
from collections import defaultdict

import numpy as np
import pandas as pd
from pyproj import Geod

from . import geometry as geo
from . import qc as QC
from . import rqi as RQ
from . import vpr as VP
from .io import read_radar

GEOD = Geod(ellps="WGS84")
AZ_GRID = np.arange(0.25, 360, 0.5)
PRODUCTS = ["LS_MP", "LS_DP", "HYZ_MP", "PHY_DP"]


def make_grid(lon_min, lon_max, lat_min, lat_max, res):
    """Cell centres; latitude decreasing (row 0 = north) like a GeoTIFF."""
    glon = np.arange(lon_min + res / 2, lon_max, res)
    glat = np.arange(lat_max - res / 2, lat_min, -res)
    return glon, glat


class PolarToGrid:
    """Nearest-bin lookup from a lon/lat grid to an (azimuth × range) polar grid. Build once per radar."""

    def __init__(self, lon0, lat0, az, r, glon, glat, max_range=240e3):
        LON, LAT = np.meshgrid(glon, glat)
        faz, _, dist = GEOD.inv(np.full(LON.shape, lon0), np.full(LON.shape, lat0), LON, LAT)
        faz = np.mod(faz, 360.0)
        az = np.asarray(az, float); order = np.argsort(az); azs = az[order]
        i1 = np.searchsorted(azs, faz) % len(azs); i0 = (i1 - 1) % len(azs)
        d1 = np.abs((azs[i1] - faz + 180) % 360 - 180); d0 = np.abs((azs[i0] - faz + 180) % 360 - 180)
        self.ia = order[np.where(d1 < d0, i1, i0)]
        gate = float(np.median(np.diff(r)))
        ir = np.round((dist - r[0]) / gate).astype(int)
        self.inside = (ir >= 0) & (ir < len(r)) & (dist <= max_range)
        self.ir = np.clip(ir, 0, len(r) - 1)
        self.dist = dist

    def __call__(self, a, fill=np.nan):
        return np.where(self.inside, a[self.ia, self.ir], fill)


def nearest_az_index(az_sweep, az_grid=AZ_GRID):
    return np.argmin(np.abs((np.asarray(az_sweep)[None, :] - az_grid[:, None] + 180) % 360 - 180), axis=1)


def radar_geometry(radar, terrain_fn, beamwidth=1.0, az_grid=AZ_GRID):
    """Geometry of every sweep on the common azimuth grid (list of dicts, see geometry.sweep_geometry)."""
    out = []
    lon0, lat0 = float(radar.longitude["data"][0]), float(radar.latitude["data"][0])
    alt = float(radar.altitude["data"][0])
    r = radar.range["data"].astype(float)
    for k in range(radar.nsweeps):
        el = float(radar.fixed_angle["data"][k])
        g = geo.sweep_geometry(lon0, lat0, alt, el, az_grid, r, terrain_fn, beamwidth)
        out.append(g)
    return out


def stack(radar, outs, keys, az_grid=AZ_GRID):
    """Put per-sweep 2-D fields onto (n_sweeps × n_az_grid × n_gates)."""
    res = {}
    for key in keys:
        a = []
        for k, o in enumerate(outs):
            az_s = radar.extract_sweeps([k]).azimuth["data"].astype(float)
            a.append(np.asarray(o[key], dtype="float32")[nearest_az_index(az_s, az_grid)])
        res[key] = np.array(a)
    return res


def process_hour(scans: pd.DataFrame, terrain_fn, calib: pd.DataFrame, glon, glat, fzl_fallback=4800.0,
                 thr_z=0.5, thr_p=0.7, clearance=150.0, max_range=240e3, verbose=True, method="weighted"):
    """Full v2 pipeline for one set of scans (normally 5 scans = one hour + the next :00 scan).

    scans : DataFrame from io.list_scans (radar, time, path)
    calib : DataFrame from qc.calibrate (index = radar; zdr_offset_db, system_phase_deg, phidp_sign)
    Returns dict with 15-min depths m15_<product> (n_times × ny × nx, mm), hourly h_<product>,
    rqi_<LS|HYZ|PHY>, melting layer, VPR profile and per-radar diagnostics.
    """
    t0 = _time.time()
    shape = (len(glat), len(glon))
    times = sorted(scans.time.unique())
    radars = sorted(scans.radar.unique())
    cache, geoms, p2g = {}, {}, {}
    hs, zs, rs = [], [], []
    # ---------------------------------------------------------------- pass 1: QC + polarimetric processing
    for t in times:
        for _, row in scans[scans.time == t].iterrows():
            rad = row.radar
            radar = read_radar(row.path)
            key = (rad, radar.nsweeps, radar.ngates)
            if key not in geoms:
                geoms[key] = radar_geometry(radar, terrain_fn)
            G = geoms[key]
            c = calib.loc[rad]
            outs = [QC.process_sweep(radar, k, c.zdr_offset_db, c.system_phase_deg, c.phidp_sign, fzl_fallback)[1]
                    for k in range(radar.nsweeps)]
            st = stack(radar, outs, ["Zc", "R_KDP", "R_A", "src", "RHOHV", "noecho", "tex"])
            H = np.array([g["h"] for g in G]); C = np.array([g["cbb"] for g in G])
            cache[(t, rad)] = (st, key)
            r2d = np.broadcast_to(G[0]["r"][None, None, :], H.shape)
            sel = np.isfinite(st["Zc"]) & (st["Zc"] > 10) & (r2d > 10e3) & (r2d < 120e3) & (C < 0.3)
            hs.append(H[sel]); zs.append(st["Zc"][sel]); rs.append(st["RHOHV"][sel])
        if verbose:
            print(f"  pass 1 {pd.Timestamp(t):%H:%M}  ({_time.time() - t0:.0f} s)")
    # ---------------------------------------------------------------- melting layer + VPR
    prof = VP.profile_from_samples(np.concatenate(hs), np.concatenate(zs), np.concatenate(rs))
    ml_bot, ml_top, ml_src = VP.melting_layer(prof, fallback_fzl=fzl_fallback)
    hg, dz, zref = VP.correction_table(prof, ml_bot)
    if verbose:
        print(f"  melting layer {ml_bot:.0f}–{ml_top:.0f} m ({ml_src}); Z_ref {zref:.1f} dBZ")
    # ---------------------------------------------------------------- pass 2: hybrid scans + mosaic
    m15 = {p: [] for p in PRODUCTS}
    rqi_max = {p: np.zeros(shape) for p in ["LS", "HYZ", "PHY"]}
    best = {p: np.full(shape, -1, int) for p in ["LS", "HYZ", "PHY"]}
    diag = []
    pick = lambda A_, k: np.take_along_axis(A_, k, 0)[0]
    for t in times:
        num = {p: np.zeros(shape) for p in PRODUCTS}; den = {p: np.zeros(shape) for p in PRODUCTS}
        qbest = {p: np.zeros(shape) for p in PRODUCTS}; vbest = {p: np.full(shape, np.nan) for p in PRODUCTS}
        for ri, rad in enumerate(radars):
            if (t, rad) not in cache:
                continue
            st, key = cache[(t, rad)]
            G = geoms[key]
            H = np.array([g["h"] for g in G]); T = np.array([g["terr"] for g in G])
            C = np.array([g["cbb"] for g in G]); CL = np.array([g["clearance"] for g in G])
            hagl = H - T
            hgt = RQ.rqi_height(hagl, h_asl_m=H, ml_bottom_m=ml_bot)
            noe = st["noecho"] > 0.5
            Zcorr = st["Zc"] + np.where(C < thr_z, -10 * np.log10(np.clip(1 - C, 0.05, 1)), np.nan) \
                + VP.apply(0.0, H, hg, dz)
            from .qpe import r_from_z
            Rz = np.where(noe, 0.0, r_from_z(Zcorr))
            Rz = np.where(noe & (H >= ml_bot), np.nan, Rz)      # no echo above the ML = unknown (overshooting)
            Rz_raw = np.where(noe, 0.0, r_from_z(st["Zc"]))
            Rp = np.where(st["src"] == 2, st["R_KDP"], np.where(st["src"] == 1, st["R_A"], np.nan))
            # LS
            q_ls = RQ.rqi_blockage(C[0]) * RQ.rqi_height(hagl[0]); q_ls[H[0] >= ml_bot] = 0.0
            dp_ls = np.where(np.isfinite(Rp[0]), Rp[0], Rz_raw[0])
            # HYZ
            okz = (C < thr_z) & (CL >= clearance)
            kz = np.where(okz.any(0), okz.argmax(0), -1); kzc = np.clip(kz, 0, None)[None]
            q_z = RQ.rqi_blockage(pick(C, kzc)) * pick(hgt, kzc); q_z[kz < 0] = 0.0
            rz = pick(Rz, kzc)
            # PHY
            okp = (C < thr_p) & (CL >= clearance) & (H < ml_bot) & np.isfinite(Rp) & (Rp > 0)
            kp = np.where(okp.any(0), okp.argmax(0), -1); kpc = np.clip(kp, 0, None)[None]
            use_p = kp >= 0
            q_p = np.where(use_p, RQ.rqi_blockage(pick(C, kpc), 0.3, thr_p) * pick(hgt, kpc)
                           * RQ.rqi_polarimetric(pick(st["tex"], kpc)), q_z)
            r_p = np.where(use_p, pick(Rp, kpc), rz)
            wet = np.isfinite(rz) & (rz > 0.5)
            diag.append(dict(time=t, radar=rad, frac_phase_bins=float(np.mean(use_p[wet])) if wet.any() else np.nan,
                             frac_hyz_higher_tilt=float(np.mean(kz > 0)), rain_bins=int(wet.sum())))
            if rad not in p2g:
                p2g[rad] = PolarToGrid(G[0]["lon0"], G[0]["lat0"], AZ_GRID, G[0]["r"], glon, glat, max_range)
            M = p2g[rad]
            Q = {"LS": M(q_ls, 0.0), "HYZ": M(q_z, 0.0), "PHY": M(q_p, 0.0)}
            V = {"LS_MP": M(Rz_raw[0]), "LS_DP": M(dp_ls), "HYZ_MP": M(rz), "PHY_DP": M(r_p)}
            for p, qk in [("LS_MP", "LS"), ("LS_DP", "LS"), ("HYZ_MP", "HYZ"), ("PHY_DP", "PHY")]:
                q = np.nan_to_num(Q[qk]); v = V[p]; ok = np.isfinite(v) & (q > 0)
                num[p][ok] += q[ok] * v[ok]; den[p][ok] += q[ok]
                b = ok & (q > qbest[p]); qbest[p][b] = q[b]; vbest[p][b] = v[b]
            for nm in ["LS", "HYZ", "PHY"]:
                q = np.nan_to_num(Q[nm]); b = q > rqi_max[nm]
                rqi_max[nm][b] = q[b]; best[nm][b] = ri
        for p in PRODUCTS:
            rate = np.where(den[p] > 0, num[p] / np.maximum(den[p], 1e-9), np.nan) if method == "weighted" else vbest[p]
            m15[p].append(rate * 0.25)                        # mm in 15 min (rate held constant)
        if verbose:
            print(f"  pass 2 {pd.Timestamp(t):%H:%M}  ({_time.time() - t0:.0f} s)")
    out = dict(glon=glon, glat=glat, times=np.array([str(pd.Timestamp(t)) for t in times]), radars=np.array(radars),
               ml_bottom=ml_bot, ml_top=ml_top, ml_source=ml_src, vpr_h=hg, vpr_dz=dz, z_ref=zref)
    for p in PRODUCTS:
        s_ = np.array(m15[p])
        out[f"m15_{p}"] = s_.astype("float32")
        out[f"h_{p}"] = hourly_from_15min(s_).astype("float32")
    for k in rqi_max:
        out[f"rqi_{k}"] = rqi_max[k].astype("float32")
        out[f"best_{k}"] = best[k]
    out["profile"] = prof.reset_index()
    out["diagnostics"] = pd.DataFrame(diag)
    if verbose:
        print(f"  done in {_time.time() - t0:.0f} s")
    return out


def hourly_from_15min(m15, min_valid=3):
    """Hourly depth from 15-min depths. With 5 scans (hh:00 … hh+1:00) the first four represent the hour.
    A cell needs ≥ min_valid of 4 values; missing ones are filled with the mean of the others."""
    s_ = np.asarray(m15)
    sh = s_[:4] if len(s_) == 5 else s_
    n = np.sum(np.isfinite(sh), 0)
    return np.where(n >= min_valid, np.nansum(sh, 0) * sh.shape[0] / np.maximum(n, 1), np.nan)


def save_npz(path, res):
    keep = {k: v for k, v in res.items() if isinstance(v, (np.ndarray, float, int, str))}
    np.savez_compressed(path, **keep)
    if "profile" in res:
        res["profile"].to_csv(str(path).replace(".npz", "_vpr_profile.csv"), index=False)
    if "diagnostics" in res:
        res["diagnostics"].to_csv(str(path).replace(".npz", "_diagnostics.csv"), index=False)


def to_geotiff(path, arr, glon, glat, nodata=np.nan):
    import rasterio
    from rasterio.transform import from_origin
    res = abs(glon[1] - glon[0])
    tr = from_origin(glon[0] - res / 2, glat[0] + res / 2, res, res)
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
                       dtype="float32", crs="EPSG:4326", transform=tr, nodata=nodata, compress="deflate") as dst:
        dst.write(np.asarray(arr, "float32"), 1)
