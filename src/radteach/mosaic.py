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
                 thr_z=0.5, thr_p=0.7, clearance=150.0, max_range=240e3, verbose=True, method="weighted",
                 keep_radars=False, z_offset_db=None, cappi_m=2000.0, cappi_zr=(250.0, 1.2), cappi_pseudo=True, keep_sweeps=False):
    """Full v2 pipeline for one set of scans (normally 5 scans = one hour + the next :00 scan).

    scans : DataFrame from io.list_scans (radar, time, path)
    calib : DataFrame from qc.calibrate (index = radar; zdr_offset_db, system_phase_deg, phidp_sign)
    Returns dict with 15-min depths m15_<product> (n_times × ny × nx, mm), hourly h_<product>,
    rqi_<LS|HYZ|PHY>, melting layer, VPR profile and per-radar diagnostics.
    keep_radars=True also returns the single-radar fields on the grid (for comparing mosaic methods, nb 09):
        radar_<product>  (n_radars × n_times × ny × nx) rain RATE (mm/h), NaN where the radar does not see the cell
        radar_rqi_<LS|HYZ|PHY> (n_radars × n_times × ny × nx) RQI of that radar,  radar_dist (n_radars × ny × nx) km
        radar_CAPPI_RT   rain rate from a CAPPI at `cappi_m` m ASL with Z = 250 R^1.2 (Rosenfeld tropical), as in
                         Mahavik et al. (2025, Spatial Information Research 33:6); quality radar_rqi_CAP = RQI_blk
                         (pseudo-CAPPI by default: the lowest/highest sweep is used below/above the beams)
    keep_sweeps=True also returns every sweep on the grid (inter-radar calibration, notebook 10):
        radar_dbz_SW (n_radars × n_times × n_sweeps × ny × nx) attenuation-corrected Z, NaN = no echo
        radar_hgt_SW, radar_cbb_SW (n_radars × n_sweeps × ny × nx) beam-centre height ASL (m) and CBB
    z_offset_db : {radar: dB} subtracted from Z before any rain estimate (inter-radar calibration, notebook 10).
                  Phase-based estimators (R(KDP)) are unaffected; R(A) from ZPHI only weakly.
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
            if z_offset_db:
                st["Zc"] = st["Zc"] - float(z_offset_db.get(rad, 0.0))
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
    if keep_radars:
        nr, nt = len(radars), len(times)
        RF = {p: np.full((nr, nt) + shape, np.nan, "float32") for p in PRODUCTS}
        RQf = {p: np.zeros((nr, nt) + shape, "float32") for p in ["LS", "HYZ", "PHY", "CAP"]}
        RD = np.full((nr,) + shape, np.nan, "float32")
        RF["CAPPI_RT"] = np.full((nr, nt) + shape, np.nan, "float32")
    if keep_sweeps:
        nr, nt, nsw = len(radars), len(times), 6
        RZS = np.full((nr, nt, nsw) + shape, np.nan, "float32")
        RHS = np.full((nr, nsw) + shape, np.nan, "float32"); RCS = np.full((nr, nsw) + shape, np.nan, "float32")
    for ti, t in enumerate(times):
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
            if keep_radars:
                for p in PRODUCTS:
                    RF[p][ri, ti] = V[p]
                for nm in Q:
                    RQf[nm][ri, ti] = np.nan_to_num(Q[nm])
                RD[ri] = np.where(M.inside, M.dist / 1e3, np.nan)
                r_cap, q_cap = cappi_rain(st["Zc"], noe, H, C, cappi_m, *cappi_zr, pseudo=cappi_pseudo)
                RF["CAPPI_RT"][ri, ti] = M(r_cap); RQf["CAP"][ri, ti] = np.nan_to_num(M(q_cap, 0.0))
            if keep_sweeps:
                for k in range(min(H.shape[0], nsw)):
                    RZS[ri, ti, k] = M(np.where(noe[k], np.nan, st["Zc"][k]))
                    if not np.isfinite(RHS[ri, k]).any():
                        RHS[ri, k] = M(H[k]); RCS[ri, k] = M(C[k])
        for p in PRODUCTS:
            rate = np.divide(num[p], den[p], out=np.full(shape, np.nan), where=den[p] > 0) if method == "weighted" else vbest[p]
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
    if keep_radars:
        for p in RF:
            out[f"radar_{p}"] = RF[p]
        for nm in RQf:
            out[f"radar_rqi_{nm}"] = RQf[nm]
        out["radar_dist"] = RD
    if keep_sweeps:
        out["radar_dbz_SW"] = RZS; out["radar_hgt_SW"] = RHS; out["radar_cbb_SW"] = RCS
    out["profile"] = prof.reset_index()
    out["diagnostics"] = pd.DataFrame(diag)
    if verbose:
        print(f"  done in {_time.time() - t0:.0f} s")
    return out


def cappi_dbz(Z, noe, H, C, h_c=2000.0, pseudo=False):
    """Reflectivity on a constant-altitude PPI at h_c m ASL: linear interpolation in dBZ between the two sweeps
    that bracket h_c. pseudo=False: 'true' CAPPI, not valid below the lowest or above the highest beam;
    pseudo=True: the lowest (highest) sweep is used where h_c is below (above) all beams.
    Arrays (n_sweeps, n_az, n_gates). Returns (dBZ, valid, CBB of the more blocked bracket, no-echo in both)."""
    pick = lambda A_, k: np.take_along_axis(A_, k, 0)[0]
    order = np.argsort(np.nanmean(H, axis=(1, 2)))
    Z, noe, H, C = Z[order], noe[order], H[order], C[order]
    k = (H <= h_c).sum(0) - 1
    valid = (k >= 0) & (k < H.shape[0] - 1)
    if H.shape[0] == 1:
        z = Z[0]; return z, np.ones_like(valid) if pseudo else np.zeros_like(valid), C[0], noe[0] > 0.5
    kc = np.clip(k, 0, H.shape[0] - 2)[None]
    h0, h1 = pick(H, kc), pick(H, kc + 1)
    z0, z1 = pick(Z, kc), pick(Z, kc + 1)
    w = np.clip((h_c - h0) / np.where(h1 > h0, h1 - h0, 1.0), 0, 1)
    if pseudo:
        valid = np.ones_like(valid)                      # w is clipped: 0 below the beams, 1 above
    zi = np.where(np.isfinite(z0) & np.isfinite(z1), z0 * (1 - w) + z1 * w, np.where(np.isfinite(z0), z0, z1))
    both_no = (pick(noe, kc) > 0.5) & (pick(noe, kc + 1) > 0.5)
    return zi, valid, np.maximum(pick(C, kc), pick(C, kc + 1)), both_no


def cappi_rain(Z, noe, H, C, h_c=2000.0, a=250.0, b=1.2, pseudo=True):
    """Rain rate (mm/h) from a CAPPI at h_c m ASL with Z = a R^b (default Rosenfeld tropical, as in
    Mahavik et al. 2025). Returns (rate, quality = RQI_blk of the more blocked bracketing sweep)."""
    from .qpe import r_from_z
    zi, valid, cbb, both_no = cappi_dbz(Z, noe, H, C, h_c, pseudo)
    R = np.where(both_no, 0.0, r_from_z(zi, a=a, b=b))
    R = np.where(valid, R, np.nan)
    return R, np.where(valid, RQ.rqi_blockage(cbb), 0.0)


COMBINE_METHODS = ["MAX", "MEAN", "NEAREST", "RQI_MAX", "RQI_WMEAN"]


def combine(R, Q=None, dist=None, method="RQI_WMEAN", k=1.0):
    """Combine single-radar fields R (n_radars × ny × nx; NaN = not seen) into one mosaic.

    MAX        max_i R_i                         (any radar that sees the cell)
    MEAN       mean_i R_i                        (any radar that sees the cell)
    NEAREST    R of the radar closest to the cell (dist, km)
    RQI_MAX    R of the radar with the highest RQI (RQI > 0)
    RQI_WMEAN  Σ RQI_i^k R_i / Σ RQI_i^k        (RQI > 0; k = 1 is MRMS, k → ∞ tends to RQI_MAX)
    The simple methods (MAX, MEAN, NEAREST) ignore quality on purpose: that is what we want to test.
    """
    R = np.asarray(R, float)
    seen = np.isfinite(R)
    nan = np.full(R.shape[1:], np.nan)
    if method == "MAX":
        with np.errstate(all="ignore"):
            return np.where(seen.any(0), np.nanmax(np.where(seen, R, -np.inf), 0), nan)
    if method == "MEAN":
        n = seen.sum(0)
        return np.where(n > 0, np.where(seen, R, 0).sum(0) / np.maximum(n, 1), nan)
    if method == "NEAREST":
        d = np.where(seen, np.nan_to_num(np.asarray(dist, float), nan=np.inf), np.inf)
        i = np.argmin(d, 0)
        return np.where(seen.any(0), np.take_along_axis(R, i[None], 0)[0], nan)
    Q = np.where(seen, np.nan_to_num(np.asarray(Q, float)), 0.0)
    if method == "RQI_MAX":
        i = np.argmax(Q, 0)
        return np.where(Q.max(0) > 0, np.take_along_axis(R, i[None], 0)[0], nan)
    if method == "RQI_WMEAN":
        w = np.where(Q > 0, Q ** k, 0.0)
        sw = w.sum(0)
        return np.divide((w * np.where(seen, R, 0)).sum(0), sw, out=nan, where=sw > 0)
    raise ValueError(f"unknown method {method}")


def combine_hour(res, product="PHY_DP", method="RQI_WMEAN", k=1.0):
    """Hourly depth (mm) of a mosaic method from process_hour(..., keep_radars=True)."""
    grp = {"LS_MP": "LS", "LS_DP": "LS", "HYZ_MP": "HYZ", "PHY_DP": "PHY", "CAPPI_RT": "CAP"}[product]
    R, Q, D = res[f"radar_{product}"], res[f"radar_rqi_{grp}"], res["radar_dist"]
    m15 = np.array([combine(R[:, t], Q[:, t], D, method, k) * 0.25 for t in range(R.shape[1])])
    return hourly_from_15min(m15)


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
