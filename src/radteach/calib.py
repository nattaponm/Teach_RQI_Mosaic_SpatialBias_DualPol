"""Relative calibration of reflectivity between neighbouring radars from their overlap (notebook 10).

Method (adapted to a common Cartesian grid):
    Gourley, J.J., Kaney, B., Maddox, R.A. (2003) Evaluating the calibrations of radars: a software approach.
        31st Conf. Radar Meteorology, AMS, P3C.1.  -> matched points at (nearly) equal distance and equal beam
        height, Z >= 15 dBZ, |dZ| <= 8 dB, and "loop closure" of the differences around radar triangles.
    Crisologo, I., Heistermann, M. (2020) Atmos. Meas. Tech. 13:645-659 -> equidistance band (10 km),
        beam-blockage quality filter.
    Hu, H. et al. (2025) Atmos. Meas. Tech. 18:5507-5526 -> 15-35 dBZ window, 3-sigma outlier removal,
        minimum sample size, per-radar bias from all neighbours.
Our addition: a weighted least-squares solution of  dZ_ij = c_i - c_j  for the whole network with sum(c) = 0,
plus a bootstrap over scans. The result is RELATIVE (who is hot / cold); the absolute level is left to the
gauge-based mean-field bias or to a spaceborne radar (GPM).
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

DEFAULTS = dict(max_dh_m=250.0, max_dd_km=30.0, max_cbb=0.2, zmin=15.0, zmax=45.0, max_absdiff=8.0,
                ml_margin_m=300.0, beamwidth_deg=1.0, min_range_km=20.0)


def overlap_pairs(res, **kw) -> pd.DataFrame:
    """Matched grid cells between every pair of radars, all scans of one process_hour result.

    Needs process_hour(..., keep_radars=True, keep_sweeps=True). Like Gourley et al. (2003) any elevation of
    one radar may be paired with any elevation of the other: for every cell the sweep pair with the smallest
    beam-height difference is chosen, so radars at very different site altitudes (TAK 943 m, PHS 72 m) can
    still be compared at the same height. Criteria (DEFAULTS): |h_i - h_j| <= max_dh_m,
    |d_i - d_j| <= max_dd_km (equidistance band = similar beam volume), both ranges >= min_range_km,
    CBB <= max_cbb for both sweeps, beam tops below the melting-layer bottom minus ml_margin_m,
    zmin <= Z <= zmax for both, |Z_i - Z_j| <= max_absdiff.
    """
    c = {**DEFAULTS, **kw}
    Z, H, C, D = res["radar_dbz_SW"], res["radar_hgt_SW"], res["radar_cbb_SW"], res["radar_dist"]
    names = [str(n) for n in res["radars"]]
    half = D[:, None] * 1e3 * np.tan(np.radians(c["beamwidth_deg"] / 2))          # (nr, 1, ny, nx)
    good = np.isfinite(H) & (C <= c["max_cbb"]) & (H + half < res["ml_bottom"] - c["ml_margin_m"])
    rows = []
    for i, j in itertools.combinations(range(len(names)), 2):
        band = (np.isfinite(D[i]) & np.isfinite(D[j]) & (np.abs(D[i] - D[j]) <= c["max_dd_km"])
                & (D[i] >= c["min_range_km"]) & (D[j] >= c["min_range_km"]))
        if not band.any():
            continue
        iy, ix = np.nonzero(band)
        Hi, Hj = H[i][:, iy, ix], H[j][:, iy, ix]                                    # (nsw, ncell)
        ok = good[i][:, None, iy, ix] & good[j][None, :, iy, ix]
        dh = np.where(ok, np.abs(Hi[:, None] - Hj[None, :]), np.inf)                 # (nsw_i, nsw_j, ncell)
        flat = dh.reshape(-1, dh.shape[-1]); best = flat.argmin(0)
        keep = flat[best, np.arange(flat.shape[1])] <= c["max_dh_m"]
        if not keep.any():
            continue
        ki, kj = np.unravel_index(best[keep], dh.shape[:2])
        iy, ix = iy[keep], ix[keep]
        for t in range(Z.shape[1]):
            zi, zj = Z[i, t, ki, iy, ix], Z[j, t, kj, iy, ix]
            k = (np.isfinite(zi) & np.isfinite(zj) & (zi >= c["zmin"]) & (zj >= c["zmin"])
                 & (zi <= c["zmax"]) & (zj <= c["zmax"]) & (np.abs(zi - zj) <= c["max_absdiff"]))
            if k.any():
                rows.append(pd.DataFrame(dict(radar_i=names[i], radar_j=names[j], scan=t, iy=iy[k], ix=ix[k],
                                              sweep_i=ki[k], sweep_j=kj[k], h_i=H[i, ki[k], iy[k], ix[k]],
                                              h_j=H[j, kj[k], iy[k], ix[k]], z_i=zi[k], z_j=zj[k], dz=zi[k] - zj[k],
                                              d_i=D[i][iy[k], ix[k]], d_j=D[j][iy[k], ix[k]])))
    cols = ["radar_i", "radar_j", "scan", "iy", "ix", "sweep_i", "sweep_j", "h_i", "h_j", "z_i", "z_j", "dz", "d_i", "d_j"]
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=cols)


def pair_stats(pairs: pd.DataFrame, min_n=30, sigma_clip=3.0) -> pd.DataFrame:
    """Median / mean / IQR of dZ = Z_i - Z_j for every radar pair (after a sigma clip, Hu et al. 2025)."""
    rows = []
    for (a, b), d in pairs.groupby(["radar_i", "radar_j"]):
        x = d.dz.values
        if sigma_clip and len(x) > 10:
            x = x[np.abs(x - x.mean()) <= sigma_clip * x.std()]
        q1, q3 = np.percentile(x, [25, 75]) if len(x) else (np.nan, np.nan)
        rows.append(dict(radar_i=a, radar_j=b, n=len(x), median=np.median(x) if len(x) else np.nan,
                         mean=x.mean() if len(x) else np.nan, iqr=q3 - q1, std=x.std() if len(x) else np.nan,
                         used=len(x) >= min_n))
    return pd.DataFrame(rows)


def components(stats: pd.DataFrame):
    """Groups of radars linked by usable pairs (offsets are only comparable inside a group)."""
    s = stats[stats.used]
    parent = {}
    find = lambda x: x if parent.setdefault(x, x) == x else find(parent[x])
    for _, r in s.iterrows():
        parent[find(r.radar_i)] = find(r.radar_j)
    groups = {}
    for x in list(parent):
        groups.setdefault(find(x), []).append(x)
    return [sorted(g) for g in groups.values()]


def solve_offsets(stats: pd.DataFrame, radars, stat="median") -> pd.Series:
    """Weighted least squares for dZ_ij = c_i - c_j (weight sqrt(n)) with sum(c) = 0 inside every connected
    group of radars. Radars without any usable pair get NaN (undetermined). Offsets of different groups are
    NOT comparable with each other: their common level must come from gauges (MFB) or a spaceborne radar."""
    s = stats[stats.used]
    radars = [str(r) for r in radars]
    out = pd.Series(np.nan, index=radars, dtype=float)
    for grp in components(stats):
        ss = s[s.radar_i.isin(grp)]
        idx = {r: k for k, r in enumerate(grp)}
        A = np.zeros((len(ss) + 1, len(grp))); y = np.zeros(len(ss) + 1); w = np.ones(len(ss) + 1)
        for m, (_, r) in enumerate(ss.iterrows()):
            A[m, idx[r.radar_i]] = 1; A[m, idx[r.radar_j]] = -1; y[m] = r[stat]; w[m] = np.sqrt(r.n)
        A[-1, :] = 1; w[-1] = 1e3                                # constraint: offsets of the group sum to zero
        out[grp] = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)[0]
    return out


def loop_closure(stats: pd.DataFrame, stat="median") -> pd.DataFrame:
    """Sum of dZ around every closed triangle of radars (should be ~0 dB, Gourley et al. 2003)."""
    s = stats[stats.used]
    d = {}
    for _, r in s.iterrows():
        d[(r.radar_i, r.radar_j)] = r[stat]; d[(r.radar_j, r.radar_i)] = -r[stat]
    names = sorted(set(s.radar_i) | set(s.radar_j))
    rows = []
    for a, b, c in itertools.combinations(names, 3):
        if (a, b) in d and (b, c) in d and (c, a) in d:
            rows.append(dict(loop=f"{a}→{b}→{c}→{a}", dZ_ab=d[(a, b)], dZ_bc=d[(b, c)], dZ_ca=d[(c, a)],
                             residual=d[(a, b)] + d[(b, c)] + d[(c, a)]))
    return pd.DataFrame(rows)


def bootstrap_offsets(pairs: pd.DataFrame, radars, n_boot=300, block="scan_id", seed=1, **kw) -> pd.DataFrame:
    """Bootstrap the network offsets by resampling whole scans (keeps the spatial correlation inside a scan)."""
    rng = np.random.default_rng(seed)
    ids = pairs[block].unique()
    out = []
    for _ in range(n_boot):
        pick = rng.choice(ids, len(ids))
        p = pd.concat([pairs[pairs[block] == i] for i in pick], ignore_index=True)
        out.append(solve_offsets(pair_stats(p, **kw), radars))
    return pd.DataFrame(out)


def rain_factor(offset_db, b=1.6):
    """Multiplicative effect on R(Z) of removing a calibration offset c (dB): R' = R * 10^(-c / (10 b))."""
    return 10 ** (-np.asarray(offset_db, float) / (10 * b))


def factor_to_db(F, b=1.6):
    """Express a multiplicative rain factor as an equivalent Z offset in dB (inverse of rain_factor)."""
    return -10 * b * np.log10(np.asarray(F, float))
