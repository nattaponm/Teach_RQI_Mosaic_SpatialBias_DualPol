"""Melting layer and apparent vertical profile of reflectivity (VPR).

Melting layer (ML): the copolar correlation ρhv drops (0.90–0.96) where snow melts into rain, because the
particles are a mix of shapes and water coatings. We search 3–6.5 km for the lowest ρhv in the median profile.

Apparent VPR: median corrected Z in 200-m height bins over all radars/tilts (10–120 km, CBB < 0.3).
Correction for a bin at height h:   ΔZ(h) = Z_ref − Z_VPR(h),  Z_ref = median Z in 0.5–2 km,
applied only above the ML bottom and capped at +12 dB (Germann & Joss 2002; Kirstetter et al. 2010 idea).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def profile_from_samples(h_m, dbz, rhohv=None, bin_m=200.0, hmax=12000.0, min_n=200):
    """Median Z (and ρhv) per height bin. h_m, dbz, rhohv: 1-D arrays of samples."""
    edges = np.arange(0, hmax + bin_m, bin_m)
    df = pd.DataFrame({"b": np.digitize(h_m, edges), "z": dbz})
    if rhohv is not None:
        df["rh"] = rhohv
    g = df.groupby("b")
    prof = pd.DataFrame({"h": edges[g.size().index - 1] + bin_m / 2, "n": g.size().values,
                         "z_med": g.z.median().values,
                         "z_p25": g.z.quantile(0.25).values, "z_p75": g.z.quantile(0.75).values})
    if rhohv is not None:
        prof["rh_med"] = g.rh.median().values
    prof = prof[(prof.n >= min_n) & (prof.h <= hmax)].set_index("h")
    return prof


def melting_layer(prof, search=(3000, 6500), dip=0.975, min_bottom=3500.0, fallback_fzl=4800.0):
    """(bottom, top, source) of the melting layer from the ρhv dip in the median profile."""
    if "rh_med" in prof:
        win = prof.loc[(prof.index >= search[0]) & (prof.index <= search[1]), "rh_med"].dropna()
        if len(win) and win.min() < 0.97:
            hs = list(win.index); i0 = int(np.argmin(win.values)); lo = hi = i0
            while lo > 0 and win.values[lo - 1] < dip:
                lo -= 1
            while hi < len(hs) - 1 and win.values[hi + 1] < dip:
                hi += 1
            return (max(float(hs[lo] - 100), min_bottom), float(hs[hi] + 100),
                    f"ρhv dip (min {win.min():.3f} at {hs[i0]:.0f} m)")
    return fallback_fzl - 500, fallback_fzl + 200, "fallback freezing level"


def correction_table(prof, ml_bottom, ref_layer=(500, 2000), max_db=12.0, top_slope_db_km=-6.0, hmax=12000):
    """ΔZ(h) on a 100-m grid. Returns (h grid, dZ, Z_ref)."""
    ref = prof.loc[(prof.index >= ref_layer[0]) & (prof.index <= ref_layer[1]), "z_med"].median()
    smooth = prof.z_med.rolling(3, center=True, min_periods=1).median()
    hg = np.arange(0, hmax + 1, 100.0)
    zp = np.interp(hg, smooth.index.values, smooth.values)
    top = smooth.index.max()
    zp[hg > top] = smooth.values[-1] + top_slope_db_km * (hg[hg > top] - top) / 1000
    dz = np.clip(ref - zp, -6, max_db)
    dz[hg < ml_bottom] = 0.0
    return hg, dz, float(ref)


def apply(dbz, h_m, hg, dz):
    return np.asarray(dbz) + np.interp(h_m, hg, dz)
