"""Point-wise evaluation of mosaic configurations with leave-one-gauge-out adjustment (notebook 10).

Every configuration is evaluated exactly like a field sampled with verify.sample(window=3): the mosaic is
computed only on the 3×3 window around each gauge, which makes per-gauge leave-one-out factors cheap.
"""
from __future__ import annotations

import numpy as np

from . import mosaic

GROUP = {"LS_MP": "LS", "LS_DP": "LS", "HYZ_MP": "HYZ", "PHY_DP": "PHY", "CAPPI_RT": "CAP"}


def _windows(glon, glat, lon, lat, w=3):
    lon = np.asarray(lon, float); lat = np.asarray(lat, float)
    ix = np.abs(glon[None, :] - lon[:, None]).argmin(1)
    iy = np.abs(glat[None, :] - lat[:, None]).argmin(1)
    off = np.arange(-(w // 2), w // 2 + 1)
    YY = np.clip(iy[:, None] + off, 0, len(glat) - 1)
    XX = np.clip(ix[:, None] + off, 0, len(glon) - 1)
    return YY[:, :, None], XX[:, None, :]                      # broadcast to (npts, w, w)


def _window_mean(a, w=3, min_frac=0.5):                       # a (npts, w, w)
    n = np.isfinite(a).sum((1, 2))
    with np.errstate(invalid="ignore"):
        m = np.nanmean(np.where(np.isfinite(a), a, np.nan).reshape(len(a), -1), 1)
    return np.where(n >= max(1, int(min_frac * w * w + 0.5)), m, np.nan)


def mosaic_at_points(res, product, method, lon, lat, scale=None, window=3, k=1.0):
    """Hourly mosaic depth (mm) at points (3×3 window mean).
    scale : None, per-radar factors (n_radars,) or per-point per-radar factors (n_points, n_radars)."""
    YY, XX = _windows(res["glon"], res["glat"], lon, lat, window)
    R = res[f"radar_{product}"][:, :, YY, XX].astype(float)    # (nr, nt, npts, w, w)
    Q = res[f"radar_rqi_{GROUP[product]}"][:, :, YY, XX]
    D = res["radar_dist"][:, YY, XX]
    if scale is not None:
        s = np.asarray(scale, float)
        R = R * (s[:, None, None, None, None] if s.ndim == 1 else s.T[:, None, :, None, None])
    nr, nt, npnt, w, _ = R.shape
    flat = lambda A: A.reshape(A.shape[:-3] + (npnt * w, w))
    m15 = np.array([mosaic.combine(flat(R[:, t]), flat(Q[:, t]), flat(D), method, k) * 0.25 for t in range(nt)])
    h = mosaic.hourly_from_15min(m15).reshape(npnt, w, w)
    return _window_mean(h, window)


def per_radar_at_points(res, product, lon, lat, window=3):
    """Hourly depth of every single radar at points: (n_points, n_radars)."""
    YY, XX = _windows(res["glon"], res["glat"], lon, lat, window)
    R = res[f"radar_{product}"][:, :, YY, XX].astype(float)
    out = [_window_mean(mosaic.hourly_from_15min(R[i] * 0.25), window) for i in range(R.shape[0])]
    return np.array(out).T


def per_radar_mfb(G, Rpts, loo=True, wet=None, min_n=1, fmin=None, fmax=None):
    """Per-radar mean-field bias F_i = sum(G) / sum(R_i) over the gauges radar i sees (Mahavik et al. 2025:
    computed once per radar and applied before the MAX mosaic). Defaults reproduce the paper (all pairs, no
    limits). A robust variant: wet=0.1 (pairs where G or R_i >= 0.1 mm), min_n=5, fmin=0.2, fmax=5.
    loo=True: the factor for gauge s is computed without gauge s -> (n_points, n_radars); else (n_radars,)."""
    G = np.broadcast_to(np.asarray(G, float)[:, None], Rpts.shape)
    ok = np.isfinite(G) & np.isfinite(Rpts)
    if wet is not None:
        ok &= (G >= wet) | (Rpts >= wet)
    g_ = np.where(ok, G, 0.0); r_ = np.where(ok, Rpts, 0.0); n_ = ok.astype(int)
    sg, sr, n = g_.sum(0), r_.sum(0), n_.sum(0)
    if loo:
        sg, sr, n = sg[None] - g_, sr[None] - r_, n[None] - n_
    F = np.where((sr > 0) & (n >= min_n), sg / np.where(sr > 0, sr, 1), 1.0)
    if fmin is not None or fmax is not None:
        F = np.clip(F, fmin if fmin is not None else -np.inf, fmax if fmax is not None else np.inf)
    return F
