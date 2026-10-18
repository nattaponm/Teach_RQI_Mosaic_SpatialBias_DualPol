"""Bias adjustment of radar rainfall with gauges.

1) Mean-field bias (MFB):   F = ΣG / ΣR  over gauges where both are wet;   R' = F · R
2) Local (spatial) bias:    at each gauge  b_i = log((G_i + ε)/(R_i + ε))   (log-ratio, ε = 0.1 mm)
                            interpolate b(x) with inverse-distance weighting (IDW) or ordinary kriging,
                            R'(x) = R(x) · exp(b(x))
   IDW:  b(x) = Σ w_i b_i / Σ w_i,   w_i = 1 / d_i^p   (p = 2), only gauges within `radius`
   Far from gauges the weights fade toward the mean-field log-bias (blend with distance).
3) Validation: leave-one-out cross-validation (LOOCV) — each gauge is removed, the adjustment is rebuilt from
   the others and evaluated at the removed gauge (honest; the gauge is not used to correct itself).
4) Spatial autocorrelation of the log-bias: Moran's I with inverse-distance weights.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

EPS = 0.1


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def mean_field_factor(g, r, wet=0.1, fmin=0.1, fmax=20.0):
    g = np.asarray(g, float); r = np.asarray(r, float)
    k = np.isfinite(g) & np.isfinite(r) & (g >= wet) & (r >= wet)
    if k.sum() < 2 or r[k].sum() <= 0:
        return 1.0, int(k.sum())
    return float(np.clip(g[k].sum() / r[k].sum(), fmin, fmax)), int(k.sum())


def log_ratio(g, r, eps=EPS):
    return np.log((np.asarray(g, float) + eps) / (np.asarray(r, float) + eps))


def idw(x_lon, x_lat, lon, lat, val, power=2.0, radius_km=50.0, background=0.0, range_km=None):
    """IDW of point values onto target points; beyond `radius_km` the value relaxes to `background`.
    range_km: if given, the IDW value is blended with the background with weight exp(−d_min/range_km)."""
    x_lon = np.atleast_1d(x_lon).ravel(); x_lat = np.atleast_1d(x_lat).ravel()
    out = np.full(x_lon.shape, background, float)
    if len(val) == 0:
        return out
    for i0 in range(0, len(x_lon), 20000):
        sl = slice(i0, i0 + 20000)
        d = haversine_km(x_lat[sl, None], x_lon[sl, None], np.asarray(lat)[None, :], np.asarray(lon)[None, :])
        d = np.maximum(d, 0.5)
        w = np.where(d <= radius_km, 1.0 / d ** power, 0.0)
        sw = w.sum(1)
        v = np.where(sw > 0, (w * np.asarray(val)[None, :]).sum(1) / np.maximum(sw, 1e-12), background)
        if range_km:
            a = np.exp(-d.min(1) / range_km)
            v = a * v + (1 - a) * background
        out[sl] = v
    return out


def empirical_variogram(lon, lat, val, bins_km=np.arange(0, 160, 10)):
    d = haversine_km(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
    g = 0.5 * (val[:, None] - val[None, :]) ** 2
    iu = np.triu_indices(len(val), 1)
    d, g = d[iu], g[iu]
    rows = []
    for a, b in zip(bins_km[:-1], bins_km[1:]):
        k = (d >= a) & (d < b)
        if k.sum() >= 5:
            rows.append(dict(h_km=(a + b) / 2, gamma=g[k].mean(), n_pairs=int(k.sum())))
    return pd.DataFrame(rows)


def exp_variogram(h, nugget, sill, range_km):
    return nugget + (sill - nugget) * (1 - np.exp(-3 * np.asarray(h) / range_km))


def fit_variogram(ev):
    from scipy.optimize import curve_fit
    p0 = [0.1 * ev.gamma.max(), ev.gamma.max(), ev.h_km.max() / 2]
    try:
        p, _ = curve_fit(exp_variogram, ev.h_km, ev.gamma, p0=p0, sigma=1 / np.sqrt(ev.n_pairs),
                         bounds=([0, 0, 5], [np.inf, np.inf, 500]))
    except Exception:
        p = p0
    return dict(nugget=float(p[0]), sill=float(p[1]), range_km=float(p[2]))


def ordinary_kriging(x_lon, x_lat, lon, lat, val, vg, max_n=16):
    """Ordinary kriging with an exponential variogram, using the max_n nearest gauges."""
    x_lon = np.atleast_1d(x_lon).ravel(); x_lat = np.atleast_1d(x_lat).ravel()
    lon = np.asarray(lon); lat = np.asarray(lat); val = np.asarray(val)
    out = np.full(x_lon.shape, np.nan); var = np.full(x_lon.shape, np.nan)
    gam = lambda h: exp_variogram(h, vg["nugget"], vg["sill"], vg["range_km"])
    for i in range(len(x_lon)):
        d = haversine_km(x_lat[i], x_lon[i], lat, lon)
        nn = np.argsort(d)[:max_n]
        dd = haversine_km(lat[nn, None], lon[nn, None], lat[None, nn], lon[None, nn])
        n = len(nn)
        A = np.ones((n + 1, n + 1)); A[:n, :n] = gam(dd); np.fill_diagonal(A[:n, :n], 0); A[n, n] = 0
        b = np.ones(n + 1); b[:n] = gam(d[nn])
        try:
            w = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            continue
        out[i] = w[:n] @ val[nn]
        var[i] = w @ b
    return out, var


def adjust_field(R, glon, glat, lon, lat, g_at, r_at, method="idw", **kw):
    """Return the adjusted field and the log-bias surface.
    method: 'mfb' (one factor) or 'idw' (local log-ratio, relaxed to the MFB far from gauges)."""
    F, _ = mean_field_factor(g_at, r_at)
    if method == "mfb":
        return R * F, np.full(R.shape, np.log(F))
    k = np.isfinite(g_at) & np.isfinite(r_at) & ((np.asarray(g_at) >= 0.1) | (np.asarray(r_at) >= 0.1))
    b = log_ratio(np.asarray(g_at)[k], np.asarray(r_at)[k])
    LON, LAT = np.meshgrid(glon, glat)
    surf = idw(LON, LAT, np.asarray(lon)[k], np.asarray(lat)[k], b, background=np.log(F), **kw).reshape(R.shape)
    return R * np.exp(surf), surf


def loocv(points: pd.DataFrame, product: str, methods=("raw", "mfb", "idw"), mfb_background=True, **kw) -> pd.DataFrame:
    """points: one row per gauge (and hour) with lon, lat, G and the radar value `product` at the gauge.
    Returns the cross-validated adjusted value for every method (column names = methods)."""
    out = points[["station_code", "lon", "lat", "G", product]].copy()
    hours = points["start_utc"] if "start_utc" in points else pd.Series(0, index=points.index)
    for m in methods:
        out[m] = np.nan
    for h in hours.unique():
        P = points[hours == h]
        for i in P.index:
            rest = P.drop(index=i)
            r_i = P.at[i, product]
            if not np.isfinite(r_i):
                continue
            F, _ = mean_field_factor(rest.G, rest[product])
            if "raw" in methods:
                out.at[i, "raw"] = r_i
            if "mfb" in methods:
                out.at[i, "mfb"] = F * r_i
            if "idw" in methods:
                k = np.isfinite(rest.G) & np.isfinite(rest[product]) & ((rest.G >= 0.1) | (rest[product] >= 0.1))
                b = idw(P.at[i, "lon"], P.at[i, "lat"], rest.lon[k], rest.lat[k],
                        log_ratio(rest.G[k], rest[product][k]), background=np.log(F) if mfb_background else 0.0, **kw)[0]
                out.at[i, "idw"] = (r_i + EPS) * np.exp(b) - EPS if r_i > 0 else 0.0
    for m in methods:
        out[m] = out[m].clip(lower=0)
    return out


def morans_i(lon, lat, val, power=1.0, max_km=100.0):
    """Global Moran's I with inverse-distance weights within max_km, and a permutation p-value."""
    lon = np.asarray(lon); lat = np.asarray(lat); x = np.asarray(val, float)
    d = haversine_km(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
    W = np.where((d > 0) & (d <= max_km), 1 / np.maximum(d, 0.5) ** power, 0.0)
    z = x - x.mean()

    def I(zz):
        return len(zz) / W.sum() * (zz @ W @ zz) / (zz @ zz)
    obs = I(z)
    rng = np.random.default_rng(0)
    sims = np.array([I(rng.permutation(z)) for _ in range(999)])
    p = (np.sum(sims >= obs) + 1) / (len(sims) + 1)
    return dict(I=float(obs), expected=-1 / (len(x) - 1), p_value=float(p))
