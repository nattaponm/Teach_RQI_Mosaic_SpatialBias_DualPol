"""Gauge–radar pairing and verification metrics.

Continuous scores (R = radar, G = gauge, n pairs):
    ME    = (1/n) Σ (R − G)                         mean error (mm)
    MAE   = (1/n) Σ |R − G|
    RMSE  = sqrt((1/n) Σ (R − G)²)
    PBIAS = 100 · Σ(R − G) / ΣG                     % (negative = radar too low)
    B     = ΣR / ΣG                                 multiplicative bias (mean-field factor = 1/B)
    r     = Pearson correlation;  ρ = Spearman rank correlation
    NSE   = 1 − Σ(R − G)² / Σ(G − Ḡ)²               Nash–Sutcliffe efficiency
    KGE   = 1 − sqrt((r − 1)² + (α − 1)² + (β − 1)²), α = σR/σG, β = μR/μG  (Gupta et al. 2009)
Categorical scores for rain/no-rain at threshold thr (contingency table hits a, false alarms b, misses c, correct
negatives d):  POD = a/(a+c), FAR = b/(a+b), CSI = a/(a+b+c),
               HSS = 2(ad − bc) / [(a+c)(c+d) + (a+b)(b+d)]
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def sample(arr, glon, glat, lon, lat, window=3, min_frac=0.5):
    """Radar value at points: nearest cell (window=1) or mean of a window×window block."""
    lon = np.asarray(lon, float); lat = np.asarray(lat, float)
    ix = np.abs(glon[None, :] - lon[:, None]).argmin(1)
    iy = np.abs(glat[None, :] - lat[:, None]).argmin(1)
    if window == 1:
        return arr[iy, ix]
    h = window // 2
    out = np.full(len(lon), np.nan)
    for k, (i, j) in enumerate(zip(iy, ix)):
        w = arr[max(i - h, 0):i + h + 1, max(j - h, 0):j + h + 1]
        if np.isfinite(w).sum() >= max(1, int(min_frac * window * window + 0.5)):
            out[k] = np.nanmean(w)
    return out


def metrics(g, r, thr=0.5):
    g = np.asarray(g, float); r = np.asarray(r, float)
    k = np.isfinite(g) & np.isfinite(r)
    g, r = g[k], r[k]
    n = len(g)
    out = {"n": n}
    if n < 3:
        return out
    d = r - g
    sg, sr = g.sum(), r.sum()
    out.update(G_mean=g.mean(), R_mean=r.mean(), ME=d.mean(), MAE=np.abs(d).mean(), RMSE=np.sqrt((d ** 2).mean()),
               PBIAS=100 * d.sum() / sg if sg > 0 else np.nan, B=sr / sg if sg > 0 else np.nan)
    if g.std() > 0 and r.std() > 0:
        rr = np.corrcoef(g, r)[0, 1]
        out.update(r=rr, rho=pd.Series(g).corr(pd.Series(r), method="spearman"),
                   NSE=1 - (d ** 2).sum() / ((g - g.mean()) ** 2).sum())
        beta = r.mean() / g.mean() if g.mean() > 0 else np.nan
        out["KGE"] = 1 - np.sqrt((rr - 1) ** 2 + (r.std() / g.std() - 1) ** 2 + (beta - 1) ** 2)
    a = int(((g >= thr) & (r >= thr)).sum()); c = int(((g >= thr) & (r < thr)).sum())
    b = int(((g < thr) & (r >= thr)).sum()); dd = int(((g < thr) & (r < thr)).sum())
    out.update(hits=a, misses=c, false_alarms=b, correct_negatives=dd,
               POD=a / (a + c) if a + c else np.nan, FAR=b / (a + b) if a + b else np.nan,
               CSI=a / (a + b + c) if a + b + c else np.nan)
    den = (a + c) * (c + dd) + (a + b) * (b + dd)
    out["HSS"] = 2 * (a * dd - b * c) / den if den else np.nan
    wet = (g >= thr) | (r >= thr)
    if wet.sum() >= 3:
        gw, rw = g[wet], r[wet]
        out.update(n_wet=int(wet.sum()), RMSE_wet=np.sqrt(((rw - gw) ** 2).mean()),
                   PBIAS_wet=100 * (rw - gw).sum() / gw.sum() if gw.sum() > 0 else np.nan,
                   r_wet=np.corrcoef(gw, rw)[0, 1] if gw.std() > 0 and rw.std() > 0 else np.nan)
    return out


def metrics_table(pairs, products, by=None, thr=0.5):
    """Metrics for every product (and every group in column `by`)."""
    rows = []
    groups = [("all", pairs)] if by is None else list(pairs.groupby(by, observed=True))
    for name, g in groups:
        for p in products:
            rows.append({"group": name, "product": p, **metrics(g.G.values, g[p].values, thr)})
    return pd.DataFrame(rows)


def pair(gauge_hour: pd.DataFrame, fields: dict, glon, glat, window=3, rqi=None):
    """gauge_hour: rows with station_code, lon, lat, G for ONE radar hour. fields: {product: 2-D array}."""
    d = gauge_hour.copy()
    for p, a in fields.items():
        d[p] = sample(a, glon, glat, d.lon.values, d.lat.values, window)
    if rqi is not None:
        for k, a in rqi.items():
            d[k] = sample(a, glon, glat, d.lon.values, d.lat.values, 1)
    return d


def lag_correlation(gauge, meta, field, glon, glat, radar_hour, lags=(0, 1, 2), window=3):
    """Correlation between the radar hour [t, t+1h) UTC and gauge values whose END stamp is t+lag (UTC)."""
    m = meta.set_index("station_code")
    rows = []
    for lag in lags:
        g = gauge[(gauge.end_utc == radar_hour + pd.Timedelta(hours=lag)) & (gauge.qc == "OK")]
        g = g[g.station_code.isin(m.index)]
        r = sample(field, glon, glat, m.loc[g.station_code, "longitude"], m.loc[g.station_code, "latitude"], window)
        k = np.isfinite(r) & np.isfinite(g.rain_qc.values)
        rows.append(dict(gauge_end_utc=radar_hour + pd.Timedelta(hours=lag),
                         gauge_end_local=radar_hour + pd.Timedelta(hours=lag + 7), n=int(k.sum()),
                         r=np.corrcoef(g.rain_qc.values[k], r[k])[0, 1] if k.sum() > 3 else np.nan))
    return pd.DataFrame(rows)
