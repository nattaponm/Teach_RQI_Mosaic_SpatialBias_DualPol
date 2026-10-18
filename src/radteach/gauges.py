"""
gauges.py — read, quality-control and time-align HII (สสน.) hourly rain-gauge data.

Data source (open data catalogue of the Hydro-Informatics Institute):
    https://tiservice.hii.or.th/opendata/data_catalog/hourly_rain/{YYYY}/{YYYYMM}/0station_metadata.csv
    https://tiservice.hii.or.th/opendata/data_catalog/hourly_rain/{YYYY}/{YYYYMM}/{STATION}.csv

Station file format (verified on 202410/BNHI.csv, MEWG.csv, ...):
    station_code,measure_datetime,rainfall_1h,quality_flag
    BNHI,2024-10-01 00:00:00,0,N
    CHM001,2024-10-03 07:00:00,null,null

Conventions used here (checked against the radar, see gauge_radar_eval.py --lagcheck):
  * measure_datetime is Thai local time (UTC+7).
  * rainfall_1h is the accumulation over the hour ENDING at measure_datetime.
    e.g. 08:00 LT = 07:00-08:00 LT = 00:00-01:00 UTC.
    Correlation with the radar hour 00-01 UTC on 3 Oct 2024 peaks at the 08:00 stamp
    (r = 0.96, n = 50) and the 10-11 UTC hour peaks at 18:00 (r = 0.83).

Quality control (each row gets one qc code; only qc == "OK" is used):
  NODATA    value is null / empty / non-numeric
  FLAG      quality_flag not in the accepted set (default {"N"}; the flag column itself
            is null when the value is missing)
  RANGE     value < 0 or > max_rate (default 150 mm/h)
  STUCK     >= stuck_frac of the hours in a centred window are non-zero but tiny and
            almost constant (typical of a leaking/vibrating tipping bucket, e.g. 0.2 or
            0.4 mm every hour for half a day). Detected per station.
  ISOLATED  (reported, not removed by default) gauge >= iso_mm while all >= 3 neighbours
            within iso_km report 0 in the same hour — a gauge-only buddy check, so the
            radar is never used to judge the gauges.

Runs on Python 3.8 (env "mossf") — only numpy, pandas, requests are needed.
"""
from __future__ import annotations

import io
import os
import time

import numpy as np
import pandas as pd

BASE_URL = "https://tiservice.hii.or.th/opendata/data_catalog/hourly_rain/"
LOCAL_TZ_OFFSET_H = 7

BASIN_EN = {
    "ปิง": "Ping", "วัง": "Wang", "ยม": "Yom", "น่าน": "Nan", "สาละวิน": "Salween",
    "โขงเหนือ": "Mekong North (Kok-Ing)", "ป่าสัก": "Pa Sak", "ชี": "Chi",
    "โขงตะวันออกเฉียงเหนือ": "Mekong NE", "แม่กลอง": "Mae Klong", "เจ้าพระยา": "Chao Phraya",
    "มูล": "Mun",
}


# --------------------------------------------------------------------------- download
def month_url(year: int, month: int, name: str) -> str:
    return f"{BASE_URL}{year:04d}/{year:04d}{month:02d}/{name}"


def _get(url, session, retries=3, timeout=60):
    for k in range(retries):
        try:
            r = session.get(url, timeout=timeout)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            r.encoding = "utf-8"
            return r.text
        except Exception as e:  # network hiccup → back off and retry
            if k == retries - 1:
                print(f"  ! failed {url}: {e}")
                return None
            time.sleep(2 * (k + 1))


def fetch_month(year, month, codes=None, bbox=None, outdir="hii_raw", session=None, pause=0.2):
    """Download metadata + station CSVs of one month to outdir/YYYYMM/ (skips files already there).

    codes : list of station codes, or None for all stations (optionally inside bbox)
    bbox  : (lon_min, lon_max, lat_min, lat_max)
    """
    import requests  # imported here so the rest of the module works without it
    session = session or requests.Session()
    d = os.path.join(outdir, f"{year:04d}{month:02d}")
    os.makedirs(d, exist_ok=True)
    meta_path = os.path.join(d, "0station_metadata.csv")
    if not os.path.exists(meta_path):
        txt = _get(month_url(year, month, "0station_metadata.csv"), session)
        if txt is None:
            raise RuntimeError(f"no metadata for {year}{month:02d}")
        open(meta_path, "w", encoding="utf-8").write(txt)
    meta = read_metadata(meta_path)
    sel = meta
    if bbox is not None:
        x0, x1, y0, y1 = bbox
        sel = sel[sel.longitude.between(x0, x1) & sel.latitude.between(y0, y1)]
    if codes is not None:
        sel = sel[sel.station_code.isin(codes)]
    n_new = 0
    for c in sel.station_code:
        p = os.path.join(d, f"{c}.csv")
        if os.path.exists(p) and os.path.getsize(p) > 0:
            continue
        txt = _get(month_url(year, month, f"{c}.csv"), session)
        if txt is None:
            continue
        open(p, "w", encoding="utf-8").write(txt)
        n_new += 1
        time.sleep(pause)  # be polite to the server
    print(f"{year}{month:02d}: {len(sel)} stations selected, {n_new} downloaded -> {d}")
    return d


# --------------------------------------------------------------------------- read
def read_metadata(path):
    m = pd.read_csv(path, encoding="utf-8")
    m.columns = [c.strip() for c in m.columns]
    m["station_code"] = m["station_code"].astype(str).str.strip()
    for c in ("latitude", "longitude"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
    m["basin_en"] = m["basin"].map(BASIN_EN).fillna(m["basin"])
    return m


def read_station_csv(path_or_text):
    """Return station_code, time_lt (end of hour, local), rain, quality_flag."""
    src = io.StringIO(path_or_text) if "\n" in str(path_or_text) else path_or_text
    df = pd.read_csv(src, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    df["time_lt"] = pd.to_datetime(df["measure_datetime"], errors="coerce")
    df["rain"] = pd.to_numeric(df["rainfall_1h"].replace({"null": "", "NULL": "", "-": ""}),
                               errors="coerce")
    df["quality_flag"] = df["quality_flag"].replace({"": "null"}).str.strip()
    return df[["station_code", "time_lt", "rain", "quality_flag"]]


def read_month_dir(d):
    """Read every <STATION>.csv in a month folder into one long table."""
    out = []
    for f in sorted(os.listdir(d)):
        if f.endswith(".csv") and not f.startswith("0station"):
            try:
                out.append(read_station_csv(os.path.join(d, f)))
            except Exception as e:
                print(f"  ! cannot read {f}: {e}")
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# --------------------------------------------------------------------------- time
def add_utc(df):
    """time_lt marks the END of the accumulation hour (local).  Adds:
       end_utc   = time_lt - 7 h
       start_utc = end_utc - 1 h   (the radar hour this value belongs to)"""
    df = df.copy()
    df["end_utc"] = df["time_lt"] - pd.Timedelta(hours=LOCAL_TZ_OFFSET_H)
    df["start_utc"] = df["end_utc"] - pd.Timedelta(hours=1)
    return df


# --------------------------------------------------------------------------- QC
def _stuck_mask(s, window=12, frac=0.9, vmax=0.6, std_max=0.15):
    """s: hourly series of one station (sorted).  True where the station looks stuck."""
    x = s.astype(float)
    nz = (x > 0) & (x <= vmax)
    roll_nz = nz.rolling(window, center=True, min_periods=window // 2).mean()
    roll_sd = x.rolling(window, center=True, min_periods=window // 2).std()
    roll_mx = x.rolling(window, center=True, min_periods=window // 2).max()
    core = ((roll_nz >= frac) & (roll_sd <= std_max) & (roll_mx <= vmax)).fillna(False)
    # dilate: every hour inside a window that triggered the test is suspect
    return core.astype(float).rolling(window, center=True, min_periods=1).max().astype(bool)


def quality_control(df, meta=None, accept_flags=("N",), max_rate=150.0,
                    stuck=True, iso_mm=10.0, iso_km=20.0, drop_isolated=False):
    """Add column 'qc' (OK / NODATA / FLAG / RANGE / STUCK / ISOLATED) and 'rain_qc'."""
    df = df.sort_values(["station_code", "time_lt"]).copy()
    qc = np.full(len(df), "OK", dtype=object)
    qc[df["rain"].isna().values] = "NODATA"
    bad_flag = ~df["quality_flag"].isin(accept_flags).values & (qc == "OK")
    qc[bad_flag] = "FLAG"
    rng = ((df["rain"] < 0) | (df["rain"] > max_rate)).values & (qc == "OK")
    qc[rng] = "RANGE"
    if stuck:
        st = df.groupby("station_code", group_keys=False)["rain"].apply(_stuck_mask)
        st = st.reindex(df.index).fillna(False).values & (qc == "OK")
        qc[st] = "STUCK"
    df["qc"] = qc

    if meta is not None and iso_mm is not None:
        df = _buddy_check(df, meta, iso_mm, iso_km)
        if drop_isolated:
            pass  # ISOLATED already written into qc
        else:  # keep the value, only report it
            df.loc[df.qc == "ISOLATED", "qc_note"] = "isolated (kept)"
            df.loc[df.qc == "ISOLATED", "qc"] = "OK"
    df["rain_qc"] = df["rain"].where(df["qc"] == "OK")
    return df


def _buddy_check(df, meta, iso_mm, iso_km):
    m = meta.set_index("station_code")[["latitude", "longitude"]]
    codes = [c for c in df.station_code.unique() if c in m.index]
    lat = np.radians(m.loc[codes, "latitude"].values)
    lon = np.radians(m.loc[codes, "longitude"].values)
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2) ** 2
    dist = 6371.0 * 2 * np.arcsin(np.sqrt(a))
    nb = {c: [codes[j] for j in np.where((dist[i] < iso_km) & (dist[i] > 0))[0]]
          for i, c in enumerate(codes)}
    ok = df[df.qc == "OK"]
    wide = ok.pivot_table(index="time_lt", columns="station_code", values="rain", aggfunc="first")
    df["qc_note"] = ""
    hits = []
    for c in wide.columns:
        if len(nb.get(c, [])) < 3:
            continue
        cols = [n for n in nb[c] if n in wide.columns]
        if len(cols) < 3:
            continue
        big = wide[c] >= iso_mm
        allzero = (wide[cols].fillna(-1) == 0).all(axis=1) & (wide[cols].notna().sum(axis=1) >= 3)
        for t in wide.index[big & allzero]:
            hits.append((c, t))
    for c, t in hits:
        df.loc[(df.station_code == c) & (df.time_lt == t), "qc"] = "ISOLATED"
    return df


def qc_report(df):
    """Counts of qc codes per station (wide) — the 'usable or not' table."""
    rep = df.pivot_table(index="station_code", columns="qc", values="time_lt",
                         aggfunc="count", fill_value=0)
    rep["n_rows"] = rep.sum(axis=1)
    rep["pct_usable"] = 100 * rep.get("OK", 0) / rep["n_rows"]
    return rep.sort_values("pct_usable")




def load_snapshot(path):
    """Read a station CSV (one or many stations) in the HII format."""
    return read_station_csv(path)


def fetch_period(start, end, codes=None, bbox=None, outdir="data/raw/gauges"):
    """Download all months between start and end (local dates) and return one long table."""
    frames, metas = [], []
    for p in pd.period_range(start, end, freq="M"):
        d = fetch_month(p.year, p.month, codes=codes, bbox=bbox, outdir=outdir)
        metas.append(read_metadata(os.path.join(d, "0station_metadata.csv")))
        frames.append(read_month_dir(d))
    meta = pd.concat(metas).drop_duplicates("station_code", keep="last")
    df = pd.concat(frames, ignore_index=True)
    t0, t1 = pd.Timestamp(start), pd.Timestamp(end) + pd.Timedelta(days=1)
    return df[(df.time_lt > t0) & (df.time_lt <= t1)].reset_index(drop=True), meta


def availability(df, freq="h"):
    """Station × time table of usable (1) / missing or rejected (0) values — for a heat map."""
    d = df.assign(ok=(df.qc == "OK").astype(int))
    return d.pivot_table(index="station_code", columns="time_lt", values="ok", aggfunc="max", fill_value=0)
