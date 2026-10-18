"""Reading radar files and building a scan inventory."""
from __future__ import annotations

import bz2
import gzip
import os
import re
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=DeprecationWarning)
NAME_RE = re.compile(r"([A-Z]{3})_?(?:240)?(\d{12})")   # LMP_202410030000.uf or LMP240202410030000.uf.bz2


def _pyart():
    import io as _io, contextlib
    with contextlib.redirect_stdout(_io.StringIO()):   # silence the Py-ART banner
        import pyart
    return pyart


def read_radar(path):
    """Read a UF/CfRadial file with Py-ART. Handles real .bz2/.gz and files that are named .bz2 but are plain UF."""
    pyart = _pyart()
    path = str(path)
    with open(path, "rb") as fh:
        magic = fh.read(3)
    if magic == b"BZh" or magic[:2] == b"\x1f\x8b":
        opener = bz2.open if magic == b"BZh" else gzip.open
        with opener(path, "rb") as src, tempfile.NamedTemporaryFile(suffix=".uf", delete=False) as dst:
            shutil.copyfileobj(src, dst)
            tmp = dst.name
        try:
            return pyart.io.read(tmp)
        finally:
            os.remove(tmp)
    return pyart.io.read(path)


def parse_name(path):
    m = NAME_RE.search(os.path.basename(str(path)))
    if not m:
        return None, None
    return m.group(1), pd.Timestamp(m.group(2)[:8] + "T" + m.group(2)[8:12])


def list_scans(folder) -> pd.DataFrame:
    """Table of radar files in a folder: radar, time (UTC), path. Sorted by time then radar."""
    rows = []
    for p in sorted(Path(folder).rglob("*")):
        if p.is_file() and re.search(r"\.uf(\.bz2|\.gz)?$", p.name, re.I):
            rad, t = parse_name(p)
            if rad:
                rows.append(dict(radar=rad, time=t, path=str(p), size_mb=p.stat().st_size / 1e6))
    return pd.DataFrame(rows).sort_values(["time", "radar"]).reset_index(drop=True)


def describe(radar) -> dict:
    """Short summary of a Py-ART radar object (used in notebook 00/01)."""
    ip = radar.instrument_parameters or {}
    return dict(
        latitude=float(radar.latitude["data"][0]), longitude=float(radar.longitude["data"][0]),
        altitude_m=float(radar.altitude["data"][0]), n_sweeps=int(radar.nsweeps),
        elevations_deg=[round(float(e), 2) for e in radar.fixed_angle["data"]],
        rays_per_sweep=[int(e - s + 1) for s, e in zip(radar.sweep_start_ray_index["data"], radar.sweep_end_ray_index["data"])],
        gate_m=float(np.median(np.diff(radar.range["data"]))), max_range_km=float(radar.range["data"][-1] / 1e3),
        frequency_ghz=float(ip["frequency"]["data"][0] / 1e9) if "frequency" in ip else None,
        beamwidth_deg=float(ip["radar_beam_width_h"]["data"][0]) if "radar_beam_width_h" in ip else None,
        fields=sorted(radar.fields),
        time_start=str(radar.time["units"]).replace("seconds since ", ""),
    )


def inventory(folder) -> pd.DataFrame:
    """Open every file in a folder and summarise it (slow for many files: ~0.5 s per file)."""
    rows = []
    for _, r in list_scans(folder).iterrows():
        try:
            d = describe(read_radar(r.path))
            d.update(radar=r.radar, time=r.time, file=os.path.basename(r.path), size_mb=round(r.size_mb, 1))
            d["fields"] = ";".join(d["fields"]); d["elevations_deg"] = ";".join(map(str, d["elevations_deg"]))
            d["rays_per_sweep"] = ";".join(map(str, d["rays_per_sweep"]))
        except Exception as e:
            d = dict(radar=r.radar, time=r.time, file=os.path.basename(r.path), error=str(e)[:120])
        rows.append(d)
    cols = ["radar", "time", "file", "size_mb", "n_sweeps", "elevations_deg", "rays_per_sweep", "gate_m",
            "max_range_km", "frequency_ghz", "beamwidth_deg", "fields"]
    df = pd.DataFrame(rows)
    return df[[c for c in cols if c in df] + [c for c in df if c not in cols]]


def field(radar, name, sweep=None):
    """Field as a float array with NaN for masked gates (optionally one sweep)."""
    a = radar.fields[name]["data"]
    if sweep is not None:
        s, e = radar.sweep_start_ray_index["data"][sweep], radar.sweep_end_ray_index["data"][sweep]
        a = a[s:e + 1]
    return np.ma.filled(a.astype(float), np.nan)
