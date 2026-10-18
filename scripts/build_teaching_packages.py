"""Instructor script: build the two release assets from the original data.

    teaching_data_v1.zip      radar UF files (plain, renamed RADAR_YYYYMMDDHHMM.uf) + SRTM tiles at 250 m
    teaching_products_v1.zip  precomputed products so that every notebook can start on its own

Usage (paths are examples):
    python scripts/build_teaching_packages.py --radar-dirs D:/uf/20241003_00 D:/uf/20241003_10 \
        --srtm D:/srtm/srtm_56_08.tif D:/srtm/srtm_57_08.tif D:/srtm/srtm_56_09.tif D:/srtm/srtm_57_09.tif \
        --out build/
Then upload both zip files to a GitHub Release tagged data-v1.0 and paste the printed sha256 values into
config/data_sources.yaml.
"""
import argparse
import bz2
import os
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from radteach import accum, colab, config, dem, io, mosaic, qc  # noqa: E402


def build_data(radar_dirs, srtm, out):
    import rasterio
    from rasterio.enums import Resampling
    base = out / "teaching_data_v1"
    for d in radar_dirs:
        sc = io.list_scans(d)
        tag = f"{sc.time.min():%Y%m%d_%H}UTC"
        dd = base / "radar" / tag
        dd.mkdir(parents=True, exist_ok=True)
        for _, r in sc.iterrows():
            raw = open(r.path, "rb").read()
            if raw[:3] == b"BZh":
                raw = bz2.decompress(raw)
            (dd / f"{r.radar}_{r.time:%Y%m%d%H%M}.uf").write_bytes(raw)
    (base / "dem").mkdir(parents=True, exist_ok=True)
    for t in srtm:
        with rasterio.open(t) as s:
            a = s.read(1, out_shape=(s.height // 3, s.width // 3), resampling=Resampling.average, masked=True)
            tr = s.transform * s.transform.scale(s.width / a.shape[1], s.height / a.shape[0])
            prof = s.profile
            prof.update(height=a.shape[0], width=a.shape[1], transform=tr, compress="deflate", predictor=2,
                        tiled=True, blockxsize=512, blockysize=512)
            with rasterio.open(base / "dem" / (Path(t).stem + "_250m.tif"), "w", **prof) as dst:
                dst.write(a.filled(-32768).astype("int16"), 1)
    return base


def build_products(data_dir, out):
    P = config.paths()
    sa, rc = config.load("study_area"), config.load("radars")
    g = sa["grid"]
    glon, glat = mosaic.make_grid(g["lon_min"], g["lon_max"], g["lat_min"], g["lat_max"], g["res_deg"])
    z, tr, _ = dem.mosaic_tiles(sorted((data_dir / "dem").glob("*.tif")),
                                bounds=(g["lon_min"], g["lat_min"], g["lon_max"], g["lat_max"]))
    prod = out / "teaching_products_v1"
    prod.mkdir(parents=True, exist_ok=True)
    dem.mosaic_tiles(sorted((data_dir / "dem").glob("*.tif")),
                     bounds=(g["lon_min"], g["lat_min"], g["lon_max"], g["lat_max"]), out_path=prod / "dem_study_area.tif")
    terr = dem.interpolator(z, tr)
    for hdir in sorted((data_dir / "radar").iterdir()):
        sc = io.list_scans(hdir)
        tag = hdir.name
        print("==", tag)
        radars = {r: [io.read_radar(p) for p in gg.path] for r, gg in sc.groupby("radar")}
        cal = qc.calibrate(radars)
        del radars
        for r in cal.index:
            cal.loc[r, "phidp_sign"] = rc["radars"][r]["phidp_sign"]
        cal.to_csv(prod / f"calibration_{tag}.csv")
        res = mosaic.process_hour(sc, terr, cal, glon, glat)
        mosaic.save_npz(prod / f"mosaic_{tag}.npz", res)
        rates = np.nan_to_num(res["m15_PHY_DP"] * 4.0)
        valid = np.all(np.isfinite(res["m15_PHY_DP"]), axis=0)
        adv = {}
        for a_ in [True, False]:
            for dt in [1, 5, 10, 15]:
                d, Vs, how = accum.accumulate(rates, 15, dt, a_)
                adv[f"{'ADV' if a_ else 'EUL'}_{dt:02d}"] = np.where(valid, d, np.nan).astype("float32")
                if a_ and dt == 5:
                    adv["V"] = Vs[0].astype("float32"); adv["motion_method"] = how
        adv["SCAN_HOLD"] = np.where(valid, accum.scan_hold(res["m15_PHY_DP"]), np.nan).astype("float32")
        np.savez_compressed(prod / f"advection_{tag}.npz", glon=glon, glat=glat, **adv)
        for p in mosaic.PRODUCTS:
            mosaic.to_geotiff(prod / f"hourly_{p}_{tag}.tif", res[f"h_{p}"], glon, glat)
        mosaic.to_geotiff(prod / f"rqi_PHY_{tag}.tif", res["rqi_PHY"], glon, glat)
    return prod


def zipdir(folder, zpath):
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f in sorted(Path(folder).rglob("*")):
            if f.is_file():
                z.write(f, f.relative_to(Path(folder).parent))
    print(zpath, f"{os.path.getsize(zpath) / 1e6:.1f} MB  sha256 {colab.sha256(zpath)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--radar-dirs", nargs="*", default=[])
    ap.add_argument("--srtm", nargs="*", default=[])
    ap.add_argument("--data", default=None, help="existing teaching_data_v1 folder (skip the data build)")
    ap.add_argument("--out", default="build")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    data = Path(a.data) if a.data else build_data(a.radar_dirs, a.srtm, out)
    prod = build_products(data, out)
    zipdir(data, out / "teaching_data_v1.zip")
    zipdir(prod, out / "teaching_products_v1.zip")
