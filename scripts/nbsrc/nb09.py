from _nb import *

F = "09_mosaic_methods_comparison.ipynb"
cells = [
header(F, "09", "เปรียบเทียบวิธีโมเสค: MAX / MEAN / NEAREST กับวิธีที่ใช้ RQI", "Comparing mosaicking methods: simple max/mean/nearest versus RQI-based",
       "🟢 ป.ตรี ปี 3–4 และ 🔵 ป.โท", "120 นาที", "04, 07", "UF 5 เรดาร์ × 2 ชั่วโมง (คำนวณสด ~3–5 นาที), ฝนสถานี สสน.",
       "โมเสค 5 วิธี × 2 ตัวประมาณ, ตาราง metrics แยกตามจำนวนเรดาร์ที่ซ้อนทับ, ดัชนีรอยต่อ, ฝนรายลุ่มน้ำ"),
md("""
## 🎯 วัตถุประสงค์การเรียนรู้
1. **อธิบาย** วิธีรวมเรดาร์หลายสถานี 5 แบบ และสมมติฐานเบื้องหลังของแต่ละแบบ
2. **สร้าง** โมเสคทั้ง 5 แบบจากค่าฝนรายเรดาร์ชุดเดียวกัน (ตรึงทุกอย่างไว้ ยกเว้นวิธีโมเสค)
3. **ประเมิน** ด้วยฝนสถานีอย่างยุติธรรม: ใช้คู่ข้อมูลชุดเดียวกัน และแยกตามจำนวนเรดาร์ที่มองเห็นพิกเซล
4. **วัด** ความต่อเนื่องเชิงพื้นที่ (รอยต่อระหว่างเรดาร์) และผลต่อฝนรายลุ่มน้ำ
5. **สรุป** ว่าวิธีใดดีกว่าในเงื่อนไขใด พร้อมความไม่แน่นอนทางสถิติ

## 1. หลักการ
ที่พิกเซล $x$ มีเรดาร์ $i = 1..n$ ที่มองเห็น ให้ค่าฝน $R_i$ และคุณภาพ $\\text{RQI}_i$ (notebook 02, 04)

| รหัส | สมการ | แนวคิด | สมมติฐานที่คาดไว้ |
|---|---|---|---|
| **MAX** | $R = \\max_i R_i$ | "เรดาร์ไหนเห็นฝนมากสุดน่าจะถูก" — ชดเชยลำคลื่นที่ลอยข้ามฝน | ฝนเกินจริง รับ clutter/ค่าสูงผิดปกติมาด้วย → FAR สูง |
| **MEAN** | $R = \\frac{1}{n}\\sum_i R_i$ | ทุกเรดาร์น่าเชื่อถือเท่ากัน | ฝนต่ำลง เพราะเฉลี่ยรวมลำคลื่นที่สูง/ถูกบัง |
| **NEAREST** | $R = R_{i^*},\\ i^*=\\arg\\min_i d_i$ | ใกล้ = ลำคลื่นต่ำ = ดี | ไม่รู้ว่าเรดาร์ใกล้อาจถูกภูเขาบัง; รอยต่อคม |
| **RQI_MAX** | $R = R_{i^*},\\ i^*=\\arg\\max_i \\text{RQI}_i$ | เลือกเรดาร์คุณภาพดีที่สุด (แนวคิด "best radar") | ดีในภูเขา แต่มีรอยต่อ |
| **RQI_WMEAN** | $R = \\dfrac{\\sum_i \\text{RQI}_i^{\\,k} R_i}{\\sum_i \\text{RQI}_i^{\\,k}}$ | MRMS (Zhang et al. 2016), $k = 1$ | ต่อเนื่อง bias และ RMSE ดีที่สุด |

- **MAX, MEAN, NEAREST จงใจไม่ใช้ RQI** — ใช้ทุกเรดาร์ที่ "เห็น" พิกเซล (อยู่ในรัศมี 240 กม. และผ่าน QC) นี่คือสิ่งที่เราต้องการทดสอบ
- **RQI_MAX และ RQI_WMEAN ใช้เฉพาะเรดาร์ที่ RQI > 0** พื้นที่ครอบคลุมจึงอาจน้อยกว่า
- 🔵 เมื่อ $k \\to \\infty$ RQI_WMEAN จะเข้าใกล้ RQI_MAX

### ทำให้การเปรียบเทียบยุติธรรม
1. **ตรึงทุกอย่าง** ยกเว้นวิธีโมเสค: QC ตัวประมาณ มุมยก กริด เวลา เหมือนกันหมด และทำซ้ำกับตัวประมาณ 2 แบบ (LS_MP, PHY_DP)
2. **คู่ข้อมูลชุดเดียวกัน**: ใช้เฉพาะสถานี-ชั่วโมงที่ *ทุกวิธี* มีค่า
3. **แยกตามจำนวนเรดาร์ที่ซ้อนทับ** (1, 2, ≥ 3) เพราะวิธีต่าง ๆ ให้ผลต่างกันเฉพาะในพื้นที่ซ้อนทับ
4. **บอกความไม่แน่นอน** ด้วย bootstrap ก่อนสรุปว่า "ดีกว่า"
"""),
setup_cell(("teaching_data", "products")),
md("## 2. คำนวณค่าฝนรายเรดาร์บนกริดเดียวกัน (2 ชั่วโมง, ~3–5 นาที)\n`process_hour(..., keep_radars=True)` เก็บค่าฝน RQI และระยะทางของเรดาร์แต่ละตัวไว้ แทนที่จะรวมทันที"),
code("""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from radteach import io, dem, mosaic, verify, gauges as gg, config, plotting as pl, mapping as mp
cfg = config.load_all(); g = cfg["study_area"]["grid"]
glon, glat = mosaic.make_grid(g["lon_min"], g["lon_max"], g["lat_min"], g["lat_max"], g["res_deg"])
z, tr, _ = dem.read(P.products / "dem_study_area.tif"); terrain = dem.interpolator(z, tr)
HOURS = {"00UTC": "2024-10-03 00:00", "10UTC": "2024-10-03 10:00"}
res = {}
for h in HOURS:
    scans = io.list_scans(P.radar / f"20241003_{h}")
    cal = pd.read_csv(P.products / f"calibration_20241003_{h}.csv", index_col=0)
    print(f"== {h}")
    res[h] = mosaic.process_hour(scans, terrain, cal, glon, glat, keep_radars=True, verbose=False)
RADARS = list(res["00UTC"]["radars"])
print("เรดาร์:", RADARS, "| ขนาดข้อมูลรายเรดาร์:", res["00UTC"]["radar_PHY_DP"].shape, "(เรดาร์, สแกน, ny, nx)")
"""),
md("## 3. สร้างโมเสค 5 วิธี × 2 ตัวประมาณ"),
code("""
METHODS = mosaic.COMBINE_METHODS
EST = ["LS_MP", "PHY_DP"]
mos = {h: {(p, m): mosaic.combine_hour(res[h], p, m) for p in EST for m in METHODS} for h in HOURS}
# จำนวนเรดาร์ที่มองเห็นพิกเซล (มีค่าในสแกนใดสแกนหนึ่ง)
nvis = {h: {p: np.isfinite(res[h][f"radar_{p}"]).any(1).sum(0) for p in EST} for h in HOURS}
# ตรวจ: RQI_WMEAN ต้องเท่ากับผลิตภัณฑ์ของ notebook 04
for h in HOURS:
    a, b = mos[h][("PHY_DP", "RQI_WMEAN")], res[h]["h_PHY_DP"]; k = np.isfinite(a) & np.isfinite(b)
    print(f"{h}: RQI_WMEAN ตรงกับ process_hour  (ผลต่างสูงสุด {np.nanmax(np.abs(a[k] - b[k])):.1e} มม.)")
cov = pd.DataFrame({(h, p): {m: 100 * np.isfinite(mos[h][(p, m)]).mean() for m in METHODS} for h in HOURS for p in EST}).round(1)
print("\\n% พื้นที่โดเมนที่มีค่า"); cov
"""),
md("## 4. แผนที่: จำนวนเรดาร์ที่ซ้อนทับ และโมเสค 5 วิธี (PHY_DP, 10 UTC)"),
code("""
from matplotlib.colors import ListedColormap, BoundaryNorm
E = cfg["study_area"]["map_extent"]["north"]
h = "10UTC"
fig, axs = pl.subplots(2, 3, width="double", height_mm=190)
cmv = ListedColormap(["#deebf7", "#9ecae1", "#4292c6", "#2171b5", "#08306b"])
nv = nvis[h]["PHY_DP"].astype(float); nv[nv == 0] = np.nan
im0 = mp.show_grid(axs[0, 0], nv, glon, glat, cmv, BoundaryNorm(np.arange(0.5, 6), 5))
mp.basemap(axs[0, 0], E, title="(a) จำนวนเรดาร์ที่มองเห็น", scale_km=None, grid_step=2)
cax = axs[0, 0].inset_axes([0.06, 0.07, 0.45, 0.03])
cb = fig.colorbar(im0, cax=cax, orientation="horizontal", ticks=range(1, 6)); cb.set_label("จำนวนเรดาร์", fontsize=7); cb.ax.xaxis.set_label_position("top"); cb.ax.tick_params(labelsize=7)
cm, nm = pl.rain_cmap()
for ax, m, lab in zip(axs.ravel()[1:], METHODS, "bcdef"):
    im = mp.show_grid(ax, mos[h][("PHY_DP", m)], glon, glat, cm, nm)
    mp.basemap(ax, E, title=f"({lab}) {m}", scale_km=None, north=False, grid_step=2)
fig.colorbar(im, ax=axs[1, :], orientation="horizontal", shrink=0.6, extend="max").set_label("ฝน PHY_DP 10–11 UTC (มม.)")
pl.save(fig, "09_mosaic_methods_maps"); plt.show()
"""),
code("""
fig, axs = pl.subplots(2, 4, width="double", height_mm=130)
dn = plt.Normalize(-3, 3); cmd = plt.get_cmap("RdBu_r", 12)
for i, hh in enumerate(HOURS):
    ref = mos[hh][("PHY_DP", "RQI_WMEAN")]
    for j, m in enumerate(["MAX", "MEAN", "NEAREST", "RQI_MAX"]):
        ax = axs[i, j]; d = mos[hh][("PHY_DP", m)] - ref
        im = mp.show_grid(ax, np.where(np.abs(d) > 0.05, d, np.nan), glon, glat, cmd, dn)
        mp.basemap(ax, E, title=m if i == 0 else None, scale_km=None, north=False, grid_step=2, radars=False)
        ax.text(0.03, 0.97, f"{hh[:2]} UTC", transform=ax.transAxes, va="top", fontsize=7.5, fontweight="bold")
        if j: ax.set_yticklabels([])
        if i == 0: ax.set_xticklabels([])
fig.colorbar(im, ax=axs, orientation="horizontal", shrink=0.5, extend="both").set_label("ผลต่างจาก RQI_WMEAN (มม.), PHY_DP")
pl.save(fig, "09_mosaic_differences"); plt.show()
"""),
md("""
**สังเกต:** ความต่างเกิดเฉพาะในพื้นที่ที่เรดาร์ซ้อนทับ — MAX สูงกว่า (สีแดง) โดยเฉพาะบริเวณที่มีเรดาร์หนึ่งให้ค่าสูงผิดปกติ
MEAN ต่ำกว่า (สีน้ำเงิน) เพราะเอาเรดาร์ที่ลำคลื่นสูง/ถูกบังมาเฉลี่ยด้วย ส่วน NEAREST กับ RQI_MAX ต่างกันเป็นแถบตามแนวที่ "เรดาร์ใกล้ไม่ใช่เรดาร์ดี"

## 5. เปรียบเทียบกับฝนสถานี (คู่ข้อมูลชุดเดียวกัน)
"""),
code("""
meta = gg.read_metadata(P.stations / "hii_station_metadata_202410.csv")
gq = gg.add_utc(gg.quality_control(gg.read_station_csv(P.stations / "hii_hourly_rain_20241003_0519LT.csv"), meta=meta))
rows = []
for h, t0 in HOURS.items():
    gh = gq[(gq.start_utc == t0) & (gq.qc == "OK")].merge(meta[["station_code", "longitude", "latitude", "basin"]], on="station_code")
    gh = gh.rename(columns={"longitude": "lon", "latitude": "lat", "rain_qc": "G"})[["station_code", "start_utc", "lon", "lat", "G", "basin"]]
    d = verify.pair(gh, {f"{p}|{m}": mos[h][(p, m)] for p in EST for m in METHODS}, glon, glat, window=3)
    for p in EST:
        d[f"nvis_{p}"] = verify.sample(nvis[h][p].astype(float), glon, glat, d.lon, d.lat, 1)
    rows.append(d)
pairs = pd.concat(rows, ignore_index=True)
pairs.to_csv(P.tables / "09_pairs.csv", index=False, encoding="utf-8-sig")
COLS = ["n", "PBIAS", "B", "RMSE", "r", "r_wet", "KGE", "POD", "FAR", "CSI", "HSS"]
def table(df, p):
    c = df.dropna(subset=[f"{p}|{m}" for m in METHODS])
    return pd.DataFrame({m: verify.metrics(c.G, c[f"{p}|{m}"]) for m in METHODS}).T.reindex(columns=COLS)
T = {p: table(pairs, p) for p in EST}
for p in EST:
    T[p].round(2).to_csv(P.tables / f"09_metrics_{p}.csv")
    print(f"\\n=== {p}  (คู่ข้อมูลชุดเดียวกัน n = {int(T[p].n.iloc[0])})"); display(T[p].round(2))
"""),
md("### 5.1 แยกตามจำนวนเรดาร์ที่ซ้อนทับ — วิธีโมเสคมีผลจริงตรงไหน"),
code("""
rows = []
for p in EST:
    c = pairs.dropna(subset=[f"{p}|{m}" for m in METHODS]).copy()
    c["ซ้อนทับ"] = pd.cut(c[f"nvis_{p}"], [0, 1.5, 2.5, 9], labels=["1", "2", "≥3"])
    for ov, cc in c.groupby("ซ้อนทับ", observed=True):
        for m in METHODS:
            rows.append({"ตัวประมาณ": p, "ซ้อนทับ": ov, "วิธี": m, **verify.metrics(cc.G, cc[f"{p}|{m}"])})
OV = pd.DataFrame(rows); OV.to_csv(P.tables / "09_metrics_by_overlap.csv", index=False, encoding="utf-8-sig")
fig, axs = pl.subplots(2, 3, width="double", height_mm=110, sharex=True)
colors = dict(zip(METHODS, ["#d62728", "#1f77b4", "#7f7f7f", "#ff7f0e", "#2ca02c"]))
for i, p in enumerate(EST):
    for j, (met, lim) in enumerate([("PBIAS", (-100, 0)), ("r_wet", (-0.2, 1)), ("FAR", (0, 1))]):
        ax = axs[i, j]; d = OV[OV["ตัวประมาณ"] == p]
        for k_, m in enumerate(METHODS):
            s = d[d["วิธี"] == m]
            ax.bar(np.arange(len(s)) + (k_ - 2) * 0.16, s[met], 0.16, color=colors[m], label=m)
        nn = d[d["วิธี"] == "RQI_WMEAN"].set_index("ซ้อนทับ").n
        ax.set_xticks(range(3)); ax.set_xticklabels([f"{c} เรดาร์\\n(n={int(nn.get(c, 0))})" for c in ["1", "2", "≥3"]]); ax.set_ylim(*lim)
        ax.set_title(f"{p}: {met}", loc="left", fontweight="bold", fontsize=8.5); ax.grid(alpha=0.3, lw=0.4, axis="y")
h_, l_ = axs[0, 0].get_legend_handles_labels()
fig.legend(h_, l_, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.06))
pl.panel_labels(axs); pl.save(fig, "09_metrics_by_overlap"); plt.show()
OV.pivot_table(index=["ตัวประมาณ", "ซ้อนทับ"], columns="วิธี", values="r_wet", observed=True)[METHODS].round(2)
"""),
md("""
### 5.2 ความได้เปรียบของ RQI มาจากไหน: "ตัดข้อมูลคุณภาพต่ำ" หรือ "เลือก/ถ่วงน้ำหนักระหว่างเรดาร์"
แยกผลสองส่วนออกจากกัน: นำ MAX และ MEAN มาใช้ **เฉพาะเรดาร์ที่ RQI > 0** (ตัดลำคลื่นที่ถูกบังเกือบหมด/สูงเกินชั้นหลอมเหลว/ไม่ผ่านเกณฑ์ ρhv ทิ้ง) แต่ยังไม่ถ่วงน้ำหนัก
- ถ้า `MAX|Q>0` ดีขึ้นจนใกล้ RQI_WMEAN → ประโยชน์หลักมาจาก **การคัดกรองคุณภาพ**
- ถ้ายังห่างกันมาก → ประโยชน์มาจาก **การถ่วงน้ำหนักตาม RQI**
"""),
code("""
def combine_masked(r, p, method):
    grp = "LS" if p == "LS_MP" else "PHY"
    R = np.where(r[f"radar_rqi_{grp}"] > 0, r[f"radar_{p}"], np.nan)
    return mosaic.hourly_from_15min(np.array([mosaic.combine(R[:, t], method=method) * 0.25 for t in range(R.shape[1])]))
rows = []
for p in EST:
    ext = {}
    for m in ["MAX", "MEAN"]:
        vals = []
        for hh, t0 in HOURS.items():
            s = pairs[pairs.start_utc == t0]
            vals.append(pd.Series(verify.sample(combine_masked(res[hh], p, m), glon, glat, s.lon, s.lat, 3), index=s.index))
        ext[f"{m}|Q>0"] = pd.concat(vals)
    c = pairs.assign(**{f"{p}|{k}": v for k, v in ext.items()}).dropna(subset=[f"{p}|{m}" for m in METHODS + list(ext)])
    for m in ["MAX", "MAX|Q>0", "MEAN", "MEAN|Q>0", "RQI_WMEAN"]:
        rows.append({"ตัวประมาณ": p, "วิธี": m, **verify.metrics(c.G, c[f"{p}|{m}"])})
DEC = pd.DataFrame(rows).set_index(["ตัวประมาณ", "วิธี"]).reindex(columns=COLS)
DEC.round(2).to_csv(P.tables / "09_quality_mask_effect.csv", encoding="utf-8-sig"); DEC.round(2)
"""),
md("### 5.3 ความไม่แน่นอน: bootstrap เทียบกับ RQI_WMEAN"),
code("""
rng = np.random.default_rng(1)
rows = []
for p in EST:
    c = pairs.dropna(subset=[f"{p}|{m}" for m in METHODS])
    by = [d for _, d in c.groupby("station_code")]
    boots = {m: [] for m in METHODS if m != "RQI_WMEAN"}
    for _ in range(500):
        s = pd.concat([by[i] for i in rng.integers(0, len(by), len(by))])
        ref_m = verify.metrics(s.G, s[f"{p}|RQI_WMEAN"])
        for m in boots:
            mm = verify.metrics(s.G, s[f"{p}|{m}"])
            boots[m].append((mm["RMSE"] - ref_m["RMSE"], mm.get("r_wet", np.nan) - ref_m.get("r_wet", np.nan)))
    for m, v in boots.items():
        v = np.array(v)
        lo, hi = np.nanpercentile(v[:, 0], [2.5, 97.5]); lo2, hi2 = np.nanpercentile(v[:, 1], [2.5, 97.5])
        rows.append({"ตัวประมาณ": p, "วิธี": m, "ΔRMSE 95% CI": f"[{lo:+.2f}, {hi:+.2f}]",
                     "Δr_wet 95% CI": f"[{lo2:+.2f}, {hi2:+.2f}]",
                     "สรุป RMSE": "แย่กว่า RQI_WMEAN" if lo > 0 else ("ดีกว่า" if hi < 0 else "ไม่ต่างอย่างมีนัยสำคัญ")})
BS = pd.DataFrame(rows); BS.to_csv(P.tables / "09_bootstrap.csv", index=False, encoding="utf-8-sig"); BS
"""),
md("""
## 6. ความต่อเนื่องเชิงพื้นที่: ดัชนีรอยต่อ (seam index)
รอยต่อเกิดเมื่อพิกเซลติดกันได้ค่าจากเรดาร์ต่างตัว นิยาม
$$\\text{Seam index} = \\frac{\\operatorname{median}|R(x) - R(x')|\\ \\text{บนแนวที่เรดาร์ที่เลือกเปลี่ยน}}{\\operatorname{median}|R(x) - R(x')|\\ \\text{ในพื้นที่ซ้อนทับอื่น ๆ}}$$
โดย $x, x'$ เป็นพิกเซลติดกันที่มีฝน ค่า ≈ 1 = ต่อเนื่อง, ค่า > 1 = มีรอยต่อชัด (แนวเปลี่ยนเรดาร์ใช้ทั้งของ NEAREST และ RQI_MAX)
"""),
code("""
def seam_index(R, boundary, overlap, wet=0.1):
    out = []
    for ax_ in (0, 1):
        a, b = (R[:-1, :], R[1:, :]) if ax_ == 0 else (R[:, :-1], R[:, 1:])
        bd = boundary[:-1, :] | boundary[1:, :] if ax_ == 0 else boundary[:, :-1] | boundary[:, 1:]
        ov = (overlap[:-1, :] & overlap[1:, :]) if ax_ == 0 else (overlap[:, :-1] & overlap[:, 1:])
        ok = np.isfinite(a) & np.isfinite(b) & ((a > wet) | (b > wet)) & ov
        out.append((np.abs(a - b)[ok & bd], np.abs(a - b)[ok & ~bd]))
    on = np.concatenate([o[0] for o in out]); off = np.concatenate([o[1] for o in out])
    return np.median(on) / np.median(off) if len(on) > 20 and len(off) > 20 else np.nan
rows = []
for hh in HOURS:
    for p in EST:
        R, Q, D = res[hh][f"radar_{p}"][:, 0], res[hh][f"radar_rqi_{'LS' if p == 'LS_MP' else 'PHY'}"][:, 0], res[hh]["radar_dist"]
        seen = np.isfinite(R)
        sel_near = np.argmin(np.where(seen, np.nan_to_num(D, nan=np.inf), np.inf), 0)
        sel_rqi = np.argmax(np.where(seen, Q, -1), 0)
        bnd = np.zeros_like(seen[0])
        for sel in (sel_near, sel_rqi):
            bnd[:-1, :] |= sel[:-1, :] != sel[1:, :]; bnd[:, :-1] |= sel[:, :-1] != sel[:, 1:]
        overlap = seen.sum(0) >= 2
        for m in METHODS:
            rows.append({"ชั่วโมง": hh, "ตัวประมาณ": p, "วิธี": m, "seam index": seam_index(mos[hh][(p, m)], bnd, overlap)})
SI = pd.DataFrame(rows).pivot_table(index=["ชั่วโมง", "ตัวประมาณ"], columns="วิธี", values="seam index")[METHODS].round(2)
SI.to_csv(P.tables / "09_seam_index.csv", encoding="utf-8-sig"); SI
"""),
md("## 7. ผลต่อฝนรายลุ่มน้ำ (ปริมาณน้ำฝนที่เข้าสู่ลุ่มน้ำ)"),
code("""
from rasterio.features import rasterize
from rasterio.transform import from_origin
mb, sb, _ = mp.boundaries()
T_ = from_origin(glon[0] - 0.005, glat[0] + 0.005, 0.01, 0.01)
zones = {r.MBASIN_E: rasterize([(r.geometry, 1)], out_shape=(len(glat), len(glon)), transform=T_).astype(bool)
         for _, r in mb[mb.MB_CODE.isin(["06", "07", "08", "09"])].iterrows()}
_hw = sb[sb.SB_CODE.isin(cfg["study_area"]["ping_headwater_subbasins"])].geometry
zones["Ping headwater"] = rasterize([((_hw.union_all() if hasattr(_hw, "union_all") else _hw.unary_union), 1)],
                                    out_shape=(len(glat), len(glon)), transform=T_).astype(bool)
rows = [{"ชั่วโมง": hh, "ลุ่มน้ำ": zn, "วิธี": m, "ฝนเฉลี่ย (มม.)": np.nanmean(mos[hh][("PHY_DP", m)][zm]),
         "ครอบคลุม (%)": 100 * np.isfinite(mos[hh][("PHY_DP", m)][zm]).mean()}
        for hh in HOURS for zn, zm in zones.items() for m in METHODS]
BM = pd.DataFrame(rows); BM.to_csv(P.tables / "09_basin_means.csv", index=False, encoding="utf-8-sig")
BM.pivot_table(index=["ชั่วโมง", "ลุ่มน้ำ"], columns="วิธี", values="ฝนเฉลี่ย (มม.)")[METHODS].round(3)
"""),
md("""
## 8. สรุปและแปลผล (อภิปรายในชั้นเรียน)
ตัวเลขด้านล่างมาจากการรันตอนสร้างบทเรียน (2 ชั่วโมง, 3 ต.ค. 2567) — ของนิสิตควรได้ใกล้เคียงกัน ให้ตรวจกับตารางข้างบน

**ตอบคำถาม "วิธีไหนดีกว่า"** — กับตัวประมาณ PHY_DP (n ≈ 409 คู่):

| | MAX | MEAN | NEAREST | RQI_MAX | RQI_WMEAN |
|---|---|---|---|---|---|
| r_wet | 0.24 | 0.38 | 0.23 | **0.59** | **0.59** |
| RMSE (มม.) | 1.67 | 1.67 | 1.68 | **1.61** | **1.61** |
| FAR | 0.24 | 0.23 | 0.24 | **0.19** | **0.19** |
| PBIAS (%) | **−60** | −72 | −62 | −73 | −73 |

1. **วิธีที่ใช้ RQI ชนะในด้านรูปแบบเชิงพื้นที่** (r, r_wet สูงกว่ามาก) และเตือนผิดน้อยกว่า (FAR) แต่ **MAX ให้ bias ติดลบน้อยที่สุด**
   — ไม่ใช่เพราะแม่นกว่า แต่เพราะหยิบค่าสูงสุดเสมอ (รวมค่าผิดปกติ) จึงได้ POD สูงและ r ต่ำ → อย่าตัดสินวิธีด้วย bias ตัวเดียว
2. **ที่มาของความได้เปรียบ (หัวข้อ 5.2)**: `MAX|Q>0` และ `MEAN|Q>0` ดีพอ ๆ กับ RQI_WMEAN (r_wet ≈ 0.58–0.59)
   → ประโยชน์เกือบทั้งหมดมาจาก **การตัดข้อมูลคุณภาพต่ำ (RQI = 0) ออก** ไม่ใช่จากสูตรถ่วงน้ำหนัก
3. **แยกตามจำนวนเรดาร์ (หัวข้อ 5.1)**: ความต่างใหญ่ที่สุดในพื้นที่ **1–2 เรดาร์** ซึ่งมักเป็นพื้นที่ภูเขาที่ลำคลื่นถูกบัง
   ส่วนพื้นที่ ≥ 3 เรดาร์ (กลางโดเมน ภูมิประเทศราบกว่า) ทุกวิธีให้ r_wet ≈ 0.6 ใกล้กัน และ n มีเพียง ~46 คู่ — สรุปไม่ได้
4. **MEAN** ให้ฝนต่ำที่สุดเกือบทุกกรณี (นำลำคลื่นสูง/ถูกบังมาเฉลี่ยด้วย) และเป็นวิธีเดียวที่ RMSE แย่กว่า RQI_WMEAN อย่างมีนัยสำคัญ (bootstrap)
5. **NEAREST** ใกล้เคียง MAX และ ~ 1 ใน 4 ของพื้นที่ซ้อนทับ "เรดาร์ใกล้สุด ≠ เรดาร์คุณภาพดีสุด" (แบบฝึกหัด 3)
6. **RQI_MAX กับ RQI_WMEAN** ให้ metrics ที่สถานีเกือบเท่ากัน แต่ RQI_WMEAN มี **seam index ต่ำกว่า** (รอยต่อน้อยกว่า) จึงเหมาะกับงานอุทกวิทยามากกว่า
7. **ตัวประมาณสำคัญกว่าวิธีโมเสค**: LS_MP (มี clutter) ทุกวิธีได้ r_wet ≈ 0 และ bootstrap ไม่พบความต่างที่มีนัยสำคัญ — โมเสคที่ดีแก้ข้อมูลต้นทางที่ไม่ดีไม่ได้
8. **ผลต่อฝนรายลุ่มน้ำ**: ในชั่วโมงฝนเบา วิธีโมเสคทำให้ฝนเฉลี่ยลุ่มต่างกันได้ 30–50 % (เช่น ลุ่มยม 00 UTC: MAX ≈ 0.70 vs RQI_WMEAN ≈ 0.46 มม.)
9. **ข้อจำกัด**: มีเพียง 2 ชั่วโมง ช่วง CI ของ MAX/NEAREST คร่อม 0 → ยังเรียกได้แค่ "มีแนวโน้มดีกว่า" ต้องทดสอบทั้งเหตุการณ์
10. **RQI ลดพื้นที่ครอบคลุม** (LS ~96 % → ~58 %) — เป็นการแลกระหว่าง "มีค่าทุกที่" กับ "ค่าที่เชื่อถือได้" ซึ่งต้องบอกผู้ใช้ข้อมูล
11. 🔵 หลังปรับแก้ MFB (แบบฝึกหัด 5) ช่องว่างยิ่งกว้างขึ้น: RQI methods RMSE ลดลง ~1.24 แต่ MAX/NEAREST กลับแย่ลง (> 2.2) เพราะรูปแบบเชิงพื้นที่ผิดถูกขยายตาม

---
## ✏️ แบบฝึกหัด
"""),
exercise(1, "ug", "เลือกพิกเซลที่มีเรดาร์มองเห็น ≥ 3 ตัวและมีฝน (โค้ดตั้งต้นให้แล้ว) คำนวณ **MAX, MEAN และ RQI_WMEAN ด้วยมือ** จากค่าฝนและ RQI ของสแกนแรก แล้วเทียบกับ `mosaic.combine`"),
todo("""
R0 = res["10UTC"]["radar_PHY_DP"][:, 0]; Q0 = res["10UTC"]["radar_rqi_PHY"][:, 0]
cand = np.argwhere((np.isfinite(R0).sum(0) >= 3) & (np.nanmax(np.where(np.isfinite(R0), R0, 0), 0) > 5))
iy, ix = cand[len(cand) // 2]
print("พิกเซล", iy, ix, "| R:", np.round(R0[:, iy, ix], 2), "| RQI:", np.round(Q0[:, iy, ix], 2))
my_max = None; my_mean = None; my_wmean = None
"""),
solution(1, """
r_, q_ = R0[:, iy, ix], Q0[:, iy, ix]; s = np.isfinite(r_)
my_max = np.max(r_[s]); my_mean = np.mean(r_[s])
w = np.where(s & (q_ > 0), q_, 0); my_wmean = np.sum(w * np.nan_to_num(r_)) / np.sum(w)
for nm_, mine in [("MAX", my_max), ("MEAN", my_mean), ("RQI_WMEAN", my_wmean)]:
    ref_ = mosaic.combine(R0[:, iy:iy + 1, ix:ix + 1], Q0[:, iy:iy + 1, ix:ix + 1], method=nm_)[0, 0]
    print(f"{nm_:10s} มือ {mine:.3f} | ฟังก์ชัน {ref_:.3f}")
    assert abs(mine - ref_) < 1e-6
"""),
exercise(2, "ug", "จากตาราง `T` วิธีใดมี **FAR สูงสุด** และ **POD สูงสุด** สำหรับ PHY_DP (ดูทั้งคอลัมน์ เพราะบางวิธีค่าใกล้กันมาก) ทั้งสองค่าบอกอะไรเกี่ยวกับวิธีที่ไม่ใช้ RQI"),
todo("""
far_max = None; pod_max = None
"""),
solution(2, """
far_max = T["PHY_DP"].FAR.idxmax(); pod_max = T["PHY_DP"].POD.idxmax()
print("FAR สูงสุด:", far_max, "| POD สูงสุด:", pod_max); display(T["PHY_DP"][["POD", "FAR", "CSI"]].round(3))
print("MAX/NEAREST/MEAN มี FAR ≈ 0.23–0.24 สูงกว่าวิธี RQI (≈ 0.19) และ MAX มี POD สูงสุด:")
print("MAX ตรวจพบฝนได้มาก แต่แลกกับการ 'เตือนผิด' มาก เพราะหยิบค่าที่สูงที่สุดเสมอ รวมถึงค่าผิดปกติจากลำคลื่นถูกบัง/clutter")
"""),
exercise(3, "ug", "คิดเป็นกี่เปอร์เซ็นต์ของพื้นที่ซ้อนทับ (≥ 2 เรดาร์) ที่ **เรดาร์ใกล้ที่สุด ≠ เรดาร์ RQI สูงสุด** (PHY, 10 UTC สแกนแรก) บริเวณเหล่านั้นอยู่ที่ไหน เพราะอะไร"),
todo("""
pct_diff = None
"""),
solution(3, """
R0 = res["10UTC"]["radar_PHY_DP"][:, 0]; Q0 = res["10UTC"]["radar_rqi_PHY"][:, 0]; D0 = res["10UTC"]["radar_dist"]
seen = np.isfinite(R0); ov = seen.sum(0) >= 2
near = np.argmin(np.where(seen, np.nan_to_num(D0, nan=np.inf), np.inf), 0); best = np.argmax(np.where(seen, Q0, -1), 0)
diff = ov & (near != best); pct_diff = 100 * diff.sum() / ov.sum()
print(f"{pct_diff:.1f}% ของพื้นที่ซ้อนทับ")
fig, ax = pl.subplots(width="one_half", height_mm=95)
mp.show_grid(ax, np.where(diff, 1.0, np.nan), glon, glat, ListedColormap(["#d62728"]), BoundaryNorm([0.5, 1.5], 1))
mp.basemap(ax, E, title="เรดาร์ใกล้สุด ≠ เรดาร์ RQI สูงสุด", scale_km=100); plt.show()
print("มักอยู่หลังภูเขาของเรดาร์ที่ใกล้ (ถูกบัง) หรือที่ระยะไกลที่ลำคลื่นเรดาร์ใกล้อยู่สูงกว่าเรดาร์อีกตัวที่มุมยกต่ำกว่า")
"""),
exercise(4, "pg", "🔵 ทดลอง **RQI^k** ด้วย k = 0, 0.5, 1, 2, 4, 8 (k = 0 คือค่าเฉลี่ยของเรดาร์ที่ RQI > 0) วาดกราฟ RMSE และ r_wet ของ PHY_DP กับ k แล้วอภิปรายว่า k ที่เหมาะสมควรเป็นเท่าใด"),
todo("""
ksweep = None
"""),
solution(4, """
rows = []
for k in [0, 0.5, 1, 2, 4, 8]:
    fl = {hh: mosaic.combine_hour(res[hh], "PHY_DP", "RQI_WMEAN", k=k) for hh in HOURS}
    vals = []
    for hh, t0 in HOURS.items():
        s = pairs[pairs.start_utc == t0]
        vals.append(pd.DataFrame({"G": s.G.values, "R": verify.sample(fl[hh], glon, glat, s.lon, s.lat, 3)}))
    v = pd.concat(vals); mm = verify.metrics(v.G, v.R)
    rows.append({"k": k, "RMSE": mm["RMSE"], "r_wet": mm.get("r_wet"), "PBIAS": mm["PBIAS"]})
ksweep = pd.DataFrame(rows); display(ksweep.round(3))
fig, ax = pl.subplots(width="single", height_mm=55)
ax.plot(ksweep.k, ksweep.RMSE, "o-", label="RMSE"); ax2 = ax.twinx(); ax2.plot(ksweep.k, ksweep.r_wet, "s--", color="#d95f02", label="r_wet")
ax.set_xlabel("k"); ax.set_ylabel("RMSE (มม.)"); ax2.set_ylabel("r_wet"); fig.legend(loc="upper center", ncol=2); plt.show()
"""),
exercise(5, "pg", "🔵 ปรับแก้ทุกวิธีด้วย **MFB แบบ LOOCV** (ใช้ `radteach.bias.loocv`) แล้วดูว่าอันดับของวิธีโมเสค (ตาม RMSE) ยังเหมือนเดิมหรือไม่ — ถ้าเปลี่ยน แปลว่าอะไรสำหรับระบบปฏิบัติการที่ปรับแก้ด้วยสถานีเสมอ"),
todo("""
after_mfb = None
"""),
solution(5, """
from radteach import bias
c = pairs.dropna(subset=[f"PHY_DP|{m}" for m in METHODS]).copy()
rows = []
for m in METHODS:
    cv = bias.loocv(c.rename(columns={f"PHY_DP|{m}": "X"}), "X", methods=("raw", "mfb"))
    rows.append({"วิธี": m, "RMSE ก่อน": verify.metrics(cv.G, cv.raw)["RMSE"], "RMSE หลัง MFB": verify.metrics(cv.G, cv.mfb)["RMSE"],
                 "r_wet หลัง MFB": verify.metrics(cv.G, cv.mfb).get("r_wet")})
after_mfb = pd.DataFrame(rows).set_index("วิธี").round(3); display(after_mfb)
print("MFB แก้ bias รวมได้ทุกวิธี จึงเหลือความต่างที่ 'รูปแบบเชิงพื้นที่' (r) — วิธีที่ใช้ RQI ยังได้เปรียบถ้า r_wet สูงกว่า")
"""),
md("""
---
## 💡 คำถามชวนคิดและต่อยอด
- ถ้าเรดาร์สองตัวปรับเทียบ Z ต่างกัน 3 dB วิธีใดได้รับผลกระทบมากที่สุด (ลองบวก 3 dB ให้ LMP แล้วรันซ้ำ) และควรปรับเทียบระหว่างเรดาร์ก่อนโมเสคอย่างไร
- JMA ใช้ค่าสูงสุดของเรดาร์ที่ซ้อนทับเพื่อชดเชยลำคลื่นที่ลอยข้ามฝน — เงื่อนไขใดที่ MAX อาจเป็นทางเลือกที่ถูกต้อง (เช่น ฝนชั้นตื้นบนภูเขา)
- สำหรับการเตือนภัยน้ำท่วมเชียงใหม่ ควรเลือกวิธีโมเสคจาก metrics ใด (POD? CSI? ปริมาตรฝนรายลุ่ม?) และต้องสื่อสารพื้นที่ที่ "ไม่มีข้อมูลเชื่อถือได้" อย่างไร
- 🔵 ทดลองใส่น้ำหนักระยะทาง $w_i = \\text{RQI}_i \\cdot e^{-d_i/d_0}$ หรือใช้ RQI ร่วมกับความแปรปรวนของความคลาดเคลื่อนที่ประมาณจากสถานี (optimal combination)

## 📚 อ้างอิง
- Zhang, J. et al. (2016). Multi-Radar Multi-Sensor (MRMS) quantitative precipitation estimation. *BAMS* 97:621–638.
- Zhang, J. et al. (2011). National Mosaic and Multi-Sensor QPE (NMQ) system. *BAMS* 92:1321–1338.
- Makihara, Y. (2000). Algorithms for precipitation nowcasting focused on detailed analysis using radar and raingauge data. *Tech. Rep. MRI* 39:63–111.
- Germann, U. et al. (2022). Weather radar in complex orography. *Remote Sensing* 14:503.
- Fornasiero, A. et al. (2005). Effects of propagation conditions on radar beam-ground interaction: impact on data quality. *Adv. Geosci.* 2:201–208.
"""),
save_cell("09"),
]
build(cells, f"notebooks/{F}")
