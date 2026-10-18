from _nb import *

F = "08_spatial_bias_correction.ipynb"
cells = [
header(F, "08", "การปรับแก้ความเอนเอียงเชิงพื้นที่ (spatial bias) ของฝนเรดาร์ด้วยฝนสถานี",
       "Spatial bias adjustment of radar rainfall with rain gauges — a gentle introduction",
       "🟢 ป.ตรี ปี 3–4 (แกนหลัก) และ 🔵 ป.โท (kriging, การถ่ายโอนข้ามลุ่มน้ำ)", "120 นาที", "07",
       "ฝนเรดาร์ PHY_DP 2 ชั่วโมง, ฝนสถานี สสน. ที่ผ่าน QC",
       "แผนที่ความเอนเอียงเชิงพื้นที่, ฝนเรดาร์หลังปรับแก้, ผล cross-validation"),
md("""
## 🎯 วัตถุประสงค์การเรียนรู้
1. **อธิบาย** ว่าทำไมความคลาดเคลื่อนของเรดาร์จึงแปรผันในพื้นที่ (spatial bias)
2. **คำนวณ** ค่าแก้แบบค่าเดียวทั้งพื้นที่ (mean-field bias, MFB)
3. **สร้าง** แผนที่ความเอนเอียงเชิงพื้นที่ด้วยการประมาณค่าในช่วงแบบ IDW และนำไปปรับแก้ฝนเรดาร์
4. **ตรวจสอบ** อย่างยุติธรรมด้วย leave-one-out cross-validation (LOOCV)
5. 🔵 **ใช้** variogram/ordinary kriging และทดสอบการถ่ายโอนค่าแก้ข้ามลุ่มน้ำ (วัง → ปิง)

## 1. หลักการ (อธิบายแบบง่าย)
### 1.1 ทำไม bias จึง "เป็นพื้นที่"
จาก notebook 07: เรดาร์ให้ฝน **ต่ำกว่า** สถานีเกือบทุกที่ แต่ "ต่ำกว่ามากน้อย" ไม่เท่ากัน
- ใกล้เรดาร์และไม่มีภูเขาบัง → ต่ำไป ~2 เท่า (การปรับเทียบ Z, Z–R)
- หลังภูเขาหรือไกลเรดาร์ (ลำคลื่นสูง) → ต่ำไป 5–10 เท่า
ถ้าความคลาดเคลื่อนมี **รูปแบบเชิงพื้นที่** (บริเวณใกล้กันผิดคล้ายกัน) เราก็ใช้สถานีที่อยู่ใกล้ ๆ ช่วยแก้ได้

### 1.2 ค่าแก้แบบค่าเดียว (Mean-Field Bias, MFB)
$$F = \\frac{\\sum_i G_i}{\\sum_i R_i}\\quad(\\text{เฉพาะสถานีที่มีฝนทั้งคู่}), \\qquad R'(x) = F\\cdot R(x)$$
ง่ายและเสถียร แต่ใช้ค่าเดียวทั้งภาคเหนือ

### 1.3 ค่าแก้เชิงพื้นที่ (local bias) ด้วย IDW
ที่สถานี $i$ คำนวณ **log-ratio** (ใช้ log เพราะความคลาดเคลื่อนเป็นแบบ "คูณ" และกันค่าติดลบ)
$$b_i = \\ln\\frac{G_i + \\varepsilon}{R_i + \\varepsilon}, \\qquad \\varepsilon = 0.1\\ \\text{มม.}$$
ประมาณค่า $b(x)$ ทุกพิกเซลด้วย **Inverse Distance Weighting**
$$b(x) = \\frac{\\sum_i w_i b_i}{\\sum_i w_i}, \\qquad w_i = \\frac{1}{d_i^{\\,p}}\\ (p = 2,\\ d_i \\le 50\\ \\text{กม.})$$
ไกลจากสถานีเกินรัศมี ให้ใช้ค่า MFB แทน แล้ว $R'(x) = R(x)\\,e^{b(x)}$

### 1.4 ตรวจสอบอย่างยุติธรรม: Leave-One-Out Cross-Validation
ถ้าใช้สถานีหนึ่งทั้ง "แก้" และ "ตรวจ" ผลจะดีเกินจริงเสมอ (ที่จุดนั้น $R' \\approx G$ โดยนิยาม)
LOOCV: เอาสถานี $k$ ออก → สร้างค่าแก้จากสถานีที่เหลือ → ประเมินที่สถานี $k$ → ทำซ้ำทุกสถานี

### 1.5 🔵 Spatial autocorrelation และ kriging
**Moran's I** วัดว่าค่าใกล้กันคล้ายกันหรือไม่ ($I > 0$ = เกาะกลุ่ม) ส่วน **ordinary kriging** ใช้ variogram
$\\gamma(h) = \\tfrac12 E[(b(x)-b(x+h))^2]$ หาน้ำหนักที่เหมาะสมที่สุดและให้ความแปรปรวนของการประมาณด้วย
"""),
setup_cell(("products",)),
md("## 2. เตรียมคู่ข้อมูล (PHY_DP กับฝนสถานี)"),
code("""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from radteach import gauges as gg, verify, bias, mosaic, config, plotting as pl, mapping as mp
cfg = config.load_all()
meta = gg.read_metadata(P.stations / "hii_station_metadata_202410.csv")
g = gg.add_utc(gg.quality_control(gg.read_station_csv(P.stations / "hii_hourly_rain_20241003_0519LT.csv"), meta=meta))
HOURS = {"00UTC": "2024-10-03 00:00", "10UTC": "2024-10-03 10:00"}
mos = {h: np.load(P.products / f"mosaic_20241003_{h}.npz") for h in HOURS}
glon, glat = mos["00UTC"]["glon"], mos["00UTC"]["glat"]
rows = []
for h, t0 in HOURS.items():
    gh = g[(g.start_utc == t0) & (g.qc == "OK")].merge(meta[["station_code", "longitude", "latitude", "basin", "basin_en"]], on="station_code")
    gh = gh.rename(columns={"longitude": "lon", "latitude": "lat", "rain_qc": "G"})
    rows.append(verify.pair(gh[["station_code", "start_utc", "lon", "lat", "G", "basin", "basin_en"]],
                            {"PHY_DP": mos[h]["h_PHY_DP"], "LS_MP": mos[h]["h_LS_MP"]}, glon, glat, window=3))
pts = pd.concat(rows, ignore_index=True).dropna(subset=["PHY_DP"])
print(f"{len(pts)} คู่ ({pts.station_code.nunique()} สถานี)")
"""),
md("## 3. MFB: ค่าแก้ค่าเดียว"),
code("""
mfb = {}
for h, t0 in HOURS.items():
    s = pts[pts.start_utc == t0]
    F, n = bias.mean_field_factor(s.G, s.PHY_DP)
    mfb[h] = F
    print(f"{h}: F = ΣG/ΣR = {F:.2f}  (จากสถานีที่มีฝนทั้งคู่ {n} สถานี) → เรดาร์ต่ำไป {F:.1f} เท่า")
"""),
md("## 4. ความเอนเอียงที่สถานี: มีรูปแบบเชิงพื้นที่หรือไม่"),
code("""
E = cfg["study_area"]["map_extent"]["north"]
cr, nr = pl.ratio_cmap(vmax_log2=4)
fig, axs = pl.subplots(1, 2, width="double", height_mm=110)
for ax, (h, t0) in zip(axs, HOURS.items()):
    s = pts[(pts.start_utc == t0) & ((pts.G >= 0.1) | (pts.PHY_DP >= 0.1))]
    b = bias.log_ratio(s.G, s.PHY_DP)
    sc = ax.scatter(s.lon, s.lat, c=b / np.log(2), cmap=cr, norm=nr, s=26, edgecolor="k", lw=0.3, zorder=13)
    mi = bias.morans_i(s.lon.values, s.lat.values, b)
    mp.basemap(ax, E, title=f"log₂(G/R) ที่สถานี  {h}", scale_km=100, north=(h == "00UTC"))
    ax.text(0.02, 0.98, f"Moran's I = {mi['I']:.2f} (p = {mi['p_value']:.3f})\\nn = {len(s)}", transform=ax.transAxes,
            fontsize=7.5, va="top", bbox=dict(fc="white", ec="0.6", lw=0.4), zorder=20)
cb = fig.colorbar(sc, ax=axs, orientation="horizontal", shrink=0.5, ticks=range(-4, 5))
cb.ax.set_xticklabels(["1/16", "1/8", "1/4", "1/2", "1", "2", "4", "8", "16"]); cb.set_label("สถานี/เรดาร์ (G/R): > 1 = เรดาร์ต่ำไป")
pl.panel_labels(axs); pl.save(fig, "08_gauge_log_ratio_map"); plt.show()
"""),
md("""
**แปลผล (ตามข้อมูลจริง):** แผนที่เห็นกลุ่มสีน้ำเงินเข้ม (เรดาร์ต่ำไปมาก) ในลุ่มยมตอนบนและลุ่มวัง แต่ **Moran's I ใกล้ 0 และ p > 0.05**
→ ยังไม่พบการเกาะกลุ่มที่มีนัยสำคัญทางสถิติในข้อมูลรายชั่วโมงชุดนี้ เพราะ (1) ฝนเบาและเป็นเซลล์เล็ก log-ratio จึงมีส่วนสุ่มสูง
(2) จำนวนสถานีที่มีฝนน้อย (37–82 สถานี) (3) bias ส่วนใหญ่เป็น "ค่าเดียวทั้งภาค" (การปรับเทียบ Z) ซึ่ง MFB จัดการได้
**ข้อคิด:** การตัดสินใจใช้ค่าแก้เชิงพื้นที่ควรทดสอบด้วยสถิติ ไม่ใช่ดูจากแผนที่อย่างเดียว — ลองทำซ้ำกับฝนสะสมทั้งเหตุการณ์

## 5. สร้างแผนที่ค่าแก้ด้วย IDW แล้วปรับแก้ฝนเรดาร์ (10 UTC)
"""),
code("""
h, t0 = "10UTC", HOURS["10UTC"]
s = pts[pts.start_utc == t0]
R = mos[h]["h_PHY_DP"]
R_mfb, _ = bias.adjust_field(R, glon, glat, s.lon, s.lat, s.G, s.PHY_DP, method="mfb")
R_idw, surf = bias.adjust_field(R, glon, glat, s.lon, s.lat, s.G, s.PHY_DP, method="idw", power=2, radius_km=50)
cm, nm = pl.rain_cmap()
fig, axs = pl.subplots(2, 2, width="double", height_mm=215)
for ax, arr, t in [(axs[0, 0], R, "(ก) เรดาร์ PHY_DP ก่อนปรับแก้"), (axs[1, 0], R_mfb, f"(ค) ปรับแก้ MFB (× {mfb[h]:.1f})"),
                   (axs[1, 1], R_idw, "(ง) ปรับแก้เชิงพื้นที่ IDW")]:
    im = mp.show_grid(ax, arr, glon, glat, cm, nm); mp.basemap(ax, E, title=t, scale_km=100, north=False, grid_step=1)
imf = mp.show_grid(axs[0, 1], np.where(np.isfinite(R), surf / np.log(2), np.nan), glon, glat, cr, nr)
axs[0, 1].scatter(s.lon, s.lat, s=5, c="k", zorder=14, label="สถานี")
mp.basemap(axs[0, 1], E, title="(ข) แผนที่ค่าแก้ e^b(x)", scale_km=100, grid_step=1); axs[0, 1].legend(loc="lower left")
fig.colorbar(im, ax=axs[1], orientation="horizontal", shrink=0.6, extend="max").set_label("ฝน 10–11 UTC (มม.)")
cb = fig.colorbar(imf, ax=axs[0, 1], orientation="horizontal", shrink=0.9, ticks=range(-4, 5))
cb.ax.set_xticklabels(["1/16", "1/8", "1/4", "1/2", "1", "2", "4", "8", "16"]); cb.set_label("ตัวคูณค่าแก้")
pl.save(fig, "08_adjustment_maps_10UTC"); plt.show()
for n_, a_ in [("ก่อนปรับแก้", R), ("MFB", R_mfb), ("IDW", R_idw)]:
    print(f"{n_:12s} ฝนเฉลี่ยทั้งโดเมน {np.nanmean(a_):.3f} มม. | สูงสุด {np.nanmax(a_):.1f} มม.")
mosaic.to_geotiff(P.rasters / "08_PHY_DP_idw_adjusted_10UTC.tif", R_idw, glon, glat)
"""),
md("""
**⚠️ ค่าสูงสุดหลังปรับแก้ (> 300 มม./ชม.) ไม่สมจริง!** เกิดจากตัวคูณขนาดใหญ่ไปคูณกับพิกเซลที่เรดาร์ให้ค่าสูงผิดปกติอยู่แล้ว
(เช่น clutter หรือลูกเห็บ) — นี่คือบทเรียนสำคัญ: **การปรับแก้ไม่ได้แก้ข้อมูลที่ผิดตั้งแต่ต้น** ในงานจริงต้อง (1) จำกัดตัวคูณ เช่น ≤ 5 เท่า
(2) ใช้ฝนสะสมช่วงยาวขึ้นในการหาค่าแก้ (3) ตรวจจุดผิดปกติก่อนปรับแก้

**ข้อสังเกต:** แผนที่ค่าแก้ (ข) สูงมาก (× 4 ถึง × 16) ในลุ่มน้ำวังและสาละวิน — บริเวณที่เรดาร์ "มองไม่ค่อยเห็น"
ส่วนบริเวณที่ไม่มีสถานีใกล้ ๆ ค่าแก้จะกลับไปเท่ากับ MFB ข้อควรระวัง: ค่าแก้ที่สูงมากจากสถานีเดียวอาจขยายความผิดพลาด (เช่นถ้าสถานีนั้นผิด)

## 6. ตรวจสอบด้วย Leave-One-Out Cross-Validation
"""),
code("""
cv = bias.loocv(pts, "PHY_DP", methods=("raw", "mfb", "idw"), power=2, radius_km=50)
tab = pd.DataFrame({m: verify.metrics(cv.G, cv[m]) for m in ["raw", "mfb", "idw"]}).T
tab = tab[["n", "PBIAS", "B", "RMSE", "MAE", "r", "r_wet", "KGE", "CSI"]].round(2)
tab.to_csv(P.tables / "08_loocv_metrics.csv"); display(tab)
fig, axs = pl.subplots(1, 3, width="double", height_mm=70)
for ax, m, t in zip(axs, ["raw", "mfb", "idw"], ["ก่อนปรับแก้", "MFB (LOOCV)", "IDW (LOOCV)"]):
    ax.scatter(np.clip(cv.G, 0.05, None), np.clip(cv[m], 0.05, None), s=12, c="#1b9e77", edgecolor="k", lw=0.3)
    pl.one_to_one(ax, 0.05, 40); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.04, 40); ax.set_ylim(0.04, 40)
    mm = verify.metrics(cv.G, cv[m]); ax.set_title(t, loc="left", fontweight="bold")
    ax.text(0.03, 0.97, f"PBIAS {mm['PBIAS']:.0f}%\\nRMSE {mm['RMSE']:.2f}\\nr(ฝน) {mm.get('r_wet', np.nan):.2f}", transform=ax.transAxes, va="top", fontsize=7,
            bbox=dict(fc="white", ec="0.6", lw=0.4))
    ax.set_xlabel(pl.L("gauge_mm"))
axs[0].set_ylabel(pl.L("radar_mm"))
pl.panel_labels(axs); pl.save(fig, "08_loocv_scatter"); plt.show()
"""),
md("""
**แปลผล:** MFB แก้ bias รวมได้เกือบหมด (PBIAS → ใกล้ 0) แต่ไม่ช่วยรูปแบบเชิงพื้นที่ ส่วน IDW ใช้ข้อมูลเพื่อนบ้านจึงลด RMSE และเพิ่ม r ได้
เมื่อสถานีใกล้กันมีความคลาดเคลื่อนคล้ายกัน — แต่ด้วยข้อมูล 2 ชั่วโมงและฝนแบบเซลล์เล็ก ผลยังผันผวนสูง

---
## ✏️ แบบฝึกหัด
"""),
exercise(1, "ug", "คำนวณ **MFB ของ 00 UTC ด้วยตนเอง** (ใช้เฉพาะคู่ที่ G ≥ 0.1 และ R ≥ 0.1) แล้วเทียบกับ `mfb['00UTC']`"),
todo("""
F00 = None
"""),
solution(1, """
s0 = pts[pts.start_utc == HOURS["00UTC"]]
k = (s0.G >= 0.1) & (s0.PHY_DP >= 0.1)
F00 = s0.G[k].sum() / s0.PHY_DP[k].sum()
print(f"F = {s0.G[k].sum():.1f} / {s0.PHY_DP[k].sum():.1f} = {F00:.2f}")
assert abs(F00 - mfb["00UTC"]) < 1e-9
"""),
exercise(2, "ug", "ทดลองเปลี่ยน **กำลัง p (1, 2, 3)** และ **รัศมี (25, 50, 100 กม.)** ของ IDW แล้วดู RMSE จาก LOOCV ค่าใดดีที่สุด และอธิบายผลของรัศมีที่เล็กเกินไป/ใหญ่เกินไป"),
todo("""
grid = None   # DataFrame power × radius → RMSE
"""),
solution(2, """
res = []
for p_ in [1, 2, 3]:
    for r_ in [25, 50, 100]:
        c_ = bias.loocv(pts, "PHY_DP", methods=("idw",), power=p_, radius_km=r_)
        res.append(dict(power=p_, radius_km=r_, RMSE=verify.metrics(c_.G, c_.idw)["RMSE"]))
grid = pd.DataFrame(res).pivot(index="power", columns="radius_km", values="RMSE").round(3)
display(grid)
print("รัศมีเล็ก: หลายพิกเซลไม่มีสถานี → ใช้ MFB; รัศมีใหญ่: ค่าแก้ถูกเกลี่ยจนเข้าใกล้ MFB; p สูง = เชื่อสถานีใกล้สุดมาก")
"""),
exercise(3, "ug", "**ถ่ายโอนค่าแก้ข้ามลุ่มน้ำ**: คำนวณ MFB จากสถานีในลุ่มน้ำวัง แล้วนำไปคูณกับเรดาร์ที่สถานีในลุ่มน้ำปิง PBIAS ของลุ่มปิงดีขึ้นหรือไม่ เทียบกับใช้ MFB ของทั้งภาค"),
todo("""
transfer = None
"""),
solution(3, """
wang = pts[pts.basin == "วัง"]; ping = pts[pts.basin == "ปิง"]
F_w, _ = bias.mean_field_factor(wang.G, wang.PHY_DP); F_all, _ = bias.mean_field_factor(pts.G, pts.PHY_DP)
transfer = pd.DataFrame({
    "ไม่ปรับแก้": verify.metrics(ping.G, ping.PHY_DP),
    f"MFB จากวัง (×{F_w:.1f})": verify.metrics(ping.G, F_w * ping.PHY_DP),
    f"MFB ทั้งภาค (×{F_all:.1f})": verify.metrics(ping.G, F_all * ping.PHY_DP)}).T[["n", "PBIAS", "RMSE", "r"]].round(2)
display(transfer)
print("ค่าแก้จากลุ่มวังสูงมาก (เรดาร์มองวังได้แย่) จึงทำให้ปิง 'เกิน' — ค่าแก้ขึ้นกับเรขาคณิตเรดาร์ ไม่ได้ขึ้นกับลุ่มน้ำอย่างเดียว")
"""),
exercise(4, "pg", "🔵 สร้าง **empirical variogram** ของ log-ratio (10 UTC) ฟิต exponential model แล้วทำ **ordinary kriging** แบบ LOOCV เปรียบเทียบ RMSE กับ IDW"),
todo("""
ok_rmse = None
"""),
solution(4, """
s = pts[(pts.start_utc == HOURS["10UTC"]) & ((pts.G >= 0.1) | (pts.PHY_DP >= 0.1))].reset_index(drop=True)
b = bias.log_ratio(s.G, s.PHY_DP)
ev = bias.empirical_variogram(s.lon.values, s.lat.values, b, bins_km=np.arange(0, 200, 15))
vg = bias.fit_variogram(ev); print("variogram:", {k: round(v, 3) for k, v in vg.items()})
pred = np.full(len(s), np.nan)
for i in range(len(s)):
    m = np.arange(len(s)) != i
    est, _ = bias.ordinary_kriging(s.lon[i], s.lat[i], s.lon[m].values, s.lat[m].values, b[m], vg)
    pred[i] = (s.PHY_DP[i] + bias.EPS) * np.exp(est[0]) - bias.EPS
ok_rmse = verify.metrics(s.G, np.clip(pred, 0, None))["RMSE"]
cv_idw = bias.loocv(s, "PHY_DP", methods=("raw", "mfb", "idw"))
print(f"RMSE (สถานีที่มีฝน 10 UTC, n = {len(s)}): ก่อนปรับแก้ {verify.metrics(s.G, s.PHY_DP)['RMSE']:.2f} | MFB {verify.metrics(cv_idw.G, cv_idw.mfb)['RMSE']:.2f} | "
      f"IDW {verify.metrics(cv_idw.G, cv_idw.idw)['RMSE']:.2f} | kriging {ok_rmse:.2f}")
print("nugget สูง (~45% ของ sill) → ความคลาดเคลื่อนระดับสถานีมีส่วนสุ่มมาก kriging ที่ไม่มีการจำกัดค่าจึงอาจแย่กว่า IDW ที่ถูกดึงเข้าหา MFB")
fig, ax = pl.subplots(width="single", height_mm=55)
hh = np.linspace(0, ev.h_km.max(), 100)
ax.plot(ev.h_km, ev.gamma, "o", label="empirical"); ax.plot(hh, bias.exp_variogram(hh, **vg), label="exponential")
ax.set_xlabel("ระยะห่าง (กม.)"); ax.set_ylabel("γ(h)"); ax.legend(); plt.show()
"""),
exercise(5, "pg", "🔵 ปรับแก้แบบ IDW กับ **LS_MP** แทน PHY_DP แล้วเทียบผล LOOCV หลังปรับแก้ ผลิตภัณฑ์เริ่มต้นที่ดีกว่ายังสำคัญหรือไม่หลังปรับแก้แล้ว"),
todo("""
compare = None
"""),
solution(5, """
pls = pts.dropna(subset=["LS_MP"])
c_ls = bias.loocv(pls, "LS_MP", methods=("raw", "idw"))
c_ph = bias.loocv(pls, "PHY_DP", methods=("raw", "idw"))
compare = pd.DataFrame({"LS_MP ก่อน": verify.metrics(c_ls.G, c_ls.raw), "LS_MP หลัง IDW": verify.metrics(c_ls.G, c_ls.idw),
                        "PHY_DP ก่อน": verify.metrics(c_ph.G, c_ph.raw), "PHY_DP หลัง IDW": verify.metrics(c_ph.G, c_ph.idw)}).T
display(compare[["n", "PBIAS", "RMSE", "r_wet", "CSI"]].round(2))
print("การปรับแก้แก้ได้แค่ bias ที่เปลี่ยนช้าในพื้นที่ แต่แก้ 'ฝนปลอม' หรือ 'ฝนที่หายไป' ระดับเซลล์ไม่ได้ → ผลิตภัณฑ์ต้นทางที่ดียังสำคัญ")
"""),
md("""
---
## 💡 คำถามชวนคิดและต่อยอด
- ปรับแก้รายชั่วโมงกับปรับแก้ด้วยฝนสะสม 24 ชั่วโมง/ทั้งเหตุการณ์ แบบใดเสถียรกว่า เพราะอะไร (ความแปรปรวนของ log-ratio)
- ควรใช้ RQI เป็นตัวแปรร่วมในการประมาณค่าแก้หรือไม่ (เช่น Kriging with External Drift, regression kriging, random forest)
- ถ้าต้นน้ำปิงมีสถานีน้อย การถ่ายโอนค่าแก้จากบริเวณที่ "เรดาร์มองเห็นคล้ายกัน" (RQI/ความสูงลำคลื่นใกล้เคียง) จะดีกว่าถ่ายโอนตามลุ่มน้ำหรือไม่
- ฝนที่ปรับแก้แล้วจะนำไปใช้กับแบบจำลองน้ำท่าหรือการเตือนภัยน้ำท่วมเชียงใหม่ได้อย่างไร ต้องบอกความไม่แน่นอนอย่างไร

## 📚 อ้างอิง
- Goudenhoofdt, E. & Delobbe, L. (2009). Evaluation of radar-gauge merging methods for quantitative precipitation estimates. *Hydrol. Earth Syst. Sci.* 13:195–203.
- Wilson, J.W. & Brandes, E.A. (1979). Radar measurement of rainfall—a summary. *Bull. Amer. Meteor. Soc.* 60:1048–1058.
- Sideris, I.V. et al. (2014). Real-time radar–rain-gauge merging using spatio-temporal co-kriging with external drift in the alpine terrain of Switzerland. *Q. J. R. Meteorol. Soc.* 140:1097–1111.
- Moran, P.A.P. (1950). Notes on continuous stochastic phenomena. *Biometrika* 37:17–23.
"""),
save_cell("08"),
]
build(cells, f"notebooks/{F}")
