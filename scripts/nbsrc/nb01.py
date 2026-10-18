from _nb import *

F = "01_radar_basics_qc.ipynb"
cells = [
header(F, "01", "ข้อมูลเรดาร์ dual-polarization และการควบคุมคุณภาพ", "Dual-pol radar variables and quality control",
       "🟢 ป.ตรี ปี 3–4 และ 🔵 ป.โท", "90–120 นาที", "00", "UF ของเรดาร์พิษณุโลก (PHS) และอีก 4 สถานี เวลา 00:00 UTC",
       "แผนที่ PPI 6 ตัวแปร, ตารางค่าปรับเทียบ ZDR/ΦDP, กราฟ ΦDP–KDP ตามแนวลำคลื่น"),
md("""
## 🎯 วัตถุประสงค์การเรียนรู้
1. **อธิบาย** ความหมายทางกายภาพของตัวแปรเรดาร์ Z, ZDR, ρhv, ΦDP, KDP, V และ W
2. **แยก** สัญญาณฝน (meteorological echo) ออกจากสิ่งรบกวน (ground clutter, AP, แมลง/นก, สัญญาณรบกวน) ด้วยเกณฑ์อย่างง่าย
3. **ปรับเทียบ** ZDR offset และ system differential phase จากข้อมูลจริง
4. **คำนวณ** KDP จาก ΦDP และอธิบายการแก้การลดทอน (attenuation) ด้วยวิธี ZPHI
5. **ประเมิน** ผลของแต่ละขั้น QC ต่อพื้นที่ฝน

## 1. หลักการ
### 1.1 สมการเรดาร์และค่าการสะท้อนกลับ
กำลังที่สะท้อนกลับจากเม็ดฝนในปริมาตรสุ่มตัวอย่างลดลงตามระยะยกกำลังสอง

$$P_r = \\frac{C\\,|K_w|^2\\,Z}{r^2}, \\qquad Z = \\int N(D)\\,D^6\\,dD \\;[\\text{mm}^6\\,\\text{m}^{-3}], \\qquad \\text{dBZ} = 10\\log_{10} Z$$

$C$ = ค่าคงที่ของเรดาร์, $|K_w|^2 \\approx 0.93$ สำหรับน้ำ, $N(D)$ = การกระจายขนาดเม็ดฝน (DSD) เพราะ $Z \\propto D^6$
เม็ดใหญ่ไม่กี่เม็ดจึงคุมค่า Z — นี่คือเหตุผลที่ความสัมพันธ์ Z–R มีความไม่แน่นอนสูง (notebook 03)

### 1.2 ตัวแปร dual-polarization
เรดาร์ส่งคลื่นทั้งแนวนอน (H) และแนวตั้ง (V) เม็ดฝนขนาดใหญ่แบนคล้ายขนมปังแฮมเบอร์เกอร์ จึงให้สัญญาณ H มากกว่า V

| ตัวแปร | นิยาม | ฝน | สิ่งรบกวน |
|---|---|---|---|
| **Z** (dBZ) | $10\\log_{10}Z_h$ | 15–55 | clutter อาจ > 50 |
| **ZDR** (dB) | $10\\log_{10}(Z_h/Z_v)$ | 0–3 (เพิ่มตามขนาดเม็ด) | แมลง/นก 3–8 |
| **ρhv** | สหสัมพันธ์ระหว่างสัญญาณ H และ V | > 0.97 | < 0.85 |
| **ΦDP** (°) | ความต่างเฟส H–V สะสมตามระยะ | เพิ่มขึ้นเรียบ ๆ | กระโดด/สุ่ม |
| **KDP** (° km⁻¹) | $\\tfrac{1}{2}\\,d\\Phi_{DP}/dr$ | 0–5 ในฝนหนัก | – |
| **V** (m s⁻¹) | ความเร็วดอปเปลอร์ตามแนวลำคลื่น | ตามลม | ≈ 0 (clutter บนพื้น) |
| **W** (m s⁻¹) | spectrum width | 1–4 | ≈ 0 (clutter) |

**ทำไม KDP สำคัญ:** KDP วัดจาก *เฟส* ไม่ใช่ *กำลัง* จึงไม่ได้รับผลจาก (1) การปรับเทียบ Z ผิด (2) การลดทอนของสัญญาณ (3) ลำคลื่นถูกบังบางส่วน
"""),
setup_cell(("teaching_data",)),
md("## 2. อ่านข้อมูลเรดาร์พิษณุโลก (PHS) เวลา 00:00 UTC"),
code("""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from radteach import io, qc, plotting as pl, mapping as mp, config
scans = io.list_scans(P.radar / "20241003_00UTC")
f = scans[(scans.radar == "PHS") & (scans.time == "2024-10-03 00:00")].path.iloc[0]
radar = io.read_radar(f)
sw = radar.extract_sweeps([0])                       # มุมยกต่ำสุด 0.48°
print(io.describe(radar)["elevations_deg"], "| ตัวแปร:", sorted(sw.fields))
"""),
md("### 2.1 แผนที่ PPI ของ 6 ตัวแปร\nแต่ละแผงใช้ช่วงสีคงที่ตามความหมายทางกายภาพ (ไม่ใช้ auto-scale) เพื่อให้เปรียบเทียบข้ามเวลาและเรดาร์ได้"),
code("""
from matplotlib.colors import Normalize, BoundaryNorm
Z, D, RH, PHI, V, W = (io.field(sw, k) for k in ["reflectivity", "corrected_differential_reflectivity",
                       "cross_correlation_ratio", "differential_phase", "velocity", "spectrum_width"])
lonlat = mp.sweep_lonlat(sw)
echo = np.isfinite(Z) & (Z > 5)                      # แสดงเฉพาะ gate ที่มีสัญญาณ
lon0, lat0 = float(sw.longitude["data"][0]), float(sw.latitude["data"][0])
E = [lon0 - 1.3, lon0 + 1.3, lat0 - 1.2, lat0 + 1.2]
panels = [("Z (dBZ)", Z, *pl.dbz_cmap()),
          ("ZDR (dB)", D, plt.get_cmap("Spectral_r", 12), Normalize(-4, 4)),
          ("ρhv", RH, plt.get_cmap("viridis", 10), Normalize(0.7, 1.0)),
          ("ΦDP (°)", PHI, plt.get_cmap("twilight", 12), Normalize(-180, 180)),
          ("V (m s⁻¹)", V, plt.get_cmap("RdBu_r", 16), Normalize(-16, 16)),
          ("W (m s⁻¹)", W, plt.get_cmap("magma_r", 10), Normalize(0, 8))]
fig, axs = pl.subplots(2, 3, width="double", height_mm=135)
for ax, (t, a, cm, nm) in zip(axs.ravel(), panels):
    im = mp.ppi(ax, sw, np.where(echo, a, np.nan), cm, nm, lonlat=lonlat)
    mp.basemap(ax, E, radars=["PHS"], rings_km=[50, 100], grid_step=1.0, scale_km=50, north=False, lw=0.7)
    fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.85, pad=0.08).set_label(t)
pl.panel_labels(axs)
fig.suptitle("เรดาร์พิษณุโลก (PHS) มุมยก 0.48°  3 ต.ค. 2567 00:00 UTC", fontweight="bold")
pl.save(fig, "01_phs_ppi_six_variables"); plt.show()
"""),
md("""
**สังเกต (ชวนอภิปราย):**
- บริเวณใกล้เรดาร์ (< 30 กม.) มี Z สูงแต่ **ρhv ต่ำ, V ≈ 0, W ต่ำ** → ลักษณะของ ground clutter / anomalous propagation (AP)
- ในฝนจริง ρhv > 0.97 และ ΦDP เพิ่มขึ้นตามระยะอย่างเรียบ
- ΦDP ไม่ได้เริ่มที่ 0° — ค่าเริ่มต้นนี้คือ **system differential phase** ของฮาร์ดแวร์ ต้องลบออกก่อนใช้
"""),
md("### 2.2 การกระจายของ ρhv และ ZDR: ฝนกับไม่ใช่ฝนต่างกันอย่างไร"),
code("""
near = (sw.range["data"][None, :] < 30e3) & echo
far = (sw.range["data"][None, :] > 50e3) & echo
fig, axs = pl.subplots(1, 2, width="double", height_mm=60)
bins = np.linspace(0.5, 1.02, 53)
axs[0].hist(RH[near], bins, alpha=0.6, density=True, label="< 30 กม. (มี clutter)", color="#d95f02")
axs[0].hist(RH[far], bins, alpha=0.6, density=True, label="> 50 กม. (ฝนส่วนใหญ่)", color="#1b9e77")
axs[0].axvline(0.85, color="k", ls="--", lw=0.8); axs[0].set_xlabel("ρhv"); axs[0].set_ylabel("ความหนาแน่น")
axs[0].legend()
b2 = np.linspace(-6, 6, 61)
axs[1].hist(D[far & (RH > 0.97)], b2, alpha=0.7, density=True, color="#1b9e77", label="ρhv > 0.97")
axs[1].set_xlabel("ZDR (dB)"); axs[1].set_ylabel("ความหนาแน่น"); axs[1].legend()
pl.panel_labels(axs); pl.save(fig, "01_rhohv_zdr_histograms"); plt.show()
print(f"มัธยฐาน ZDR ในฝน = {np.nanmedian(D[far & (RH > 0.97)]):.2f} dB  ← ฝนเบาควรมี ZDR ≈ 0.2–0.4 dB")
"""),
md("""
มัธยฐาน ZDR ของ PHS ติดลบประมาณ −2 dB ซึ่งเป็นไปไม่ได้ทางกายภาพสำหรับฝน (เม็ดฝนแบนเสมอ ZDR ≥ 0) แปลว่า
**ช่อง H และ V ของเรดาร์มีความต่างของการขยายสัญญาณ (ZDR bias)** ต้องปรับเทียบก่อนใช้

## 3. การปรับเทียบ (calibration)
$$\\text{ZDR}_{\\text{offset}} = \\operatorname{median}\\left(\\text{ZDR}\\;\\middle|\\; 20 \\le Z \\le 28\\ \\text{dBZ},\\ \\rho_{hv} > 0.97\\right) - 0.25\\ \\text{dB}$$

ฝนเบา (20–28 dBZ) มีเม็ดเล็กเกือบกลม ZDR ตามธรรมชาติประมาณ 0.2–0.3 dB (Gourley et al. 2009)
**system phase** ประมาณจากค่า ΦDP เฉลี่ยแบบวงกลมของ 5 gate แรกที่เป็นฝนในแต่ละ ray และ **ทิศการสะสมของ ΦDP**
(บวก/ลบ ขึ้นกับผู้ผลิต) ดูจากความชันของ ΦDP ตามระยะ
"""),
code("""
radars = {r: [io.read_radar(p) for p in g.path] for r, g in scans.groupby("radar")}
cal = qc.calibrate(radars)
cal.round(2)
"""),
code("""
fig, ax = pl.subplots(width="double", height_mm=70)
for name, lst in radars.items():
    vals = []
    for r in lst:
        s = r.extract_sweeps([0])
        z, d, rh = io.field(s, "reflectivity"), io.field(s, "corrected_differential_reflectivity"), io.field(s, "cross_correlation_ratio")
        vals.append(d[(z >= 20) & (z <= 28) & (rh > 0.97)])
    v = np.concatenate(vals)
    if v.size > 50:
        ax.hist(v, np.linspace(-6, 4, 81), histtype="step", lw=1.3, density=True, color=mp.RADAR_COLORS[name],
                label=f"{name} (มัธยฐาน {np.median(v):+.2f} dB)")
ax.axvline(0.25, color="k", ls="--", lw=0.8); ax.text(0.35, ax.get_ylim()[1] * 0.9, "ค่าคาดหวัง 0.25 dB", fontsize=7)
ax.set_xlabel("ZDR ในฝนเบา (dB)"); ax.set_ylabel("ความหนาแน่น"); ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1))
pl.save(fig, "01_zdr_calibration"); plt.show()
cal.to_csv(P.tables / "01_radar_calibration_00UTC.csv")
"""),
md("""
**แปลผล:** PHS และ TAK มี ZDR bias −2 ถึง −3 dB ขณะที่ CRI ใกล้ศูนย์ — เรดาร์แต่ละตัวต้องปรับเทียบแยกกัน
ถ้าไม่แก้ ZDR การจำแนกชนิดอนุภาคและสูตรฝนที่ใช้ ZDR จะผิดทั้งหมด (MHS มีจุดฝนเบาน้อยเกินไป ค่า offset จึงไม่น่าเชื่อถือ)

## 4. ขั้นตอน QC ทีละขั้น
"""),
code("""
c = cal.loc["PHS"]
s, out = qc.process_sweep(radar, 0, c.zdr_offset_db, c.system_phase_deg, -1.0)
r_ = sw.range["data"][None, :]
steps = {
    "1. มีสัญญาณ (Z ≥ 5 dBZ)": ~out["noecho"],
    "2. + ρhv ≥ 0.85": ~out["noecho"] & (out["RHOHV"] >= 0.85),
    "3. + ระยะ > 5 กม.": ~out["noecho"] & (out["RHOHV"] >= 0.85) & (r_ > 5e3),
    "4. + speckle + clutter (= met)": out["met"],
}
cnt = pd.Series({k: int(v.sum()) for k, v in steps.items()})
display(cnt.to_frame("จำนวน gate"))
fig, axs = pl.subplots(1, 2, width="double", height_mm=105)
cm, nm = pl.rain_cmap()
raw = qc.qpe.r_from_z(Z)
for ax, a, t in [(axs[0], np.where(np.isfinite(Z) & (Z > 5), raw, np.nan), "ก่อน QC: R(Z) จากทุก gate"),
                 (axs[1], np.where(out["met"], out["R_MP"], np.nan), "หลัง QC: เฉพาะสัญญาณฝน")]:
    im = mp.ppi(ax, sw, a, cm, nm, lonlat=lonlat)
    mp.basemap(ax, E, radars=["PHS"], rings_km=[50, 100], scale_km=50, north=False, title=t, lw=0.7)
fig.colorbar(im, ax=axs, orientation="horizontal", shrink=0.6, extend="max").set_label(pl.L("rain_mmh"))
pl.save(fig, "01_qc_before_after"); plt.show()
"""),
md("""
## 5. ΦDP → KDP และการแก้การลดทอน (ZPHI)
ที่ย่าน C (λ ≈ 5.3 ซม.) ฝนหนักทำให้สัญญาณอ่อนลงระหว่างทาง (attenuation) ค่า Z ที่อยู่หลังเซลล์ฝนจึงต่ำกว่าความจริง
การลดทอนรวมตามทาง (path-integrated attenuation, PIA) แปรผันตาม ΦDP:

$$\\text{PIA}(r) \\approx \\alpha\\,[\\Phi_{DP}(r) - \\Phi_{DP}(r_0)], \\qquad \\alpha \\approx 0.06\\text{–}0.10\\ \\text{dB}/^{\\circ}\\ (\\text{C band})$$

วิธี **ZPHI** (Testud et al. 2000) ใช้ ΔΦDP ของทั้ง ray เป็นข้อจำกัด แล้วกระจาย PIA ตามรูปร่างของ Z เพื่อหา
ค่าการลดทอนจำเพาะ $A(r)$ (dB km⁻¹) และ $Z_c = Z + 2\\int_0^r A\\,ds$ ซึ่ง $A$ เองก็ใช้ประมาณฝนได้ (R(A), notebook 03)

KDP คำนวณด้วยวิธี Vulpiani et al. (2012) ใน Py-ART: กรองสัญญาณรบกวนแล้วหาอนุพันธ์ของ ΦDP ในหน้าต่าง ~6 กม.
"""),
code("""
span = out["phidp_span"]
iray = int(np.argmax(span))
rk = sw.range["data"] / 1e3
fig, axs = pl.subplots(3, 1, width="one_half", height_mm=120, sharex=True)
axs[0].plot(rk, out["PHIDP"][iray], ".", ms=2, color="0.6", label="ΦDP หลังลบ system phase")
axs[0].plot(rk, out["PHIDP_F"][iray], "-", color="#1b9e77", lw=1.2, label="ΦDP กรองแล้ว (Vulpiani)")
axs[0].set_ylabel("ΦDP (°)"); axs[0].legend(loc="upper left")
axs[1].plot(rk, out["KDP"][iray], color="#7570b3"); axs[1].set_ylabel("KDP (° km⁻¹)")
axs[2].plot(rk, out["Z"][iray], color="0.5", lw=1, label="Z วัดได้")
axs[2].plot(rk, out["Zc"][iray], color="#d95f02", lw=1.2, label="Z แก้การลดทอน (ZPHI)")
axs[2].set_ylabel("Z (dBZ)"); axs[2].set_xlabel(pl.L("range_km")); axs[2].legend(loc="upper right")
az = sw.azimuth["data"][iray]
axs[0].set_title(f"ray ที่ azimuth {az:.0f}° (ΔΦDP สูงสุด {span[iray]:.0f}°)", loc="left", fontweight="bold")
for a in axs: a.grid(alpha=0.3, lw=0.4)
pl.panel_labels(axs); pl.save(fig, "01_phidp_kdp_zphi_ray"); plt.show()
print(f"สัดส่วน gate ฝนที่ใช้ R(KDP) = {out['frac_kdp']:.3f}, ใช้ R(A) = {out['frac_a']:.3f}")
"""),
md("""
**แปลผล:** ในชั่วโมงนี้ฝนส่วนใหญ่เป็นฝนเบา ΔΦDP จึงน้อย (ไม่กี่องศา) การแก้ attenuation ทำให้ Z เพิ่มเพียงเล็กน้อย
และ KDP มีสัญญาณน้อย — dual-pol จะแสดงพลังชัดในฝนหนัก (notebook 03 ใช้เวลา 10 UTC)

---
## ✏️ แบบฝึกหัด
"""),
exercise(1, "ug", "นับว่า **แต่ละขั้น QC ตัด gate ออกไปกี่ gate** (ใช้ dict `steps`) และคิดเป็นร้อยละของ gate ที่มีสัญญาณ ขั้นใดตัดมากที่สุด"),
todo("""
removed = None   # pd.Series: step -> number of gates removed by that step
"""),
solution(1, """
v = cnt.values
removed = pd.Series(np.r_[0, -np.diff(v)], index=cnt.index)
display(pd.DataFrame({"เหลือ": cnt, "ตัดออก": removed, "% ของ gate ที่มีสัญญาณ": (100 * removed / cnt.iloc[0]).round(1)}))
print("ขั้นที่ตัดมากที่สุด:", removed.idxmax())
"""),
exercise(2, "ug", "คำนวณ **ZDR offset ของ TAK** ด้วยตนเองจาก sweep ต่ำสุดของทุกไฟล์ (มัธยฐานของ ZDR ที่ 20 ≤ Z ≤ 28 dBZ และ ρhv > 0.97 ลบ 0.25 dB) แล้วเทียบกับตาราง `cal`"),
todo("""
zdr_off_tak = None
"""),
solution(2, """
vals = []
for r in radars["TAK"]:
    s0 = r.extract_sweeps([0])
    z, d, rh = (io.field(s0, k) for k in ["reflectivity", "corrected_differential_reflectivity", "cross_correlation_ratio"])
    vals.append(d[(z >= 20) & (z <= 28) & (rh > 0.97)])
zdr_off_tak = np.median(np.concatenate(vals)) - 0.25
print(f"ZDR offset TAK = {zdr_off_tak:.2f} dB (ตาราง: {cal.loc['TAK', 'zdr_offset_db']:.2f} dB)")
assert abs(zdr_off_tak - cal.loc["TAK", "zdr_offset_db"]) < 1e-6
"""),
exercise(3, "ug", "เปลี่ยนเกณฑ์ `min_rhohv` จาก 0.85 เป็น **0.95** (`qc.process_sweep(..., qc={'min_rhohv': 0.95})`) พื้นที่ฝน (จำนวน gate ที่ R ≥ 0.5 มม./ชม.) ลดลงกี่เปอร์เซ็นต์ ข้อดีข้อเสียของเกณฑ์เข้มคืออะไร"),
todo("""
change_pct = None
"""),
solution(3, """
_, o95 = qc.process_sweep(radar, 0, c.zdr_offset_db, c.system_phase_deg, -1.0, qc={"min_rhohv": 0.95})
n85 = int(np.nansum(out["R_MP"] >= 0.5)); n95 = int(np.nansum(o95["R_MP"] >= 0.5))
change_pct = 100 * (n95 - n85) / n85
print(f"gate ฝน: 0.85 → {n85}, 0.95 → {n95}  ({change_pct:+.1f}%)")
print("เกณฑ์เข้มตัด clutter ได้มากขึ้น แต่ตัดฝนจริงที่ ρhv ต่ำ (ขอบเซลล์, ฝนปนลูกเห็บ, ระยะไกลที่ SNR ต่ำ) ด้วย")
"""),
exercise(4, "pg", "ประมาณ **KDP ด้วยการถดถอยเชิงเส้น** ของ ΦDP (กรองแล้ว) ในหน้าต่างเลื่อน 6 กม. บน ray `iray` แล้วเปรียบเทียบกับ KDP ของ Vulpiani ทำไม KDP = ½ × ความชัน"),
todo("""
kdp_ls = None   # array along the ray
"""),
solution(4, """
phi = out["PHIDP_F"][iray]; gate = float(np.median(np.diff(sw.range["data"]))) / 1e3
n = int(round(6 / gate)) | 1; h = n // 2
kdp_ls = np.full(phi.shape, np.nan)
for i in range(h, len(phi) - h):
    y = phi[i - h:i + h + 1]; x = rk[i - h:i + h + 1]; k = np.isfinite(y)
    if k.sum() >= n * 0.7:
        kdp_ls[i] = 0.5 * np.polyfit(x[k], y[k], 1)[0]
fig, ax = pl.subplots(width="one_half", height_mm=50)
ax.plot(rk, out["KDP"][iray], label="Vulpiani (Py-ART)"); ax.plot(rk, kdp_ls, label="ถดถอยเชิงเส้น 6 กม.")
ax.set_xlabel(pl.L("range_km")); ax.set_ylabel("KDP (° km⁻¹)"); ax.legend(); plt.show()
print("ΦDP สะสมทั้งขาไปและขากลับ (two-way) จึงต้องหาร 2 เพื่อให้ได้ค่าต่อทางเดียว")
"""),
md("""
---
## 💡 คำถามชวนคิดและต่อยอด
- ทำไม ρhv จึงลดลงในชั้นหลอมเหลว (melting layer) — เราจะใช้คุณสมบัตินี้หาความสูงชั้นหลอมเหลวใน notebook 05
- ถ้า ZDR bias ของเรดาร์เปลี่ยนตามเวลา (อุณหภูมิอุปกรณ์) ควรปรับเทียบรายชั่วโมง รายวัน หรือรายเหตุการณ์?
- ทิศการสะสมของ ΦDP ที่ติดลบใน CRI/LMP/PHS มาจากอะไร? ควรยืนยันกับผู้ผลิตหรือกรมอุตุฯ อย่างไร
- 🔵 ลองออกแบบตัวจำแนกชนิดอนุภาคอย่างง่าย (ฝน/clutter/แมลง) จาก Z, ZDR, ρhv, texture ของ ΦDP แบบ fuzzy logic

## 📚 อ้างอิง
- Bringi, V.N. & Chandrasekar, V. (2001). *Polarimetric Doppler Weather Radar*. Cambridge University Press.
- Ryzhkov, A.V. & Zrnić, D.S. (2019). *Radar Polarimetry for Weather Observations*. Springer.
- Testud, J. et al. (2000). The rain profiling algorithm applied to polarimetric weather radar. *J. Atmos. Oceanic Technol.* 17:332–356.
- Vulpiani, G. et al. (2012). On the use of dual-polarized C-band radar for operational rainfall retrieval in mountainous areas. *J. Appl. Meteor. Climatol.* 51:405–425.
- Gourley, J.J. et al. (2009). Absolute calibration of radar reflectivity using redundancy of the polarization observations. *J. Atmos. Oceanic Technol.* 26:689–703.
"""),
save_cell("01"),
]
build(cells, f"notebooks/{F}")
