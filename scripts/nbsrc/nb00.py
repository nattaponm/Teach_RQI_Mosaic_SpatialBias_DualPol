from _nb import *

F = "00_setup_and_data.ipynb"
cells = [
header(F, "00", "เตรียมสภาพแวดล้อมและข้อมูล", "Setting up Colab, the repository and the teaching data",
       "🟢 ป.ตรี ปี 3–4 และ 🔵 ป.โท", "45–60 นาที", "–",
       "ข้อมูลเรดาร์ UF 5 สถานี, DEM, ขอบลุ่มน้ำ, ฝนสถานี สสน.", "ตารางรายการไฟล์, แผนที่พื้นที่ศึกษา, token สำหรับบันทึกผล"),
md("""
## 🎯 วัตถุประสงค์การเรียนรู้
เมื่อจบ notebook นี้ นิสิตจะสามารถ
1. **อธิบาย** โครงสร้างของ repository รายวิชา (data, notebooks, config, src, mapping, outputs) และหน้าที่ของแต่ละส่วน
2. **ดาวน์โหลดและตรวจสอบ** ข้อมูลชุดสอนด้วย checksum (SHA-256)
3. **อธิบาย** โครงสร้างข้อมูลเรดาร์แบบ volume → sweep → ray → gate และอ่านไฟล์ UF ด้วย Py-ART
4. **สร้าง** ตารางรายการข้อมูล (inventory) และแผนที่พื้นที่ศึกษาตามมาตรฐานวารสาร
5. **ตั้งค่า** fork และ token เพื่อบันทึกผลลัพธ์ขึ้น GitHub ของตนเอง

## 🌧️ กรณีศึกษา: น้ำท่วมเชียงใหม่ ตุลาคม 2567
ปลายเดือนกันยายน–ต้นตุลาคม 2567 ร่องความกดอากาศต่ำพาดผ่านภาคเหนือตอนบน ทำให้ฝนตกหนักต่อเนื่องบนพื้นที่ต้นน้ำปิง
(แม่แตง เชียงดาว แม่งัด) ระดับน้ำที่สถานี P.1 สะพานนวรัฐ สูงเป็นประวัติการณ์ราว 5.3 ม. ในคืนวันที่ 4–5 ต.ค. 2567
รายวิชานี้ใช้ข้อมูล **วันที่ 3 ต.ค. 2567 สองช่วงเวลา** คือ 00–01 UTC (07–08 น. เวลาไทย, ฝนแบบ stratiform) และ
10–11 UTC (17–18 น., ฝนแบบพาความร้อน) จากเรดาร์ dual-polarization ย่าน C ของกรมอุตุนิยมวิทยา 5 สถานี
"""),
md("""
## 1. หลักการ: ข้อมูลเรดาร์ตรวจอากาศมีหน้าตาอย่างไร

เรดาร์หมุนรอบตัวเองที่มุมยก (elevation) หลายค่า ข้อมูลหนึ่งรอบการสแกนเรียกว่า **volume** (กรมอุตุฯ ทุก 15 นาที)

| ระดับ | ความหมาย | ตัวอย่างในข้อมูลนี้ |
|---|---|---|
| volume | การสแกนครบทุกมุมยก 1 ครั้ง | 1 ไฟล์ UF |
| sweep | การหมุน 360° ที่มุมยกเดียว | LMP มี 4 sweep: 1.59°, 2.39°, 3.39°, 4.30° |
| ray | ลำคลื่น 1 ทิศ (azimuth) | ประมาณ 360 ray ต่อ sweep |
| gate (bin) | ช่วงระยะตามลำคลื่น | ทุก 0.5–1 กม. ไกลสุด 240 กม. |

ตำแหน่งของ gate ในอวกาศขึ้นกับระยะ $r$ มุมยก $\\theta$ และความโค้งของโลก (เรียนละเอียดใน notebook 02)

$$h(r,\\theta) = \\sqrt{r^2 + (k_e a)^2 + 2 r k_e a \\sin\\theta} - k_e a + h_0, \\qquad k_e = 4/3,\\; a = 6371\\ \\text{km}$$

**รูปแบบไฟล์ UF (Universal Format)** เป็นรูปแบบเก่าที่ยังใช้แลกเปลี่ยนข้อมูลเรดาร์ ภายในเก็บตัวแปรหลายตัวต่อ gate
(Z, ZDR, ρhv, ΦDP, V, W) Py-ART อ่านแล้วแปลงเป็นวัตถุ `Radar` ที่มี `fields`, `range`, `azimuth`, `fixed_angle`

> ⏰ **เวลา**: ชื่อไฟล์เรดาร์เป็นเวลา **UTC** ส่วนฝนสถานี สสน. เป็นเวลา **ไทย (UTC+7)** — ต้องแปลงให้ตรงกันเสมอ (notebook 06)
"""),
md("## 2. เตรียมสภาพแวดล้อม\nเซลล์ต่อไปนี้จะ clone repository (ถ้ารันใน Colab) ติดตั้ง Py-ART และดาวน์โหลดข้อมูลชุดสอน (~180 MB, ครั้งแรกใช้เวลา 1–3 นาที)"),
setup_cell(("teaching_data", "products")),
md("""
### 2.1 โครงสร้าง repository
ทุกพารามิเตอร์อยู่ใน `config/*.yaml` เพื่อให้ทดลองเปลี่ยนค่าได้โดยไม่ต้องแก้โค้ด ส่วนฟังก์ชันคำนวณอยู่ใน `src/radteach/`
"""),
code("""
from pathlib import Path
def tree(p, depth=2, prefix=""):
    p = Path(p)
    items = sorted([x for x in p.iterdir() if not x.name.startswith(".") and x.name not in ("__pycache__",)])
    for i, x in enumerate(items):
        last = i == len(items) - 1
        print(prefix + ("└── " if last else "├── ") + x.name + ("/" if x.is_dir() else ""))
        if x.is_dir() and depth > 1 and x.name not in ("teaching_data_v1", "teaching_products_v1", "fonts"):
            tree(x, depth - 1, prefix + ("    " if last else "│   "))
tree(P.root, depth=2)
"""),
code("""
from radteach import config
cfg = config.load_all()
print("ไฟล์ config:", list(cfg))
print("กริดวิเคราะห์:", cfg["study_area"]["grid"])
print("สัมประสิทธิ์ Z–R (Marshall–Palmer):", cfg["qpe"]["estimators"]["zr_marshall_palmer"])
"""),
md("""
### 2.2 ตรวจสอบความถูกต้องของไฟล์ด้วย checksum
ข้อมูลวิทยาศาสตร์ต้อง **ทำซ้ำได้ (reproducible)** ค่า SHA-256 เป็น "ลายนิ้วมือ" ของไฟล์ ถ้าไฟล์เปลี่ยนแม้ 1 ไบต์ ค่าจะเปลี่ยนทั้งหมด
`colab.setup()` ตรวจค่านี้กับ `config/data_sources.yaml` ให้อัตโนมัติแล้ว ลองคำนวณเองกับไฟล์หนึ่งไฟล์
"""),
code("""
from radteach import io
scans = io.list_scans(P.radar)
print(f"พบไฟล์เรดาร์ {len(scans)} ไฟล์, รวม {scans.size_mb.sum():.0f} MB")
f0 = scans.path.iloc[0]
print(Path(f0).name, "SHA-256 =", colab.sha256(f0)[:16], "…")
scans.head(10)
"""),
md("## 3. อ่านไฟล์เรดาร์ 1 ไฟล์ด้วย Py-ART"),
code("""
f = [p for p in scans.path if "LMP_202410030000" in p][0]
radar = io.read_radar(f)
info = io.describe(radar)
for k, v in info.items():
    print(f"{k:16s}: {v}")
"""),
md("""
**สังเกต:** LMP เริ่มสแกนที่ 1.59° (สูงกว่าเรดาร์อื่น) และมี gate ทุก 1 กม. ส่วน PHS เริ่มที่ 0.48° — มุมยกต่ำสุดมีผลมาก
ต่อความสูงของลำคลื่นและการถูกภูเขาบัง ซึ่งเป็นหัวใจของ RQI ใน notebook 02
"""),
md("## 4. ตารางรายการข้อมูล (inventory) ทั้ง 2 ชั่วโมง"),
code("""
inv = io.inventory(P.radar)          # เปิดทุกไฟล์ ~0.3 วินาที/ไฟล์
inv.to_csv(P.tables / "00_radar_inventory.csv", index=False)
summary = (inv.groupby("radar")
             .agg(n_files=("file", "size"), n_sweeps=("n_sweeps", "first"), elevations=("elevations_deg", "first"),
                  gate_m=("gate_m", "first"), max_range_km=("max_range_km", "first"), size_mb=("size_mb", "sum")))
summary
"""),
md("### 4.1 กราฟความพร้อมของข้อมูล (data availability)\nกราฟมาตรฐาน: กว้าง 140 มม., มีหน่วย, ระบุเวลาเป็น UTC"),
code("""
import numpy as np, matplotlib.pyplot as plt
from radteach import plotting as pl
avail = inv.assign(t=inv.time.dt.strftime("%H:%M")).pivot_table(index="radar", columns="t", values="n_sweeps", aggfunc="first")
fig, ax = pl.subplots(width="one_half", height_mm=55)
im = ax.imshow(avail.values, cmap="Blues", vmin=0, vmax=6, aspect="auto")
ax.set_xticks(range(avail.shape[1])); ax.set_xticklabels(avail.columns, rotation=90)
ax.set_yticks(range(avail.shape[0])); ax.set_yticklabels(avail.index)
for (i, j), v in np.ndenumerate(avail.values):
    ax.text(j, i, "–" if np.isnan(v) else int(v), ha="center", va="center", fontsize=7)
ax.set_xlabel("เวลา (UTC) วันที่ 3 ต.ค. 2567"); ax.set_title("จำนวน sweep ต่อไฟล์", loc="left", fontweight="bold")
cb = fig.colorbar(im, ax=ax, shrink=0.9); cb.ax.set_title("sweep", fontsize=8)
pl.save(fig, "00_data_availability")
plt.show()
"""),
md("""
## 5. แผนที่พื้นที่ศึกษา (มาตรฐานวารสาร)
องค์ประกอบที่แผนที่วิทยาศาสตร์ต้องมี: **ระบบพิกัด/กริดพิกัด, มาตราส่วน, ทิศเหนือ, คำอธิบายสัญลักษณ์, แผนที่บอกตำแหน่ง (inset), แหล่งข้อมูล**
ฟังก์ชัน `radteach.mapping.basemap()` วาดองค์ประกอบเหล่านี้ให้เหมือนกันทุกแผนที่ในรายวิชา
"""),
code("""
from radteach import dem, mapping as mp
from matplotlib.lines import Line2D
z, tr, _ = dem.read(P.products / "dem_study_area.tif")
E = cfg["study_area"]["map_extent"]["north"]
fig = pl.figure("double", 150)
ax = fig.add_axes([0.06, 0.08, 0.62, 0.86])
mp.hillshade(ax, z, dem.extent(z, tr), alpha=0.45)
ax.imshow(z, extent=dem.extent(z, tr), cmap="terrain", vmin=0, vmax=2500, alpha=0.35, zorder=0.5)
mp.basemap(ax, E, rings_km=240, subbasins=cfg["study_area"]["ping_headwater_subbasins"], scale_km=100)
st = __import__("pandas").read_csv(P.stations / "hii_station_metadata_202410.csv")
st = st[(st.longitude.between(E[0], E[1])) & (st.latitude.between(E[2], E[3]))]
ax.scatter(st.longitude, st.latitude, s=6, c="#1f3a93", lw=0, zorder=12)
rid = __import__("pandas").read_csv(P.stations / "rid_stations.csv")
ax.plot(rid.lon, rid.lat, "s", mfc="yellow", mec="k", ms=5, zorder=13)
mb, sb, _ = mp.boundaries()
for code_, row in mb[mb.MB_CODE.isin(cfg["study_area"]["basins_of_interest"])].iterrows():
    p = row.geometry.representative_point()
    ax.text(p.x, p.y + 0.12, row.MBASIN_T, ha="center", fontsize=9, color="0.15", zorder=12,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))
mp.source_note(ax, "DEM: SRTM 250 m · ลุ่มน้ำ: สทนช. · สถานี: สสน. · เรดาร์: กรมอุตุนิยมวิทยา")
mp.inset_thailand(fig, [0.70, 0.62, 0.17, 0.32], E)
handles = [Line2D([], [], marker="^", ls="", color="#c0392b", mec="k", label="เรดาร์ (กรมอุตุฯ)"),
           Line2D([], [], color="k", ls="--", lw=0.7, label="รัศมี 240 กม."),
           Line2D([], [], color="k", lw=0.9, label="ลุ่มน้ำหลัก 6 ลุ่ม"),
           Line2D([], [], color="k", lw=0.7, ls="--", label="ต้นน้ำปิงเหนือเชียงใหม่"),
           Line2D([], [], marker="o", ls="", mfc="#1f3a93", mec="none", ms=3, label=f"สถานีฝน สสน. ({len(st)})"),
           Line2D([], [], marker="s", ls="", mfc="yellow", mec="k", ms=5, label="สถานีวัดน้ำ กรมชลประทาน")]
fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.69, 0.10), frameon=False)
pl.save(fig, "00_study_area_map")
plt.show()
"""),
md("""
**อ่านแผนที่:** เรดาร์ 5 สถานีครอบคลุมภาคเหนือเกือบทั้งหมดภายในรัศมี 240 กม. แต่พื้นที่ต้นน้ำปิง (เส้นประ) เป็นเทือกเขาสูง
1,500–2,500 ม. — ลำคลื่นจาก LMP และ CRI ต้องผ่านภูเขา คำถามหลักของรายวิชาคือ **"เรดาร์เห็นฝนบนภูเขาได้ดีเพียงใด และจะรวมเรดาร์หลายตัวอย่างไรให้ดีที่สุด"**
"""),
md("""
## 6. ตั้งค่าการบันทึกผลขึ้น GitHub ของตนเอง (ทำครั้งเดียว)

repository ของอาจารย์เป็น **อ่านอย่างเดียว** สำหรับนิสิต นิสิตเก็บผลงานใน **fork** ของตนเอง:

1. ไปที่ https://github.com/nattaponm/Teach_RQI_Mosaic_SpatialBias_DualPol แล้วกด **Fork** (มุมขวาบน)
2. สร้าง token: GitHub ▸ Settings ▸ Developer settings ▸ *Fine-grained tokens* ▸ Generate new token
   - Repository access: **Only select repositories** → เลือก fork ของตนเองเท่านั้น
   - Permissions ▸ Repository permissions ▸ **Contents: Read and write**
   - อายุ token: จบภาคการศึกษา
3. ใน Colab คลิกไอคอน 🔑 (Secrets) แถบซ้าย เพิ่ม 2 รายการ และเปิด *Notebook access*
   - `GITHUB_USER` = ชื่อบัญชี GitHub ของคุณ
   - `GITHUB_TOKEN` = token จากข้อ 2
4. ท้ายทุก notebook มีเซลล์ 💾 ใส่รหัสนิสิตแล้วรัน ผลใน `outputs/` จะถูก push ไปที่ branch `results/<รหัสนิสิต>` ของ fork

> 🔒 **ห้าม** พิมพ์ token ลงใน notebook โดยตรง และห้ามแชร์ notebook ที่มี token ฟังก์ชัน `save_to_fork` จะไม่แสดง token ในข้อความผิดพลาด
> และปฏิเสธการ push ไปยังบัญชีของอาจารย์
>
> ต้องการบันทึก *ตัว notebook* ด้วย: เมนู **File ▸ Save a copy in GitHub** แล้วเลือก fork ของตนเอง
"""),
md("---\n## ✏️ แบบฝึกหัด"),
exercise(1, "ug", "ใช้ตาราง `inv` หา **จำนวนไฟล์และขนาดรวม (MB) ของแต่ละเรดาร์** และบอกว่าเรดาร์ใดมีไฟล์ใหญ่ที่สุด เพราะอะไร (คำใบ้: จำนวน sweep และขนาด gate)"),
todo("""
ans1 = None   # DataFrame: index = radar, columns = n_files, total_mb
"""),
solution(1, """
ans1 = inv.groupby("radar").agg(n_files=("file", "size"), total_mb=("size_mb", "sum"),
                                n_sweeps=("n_sweeps", "first"), gate_m=("gate_m", "first"))
display(ans1.sort_values("total_mb", ascending=False))
print("MHS ใหญ่ที่สุด: มี 5 sweep และ gate 500 ม. (จำนวน gate ต่อ ray มากกว่า LMP ที่ใช้ 1 กม.)")
"""),
exercise(2, "ug", "เรดาร์ใด **มีมุมยกต่ำสุดต่ำที่สุด** และใด **สูงที่สุด**? ข้อนี้บอกอะไรเกี่ยวกับความสามารถในการมองฝนใกล้พื้นที่ระยะไกล"),
todo("""
lowest = None   # dict radar -> lowest elevation (deg)
"""),
solution(2, """
lowest = {r: float(e.split(";")[0]) for r, e in inv.groupby("radar").elevations_deg.first().items()}
print(dict(sorted(lowest.items(), key=lambda kv: kv[1])))
print("ต่ำสุด PHS 0.48°, สูงสุด MHS 3.0° → MHS มองฝนใกล้พื้นได้แค่ระยะใกล้ (ลำคลื่นสูงเร็ว)")
assert min(lowest, key=lowest.get) == "PHS"
"""),
exercise(3, "pg", "ใช้สมการความสูงลำคลื่น คำนวณว่าที่ระยะ 150 กม. ลำคลื่นของ PHS (0.48°, ฐานสูง 72 ม.) และ MHS (3.0°, 280 ม.) สูงจากระดับน้ำทะเลเท่าใด เทียบกับชั้นหลอมเหลว (~4.8 กม.)"),
todo("""
h_phs = None; h_mhs = None   # km
"""),
solution(3, """
from radteach import geometry as g
h_phs = g.beam_height(150e3, 0.48, 72) / 1e3
h_mhs = g.beam_height(150e3, 3.0, 280) / 1e3
print(f"PHS {h_phs:.2f} km, MHS {h_mhs:.2f} km → MHS อยู่เหนือชั้นหลอมเหลว จึงวัดหิมะ/ชั้นสว่าง ไม่ใช่ฝน")
assert 2.5 < h_phs < 3.2 and 9.0 < h_mhs < 10.5
"""),
md("""
---
## 💡 คำถามชวนคิดและต่อยอด
- ถ้าข้อมูลเรดาร์มีทุก 6 นาที (อย่าง NEXRAD) แทน 15 นาที จะมีผลต่อการสะสมฝนรายชั่วโมงอย่างไร (ดู notebook 05)
- ทำไมงานวิจัยที่ดีต้องเผยแพร่ทั้ง **ข้อมูล โค้ด และ checksum**? สิ่งนี้เกี่ยวกับหลัก FAIR data อย่างไร
- ลองหาว่าประเทศไทยมีเรดาร์ของหน่วยงานใดอีกบ้าง (กรมฝนหลวงฯ, กทม.) และถ้านำมารวม จะต้องปรับมาตรฐานข้อมูลอะไรบ้าง

## 📚 อ้างอิง
- Helmus, J.J. & Collis, S.M. (2016). The Python ARM Radar Toolkit (Py-ART). *J. Open Research Software* 4(1):e25. doi:10.5334/jors.119
- Doviak, R.J. & Zrnić, D.S. (1993). *Doppler Radar and Weather Observations*, 2nd ed. Academic Press.
- Wilkinson, M.D. et al. (2016). The FAIR Guiding Principles. *Scientific Data* 3:160018. doi:10.1038/sdata.2016.18
"""),
save_cell("00"),
]
build(cells, f"notebooks/{F}")
