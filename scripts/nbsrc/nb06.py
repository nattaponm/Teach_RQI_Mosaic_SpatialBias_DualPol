from _nb import *

F = "06_rain_gauges_hii.ipynb"
cells = [
header(F, "06", "ฝนสถานีโทรมาตร สสน.: ดึงข้อมูล ตรวจคุณภาพ และจัดเวลา", "HII hourly rain gauges: download, quality control and time alignment",
       "🟢 ป.ตรี ปี 3–4 และ 🔵 ป.โท", "90 นาที", "00", "Open data สสน. (ออนไลน์) หรือไฟล์ snapshot 221 สถานี วันที่ 3 ต.ค. 2567",
       "ตารางฝนรายชั่วโมงที่ผ่าน QC (เวลา UTC), รายงานคุณภาพรายสถานี, แผนที่ฝนสถานี, กราฟความพร้อมของข้อมูล"),
md("""
## 🎯 วัตถุประสงค์การเรียนรู้
1. **ดึง** ข้อมูลฝนรายชั่วโมงจากคลังข้อมูลเปิดของ สสน. (HII open data catalogue) ด้วย Python
2. **ใช้** `quality_flag` และการตรวจเพิ่มเติม (ช่วงค่า, ค่าค้าง, เทียบเพื่อนบ้าน) เพื่อคัดข้อมูลที่ใช้ได้
3. **แปลง** เวลาไทยเป็น UTC และจับคู่ "ชั่วโมงฝน" ของสถานีกับ "ชั่วโมงเรดาร์" ได้ถูกต้อง
4. **พิสูจน์** ข้อตกลงด้านเวลา (ต้นชั่วโมง/ปลายชั่วโมง) ด้วยข้อมูลจริง
5. **นำเสนอ** ความพร้อมของข้อมูลและฝนสถานีด้วยกราฟและแผนที่มาตรฐาน

## 1. หลักการ
### 1.1 สถานีวัดฝนอัตโนมัติ
สถานีโทรมาตรของ สสน. ส่วนใหญ่ใช้ **กระบอกวัดฝนแบบกระดกน้ำ (tipping bucket)** 1 ครั้งที่กระดก = 0.2 มม.
ข้อดี: วัดเวลาละเอียดและส่งข้อมูลทันที ข้อจำกัด: วัดต่ำกว่าจริงในฝนหนักมาก/ลมแรง และอุดตันหรือเสียได้ → **ต้อง QC เสมอ**

### 1.2 โครงสร้างคลังข้อมูลเปิด
```
https://tiservice.hii.or.th/opendata/data_catalog/hourly_rain/{YYYY}/{YYYYMM}/0station_metadata.csv
https://tiservice.hii.or.th/opendata/data_catalog/hourly_rain/{YYYY}/{YYYYMM}/{STATION}.csv
```
| คอลัมน์ | ความหมาย |
|---|---|
| `station_code` | รหัสสถานี เช่น `BNHI` |
| `measure_datetime` | เวลา **ไทย (UTC+7)** ของ **ปลายชั่วโมง** |
| `rainfall_1h` | ฝนสะสม 1 ชั่วโมง (มม.) หรือ `null` |
| `quality_flag` | `N` = ปกติ, `null` = ไม่มีข้อมูล |

### 1.3 เวลา: ต้องตรงกับเรดาร์
$$t_{UTC} = t_{LT} - 7\\ \\text{ชม.}, \\qquad \\text{ฝนที่ประทับเวลา } t_{LT} = \\text{ฝนช่วง } [t_{LT}-1\\text{ ชม.},\\ t_{LT})$$
เรดาร์ชั่วโมง **00–01 UTC** จึงคู่กับค่าที่ประทับเวลา **08:00 น.** (ไม่ใช่ 07:00 น.) — จะพิสูจน์ในหัวข้อ 6

### 1.4 การตรวจคุณภาพ (QC)
| รหัส | เงื่อนไข | เหตุผล |
|---|---|---|
| NODATA | ค่าว่าง/`null` | สถานีไม่ส่งข้อมูล |
| FLAG | `quality_flag` ไม่ใช่ `N` | สสน. ระบุว่าใช้ไม่ได้ |
| RANGE | < 0 หรือ > 150 มม./ชม. | เกินขีดจำกัดทางกายภาพ |
| STUCK | มีค่าเล็ก ๆ คงที่ (≤ 0.6 มม.) เกือบทุกชั่วโมง | กระบอกรั่ว/สั่น/ส่งค่าค้าง |
| ISOLATED | ≥ 10 มม. แต่เพื่อนบ้าน ≥ 3 สถานีในรัศมี 20 กม. เป็น 0 | รายงานไว้ ไม่ตัดทิ้ง (ฝนพาความร้อนอาจเฉพาะจุดจริง) |

> ⚠️ **ห้ามใช้เรดาร์ตัดสินว่าฝนสถานีผิด** เพราะสถานีคือ "ค่าอ้างอิง" ที่จะใช้ประเมินเรดาร์ (หลีกเลี่ยงการให้เหตุผลวนกลับ)
"""),
setup_cell(("products",)),
md("## 2. ข้อมูลสถานี (metadata)"),
code("""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from radteach import gauges as gg, verify, config, plotting as pl, mapping as mp
cfg = config.load_all()
meta = gg.read_metadata(P.stations / "hii_station_metadata_202410.csv")
print(f"สถานีทั้งประเทศ {len(meta)} สถานี:", meta.station_type_telem.value_counts().to_dict(), "(R = ฝน, W = ระดับน้ำ+ฝน)")
meta.head()
"""),
md("""
## 3. ดึงข้อมูลจากเซิร์ฟเวอร์ สสน. (ออนไลน์)
เซลล์นี้ดาวน์โหลดข้อมูลเดือน ต.ค. 2567 ของ 3 สถานีตัวอย่าง ถ้าเครือข่ายใช้ไม่ได้จะข้ามไปใช้ไฟล์ snapshot ที่เตรียมไว้
(นำไปใช้กับทั้งเหตุการณ์ได้ด้วย `gg.fetch_period("2024-09-20", "2024-10-07", bbox=(96.5, 102, 15.3, 21.2))`)
"""),
code("""
LIVE_FETCH = True   #@param {type:"boolean"}
live = None
if LIVE_FETCH:
    try:
        d = gg.fetch_month(2024, 10, codes=["MEWG", "TGLG", "MMMO"], outdir=str(P.gauges))
        live = gg.read_month_dir(d)
        print(live[live.time_lt.dt.day == 3].head(12))
    except Exception as e:
        print("ดึงข้อมูลออนไลน์ไม่สำเร็จ:", type(e).__name__, "→ ใช้ snapshot แทน")
"""),
md("## 4. อ่าน snapshot 221 สถานี (05:00–19:00 น. วันที่ 3 ต.ค. 2567) และทำ QC"),
code("""
raw = gg.read_station_csv(P.stations / "hii_hourly_rain_20241003_0519LT.csv")
print("แถวทั้งหมด:", len(raw), "| สถานี:", raw.station_code.nunique())
print("ค่า quality_flag ที่พบ:", raw.quality_flag.value_counts().to_dict())
gc = cfg["gauges"]
qcd = gg.quality_control(raw, meta=meta, accept_flags=tuple(gc["accept_flags"]), max_rate=gc["max_rate_mm_h"],
                         iso_mm=gc["buddy"]["min_mm"], iso_km=gc["buddy"]["radius_km"])
qcd = gg.add_utc(qcd)
print(qcd.qc.value_counts().to_string())
rep = gg.qc_report(qcd); rep.to_csv(P.tables / "06_gauge_qc_report.csv")
rep.head(8)
"""),
md("### 4.1 ตัวอย่างสถานีค่าค้าง (STUCK) เทียบกับสถานีปกติ"),
code("""
stuck = qcd[qcd.qc == "STUCK"].station_code.unique()
print("สถานีค่าค้าง:", list(stuck))
fig, ax = pl.subplots(width="one_half", height_mm=60)
for c, col in zip(list(stuck)[:3], ["#d95f02", "#7570b3", "#e7298a"]):
    s = qcd[qcd.station_code == c]
    ax.step(s.time_lt, s.rain, where="post", color=col, lw=1.3, label=f"{c} (STUCK)")
s = qcd[qcd.station_code == "MEWG"]
ax.step(s.time_lt, s.rain, where="post", color="0.3", lw=1, label="MEWG (ปกติ)")
ax.set_ylabel("ฝน 1 ชม. (มม.)"); ax.set_xlabel("เวลาไทย 3 ต.ค. 2567"); ax.set_yscale("symlog", linthresh=1)
ax.legend(ncol=2); ax.grid(alpha=0.3, lw=0.4)
import matplotlib.dates as mdates; ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
pl.save(fig, "06_stuck_gauges"); plt.show()
"""),
md("## 5. ความพร้อมของข้อมูล (ลุ่มน้ำปิงและวัง)"),
code("""
sub = qcd.merge(meta[["station_code", "basin", "basin_en"]], on="station_code")
sub = sub[sub.basin.isin(["ปิง", "วัง"])]
code_map = {"OK": 0, "ISOLATED": 1, "STUCK": 2, "FLAG": 3, "RANGE": 3, "NODATA": 4}
mat = sub.assign(c=sub.qc.map(code_map)).pivot_table(index=["basin_en", "station_code"], columns="time_lt", values="c", aggfunc="first")
from matplotlib.colors import ListedColormap, BoundaryNorm
cmq = ListedColormap(["#1b9e77", "#e6ab02", "#d95f02", "#7570b3", "#bdbdbd"])
fig, ax = pl.subplots(width="one_half", height_mm=170)
im = ax.imshow(mat.values, aspect="auto", cmap=cmq, norm=BoundaryNorm(np.arange(-0.5, 5), 5), interpolation="nearest")
ax.set_yticks(range(len(mat))); ax.set_yticklabels([f"{b[:4]} {c}" for b, c in mat.index], fontsize=5.5)
ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels([t.strftime("%H") for t in mat.columns])
ax.set_xlabel("ชั่วโมง (เวลาไทย, ปลายชั่วโมง)")
cb = fig.colorbar(im, ax=ax, ticks=range(5), shrink=0.4); cb.ax.set_yticklabels(["ใช้ได้", "โดดเดี่ยว (เก็บไว้)", "ค่าค้าง", "flag/ช่วงค่า", "ไม่มีข้อมูล"])
pl.save(fig, "06_availability_ping_wang"); plt.show()
print(f"ใช้ได้ {100*(sub.qc=='OK').mean():.1f}% ของสถานี-ชั่วโมง ในลุ่มปิง+วัง")
"""),
md("## 6. พิสูจน์ข้อตกลงด้านเวลาด้วยเรดาร์\nถ้า `measure_datetime` คือปลายชั่วโมง ค่าที่ประทับ 08:00 น. ต้องสัมพันธ์กับเรดาร์ 00–01 UTC สูงสุด"),
code("""
rows = []
for tag, t0 in [("00UTC", "2024-10-03 00:00"), ("10UTC", "2024-10-03 10:00")]:
    d = np.load(P.products / f"mosaic_20241003_{tag}.npz")
    lc = verify.lag_correlation(qcd, meta, d["h_PHY_DP"], d["glon"], d["glat"], pd.Timestamp(t0))
    lc["radar_hour"] = tag; rows.append(lc)
lag = pd.concat(rows); lag.to_csv(P.tables / "06_lag_check.csv", index=False)
fig, ax = pl.subplots(width="single", height_mm=60)
for tag, g_ in lag.groupby("radar_hour"):
    ax.plot(g_.gauge_end_local.dt.strftime("%H:%M"), g_.r, "o-", label=f"เรดาร์ {tag}")
ax.set_xlabel("เวลาประทับของฝนสถานี (เวลาไทย)"); ax.set_ylabel("สหสัมพันธ์ r"); ax.set_ylim(0, 1); ax.grid(alpha=0.3, lw=0.4); ax.legend()
pl.save(fig, "06_time_alignment_check"); plt.show()
lag
"""),
md("""
**แปลผล:** สหสัมพันธ์สูงสุดที่ 08:00 น. (กับ 00–01 UTC) และ 18:00 น. (กับ 10–11 UTC) ยืนยันว่า **เวลาประทับคือปลายชั่วโมง**
ถ้าใช้ผิด (07:00 หรือ 17:00) ค่า r จะลดลงอย่างชัดเจน — ความผิดพลาดเรื่องเวลาเป็นสาเหตุอันดับต้นของการประเมินเรดาร์ผิด

## 7. แผนที่ฝนสถานี 18:00 น. (= 10–11 UTC)
"""),
code("""
h18 = qcd[(qcd.time_lt == "2024-10-03 18:00") & (qcd.qc == "OK")].merge(meta, on="station_code")
E = cfg["study_area"]["map_extent"]["north"]
fig, ax = pl.subplots(width="one_half", height_mm=150)
cm, nm = pl.rain_cmap()
dry = h18[h18.rain_qc < 0.1]; wet = h18[h18.rain_qc >= 0.1].sort_values("rain_qc")
ax.scatter(dry.longitude, dry.latitude, s=8, c="white", edgecolors="0.4", lw=0.5, zorder=12, label="< 0.1 มม.")
sc = ax.scatter(wet.longitude, wet.latitude, s=12 + 5 * wet.rain_qc, c=wet.rain_qc, cmap=cm, norm=nm, edgecolors="k", lw=0.4, zorder=13)
mp.basemap(ax, E, title="ฝนสถานี สสน. 17:00–18:00 น. (10–11 UTC)", scale_km=100)
for _, r in wet[wet.rain_qc >= 10].iterrows():
    ax.annotate(f"{r.station_code} {r.rain_qc:.0f}", (r.longitude, r.latitude), xytext=(4, -8), textcoords="offset points", fontsize=6.5, zorder=14)
fig.colorbar(sc, ax=ax, orientation="horizontal", shrink=0.8, extend="max").set_label("ฝน 1 ชั่วโมง (มม.)")
ax.legend(loc="upper left"); mp.source_note(ax, "ข้อมูล: สสน. (HII) open data")
pl.save(fig, "06_gauge_map_18LT"); plt.show()
"""),
md("## 8. บันทึกตารางฝนสถานีที่สะอาดแล้ว (ใช้ต่อใน notebook 07–08)"),
code("""
clean = qcd.merge(meta[["station_code", "station_name", "latitude", "longitude", "basin", "basin_en", "sub_basin", "province"]], on="station_code")
clean.to_csv(P.tables / "06_gauge_hourly_qc.csv", index=False, encoding="utf-8-sig")
print(clean.shape); clean.head()
"""),
md("---\n## ✏️ แบบฝึกหัด"),
exercise(1, "ug", "ค่าฝนของสถานีที่ประทับเวลา **2024-10-03 18:00 (เวลาไทย)** เป็นฝนช่วงเวลา UTC ใด และตรงกับไฟล์เรดาร์ช่วงใด"),
todo("""
start_utc = None; end_utc = None
"""),
solution(1, """
t = pd.Timestamp("2024-10-03 18:00")
end_utc = t - pd.Timedelta(hours=7); start_utc = end_utc - pd.Timedelta(hours=1)
print(f"ฝนช่วง {start_utc} ถึง {end_utc} UTC → ใช้ไฟล์เรดาร์ 10:00, 10:15, 10:30, 10:45 UTC")
assert start_utc == pd.Timestamp("2024-10-03 10:00")
"""),
exercise(2, "ug", "สถานีใดถูกตัดเป็น **STUCK** และมีค่ารวม 15 ชั่วโมงเท่าใด ถ้าไม่ตัดทิ้ง จะส่งผลต่อการประเมินเรดาร์อย่างไร (คิดถึงเรดาร์ที่บอกว่าฝนเป็น 0)"),
todo("""
stuck_sum = None
"""),
solution(2, """
stuck_sum = qcd[qcd.station_code.isin(stuck)].groupby("station_code").rain.sum()
print(stuck_sum)
print("ทำให้นับเป็น 'เรดาร์พลาดฝน' (miss) ทั้งที่จริงอาจไม่มีฝน → POD ต่ำลง และ bias ติดลบเกินจริง")
"""),
exercise(3, "ug", "ช่วง **17:00 และ 18:00 น.** มีกี่สถานีที่ฝน ≥ 10 มม./ชม. อยู่ลุ่มน้ำใดบ้าง"),
todo("""
heavy = None
"""),
solution(3, """
heavy = clean[(clean.time_lt.isin(pd.to_datetime(["2024-10-03 17:00", "2024-10-03 18:00"]))) & (clean.rain_qc >= 10)]
display(heavy[["time_lt", "station_code", "station_name", "basin_en", "rain_qc"]].sort_values("rain_qc", ascending=False))
"""),
exercise(4, "pg", "สร้าง **แผนที่ฝนจากสถานีอย่างเดียว** ด้วย IDW (กำลัง 2, รัศมี 50 กม.) ที่ 18:00 น. บนกริด 0.02° (ใช้ `radteach.bias.idw`) แล้วเทียบกับแผนที่เรดาร์ PHY_DP 10 UTC — อะไรที่สถานีมองไม่เห็น"),
todo("""
idw_map = None
"""),
solution(4, """
from radteach import bias
lo = np.arange(97.0, 101.6, 0.02); la = np.arange(21.0, 15.6, -0.02)
LO, LA = np.meshgrid(lo, la)
idw_map = bias.idw(LO, LA, h18.longitude.values, h18.latitude.values, h18.rain_qc.values, power=2, radius_km=50, background=np.nan).reshape(LO.shape)
d = np.load(P.products / "mosaic_20241003_10UTC.npz")
fig, axs = pl.subplots(1, 2, width="double", height_mm=105)
mp.show_grid(axs[0], idw_map, lo, la, cm, nm); mp.basemap(axs[0], E, title="IDW จากสถานี 18:00 น.", scale_km=100)
im = mp.show_grid(axs[1], d["h_PHY_DP"], d["glon"], d["glat"], cm, nm); mp.basemap(axs[1], E, title="เรดาร์ PHY_DP 10–11 UTC", scale_km=100, north=False)
fig.colorbar(im, ax=axs, orientation="horizontal", shrink=0.5, extend="max").set_label("มม."); plt.show()
print("สถานีห่างกัน 10–30 กม. มองไม่เห็นเซลล์พายุขนาดเล็ก (< 10 กม.) ที่ตกระหว่างสถานี แต่เป็นค่าจริงที่พื้น")
"""),
md("""
---
## 💡 คำถามชวนคิดและต่อยอด
- ถ้าจะประเมินเรดาร์ทั้งเหตุการณ์ (20 ก.ย.–7 ต.ค. 2567) ควรดึงข้อมูลอย่างไรให้สุภาพต่อเซิร์ฟเวอร์ (ความถี่การเรียก, cache, checksum)
- สถานีของกรมอุตุนิยมวิทยาและกรมชลประทานมีมาตรฐานเวลาและการ QC ต่างจาก สสน. อย่างไร ถ้าจะรวมหลายเครือข่ายต้องระวังอะไร
- ความหนาแน่นของสถานีในต้นน้ำปิงเพียง ~1 สถานี/840 กม.² เพียงพอหรือไม่สำหรับการปรับแก้เรดาร์ (notebook 08)
- 🔵 ออกแบบ QC แบบ spatial–temporal ที่ใช้ทั้งเพื่อนบ้านและประวัติของสถานีเอง (เช่น median ของ 30 วันที่ผ่านมา)

## 📚 อ้างอิง
- สถาบันสารสนเทศทรัพยากรน้ำ (องค์การมหาชน). คลังข้อมูลเปิด hourly_rain. https://tiservice.hii.or.th/opendata/
- Habib, E. et al. (2001). Sampling errors of tipping-bucket rain gauge measurements. *J. Hydrol. Eng.* 6:159–166.
- Sevruk, B. (2005). Rainfall measurement: gauges. In *Encyclopedia of Hydrological Sciences*. Wiley.
- Lewis, E. et al. (2021). Quality control of a global hourly rainfall dataset. *Environ. Model. Softw.* 144:105169.
"""),
save_cell("06"),
]
build(cells, f"notebooks/{F}")
