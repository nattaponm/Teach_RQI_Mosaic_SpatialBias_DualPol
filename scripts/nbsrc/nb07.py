from _nb import *

F = "07_verification_interpretation.ipynb"
cells = [
header(F, "07", "การตรวจสอบความถูกต้องของฝนเรดาร์และการแปลผล", "Verifying radar rainfall against gauges and interpreting the results",
       "🟢 ป.ตรี ปี 3–4 และ 🔵 ป.โท", "120 นาที", "04, 05, 06", "ผลสำเร็จรูป (โมเสค + advection 2 ชั่วโมง), ฝนสถานี สสน. 221 สถานี",
       "ตาราง metrics รายผลิตภัณฑ์/ลุ่มน้ำ/ชั้น RQI, scatter plot, แผนที่ bias"),
md("""
## 🎯 วัตถุประสงค์การเรียนรู้
1. **จับคู่** ฝนสถานีกับฝนเรดาร์ให้ถูกทั้งตำแหน่งและเวลา และเข้าใจปัญหาจุด (point) เทียบพื้นที่ (area)
2. **คำนวณและตีความ** ดัชนีความถูกต้องแบบต่อเนื่อง (PBIAS, RMSE, r, NSE, KGE) และแบบหมวดหมู่ (POD, FAR, CSI, HSS)
3. **เปรียบเทียบ** ผลิตภัณฑ์ทุกวิธีแยกตามลุ่มน้ำและตามชั้นคุณภาพ RQI
4. **อธิบาย** สาเหตุทางกายภาพของความคลาดเคลื่อน (clutter, การบัง, ลำคลื่นสูง, การปรับเทียบ, Z–R)
5. **วิพากษ์** ข้อจำกัดของการประเมินที่มีตัวอย่างน้อย

## 1. หลักการ
### 1.1 จุดเทียบกับพื้นที่
สถานีวัดฝนเก็บน้ำจากปากกระบอก ~0.03 ม.² ขณะที่พิกเซลเรดาร์ ~1 กม.² และลำคลื่นอยู่สูงจากพื้นหลายร้อยเมตรถึงหลายกิโลเมตร
→ แม้เรดาร์สมบูรณ์แบบก็ไม่ตรงกับสถานี 100% (**representativeness error**) จึงนิยมใช้ค่าเฉลี่ย 3×3 พิกเซล (~3 กม.) ลดผลการเลื่อนตำแหน่ง

### 1.2 ดัชนีแบบต่อเนื่อง (R = เรดาร์, G = สถานี)
$$\\text{PBIAS} = 100\\frac{\\sum(R-G)}{\\sum G},\\quad B = \\frac{\\sum R}{\\sum G},\\quad \\text{RMSE}=\\sqrt{\\tfrac1n\\sum(R-G)^2},\\quad \\text{NSE}=1-\\frac{\\sum(R-G)^2}{\\sum(G-\\bar G)^2}$$
$$\\text{KGE} = 1-\\sqrt{(r-1)^2+(\\alpha-1)^2+(\\beta-1)^2},\\qquad \\alpha=\\sigma_R/\\sigma_G,\\ \\beta=\\mu_R/\\mu_G\\quad(\\text{Gupta et al. 2009})$$

### 1.3 ดัชนีแบบหมวดหมู่ (ฝน/ไม่ฝน ที่เกณฑ์ 0.5 มม./ชม.)
| | สถานี: ฝน | สถานี: ไม่ฝน |
|---|---|---|
| **เรดาร์: ฝน** | hits *a* | false alarms *b* |
| **เรดาร์: ไม่ฝน** | misses *c* | correct negatives *d* |

$$\\text{POD}=\\frac{a}{a+c},\\quad \\text{FAR}=\\frac{b}{a+b},\\quad \\text{CSI}=\\frac{a}{a+b+c},\\quad \\text{HSS}=\\frac{2(ad-bc)}{(a+c)(c+d)+(a+b)(b+d)}$$

> ⚠️ ชั่วโมงส่วนใหญ่ "แห้งทั้งคู่" ทำให้ r และ NSE ดูดีเกินจริง จึงรายงาน **r เฉพาะชั่วโมงที่มีฝน (r_wet)** ด้วย
"""),
setup_cell(("products",)),
md("## 2. เตรียมฝนสถานี (QC + UTC) และสร้างคู่ข้อมูล"),
code("""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from radteach import gauges as gg, verify, mosaic, config, plotting as pl, mapping as mp
cfg = config.load_all()
meta = gg.read_metadata(P.stations / "hii_station_metadata_202410.csv")
raw = gg.read_station_csv(P.stations / "hii_hourly_rain_20241003_0519LT.csv")
g = gg.add_utc(gg.quality_control(raw, meta=meta))
HEAD = ["แม่น้ำปิงตอนบนส่วนที่ 1", "น้ำแม่งัด", "แม่น้ำปิงตอนบนส่วนที่ 2", "แม่น้ำแม่แตง", "แม่น้ำปิงส่วนที่ 2/1", "น้ำแม่ริม"]
PRODUCTS = ["LS_MP", "LS_DP", "HYZ_MP", "PHY_DP", "SCAN_HOLD", "ADV_15", "ADV_05"]
rows = []
for tag, t0 in [("00UTC", "2024-10-03 00:00"), ("10UTC", "2024-10-03 10:00")]:
    mo = np.load(P.products / f"mosaic_20241003_{tag}.npz"); ad = np.load(P.products / f"advection_20241003_{tag}.npz", allow_pickle=True)
    fields = {p: mo[f"h_{p}"] for p in mosaic.PRODUCTS}
    fields.update({p: ad[p] for p in ["SCAN_HOLD", "ADV_15", "ADV_05"]})
    gh = g[(g.start_utc == t0) & (g.qc == "OK")].merge(meta[["station_code", "longitude", "latitude"]], on="station_code")
    gh = gh.rename(columns={"longitude": "lon", "latitude": "lat", "rain_qc": "G"})[["station_code", "start_utc", "lon", "lat", "G"]]
    rows.append(verify.pair(gh, fields, mo["glon"], mo["glat"], window=3,
                            rqi={"rqi_LS": mo["rqi_LS"], "rqi_HYZ": mo["rqi_HYZ"], "rqi_PHY": mo["rqi_PHY"]}))
pairs = pd.concat(rows, ignore_index=True).merge(meta[["station_code", "station_name", "basin", "basin_en", "sub_basin"]], on="station_code")
pairs["headwater"] = pairs.sub_basin.isin(HEAD)
pairs.to_csv(P.tables / "07_pairs.csv", index=False, encoding="utf-8-sig")
print(f"คู่ข้อมูล {len(pairs)} (สถานี-ชั่วโมง) จาก {pairs.station_code.nunique()} สถานี")
pairs.sort_values("G", ascending=False).head(8)[["station_code", "start_utc", "basin_en", "G", "LS_MP", "PHY_DP", "ADV_05", "rqi_PHY"]].round(2)
"""),
md("## 3. ภาพรวม: ทุกผลิตภัณฑ์ทุกสถานี"),
code("""
M = verify.metrics_table(pairs, PRODUCTS, thr=0.5)
show = ["product", "n", "G_mean", "R_mean", "PBIAS", "B", "RMSE", "r", "r_wet", "KGE", "POD", "FAR", "CSI", "HSS"]
M[show].round(2)
"""),
code("""
FLOOR, TOP = 0.05, 40
show_p = ["LS_MP", "LS_DP", "HYZ_MP", "PHY_DP", "SCAN_HOLD", "ADV_05"]
title = {"LS_MP": "Z–R มุมต่ำสุด", "LS_DP": "Dual-pol มุมต่ำสุด", "HYZ_MP": "Z–R hybrid + VPR",
         "PHY_DP": "Dual-pol hybrid + VPR", "SCAN_HOLD": "PHY_DP คงค่า 15 นาที", "ADV_05": "PHY_DP advection 5 นาที"}
rqi_of = {"LS_MP": "rqi_LS", "LS_DP": "rqi_LS", "HYZ_MP": "rqi_HYZ"}
cq, nq = pl.rqi_cmap((0, 0.3, 0.6, 0.9, 1.0))
fig, axs = pl.subplots(2, 3, width="double", height_mm=140)
for ax, p in zip(axs.ravel(), show_p):
    for h, mk in zip(sorted(pairs.start_utc.unique()), ["o", "^"]):
        s = pairs[pairs.start_utc == h]
        sc = ax.scatter(np.clip(s.G, FLOOR, None), np.clip(s[p], FLOOR, None), c=s[rqi_of.get(p, "rqi_PHY")], cmap=cq, norm=nq,
                        s=14, marker=mk, edgecolor="k", lw=0.3, label=f"{pd.Timestamp(h):%H}–{pd.Timestamp(h) + pd.Timedelta(hours=1):%H} UTC")
    pl.one_to_one(ax, FLOOR, TOP)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(FLOOR * 0.8, TOP); ax.set_ylim(FLOOR * 0.8, TOP)
    m = verify.metrics(pairs.G, pairs[p])
    ax.text(0.03, 0.97, f"PBIAS {m['PBIAS']:.0f}%\\nr(ฝน) {m.get('r_wet', np.nan):.2f}\\nCSI {m['CSI']:.2f}", transform=ax.transAxes,
            va="top", fontsize=7, bbox=dict(fc="white", ec="0.6", lw=0.4, alpha=0.8))
    ax.set_title(title[p], loc="left", fontweight="bold", fontsize=8.5)
for ax in axs[1]: ax.set_xlabel(pl.L("gauge_mm"))
for ax in axs[:, 0]: ax.set_ylabel(pl.L("radar_mm"))
axs[0, 0].legend(loc="lower right", fontsize=6.5)
fig.colorbar(sc, ax=axs, orientation="horizontal", shrink=0.5).set_label("RQI ที่ตำแหน่งสถานี  (ค่า < 0.05 มม. วางไว้ที่ขอบ)")
pl.panel_labels(axs); pl.save(fig, "07_scatter_products"); plt.show()
"""),
md("## 4. แยกตามลุ่มน้ำ"),
code("""
groups = {"ทั้งพื้นที่": pairs, "ลุ่มน้ำปิง": pairs[pairs.basin == "ปิง"], "ต้นน้ำปิง": pairs[pairs.headwater],
          "ลุ่มน้ำวัง": pairs[pairs.basin == "วัง"], "ลุ่มน้ำยม": pairs[pairs.basin == "ยม"], "ลุ่มน้ำน่าน": pairs[pairs.basin == "น่าน"],
          "สาละวิน": pairs[pairs.basin == "สาละวิน"], "โขงเหนือ": pairs[pairs.basin == "โขงเหนือ"]}
MB = pd.concat([verify.metrics_table(d, PRODUCTS).assign(group=k) for k, d in groups.items()])
MB.to_csv(P.tables / "07_metrics_by_basin.csv", index=False, encoding="utf-8-sig")
fig, axs = pl.subplots(1, 3, width="double", height_mm=70)
x = np.arange(len(PRODUCTS)); cols = plt.get_cmap("tab10")
for i, (k, d) in enumerate(MB.groupby("group", sort=False)):
    d = d.set_index("product").reindex(PRODUCTS)
    for ax, col_, lim in [(axs[0], "PBIAS", (-100, 40)), (axs[1], "r_wet", (-0.2, 1)), (axs[2], "CSI", (0, 1))]:
        ax.plot(x, d[col_], "-o", ms=3, lw=0.9, color=cols(i), label=k); ax.set_ylim(*lim)
for ax, t in zip(axs, ["PBIAS (%)", "r (ชั่วโมงที่มีฝน)", "CSI (≥ 0.5 มม./ชม.)"]):
    ax.set_xticks(x); ax.set_xticklabels(PRODUCTS, rotation=60, ha="right", fontsize=6.5); ax.set_title(t, loc="left", fontweight="bold"); ax.grid(alpha=0.3, lw=0.4)
axs[0].axhline(0, color="k", lw=0.6)
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="upper center", ncol=8, bbox_to_anchor=(0.5, 1.10), fontsize=7)
pl.panel_labels(axs); pl.save(fig, "07_metrics_by_basin"); plt.show()
MB.pivot_table(index="group", columns="product", values="PBIAS", sort=False)[PRODUCTS].round(0)
"""),
md("## 5. แยกตามชั้นคุณภาพ RQI — RQI บอกความน่าเชื่อถือได้จริงหรือไม่"),
code("""
bins = [0, 0.3, 0.6, 0.9, 1.01]; labs = ["< 0.3", "0.3–0.6", "0.6–0.9", "≥ 0.9"]
rq = []
for p in ["LS_MP", "HYZ_MP", "PHY_DP"]:
    cls = pd.cut(pairs[rqi_of.get(p, "rqi_PHY")], bins, labels=labs, right=False)
    for c in labs:
        s = pairs[cls == c]; rq.append({"product": p, "RQI": c, **verify.metrics(s.G, s[p])})
RQ = pd.DataFrame(rq); RQ.to_csv(P.tables / "07_metrics_by_rqi.csv", index=False)
fig, axs = pl.subplots(1, 2, width="one_half", height_mm=60)
for i, p in enumerate(["LS_MP", "HYZ_MP", "PHY_DP"]):
    d = RQ[RQ["product"] == p]
    axs[0].bar(np.arange(4) + (i - 1) * 0.27, d.PBIAS, 0.27, label=p)
    axs[1].bar(np.arange(4) + (i - 1) * 0.27, d.r, 0.27, label=p)
for ax, t in zip(axs, ["PBIAS (%)", "r"]):
    ax.set_xticks(range(4)); ax.set_xticklabels(labs); ax.set_xlabel("ชั้น RQI"); ax.set_title(t, loc="left", fontweight="bold"); ax.grid(alpha=0.3, lw=0.4, axis="y")
axs[0].axhline(0, color="k", lw=0.6); axs[1].legend(fontsize=6.5)
pl.panel_labels(axs); pl.save(fig, "07_metrics_by_rqi"); plt.show()
RQ[["product", "RQI", "n", "G_mean", "R_mean", "PBIAS", "r", "RMSE"]].round(2)
"""),
md("## 6. แผนที่อัตราส่วนเรดาร์/สถานี"),
code("""
st = pairs.groupby("station_code").agg(lon=("lon", "first"), lat=("lat", "first"), G=("G", "sum"), LS_MP=("LS_MP", "sum"), PHY_DP=("PHY_DP", "sum"))
cr, nr = pl.ratio_cmap()
E = cfg["study_area"]["map_extent"]["north"]
fig, axs = pl.subplots(1, 2, width="double", height_mm=110)
for ax, p, t in [(axs[0], "LS_MP", "Z–R มุมต่ำสุด"), (axs[1], "PHY_DP", "Dual-pol hybrid + VPR")]:
    wet = (st.G >= 1) | (st[p] >= 1)
    ax.scatter(st.lon[~wet], st.lat[~wet], s=5, c="0.7", lw=0, zorder=12, label="แห้งทั้งคู่")
    sc = ax.scatter(st.lon[wet], st.lat[wet], c=np.log2((st[p][wet] + 0.2) / (st.G[wet] + 0.2)), cmap=cr, norm=nr, s=28, edgecolor="k", lw=0.3, zorder=13)
    mp.basemap(ax, E, title=t, scale_km=100, north=(p == "LS_MP"))
axs[0].legend(loc="upper left")
cb = fig.colorbar(sc, ax=axs, orientation="horizontal", shrink=0.5, ticks=[-3, -2, -1, 0, 1, 2, 3])
cb.ax.set_xticklabels(["1/8", "1/4", "1/2", "1", "2", "4", "8"]); cb.set_label("ฝนเรดาร์ / ฝนสถานี (รวม 2 ชั่วโมง)")
pl.panel_labels(axs); pl.save(fig, "07_bias_ratio_map"); plt.show()
"""),
md("""
## 7. แปลผล (อภิปรายในชั้นเรียน)
1. **เรดาร์ต่ำกว่าสถานีมากทุกวิธี** (ΣR/ΣG ≈ 0.2–0.4) แม้ที่ RQI ≥ 0.9 → ไม่ใช่แค่ปัญหาภูเขา แต่มี **bias ของการปรับเทียบ Z และ/หรือ Z–R** ร่วมด้วย
   → ขั้นการปรับแก้ด้วยฝนสถานี (notebook 08) จำเป็น
2. **มุมต่ำสุด (LS) ให้ PBIAS ดีที่สุดแต่ r ต่ำมาก** เพราะมี "ฝนปลอม" จาก clutter/AP รอบเรดาร์พิษณุโลก (สถานีแห้งแต่เรดาร์ 1–8 มม.)
   hybrid scan ตัดออกได้ — **ค่าเฉลี่ยดีไม่ได้แปลว่าผลิตภัณฑ์ดี** ต้องดูหลายดัชนีพร้อมกัน
3. **Dual-pol ดีกว่า Z–R** ทั้งแบบ LS และ hybrid (PBIAS ดีขึ้น 5–7 จุด, r_wet สูงขึ้น 0.1–0.15)
4. **RQI มีความหมายจริง**: ชั้น RQI ต่ำ bias รุนแรงและ r ≈ 0
5. **Advection ไม่ช่วยในระดับชั่วโมง-จุด** สำหรับ 2 ชั่วโมงนี้
6. **ข้อจำกัด**: 2 ชั่วโมง ฝนเบา คู่ที่มีฝน ~60 คู่ ถูกครอบงำโดยไม่กี่สถานี → ต้องขยายเป็นทั้งเหตุการณ์ก่อนสรุปลงงานวิจัย

---
## ✏️ แบบฝึกหัด
"""),
exercise(1, "ug", "คำนวณ **PBIAS และ RMSE ของ PHY_DP ด้วยตนเอง** (ใช้ numpy) แล้วเทียบกับ `verify.metrics`"),
todo("""
pbias = None; rmse = None
"""),
solution(1, """
k = np.isfinite(pairs.G) & np.isfinite(pairs.PHY_DP)
G_, R_ = pairs.G[k].values, pairs.PHY_DP[k].values
pbias = 100 * (R_ - G_).sum() / G_.sum(); rmse = np.sqrt(np.mean((R_ - G_) ** 2))
m = verify.metrics(pairs.G, pairs.PHY_DP)
print(f"PBIAS {pbias:.2f}% (ฟังก์ชัน {m['PBIAS']:.2f}%),  RMSE {rmse:.3f} มม. (ฟังก์ชัน {m['RMSE']:.3f})")
assert abs(pbias - m["PBIAS"]) < 1e-9
"""),
exercise(2, "ug", "สร้าง **ตาราง contingency ที่เกณฑ์ 1 มม./ชม.** สำหรับ LS_MP และ PHY_DP แล้วคำนวณ POD, FAR, CSI ผลิตภัณฑ์ใดมี false alarm มากกว่า เพราะอะไร"),
todo("""
ct = None
"""),
solution(2, """
ct = {}
for p in ["LS_MP", "PHY_DP"]:
    k = np.isfinite(pairs.G) & np.isfinite(pairs[p]); G_, R_ = pairs.G[k], pairs[p][k]
    a = int(((G_ >= 1) & (R_ >= 1)).sum()); b = int(((G_ < 1) & (R_ >= 1)).sum()); c = int(((G_ >= 1) & (R_ < 1)).sum())
    ct[p] = dict(hits=a, false_alarms=b, misses=c, POD=a / (a + c), FAR=b / max(a + b, 1), CSI=a / (a + b + c))
display(pd.DataFrame(ct).T.round(2))
print("LS_MP มี false alarm มากกว่า (clutter/AP ใกล้ PHS ที่มุม 0.48°)")
"""),
exercise(3, "ug", "เรียงลำดับผลิตภัณฑ์ตาม **KGE** และตาม **r_wet** ลำดับเหมือนกันหรือไม่ ถ้าต้องเลือก 1 ผลิตภัณฑ์ไปใช้เตือนภัยน้ำท่วม จะเลือกอะไร เพราะอะไร"),
todo("""
rank = None
"""),
solution(3, """
rank = M.set_index("product")[["KGE", "r_wet", "PBIAS", "CSI"]].sort_values("KGE", ascending=False).round(2)
display(rank)
print("KGE ลงโทษ bias หนัก ส่วน r_wet ดูรูปแบบ — สำหรับเตือนภัย รูปแบบเชิงพื้นที่ (r, CSI) สำคัญ และแก้ bias ได้ด้วยสถานี (nb 08) → PHY_DP")
"""),
exercise(4, "pg", "ประเมิน **ความไม่แน่นอน** ของผลต่าง PBIAS ระหว่าง PHY_DP กับ LS_MP ด้วย bootstrap (สุ่มสถานีแบบใส่คืน 1,000 ครั้ง) ช่วงความเชื่อมั่น 95% ครอบคลุม 0 หรือไม่ หมายความว่าอะไร"),
todo("""
ci = None
"""),
solution(4, """
rng = np.random.default_rng(42); codes = pairs.station_code.unique(); diffs = []
by = {c: d for c, d in pairs.groupby("station_code")}
for _ in range(1000):
    s = pd.concat([by[c] for c in rng.choice(codes, len(codes))])
    diffs.append(verify.metrics(s.G, s.PHY_DP)["PBIAS"] - verify.metrics(s.G, s.LS_MP)["PBIAS"])
ci = np.nanpercentile(diffs, [2.5, 97.5])
print(f"ΔPBIAS (PHY − LS) 95% CI = [{ci[0]:.1f}, {ci[1]:.1f}] จุดร้อยละ")
print("ถ้าครอบคลุม 0 แปลว่าตัวอย่างน้อยเกินกว่าจะสรุปว่าวิธีหนึ่งดีกว่าอีกวิธีในด้าน bias")
"""),
exercise(5, "pg", "ทดสอบความไวต่อ **ขนาดหน้าต่างการสุ่ม** (1×1, 3×3, 5×5) สำหรับ PHY_DP ดัชนีใดเปลี่ยนมากที่สุด เพราะอะไร"),
todo("""
sens = None
"""),
solution(5, """
res = []
for w in [1, 3, 5]:
    rr = []
    for tag, t0 in [("00UTC", "2024-10-03 00:00"), ("10UTC", "2024-10-03 10:00")]:
        mo = np.load(P.products / f"mosaic_20241003_{tag}.npz"); s = pairs[pairs.start_utc == t0]
        rr.append(pd.DataFrame({"G": s.G.values, "R": verify.sample(mo["h_PHY_DP"], mo["glon"], mo["glat"], s.lon, s.lat, w)}))
    d = pd.concat(rr); res.append({"window": f"{w}x{w}", **verify.metrics(d.G, d.R)})
sens = pd.DataFrame(res).set_index("window")[["PBIAS", "RMSE", "r", "r_wet", "POD", "FAR", "CSI"]].round(2)
display(sens)
print("หน้าต่างใหญ่ขึ้น: PBIAS ติดลบน้อยลง (ดึงฝนจากพิกเซลรอบ ๆ เข้ามา) แต่ r ลดลงเพราะค่าถูกเกลี่ย และ FAR เพิ่มเพราะฝนข้างเคียง 'ไหล' เข้าจุดที่แห้ง")
print("→ ไม่มีขนาดหน้าต่างที่ถูกเสมอ ต้องรายงานว่าใช้ขนาดใดและทดสอบความไว")
"""),
md("""
---
## 💡 คำถามชวนคิดและต่อยอด
- ถ้าประเมินด้วยฝนสะสมทั้งเหตุการณ์ (เช่น 72 ชั่วโมง) แทนรายชั่วโมง ดัชนีจะเปลี่ยนอย่างไร และเหมาะกับคำถามแบบไหน
- ควรรายงานผลแยกตามความสูงภูมิประเทศของสถานี (ที่ราบ/เชิงเขา/ภูเขา) หรือไม่ เกี่ยวกับ "ฝนบนภูเขา" อย่างไร
- การตรวจสอบนี้ใช้สถานีทั้งหมด แต่ถ้าใช้สถานีไปปรับแก้เรดาร์แล้ว จะตรวจสอบอย่างไรให้ "ยุติธรรม" (notebook 08: cross-validation)
- 🔵 ลองดัชนีเชิงพื้นที่ เช่น Fractions Skill Score (FSS) ที่ไม่ลงโทษการเลื่อนตำแหน่งเล็กน้อย

## 📚 อ้างอิง
- Gupta, H.V. et al. (2009). Decomposition of the mean squared error and NSE performance criteria. *J. Hydrol.* 377:80–91.
- Wilks, D.S. (2019). *Statistical Methods in the Atmospheric Sciences*, 4th ed., ch. 9. Elsevier.
- Ciach, G.J. & Krajewski, W.F. (1999). On the estimation of radar rainfall error variance. *Adv. Water Resour.* 22:585–595.
- Roberts, N.M. & Lean, H.W. (2008). Scale-selective verification of rainfall accumulations. *Mon. Wea. Rev.* 136:78–97.
"""),
save_cell("07"),
]
build(cells, f"notebooks/{F}")
