# ข้อมูลและเงื่อนไขการใช้ / Data sources and terms

| ข้อมูล | ที่อยู่ใน repo | แหล่ง | เงื่อนไข / การอ้างอิง |
|---|---|---|---|
| เรดาร์ UF 5 สถานี (3 ต.ค. 2567, 00–01 และ 10–11 UTC) | Release `data-v1.0` → `teaching_data_v1.zip` | กรมอุตุนิยมวิทยา (Thai Meteorological Department) | เผยแพร่เพื่อการศึกษาโดยได้รับความเห็นชอบจากผู้ดูแลรายวิชา; อ้างอิง "Thai Meteorological Department" ทุกครั้งที่ใช้ |
| ผลลัพธ์สำเร็จรูป (โมเสค, RQI, VPR, advection) | Release → `teaching_products_v1.zip` | สร้างจากข้อมูลข้างต้นด้วยโค้ดใน repo นี้ | CC BY 4.0 + อ้างอิงกรมอุตุนิยมวิทยา |
| DEM SRTM v4.1 ลดความละเอียด 250 ม. | Release → `teaching_data_v1.zip/dem` | CGIAR-CSI (Jarvis et al. 2008) จาก NASA SRTM | อ้างอิง Jarvis, A., Reuter, H.I., Nelson, A., Guevara, E. (2008) Hole-filled SRTM for the globe Version 4, CGIAR-CSI SRTM 90m Database. เงื่อนไข CIAT: ห้ามแจกจ่ายต่อ (ทั้งเชิงพาณิชย์และไม่คิดค่าใช้จ่าย) โดยไม่ได้รับอนุญาตเป็นลายลักษณ์อักษรจาก CIAT, ต้องอ้างอิง CIAT และแจ้ง CIAT เมื่อมีการแจกจ่าย — **สถานะ: อยู่ระหว่างขออนุญาต / พิจารณาเปลี่ยนเป็น DEM ที่แจกจ่ายต่อได้ (NASA SRTMGL3 v003 หรือ Copernicus GLO-90)** |
| ฝนสถานีรายชั่วโมง (snapshot 221 สถานี, metadata) | `data/static/stations/` | สถาบันสารสนเทศทรัพยากรน้ำ (สสน.) open data: https://tiservice.hii.or.th/opendata/ | อ้างอิง สสน. |
| ขอบเขตลุ่มน้ำหลัก/ลุ่มน้ำย่อย (simplify 0.002°) | `data/static/boundaries/` | สำนักงานทรัพยากรน้ำแห่งชาติ (สทนช.) | อ้างอิง สทนช. |
| ตำแหน่งสถานีวัดน้ำ P.1, P.67, W.1C | `data/static/stations/rid_stations.csv` | กรมชลประทาน | อ้างอิง กรมชลประทาน |
| ฟอนต์ Sarabun | `mapping/fonts/` | Cadson Demak / Google Fonts | SIL Open Font License 1.1 (`mapping/fonts/OFL.txt`) |

ข้อมูลทุกชุดใช้เพื่อการเรียนการสอนและการวิจัย ไม่ใช่ข้อมูลทางการสำหรับการเตือนภัยหรือการตัดสินใจ
Data are provided for teaching and research only and are not an official product for warnings or decisions.
