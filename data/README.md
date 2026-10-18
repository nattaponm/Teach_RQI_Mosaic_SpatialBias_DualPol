# data/

| โฟลเดอร์ | อยู่ใน git | เนื้อหา |
|---|---|---|
| `static/boundaries/` | ✔ | `onwr_main_basins.gpkg`, `onwr_sub_basins.gpkg` (สทนช., simplify 0.002°), `thailand_outline.gpkg` |
| `static/stations/` | ✔ | `hii_station_metadata_202410.csv` (สสน.), `hii_hourly_rain_20241003_0519LT.csv` (snapshot 221 สถานี 05–19 น. เวลาไทย, รูปแบบเดียวกับ open data), `rid_stations.csv` |
| `static/radar_sites.csv` | ✔ | ตำแหน่ง ความสูง มุมยกของเรดาร์ 5 สถานี |
| `raw/teaching_data_v1/` | ✘ | ดาวน์โหลดจาก Release: `radar/20241003_00UTC/*.uf`, `radar/20241003_10UTC/*.uf`, `dem/srtm_*_250m.tif` |
| `raw/teaching_products_v1/` | ✘ | ผลสำเร็จรูป: `mosaic_*.npz`, `advection_*.npz`, `calibration_*.csv`, `dem_study_area.tif`, GeoTIFF รายชั่วโมง |
| `raw/gauges/` | ✘ | ข้อมูล สสน. ที่ดึงสดใน notebook 06 |
| `interim/`, `processed/` | ✘ | ผลกลาง/ผลที่นิสิตคำนวณเอง |

ชื่อไฟล์เรดาร์: `<RADAR>_<YYYYMMDDHHMM>.uf` เวลา **UTC** · ฝนสถานี: `measure_datetime` เวลา **ไทย** ปลายชั่วโมง
