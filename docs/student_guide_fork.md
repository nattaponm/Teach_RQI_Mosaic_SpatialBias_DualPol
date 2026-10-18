# คู่มือนิสิต: Fork → Colab → บันทึกผลขึ้น GitHub ของตนเอง

repository ของอาจารย์ (`nattaponm/Teach_RQI_Mosaic_SpatialBias_DualPol`) เป็น **อ่านอย่างเดียว** สำหรับนิสิต
ผลงานของนิสิตอยู่ใน **fork ของตนเอง** ซึ่งเป็นสำเนาส่วนตัวบนบัญชี GitHub ของนิสิต

## ขั้นที่ 1 — Fork (ทำครั้งเดียว)
1. เข้าสู่ระบบ GitHub แล้วเปิด https://github.com/nattaponm/Teach_RQI_Mosaic_SpatialBias_DualPol
2. กดปุ่ม **Fork** (มุมขวาบน) → *Create fork*
3. จะได้ `https://github.com/<ชื่อบัญชีของคุณ>/Teach_RQI_Mosaic_SpatialBias_DualPol`

## ขั้นที่ 2 — สร้าง token ที่ใช้ได้กับ fork เท่านั้น (ทำครั้งเดียว)
1. GitHub ▸ รูปโปรไฟล์ ▸ **Settings** ▸ **Developer settings** ▸ **Personal access tokens** ▸ **Fine-grained tokens** ▸ *Generate new token*
2. Token name: `colab-radteach` · Expiration: วันสิ้นภาคการศึกษา
3. **Repository access**: *Only select repositories* → เลือก fork ของคุณ **อันเดียว**
4. **Permissions ▸ Repository permissions ▸ Contents: Read and write** (อย่างอื่นไม่ต้องเปิด)
5. กด *Generate token* แล้วคัดลอกเก็บไว้ (แสดงครั้งเดียว)

## ขั้นที่ 3 — ใส่ token ใน Colab Secrets (ทำครั้งเดียวต่อบัญชี Google)
1. เปิด notebook ใดก็ได้ใน Colab → แถบซ้ายคลิกไอคอน 🔑 **Secrets**
2. *Add new secret*: `GITHUB_USER` = ชื่อบัญชี GitHub ของคุณ
3. *Add new secret*: `GITHUB_TOKEN` = token จากขั้นที่ 2
4. เปิดสวิตช์ **Notebook access** ของทั้งสองรายการ

> 🔒 ห้ามพิมพ์ token ลงในเซลล์โค้ด ห้ามแชร์ภาพหน้าจอที่มี token ถ้าเผลอเปิดเผย ให้ลบ token ทิ้งที่ GitHub แล้วสร้างใหม่

## ขั้นที่ 4 — บันทึกผลหลังทำ notebook เสร็จ
- เลื่อนไปเซลล์สุดท้าย 💾 **บันทึกผลลัพธ์ขึ้น GitHub** ใส่ `STUDENT_ID` แล้วรัน
- ไฟล์ทั้งหมดใน `outputs/` (ภาพ ตาราง GeoTIFF ≤ 50 MB/ไฟล์) จะไปอยู่ที่
  `https://github.com/<คุณ>/Teach_RQI_Mosaic_SpatialBias_DualPol/tree/results/<รหัสนิสิต>/outputs/<รหัสนิสิต>`
- รันซ้ำได้ ไฟล์จะถูกอัปเดตใน branch เดิม

## ขั้นที่ 5 — บันทึกตัว notebook ที่ทำแบบฝึกหัดแล้ว
*File ▸ Save a copy in GitHub* → Repository: fork ของคุณ · Branch: `results/<รหัสนิสิต>` · File path: `notebooks/<ชื่อไฟล์>`

## ส่งงาน
ส่ง **ลิงก์ branch** `results/<รหัสนิสิต>` ใน fork ของคุณผ่านระบบที่ผู้สอนกำหนด (ไม่ต้องเปิด Pull Request — PR จะถูกปิดอัตโนมัติ)

## ปัญหาที่พบบ่อย
| อาการ | สาเหตุ/วิธีแก้ |
|---|---|
| `Set the Colab secrets GITHUB_USER and GITHUB_TOKEN first` | ยังไม่ได้เพิ่ม secret หรือยังไม่เปิด *Notebook access* |
| `git failed: ... 403` | token ไม่มีสิทธิ์ Contents: Read and write หรือเลือก repo ไม่ถูก |
| `git failed: ... not found` | ยังไม่ได้ fork หรือ `GITHUB_USER` สะกดผิด |
| `GITHUB_USER is the instructor account` | ใส่ชื่อบัญชีของอาจารย์ — ต้องเป็นบัญชีของนิสิตเอง |
