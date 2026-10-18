# สรุปสมการหลัก

| บท | สมการ | ความหมาย |
|---|---|---|
| 01 | $\text{dBZ} = 10\log_{10}Z$, $Z=\int N(D)D^6dD$ | ค่าการสะท้อนกลับ |
| 01 | $\text{ZDR}=10\log_{10}(Z_h/Z_v)$ | ความแบนของเม็ดฝน |
| 01 | $K_{DP}=\tfrac12\,d\Phi_{DP}/dr$ | อัตราเปลี่ยนเฟส |
| 01 | $\text{PIA}\approx\alpha\,\Delta\Phi_{DP}$ | การลดทอนรวมตามทาง (ZPHI) |
| 02 | $h=\sqrt{r^2+(k_ea)^2+2rk_ea\sin\theta}-k_ea+h_0$ | ความสูงลำคลื่น |
| 02 | $\text{PBB}=\dfrac{a\sqrt{1-a^2}+\arcsin a+\pi/2}{\pi}$, $a=(z_{terr}-h)/w$ | การบังลำคลื่นบางส่วน |
| 02 | $\text{RQI}=\text{RQI}_{blk}\cdot\text{RQI}_{hgt}$ | ดัชนีคุณภาพ |
| 03 | $Z=aR^b$, $R=cK_{DP}^d$, $R=cA^d$ | ตัวประมาณฝน |
| 04 | $R_{mos}=\sum_i \text{RQI}_iR_i/\sum_i\text{RQI}_i$ | โมเสคถ่วงน้ำหนัก |
| 04 | $\Delta Z_{blk}=-10\log_{10}(1-\text{CBB})$ | แก้การบังบางส่วน |
| 05 | $\Delta Z(h)=Z_{ref}-Z_{VPR}(h)$ | แก้ VPR |
| 05 | $R(\mathbf x,\alpha)=(1-\alpha)R_0(\mathbf x-\alpha\mathbf V)+\alpha R_1(\mathbf x+(1-\alpha)\mathbf V)$ | advection |
| 06 | $t_{UTC}=t_{LT}-7$ ชม. (ปลายชั่วโมง) | จัดเวลา |
| 07 | PBIAS, RMSE, NSE, KGE, POD, FAR, CSI, HSS | ดัชนีตรวจสอบ |
| 08 | $F=\sum G/\sum R$; $b_i=\ln\frac{G_i+\varepsilon}{R_i+\varepsilon}$; $b(x)=\sum w_ib_i/\sum w_i$ | ปรับแก้ bias |
| 08 | $I=\dfrac{n}{\sum w_{ij}}\dfrac{\sum_{ij}w_{ij}z_iz_j}{\sum z_i^2}$ | Moran's I |
| 09 | MAX: $\max_i R_i$ · MEAN: $\tfrac1n\sum_i R_i$ · NEAREST: $R_{\arg\min_i d_i}$ | โมเสคแบบไม่ใช้คุณภาพ |
| 09 | RQI_MAX: $R_{\arg\max_i \text{RQI}_i}$ · RQI_WMEAN: $\sum_i \text{RQI}_i^kR_i/\sum_i\text{RQI}_i^k$ | โมเสคตามคุณภาพ |
| 09 | Seam index $=\operatorname{med}\lvert\Delta R\rvert_{\text{แนวเปลี่ยนเรดาร์}}/\operatorname{med}\lvert\Delta R\rvert_{\text{อื่น ๆ}}$ | ความต่อเนื่องเชิงพื้นที่ |
| 10 | $\Delta Z_{ij}=c_i-c_j$, $\min\sum n_{ij}(\Delta Z_{ij}-c_i+c_j)^2$ s.t. $\sum c=0$ | ปรับเทียบ Z ระหว่างเรดาร์ (Gourley et al. 2003) |
| 10 | $\Delta Z_{ab}+\Delta Z_{bc}+\Delta Z_{ca}\approx0$ | loop closure |
| 10 | $R'=R\,10^{-c/(10b)}$ | ผลของ offset ต่อฝน Z–R |
| 10 | $F_i=\overline{G}/\overline{R_i}$ (รายเรดาร์) vs $F=\sum G/\sum R$ (ทั้งโดเมน) | MFB |
