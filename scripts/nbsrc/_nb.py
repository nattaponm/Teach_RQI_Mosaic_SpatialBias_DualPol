"""Tiny DSL to write the course notebooks as Python files (easier to review and diff than .ipynb JSON).

    md("text")           markdown cell
    code("python")       code cell
    exercise(n, lvl, t)  exercise statement (lvl 'ug' = undergraduate core, 'pg' = graduate extension)
    todo("python")       student cell (runs without error even if not completed)
    solution(n, "py")    collapsed solution cell (Colab form: code hidden until double-click)
"""
import textwrap

import nbformat as nbf

OWNER = "nattaponm"
REPO = "Teach_RQI_Mosaic_SpatialBias_DualPol"
LEVEL = {"ug": "🟢 ป.ตรี", "pg": "🔵 ป.โท / ส่วนขยาย"}


def _clean(s):
    return textwrap.dedent(s).strip("\n")


def md(s):
    return nbf.v4.new_markdown_cell(_clean(s))


def code(s, hidden=False, title=None):
    src = _clean(s)
    if title:
        src = f"#@title {title}\n" + src
    c = nbf.v4.new_code_cell(src)
    if hidden:
        c.metadata["cellView"] = "form"
        c.metadata["jupyter"] = {"source_hidden": True}
    return c


def exercise(n, level, text):
    return md(f"### แบบฝึกหัด {n} · {LEVEL[level]}\n\n{_clean(text)}")


def todo(s):
    return code("# ✏️ เขียนคำตอบของคุณในเซลล์นี้\n" + _clean(s))


def solution(n, s):
    c = code(f'#@title 🔑 เฉลยข้อ {n} — ลองทำเองก่อน แล้วดับเบิลคลิกเพื่อดูโค้ด {{ display-mode: "form" }}\n' + _clean(s))
    c.metadata["cellView"] = "form"
    c.metadata["jupyter"] = {"source_hidden": True}
    return c


def badge(fname):
    url = f"https://colab.research.google.com/github/{OWNER}/{REPO}/blob/main/notebooks/{fname}"
    return f"[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({url})"


def header(fname, num, title, title_en, level, minutes, prereq, inputs, outputs):
    return md(f"""
    {badge(fname)}

    # {num} · {title}
    *{title_en}*

    | | |
    |---|---|
    | **ระดับ** | {level} |
    | **เวลาโดยประมาณ** | {minutes} |
    | **ควรผ่านมาก่อน** | {prereq} |
    | **ข้อมูลเข้า** | {inputs} |
    | **ผลลัพธ์** | {outputs} |

    > 📌 **วิธีใช้:** กด *Runtime ▸ Run all* ได้ทันที หรือรันทีละเซลล์ (Shift+Enter) พร้อมอ่านคำอธิบาย
    > เซลล์ที่มีหัวข้อ 🔑 เป็นเฉลยที่ซ่อนไว้ ลองทำแบบฝึกหัดก่อนแล้วค่อยเปิดดู
    """)


SETUP = '''
#@title ⚙️ เตรียมสภาพแวดล้อม (รันเซลล์นี้ก่อนทุกครั้ง) { display-mode: "form" }
import os, sys, subprocess
REPO, OWNER = "Teach_RQI_Mosaic_SpatialBias_DualPol", "nattaponm"
if "google.colab" in sys.modules:
    if not os.path.exists(f"/content/{REPO}"):
        subprocess.run(["git", "clone", "-q", "--depth", "1",
                        f"https://github.com/{OWNER}/{REPO}.git", f"/content/{REPO}"], check=True)
    os.chdir(f"/content/{REPO}")
else:  # running on your own computer: move to the repository root
    while not os.path.exists("config/study_area.yaml") and os.getcwd() != os.path.dirname(os.getcwd()):
        os.chdir("..")
sys.path.insert(0, os.path.abspath("src"))
from radteach import colab
P = colab.setup(packages=PACKAGES)
'''


def setup_cell(packages=("teaching_data", "products")):
    return code(f"PACKAGES = {tuple(packages)!r}\n" + _clean(SETUP), hidden=True)


SAVE = '''
#@title 💾 บันทึกผลลัพธ์ขึ้น GitHub (fork ของคุณเอง) { display-mode: "form" }
#@markdown ใส่รหัสนิสิต แล้วรันเซลล์นี้ (ต้องตั้ง Colab Secrets `GITHUB_USER` และ `GITHUB_TOKEN` ก่อน — ดูขั้นตอนใน notebook 00)
STUDENT_ID = ""  #@param {type:"string"}
if STUDENT_ID:
    colab.save_to_fork(STUDENT_ID, message="NB{num}: results")
else:
    print("ยังไม่ได้ใส่รหัสนิสิต — ข้ามการบันทึก (ผลยังอยู่ในโฟลเดอร์ outputs/ ของ Colab)")
'''


def save_cell(num):
    c = code(_clean(SAVE).replace("{num}", num), hidden=True)
    return c


def build(cells, path):
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
        "colab": {"provenance": [], "toc_visible": True},
    }
    nbf.write(nb, path)
