"""Build notebooks/*.ipynb from the sources in scripts/nbsrc/nbXX.py (instructor tool).

    python scripts/build_notebooks.py            # all notebooks
    python scripts/build_notebooks.py 03 08      # selected ones
The generated notebooks contain no outputs (solution outputs would reveal the answers).
"""
import os
import runpy
import shutil
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
src = root / "scripts" / "nbsrc"
sel = sys.argv[1:] or sorted(p.stem[2:] for p in src.glob("nb[0-9][0-9].py"))
os.chdir(src)
sys.path.insert(0, str(src))
(src / "notebooks").mkdir(exist_ok=True)
for n in sel:
    runpy.run_path(str(src / f"nb{n}.py"))
    for f in (src / "notebooks").glob("*.ipynb"):
        shutil.move(str(f), root / "notebooks" / f.name)
        print("built", f.name)
