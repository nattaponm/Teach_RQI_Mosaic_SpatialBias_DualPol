"""Execute notebooks headlessly (CI and instructor checks). Outputs are written to a scratch copy, not the repo.

    python scripts/execute_notebooks.py                 # all
    python scripts/execute_notebooks.py 00 07           # selected
Environment: RADTEACH_DATA_URL=file:///path/to/zips/  to use local release zips.
"""
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

root = Path(__file__).resolve().parents[1]
sel = sys.argv[1:]
nbs = sorted((root / "notebooks").glob("*.ipynb"))
if sel:
    nbs = [p for p in nbs if p.name[:2] in sel]
failed = []
for p in nbs:
    nb = nbformat.read(p, 4)
    t = time.time()
    try:
        NotebookClient(nb, timeout=1800, kernel_name="python3", resources={"metadata": {"path": str(root / "notebooks")}}).execute()
        print(f"OK    {p.name}  {time.time() - t:.0f} s")
    except Exception as e:
        failed.append(p.name)
        print(f"FAIL  {p.name}: {str(e)[-800:]}")
sys.exit(1 if failed else 0)
