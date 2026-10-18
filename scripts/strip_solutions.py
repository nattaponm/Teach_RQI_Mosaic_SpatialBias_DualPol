"""Create notebooks_no_solutions/ with every 🔑 solution cell removed (for exams). Do not commit the output.

    python scripts/strip_solutions.py
"""
from pathlib import Path

import nbformat

root = Path(__file__).resolve().parents[1]
out = root / "notebooks_no_solutions"
out.mkdir(exist_ok=True)
for p in sorted((root / "notebooks").glob("*.ipynb")):
    nb = nbformat.read(p, 4)
    keep = [c for c in nb.cells if not (c.cell_type == "code" and c.source.lstrip().startswith("#@title 🔑"))]
    print(f"{p.name}: removed {len(nb.cells) - len(keep)} solution cells")
    nb.cells = keep
    nbformat.write(nb, out / p.name)
