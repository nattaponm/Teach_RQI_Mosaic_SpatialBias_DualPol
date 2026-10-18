"""Journal-style figures: sizes in millimetres, Thai font, colour maps, saving.

Standards used in every notebook
  * width 90 mm (1 column), 140 mm, 190 mm (2 columns); font >= 8 pt at print size
  * PNG at 600 dpi for slides/reports + PDF with embedded TrueType fonts for journals
  * discrete rain levels, colour-blind-safe sequential maps, diverging maps centred on 1 (ratios)
  * every axis has a quantity and a unit; panels are labelled (a), (b), ...
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm
from matplotlib.colors import BoundaryNorm, ListedColormap

MM = 1 / 25.4
WIDTH = {"single": 90, "one_half": 140, "double": 190}
LANG = "th"

# bilingual labels used across notebooks
LABELS = {
    "rain_mm": {"th": "ฝน (มม.)", "en": "Rainfall (mm)"},
    "rain_mmh": {"th": "ความเข้มฝน (มม./ชม.)", "en": "Rain rate (mm h$^{-1}$)"},
    "gauge_mm": {"th": "ฝนสถานีวัด (มม./ชม.)", "en": "Gauge rainfall (mm h$^{-1}$)"},
    "radar_mm": {"th": "ฝนเรดาร์ (มม./ชม.)", "en": "Radar rainfall (mm h$^{-1}$)"},
    "rqi": {"th": "ดัชนีคุณภาพเรดาร์ (RQI)", "en": "Radar quality index (RQI)"},
    "height_km": {"th": "ความสูง (กม.)", "en": "Height (km)"},
    "range_km": {"th": "ระยะจากเรดาร์ (กม.)", "en": "Range (km)"},
    "dbz": {"th": "การสะท้อนกลับ (dBZ)", "en": "Reflectivity (dBZ)"},
    "elev_m": {"th": "ความสูงภูมิประเทศ (ม. รทก.)", "en": "Elevation (m a.s.l.)"},
    "cbb": {"th": "สัดส่วนลำคลื่นถูกบัง (CBB)", "en": "Cumulative beam blockage"},
}


def L(key: str, lang: str | None = None) -> str:
    """Bilingual label: L('rain_mm') -> 'ฝน (มม.)' when LANG == 'th'."""
    return LABELS.get(key, {}).get(lang or LANG, key)


def _font_dir():
    here = Path(__file__).resolve()
    for p in here.parents:
        if (p / "mapping" / "fonts").exists():
            return p / "mapping" / "fonts"
    return None


def setup(lang: str = "th", size: float = 9.0):
    """Apply the course style. Thai text uses the bundled Sarabun (SIL OFL); symbols fall back to DejaVu Sans."""
    global LANG
    LANG = lang
    fd = _font_dir()
    if fd is not None:
        for f in fd.glob("*.ttf"):
            try:
                fm.fontManager.addfont(str(f))
            except Exception:
                pass
    have = {f.name for f in fm.fontManager.ttflist}
    thai = next((n for n in ["Sarabun", "TH Sarabun New", "Leelawadee UI", "Tahoma", "Noto Sans Thai", "Loma"]
                 if n in have), None)
    family = [thai, "DejaVu Sans"] if thai else ["DejaVu Sans"]
    scale = 1.2 if thai and "Sarabun" in thai else 1.0     # Sarabun is optically small
    mpl.rcParams.update({
        "font.family": family,
        "font.size": size * scale,
        "axes.titlesize": size * scale * 1.05,
        "axes.labelsize": size * scale,
        "xtick.labelsize": size * scale * 0.92,
        "ytick.labelsize": size * scale * 0.92,
        "legend.fontsize": size * scale * 0.9,
        "axes.unicode_minus": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": True, "axes.spines.right": True,
        "legend.frameon": False,
        "figure.dpi": 110,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "image.interpolation": "nearest",
    })
    return thai or "DejaVu Sans"


def figure(width="double", height_mm=100.0, **kw):
    """New figure with width 'single'/'one_half'/'double' or a number in mm."""
    w = WIDTH.get(width, width) if isinstance(width, str) else width
    return plt.figure(figsize=(w * MM, height_mm * MM), **kw)


def subplots(nrows=1, ncols=1, width="double", height_mm=100.0, **kw):
    w = WIDTH.get(width, width) if isinstance(width, str) else width
    kw.setdefault("constrained_layout", True)
    return plt.subplots(nrows, ncols, figsize=(w * MM, height_mm * MM), **kw)


def panel_labels(axes, x=-0.02, y=1.02, start="a", **kw):
    axes = np.atleast_1d(axes).ravel()
    for i, ax in enumerate(axes):
        ax.text(x, y, f"({chr(ord(start) + i)})", transform=ax.transAxes, fontweight="bold",
                ha="right", va="bottom", **kw)


def save(fig, name: str, outdir: str | os.PathLike | None = None, formats=("png", "pdf"), dpi=600):
    """Save as outputs/figures/<name>.png and .pdf (fonts embedded). Returns the list of paths."""
    if outdir is None:
        from .config import paths
        outdir = paths().figures
    Path(outdir).mkdir(parents=True, exist_ok=True)
    out = []
    for f in formats:
        p = Path(outdir) / f"{name}.{f}"
        fig.savefig(p, dpi=dpi)
        out.append(p)
    return out


# --------------------------------------------------------------------------- colour maps
def rain_cmap(levels=(0.1, 0.5, 1, 2, 5, 10, 20, 35, 50)):
    """Discrete rain colours (light blue -> dark purple). Values below the first level are transparent."""
    base = ["#c6e9f7", "#7dc8ee", "#3a96d4", "#1f6fb2", "#2ca25f", "#a6d854", "#ffd92f",
            "#fc8d59", "#d7301f", "#7a0177"]
    n = len(levels) - 1
    cols = [base[int(round(i))] for i in np.linspace(0, len(base) - 2, n)]
    cmap = ListedColormap(cols)
    cmap.set_over(base[-1])
    cmap.set_under((1, 1, 1, 0))
    cmap.set_bad((1, 1, 1, 0))
    return cmap, BoundaryNorm(levels, cmap.N)


def rqi_cmap(levels=(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)):
    cmap = plt.get_cmap("cividis", len(levels) - 1).copy()
    cmap.set_bad("0.85")
    return cmap, BoundaryNorm(levels, cmap.N)


def ratio_cmap(vmax_log2=3, n=12):
    """Diverging map for multiplicative ratios R/G on a log2 scale (blue = radar too high)."""
    from matplotlib.colors import TwoSlopeNorm
    cmap = plt.get_cmap("RdBu", n).copy()
    return cmap, TwoSlopeNorm(vmin=-vmax_log2, vcenter=0, vmax=vmax_log2)


def dbz_cmap():
    """Reflectivity colours (NWS-like, 5 dB steps from 5 to 70 dBZ)."""
    lv = np.arange(5, 75, 5)
    cols = ["#04e9e7", "#019ff4", "#0300f4", "#02fd02", "#01c501", "#008e00", "#fdf802",
            "#e5bc00", "#fd9500", "#fd0000", "#d40000", "#bc0000", "#f800fd"]
    cmap = ListedColormap(cols[:len(lv) - 1])
    cmap.set_over(cols[-1]); cmap.set_under((1, 1, 1, 0)); cmap.set_bad((1, 1, 1, 0))
    return cmap, BoundaryNorm(lv, cmap.N)


def one_to_one(ax, lo, hi, factors=(0.5, 2.0)):
    """1:1 line plus factor-of-two guides (dotted)."""
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, zorder=1)
    for f in factors:
        ax.plot([lo, hi], [lo * f, hi * f], ":", color="0.5", lw=0.6, zorder=1)
