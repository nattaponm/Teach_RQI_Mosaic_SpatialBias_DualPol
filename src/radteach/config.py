"""Paths and YAML configuration shared by all notebooks."""
from __future__ import annotations

import os
from pathlib import Path

import yaml


def find_root(start: str | os.PathLike | None = None) -> Path:
    """Return the repository root (the folder that contains config/ and src/)."""
    p = Path(start or os.getcwd()).resolve()
    for q in [p, *p.parents]:
        if (q / "config" / "study_area.yaml").exists() and (q / "src" / "radteach").exists():
            return q
    env = os.environ.get("RADTEACH_ROOT")
    if env:
        return Path(env)
    raise FileNotFoundError("repository root not found — run notebook 00 or set RADTEACH_ROOT")


ROOT = None


class Paths:
    """All folders used by the course, relative to the repository root."""

    def __init__(self, root: Path):
        self.root = Path(root)
        self.config = self.root / "config"
        self.static = self.root / "data" / "static"
        self.raw = self.root / "data" / "raw"
        self.interim = self.root / "data" / "interim"
        self.processed = self.root / "data" / "processed"
        self.radar = self.raw / "teaching_data_v1" / "radar"
        self.dem = self.raw / "teaching_data_v1" / "dem"
        self.products = self.raw / "teaching_products_v1"
        self.gauges = self.raw / "gauges"
        self.outputs = self.root / "outputs"
        self.figures = self.outputs / "figures"
        self.tables = self.outputs / "tables"
        self.rasters = self.outputs / "rasters"
        self.fonts = self.root / "mapping" / "fonts"
        self.boundaries = self.static / "boundaries"
        self.stations = self.static / "stations"

    def make(self):
        for p in [self.raw, self.interim, self.processed, self.gauges, self.figures, self.tables, self.rasters]:
            p.mkdir(parents=True, exist_ok=True)
        return self

    def __repr__(self):
        return f"Paths(root={self.root})"


def paths(root=None) -> Paths:
    return Paths(root or find_root()).make()


def load(name: str, root=None) -> dict:
    """Load config/<name>.yaml as a dict, e.g. load('qpe')."""
    root = Path(root or find_root())
    with open(root / "config" / f"{name}.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_all(root=None) -> dict:
    root = Path(root or find_root())
    return {p.stem: load(p.stem, root) for p in sorted((root / "config").glob("*.yaml"))}
