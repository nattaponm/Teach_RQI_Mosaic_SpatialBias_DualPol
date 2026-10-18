"""Radar Quality Index (after Zhang et al. 2012, IAHS Publ. 351:388–393; MRMS).

    RQI = RQI_blk × RQI_hgt (× RQI_pol for phase-based estimates)

    RQI_blk = 1                                   if CBB ≤ c₁ (0.1)
            = (c₂ − CBB) / (c₂ − c₁)              if c₁ < CBB < c₂
            = 0                                   if CBB ≥ c₂ (0.5 for Z; 0.7 for KDP/A)
    RQI_hgt = exp(−[max(0, h_AGL − h_ref) / h_s]²)   × exp(−(h − ML_bottom)/1 km) above the melting layer
    RQI_pol = clip((30° − texture(ΦDP)) / 20°, 0, 1)
"""
from __future__ import annotations

import numpy as np


def rqi_blockage(cbb, c1=0.1, c2=0.5):
    return np.clip((c2 - np.asarray(cbb, float)) / (c2 - c1), 0, 1)


def rqi_height(hagl_m, h_ref=1000.0, h_scale=1500.0, h_asl_m=None, ml_bottom_m=None, ml_decay_m=1000.0):
    q = np.exp(-(np.maximum(0, np.asarray(hagl_m, float) - h_ref) / h_scale) ** 2)
    if ml_bottom_m is not None and h_asl_m is not None:
        q = q * np.where(h_asl_m > ml_bottom_m, np.exp(-(h_asl_m - ml_bottom_m) / ml_decay_m), 1.0)
    return q


def rqi_polarimetric(texture_deg, good=10.0, bad=30.0):
    return np.clip((bad - np.asarray(texture_deg, float)) / (bad - good), 0, 1)


def rqi(cbb, hagl_m, c1=0.1, c2=0.5, h_ref=1000.0, h_scale=1500.0, h_asl_m=None, ml_bottom_m=None):
    return rqi_blockage(cbb, c1, c2) * rqi_height(hagl_m, h_ref, h_scale, h_asl_m, ml_bottom_m)
