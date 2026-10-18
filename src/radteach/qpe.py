"""Rain-rate estimators (C band).

    Z–R power law          Z = a Rᵇ            → R = (Z/a)^(1/b),   Z in mm⁶ m⁻³ = 10^(dBZ/10)
    Specific diff. phase   R = c · KDPᵈ        (KDP in ° km⁻¹)
    Specific attenuation   R = c · Aᵈ          (A in dB km⁻¹)
    Blended (decision tree, after Ryzhkov et al. 2014 / Giangrande & Ryzhkov 2008 idea):
        R(KDP)  if Z ≥ 40 dBZ and KDP ≥ 0.3 ° km⁻¹   (heavy rain: phase is robust, immune to calibration/blockage)
        R(A)    if ΔΦDP along the ray ≥ 5° and Z ≥ 25 dBZ (moderate rain with enough phase signal)
        R(Z)    otherwise (light rain)
"""
from __future__ import annotations

import numpy as np

DEFAULTS = {
    "zr_marshall_palmer": {"a": 200.0, "b": 1.6},
    "zr_rosenfeld": {"a": 250.0, "b": 1.2},
    "r_kdp": {"c": 29.7, "d": 0.85},
    "r_a": {"c": 250.0, "d": 0.91},
}


def dbz_to_z(dbz):
    return 10.0 ** (np.asarray(dbz, float) / 10.0)


def z_to_dbz(z):
    return 10.0 * np.log10(np.asarray(z, float))


def r_from_z(dbz, a=200.0, b=1.6, cap_dbz=55.0):
    """R (mm/h) from reflectivity (dBZ); dBZ capped at cap_dbz to limit hail contamination."""
    dbz = np.clip(np.asarray(dbz, float), -10, cap_dbz)
    return (dbz_to_z(dbz) / a) ** (1.0 / b)


def z_from_r(r, a=200.0, b=1.6):
    """dBZ that corresponds to a rain rate (inverse Z–R)."""
    return z_to_dbz(a * np.asarray(r, float) ** b)


def r_from_kdp(kdp, c=29.7, d=0.85):
    kdp = np.asarray(kdp, float)
    return np.where(kdp > 0, c * np.abs(kdp) ** d, 0.0)


def r_from_a(a_spec, c=250.0, d=0.91):
    a_spec = np.asarray(a_spec, float)
    return np.where(a_spec > 0, c * a_spec ** d, 0.0)


def blended(dbz, kdp, a_spec, phidp_span, cfg=None, coef=None):
    """Dual-pol blended rain rate. Returns (R, source) with source 0 = R(Z), 1 = R(A), 2 = R(KDP)."""
    cfg = cfg or {"kdp_min_dbz": 40, "kdp_min": 0.3, "a_min_dbz": 25, "a_min_phidp_span": 5.0}
    coef = coef or DEFAULTS
    rz = r_from_z(dbz, **coef["zr_marshall_palmer"])
    rk = r_from_kdp(kdp, **coef["r_kdp"])
    ra = r_from_a(a_spec, **coef["r_a"])
    span = np.broadcast_to(np.asarray(phidp_span, float), np.shape(dbz))
    use_k = (dbz >= cfg["kdp_min_dbz"]) & (np.nan_to_num(kdp) >= cfg["kdp_min"])
    use_a = (np.nan_to_num(a_spec) > 0) & ~use_k & (span >= cfg["a_min_phidp_span"]) & (dbz >= cfg["a_min_dbz"])
    src = np.where(use_k, 2, np.where(use_a, 1, 0))
    return np.where(use_k, rk, np.where(use_a, ra, rz)), src
