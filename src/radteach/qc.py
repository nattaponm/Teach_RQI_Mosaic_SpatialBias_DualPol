"""Quality control and polarimetric processing of one sweep.

Steps (process_sweep):
  1. no-echo mask: Z < 5 dBZ or missing  → rain = 0 mm (known dry)
  2. meteorological mask: ρhv ≥ 0.85, range > 5 km, speckle filter (≥5 of 9 neighbours), clutter test
     (|V| < 0.3 m/s and W < 0.6 m/s and Z > 35 dBZ = stationary target)
  3. ZDR offset removed (estimated from light rain, see calibrate)
  4. ΦDP: subtract system phase, apply sign convention, keep only smooth gates (texture < 20°)
  5. KDP (Vulpiani et al. 2012, Py-ART) = ½ dΦDP/dr over a ~6 km window
  6. attenuation correction ZPHI (Testud et al. 2000, Py-ART): specific attenuation A and corrected Z
  7. rain rates R(Z) MP, R(Z) Rosenfeld, R(KDP), R(A), blended
"""
from __future__ import annotations

import contextlib
import io as _io
import warnings

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter, uniform_filter1d

from . import qpe

warnings.filterwarnings("ignore")
QC_DEFAULT = {"min_dbz": 5.0, "min_rhohv": 0.85, "min_range_km": 5.0, "speckle_fraction": 5 / 9,
              "clutter": {"max_abs_velocity": 0.3, "max_spectrum_width": 0.6, "min_dbz": 35},
              "phase_texture_max_deg": 20}
F = {"Z": "reflectivity", "ZDR": "corrected_differential_reflectivity", "RHO": "cross_correlation_ratio",
     "PHI": "differential_phase", "V": "velocity", "W": "spectrum_width"}


def _pyart():
    with contextlib.redirect_stdout(_io.StringIO()):
        import pyart
    return pyart


def arr(radar_or_sweep, key):
    return np.ma.filled(radar_or_sweep.fields[F.get(key, key)]["data"].astype(float), np.nan)


def phase_texture(phi_deg, n=5):
    """Circular standard deviation (°) of ΦDP over n gates — small in rain, large in noise/clutter."""
    c = np.exp(1j * np.radians(phi_deg))
    c = np.where(np.isfinite(phi_deg), c, np.nan)
    re_ = uniform_filter1d(np.nan_to_num(c.real), n, axis=1)
    im_ = uniform_filter1d(np.nan_to_num(c.imag), n, axis=1)
    R = np.clip(np.hypot(re_, im_), 1e-6, 1)
    return np.degrees(np.sqrt(-2 * np.log(R)))


def calibrate(radars: dict, light_rain_dbz=(20, 28), min_rhohv=0.97, zdr_expected=0.25) -> pd.DataFrame:
    """Per-radar calibration constants from the lowest sweep of many scans.

    radars: {name: [pyart radar, ...]}
    zdr_off      median ZDR in light rain − expected intrinsic value (0.25 dB)
    sysphase     circular mean of ΦDP over the first rain gates of each ray (system differential phase)
    sign         +1 or −1: direction in which ΦDP accumulates with range (vendor convention)
    """
    rows = {}
    for name, lst in radars.items():
        zd, sp, slopes = [], [], []
        for r in lst:
            s = r.extract_sweeps([0])
            if F["PHI"] not in s.fields:
                continue
            Z, RH, D, P = arr(s, "Z"), arr(s, "RHO"), arr(s, "ZDR"), arr(s, "PHI")
            zd.append(D[(Z >= light_rain_dbz[0]) & (Z <= light_rain_dbz[1]) & (RH > min_rhohv)])
            rng = s.range["data"]
            for i in range(P.shape[0]):
                m = (Z[i] > 20) & (RH[i] > 0.95) & (rng > 5e3)
                if m.sum() >= 5:
                    sp.append(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(P[i][m][:5]))))))
                if m.sum() >= 30:
                    ph = np.unwrap(np.radians(P[i][m]))
                    slopes.append(np.polyfit(rng[m] / 1e3, np.degrees(ph), 1)[0])
        zz = np.concatenate(zd) if zd else np.array([])
        med = float(np.nanmedian(zz)) if zz.size > 50 else np.nan
        rows[name] = dict(zdr_light_rain_db=med, zdr_offset_db=(med - zdr_expected) if np.isfinite(med) else 0.0,
                          n_zdr_gates=int(zz.size),
                          system_phase_deg=float(np.degrees(np.angle(np.mean(np.exp(1j * np.radians(sp)))))) if sp else 0.0,
                          n_phase_rays=len(sp),
                          phidp_slope_deg_km=float(np.median(slopes)) if slopes else np.nan,
                          phidp_sign=(-1.0 if (slopes and np.median(slopes) < 0) else 1.0))
    return pd.DataFrame(rows).T


def process_sweep(radar, sweep=0, zdr_offset=0.0, system_phase=0.0, phidp_sign=1.0, fzl_m=4800.0,
                  qc=None, coef=None, blend_cfg=None, do_attenuation=True):
    """Quality-control and rain rates for one sweep. Returns (pyart sweep object, dict of 2-D arrays)."""
    pyart = _pyart()
    qc = {**QC_DEFAULT, **(qc or {})}
    s = radar.extract_sweeps([sweep])
    s.instrument_parameters = s.instrument_parameters or {}
    s.instrument_parameters["frequency"] = {"data": np.array([5.6e9])}
    Zm = s.fields[F["Z"]]["data"]
    Z = np.ma.filled(Zm.astype(float), np.nan)
    noecho = np.ma.getmaskarray(Zm) | ~np.isfinite(Z) | (Z < qc["min_dbz"])
    rng = s.range["data"].astype(float)
    nan = np.full_like(Z, np.nan)
    out = dict(Z=Z, noecho=noecho, range=rng, az=s.azimuth["data"].astype(float),
               elev=float(s.fixed_angle["data"][0]))
    dual = all(F[k] in s.fields for k in ["PHI", "RHO", "ZDR"])
    span = np.zeros(Z.shape[0])
    K, A, Zc, tex, RH, D = nan, nan, Z.copy(), nan, nan, nan
    if dual:
        RH, P = arr(s, "RHO"), arr(s, "PHI")
        D = arr(s, "ZDR") - zdr_offset
        tex = phase_texture(np.where(noecho, np.nan, P))
        echo = (~noecho) & (RH >= qc["min_rhohv"])
        met = echo & (rng[None, :] > qc["min_range_km"] * 1e3)
        met &= uniform_filter(echo.astype(float), size=3, mode="wrap") >= qc["speckle_fraction"]
        if F["V"] in s.fields and F["W"] in s.fields:
            c = qc["clutter"]
            V, W = arr(s, "V"), arr(s, "W")
            met &= ~((Z > c["min_dbz"]) & (np.abs(V) < c["max_abs_velocity"]) & (W < c["max_spectrum_width"]))
        phase_ok = met & (tex < qc["phase_texture_max_deg"])
        Pc = (phidp_sign * (P - system_phase) + 420) % 360 - 60
        s.add_field_like(F["PHI"], "PSIDP", np.ma.masked_where(~phase_ok, Pc), replace_existing=True)
        s.add_field_like(F["ZDR"], "ZDRC", np.ma.masked_invalid(D), replace_existing=True)
        gf = pyart.filters.GateFilter(s); gf.exclude_masked("PSIDP")
        gate = float(np.median(np.diff(rng)))
        win = int(max(4, round(6000 / gate / 2) * 2))
        with contextlib.redirect_stdout(_io.StringIO()):
            kdp, phidp = pyart.retrieve.kdp_vulpiani(s, gatefilter=gf, psidp_field="PSIDP", windsize=win,
                                                     band="C", n_iter=10)
        s.add_field("KDP", kdp, replace_existing=True)
        s.add_field("PHIDP_F", phidp, replace_existing=True)
        K = np.ma.filled(kdp["data"].astype(float), np.nan)
        PF = np.where(phase_ok, np.ma.filled(phidp["data"].astype(float), np.nan), np.nan)
        with np.errstate(all="ignore"):
            span = np.nan_to_num(np.nanmax(PF, axis=1) - np.nanmin(PF, axis=1)) if np.isfinite(PF).any() else span
        if do_attenuation:
            try:
                with contextlib.redirect_stdout(_io.StringIO()):
                    spec_at, _, cor_z, *_ = pyart.correct.calculate_attenuation_zphi(
                        s, doc=0, fzl=fzl_m, gatefilter=gf, refl_field=F["Z"], phidp_field="PHIDP_F",
                        zdr_field="ZDRC", temp_ref="fixed_fzl")
                A = np.ma.filled(spec_at["data"].astype(float), np.nan)
                Zc = np.ma.filled(cor_z["data"].astype(float), np.nan)
                Zc = np.where(np.isfinite(Zc), Zc, Z)
            except Exception as exc:  # ZPHI can fail on rays without enough phase
                out["zphi_error"] = str(exc)
        K = np.where(phase_ok, K, np.nan); A = np.where(phase_ok, A, np.nan)
        out.update(PHIDP=Pc, PHIDP_F=PF, phase_ok=phase_ok)
    else:
        met = ~noecho
    coef = coef or qpe.DEFAULTS
    Zc = np.clip(Zc, -10, 55)
    R_MP = qpe.r_from_z(Zc, **coef["zr_marshall_palmer"])
    R_ROS = qpe.r_from_z(Zc, **coef["zr_rosenfeld"])
    R_KDP = qpe.r_from_kdp(K, **coef["r_kdp"])
    R_A = qpe.r_from_a(A, **coef["r_a"])
    R_DP, src = qpe.blended(Zc, K, A, span[:, None], blend_cfg, coef)
    if not dual:
        R_KDP = R_A = nan; R_DP = R_MP.copy(); src = np.zeros_like(Z, dtype=int)
    for k, v in dict(R_MP=R_MP, R_ROS=R_ROS, R_KDP=R_KDP, R_A=R_A, R_DP=R_DP).items():
        v = np.clip(np.where(met, v, np.nan), 0, 200)
        v[noecho] = 0.0
        out[k] = v
    out.update(met=met, dual=dual, Zc=np.where(met, Zc, np.nan), KDP=K, A=A, ZDR=D, RHOHV=RH, tex=tex,
               src=np.where(met, src, -1), phidp_span=span,
               frac_kdp=float(np.mean(src[met] == 2)) if met.any() else np.nan,
               frac_a=float(np.mean(src[met] == 1)) if met.any() else np.nan)
    return s, out
