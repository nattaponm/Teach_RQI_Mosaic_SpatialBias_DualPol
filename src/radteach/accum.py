"""From 15-min scans to accumulations: scan-hold vs advection correction.

Scan-hold (operational practice):  D = Σₖ Rₖ · Δt   (each scan's rate held for 15 min)
Advection correction (Anagnostou & Krajewski 1999; pysteps example `advection_correction.py`):
    motion V from optical flow (Lucas–Kanade) between consecutive scans, then for a time fraction α ∈ [0,1]
    R(x, α) = (1 − α) · R₀(x − αV) + α · R₁(x + (1 − α)V)
    D = ∫ R dt  (trapezoidal, sampling every dt minutes)
Optical flow: pysteps (LK) if installed, else OpenCV Farnebäck, else zero motion.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates


def motion_field(r0, r1):
    """Displacement (2, ny, nx) in pixels per interval, from rain rates r0 → r1."""
    R = np.array([np.nan_to_num(r0), np.nan_to_num(r1)])
    try:
        from pysteps import motion
        return motion.get_method("LK")(np.log(np.maximum(R, 0.05)), fd_kwargs={"buffer_mask": 10}), "pysteps-LK"
    except Exception:
        pass
    try:
        import cv2
        def to8(a):
            return np.clip(np.log(np.maximum(a, 0.05)) - np.log(0.05), 0, None) / 8.0 * 255
        flow = cv2.calcOpticalFlowFarneback(to8(R[0]).astype(np.uint8), to8(R[1]).astype(np.uint8), None,
                                            0.5, 4, 25, 5, 7, 1.5, 0)
        return np.array([flow[..., 0], flow[..., 1]]), "opencv-farneback"
    except Exception:
        return np.zeros((2,) + R.shape[1:]), "none"


def interp(r0, r1, V, alpha, x, y):
    if alpha <= 0:
        return r0
    if alpha >= 1:
        return r1
    p1 = (y - alpha * V[1], x - alpha * V[0])
    p2 = (y + (1 - alpha) * V[1], x + (1 - alpha) * V[0])
    return (1 - alpha) * map_coordinates(r0, p1, order=1) + alpha * map_coordinates(r1, p2, order=1)


def accumulate(rates, T=15.0, dt=5.0, advect=True):
    """rates: (n, ny, nx) mm/h at 0, T, 2T, … minutes. Returns (depth mm over (n−1)·T, list of motion fields, method)."""
    rates = np.nan_to_num(np.asarray(rates, float))
    ny, nx = rates.shape[1:]
    x, y = np.meshgrid(np.arange(nx, dtype=float), np.arange(ny, dtype=float))
    tot = (len(rates) - 1) * T
    ts = np.arange(0, tot + 1e-9, dt)
    if ts[-1] < tot:
        ts = np.append(ts, tot)
    Vs, how = [], "eulerian"
    for k in range(len(rates) - 1):
        if advect:
            v, how = motion_field(rates[k], rates[k + 1])
        else:
            v = np.zeros((2, ny, nx))
        Vs.append(v)
    samples = np.array([interp(rates[min(int(t // T), len(rates) - 2)], rates[min(int(t // T), len(rates) - 2) + 1],
                               Vs[min(int(t // T), len(rates) - 2)],
                               (t - min(int(t // T), len(rates) - 2) * T) / T, x, y) for t in ts])
    trapz = getattr(np, "trapezoid", None) or np.trapz
    return trapz(samples, ts, axis=0) / 60.0, Vs, how


def scan_hold(m15):
    """Hourly depth = sum of the first four 15-min depths (needs all four)."""
    return np.sum(np.asarray(m15)[:4], axis=0)
