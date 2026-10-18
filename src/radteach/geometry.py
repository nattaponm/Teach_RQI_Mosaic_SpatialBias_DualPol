"""Radar beam geometry and beam blockage.

Beam-centre height with the effective earth radius model (Doviak & Zrnić 1993, eq. 2.28b):
    h(r, θ) = sqrt(r² + (kₑa)² + 2 r kₑa sin θ) − kₑa + h₀,     kₑ = 4/3, a = 6371 km
Half-power beam radius:  w(r) = r · tan(θ_3dB / 2)
Partial beam blockage of a circular beam cross-section (Bech et al. 2003, J. Atmos. Ocean. Tech. 20:845–855):
    y = terrain − h,  a = y / w
    PBB = [a·sqrt(1−a²) + arcsin(a) + π/2] / π        (0 if y < −w, 1 if y > w)
Cumulative beam blockage CBB = running maximum of PBB along each ray (what is blocked stays blocked).
These are the same formulas as wradlib.qual.beam_block_frac / cum_beam_block_frac, written out so that
students can read them.
"""
from __future__ import annotations

import numpy as np
from pyproj import Geod

EARTH_RADIUS = 6371_000.0
GEOD = Geod(ellps="WGS84")


def beam_height(r_m, elev_deg, alt_m=0.0, ke=4 / 3):
    """Height (m above sea level) of the beam centre at slant range r_m."""
    r = np.asarray(r_m, float)
    a = ke * EARTH_RADIUS
    return np.sqrt(r ** 2 + a ** 2 + 2 * r * a * np.sin(np.radians(elev_deg))) - a + alt_m


def ground_distance(r_m, elev_deg, ke=4 / 3):
    """Distance along the earth surface below the beam (m)."""
    r = np.asarray(r_m, float)
    a = ke * EARTH_RADIUS
    h = beam_height(r, elev_deg, 0.0, ke)
    return a * np.arcsin(r * np.cos(np.radians(elev_deg)) / (a + h))


def range_for_height(h_m, elev_deg, alt_m=0.0, rmax=400e3):
    """Slant range at which the beam centre reaches height h_m (m)."""
    r = np.linspace(0, rmax, 40001)
    hh = beam_height(r, elev_deg, alt_m)
    return float(r[np.argmax(hh >= h_m)]) if hh[-1] >= h_m else np.nan


def beam_radius(r_m, beamwidth_deg=1.0):
    """Half-power beam radius (m)."""
    return np.asarray(r_m, float) * np.tan(np.radians(beamwidth_deg / 2))


def partial_blockage(terrain_m, beam_h_m, radius_m):
    """Fraction of the beam cross-section below the terrain (Bech et al. 2003)."""
    y = np.asarray(terrain_m, float) - np.asarray(beam_h_m, float)
    a = np.clip(y / np.maximum(radius_m, 1e-6), -1, 1)
    return (a * np.sqrt(1 - a ** 2) + np.arcsin(a) + np.pi / 2) / np.pi


def cumulative_blockage(pbb):
    """Running maximum along range (last axis)."""
    return np.maximum.accumulate(np.nan_to_num(pbb, nan=0.0), axis=-1)


def polar_lonlat(lon0, lat0, az_deg, r_m, elev_deg=0.0):
    """Longitude/latitude of every (azimuth, range) bin, using the ground distance under the beam."""
    AZ, R = np.meshgrid(np.asarray(az_deg, float), np.asarray(r_m, float), indexing="ij")
    s = ground_distance(R, elev_deg)
    lon, lat, _ = GEOD.fwd(np.full(AZ.shape, lon0), np.full(AZ.shape, lat0), AZ, s)
    return lon, lat


def sweep_geometry(lon0, lat0, alt_m, elev_deg, az_deg, r_m, terrain_fn, beamwidth_deg=1.0):
    """Geometry of one sweep on an (azimuth × range) grid.

    Returns dict with lon, lat, h (beam centre, m asl), terr (m asl), hagl (m above ground),
    radius (m), pbb, cbb, clearance (beam bottom above terrain, m).
    """
    lon, lat = polar_lonlat(lon0, lat0, az_deg, r_m, elev_deg)
    h = beam_height(r_m, elev_deg, alt_m)[None, :] * np.ones((len(az_deg), 1))
    terr = terrain_fn(lat, lon)
    rad = beam_radius(r_m, beamwidth_deg)[None, :]
    pbb = partial_blockage(terr, h, rad)
    cbb = cumulative_blockage(pbb)
    return dict(lon0=lon0, lat0=lat0, alt=alt_m, el=elev_deg, az=np.asarray(az_deg, float), r=np.asarray(r_m, float),
                lon=lon, lat=lat, h=h, terr=terr, hagl=h - terr, radius=np.broadcast_to(rad, h.shape),
                pbb=pbb, cbb=cbb, clearance=h - rad - terr)
