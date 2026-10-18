import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from radteach import bias, geometry, gauges, qpe, rqi, verify  # noqa: E402


def test_beam_height():
    assert abs(geometry.beam_height(100e3, 0.5) / 1e3 - 1.46) < 0.05
    assert geometry.range_for_height(4800, 1.59, 328) > 100e3


def test_partial_blockage_limits():
    w = np.array([100.0])
    assert geometry.partial_blockage(np.array([0.0]), np.array([1000.0]), w)[0] == 0
    assert geometry.partial_blockage(np.array([2000.0]), np.array([1000.0]), w)[0] == 1
    assert abs(geometry.partial_blockage(np.array([1000.0]), np.array([1000.0]), w)[0] - 0.5) < 1e-9


def test_zr_roundtrip():
    r = qpe.r_from_z(qpe.z_from_r(10.0))
    assert abs(r - 10.0) < 1e-6


def test_rqi_bounds():
    assert rqi.rqi_blockage(0.05) == 1 and rqi.rqi_blockage(0.6) == 0
    assert rqi.rqi_height(500) == 1


def test_metrics_perfect():
    g = np.array([0, 1, 2, 5.0]); m = verify.metrics(g, g)
    assert m["PBIAS"] == 0 and abs(m["r"] - 1) < 1e-12 and m["CSI"] == 1


def test_mfb():
    F, n = bias.mean_field_factor(np.array([2, 4.0]), np.array([1, 2.0]))
    assert F == 2 and n == 2


def test_gauge_time_convention():
    df = pd.DataFrame({"station_code": ["X"], "time_lt": [pd.Timestamp("2024-10-03 08:00")], "rain": [1.0], "quality_flag": ["N"]})
    d = gauges.add_utc(df)
    assert d.start_utc.iloc[0] == pd.Timestamp("2024-10-03 00:00")


def test_combine_methods():
    from radteach import mosaic
    R = np.array([[[1.0]], [[3.0]], [[np.nan]]])
    Q = np.array([[[0.2]], [[0.6]], [[0.9]]])
    D = np.array([[[50.0]], [[10.0]], [[5.0]]])
    assert mosaic.combine(R, method="MAX")[0, 0] == 3
    assert mosaic.combine(R, method="MEAN")[0, 0] == 2
    assert mosaic.combine(R, dist=D, method="NEAREST")[0, 0] == 3       # radar 3 is closer but does not see the cell
    assert mosaic.combine(R, Q, method="RQI_MAX")[0, 0] == 3
    assert abs(mosaic.combine(R, Q, method="RQI_WMEAN")[0, 0] - (0.2 * 1 + 0.6 * 3) / 0.8) < 1e-12


def test_calib_offsets_and_loop():
    import pandas as pd
    from radteach import calib
    true = {"A": 1.0, "B": -0.5, "C": -0.5, "D": 0.0}
    rows = [dict(radar_i=a, radar_j=b, n=100, median=true[a] - true[b], used=True) for a, b in [("A", "B"), ("B", "C"), ("A", "C")]]
    st = pd.DataFrame(rows)
    off = calib.solve_offsets(st, list(true))
    assert abs(off["A"] - off["B"] - 1.5) < 1e-6 and abs(off[["A", "B", "C"]].sum()) < 1e-3 and np.isnan(off["D"])
    assert abs(calib.loop_closure(st).residual.iloc[0]) < 1e-9
    assert abs(calib.rain_factor(calib.factor_to_db(2.0, 1.2), 1.2) - 2.0) < 1e-9


def test_per_radar_mfb_loo():
    from radteach import pipeline
    G = np.array([1.0, 2.0, 3.0]); R = np.array([[0.5, np.nan], [1.0, 1.0], [1.5, 3.0]])
    F = pipeline.per_radar_mfb(G, R, loo=False)
    assert np.allclose(F, [2.0, 5 / 4])
    Fl = pipeline.per_radar_mfb(G, R, loo=True)
    assert np.isclose(Fl[0, 0], 5 / 2.5) and np.isclose(Fl[0, 1], 5 / 4)
