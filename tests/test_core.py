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
