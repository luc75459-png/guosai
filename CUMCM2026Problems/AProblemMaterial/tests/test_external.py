"""L1 外部实测对照（validation）的单元测试。"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from aproblem.external import (
    DATA_RELATIVE,
    CaseSettings,
    CaseTrajectory,
    compute_metrics,
    fit_diffusivity,
    load_conditions,
    load_measured,
    newman_slab_mean_factor,
    run_case,
    sample_at,
)


DATA_DIR = Path(__file__).resolve().parents[1] / DATA_RELATIVE


def test_conditions_and_measured_data_are_consistent() -> None:
    """数据契约：两工况各 9 个采样点，且全部通过数字化质量门。"""

    conditions = load_conditions(DATA_DIR / "conditions.json")
    assert set(conditions["cases"]) == {"60C", "70C"}
    expected_times = np.asarray(conditions["sampling_times_min"], dtype=float)
    for key in ("60C", "70C"):
        measured = load_measured(DATA_DIR / f"figure4_{key}.csv")
        assert measured["time_min"].shape == (9,)
        assert np.allclose(measured["time_min"], expected_times)
        assert measured["included"].all()
        assert np.all(measured["x_db"] > 0.0)
        assert np.all(measured["sigma"] > 0.0)


def test_metrics_match_hand_computed_values() -> None:
    """RMSE / R² / 偏差的手算校验。"""

    measured = np.array([1.0, 2.0, 3.0])
    predicted = np.array([1.5, 2.0, 2.5])
    metrics = compute_metrics(measured, predicted)
    assert np.isclose(metrics.rmse, np.sqrt(1.0 / 6.0))
    assert np.isclose(metrics.bias, 0.0)
    assert np.isclose(metrics.max_abs, 0.5)
    assert np.isclose(metrics.r2, 1.0 - 0.5 / 2.0)


def _short_settings(diffusivity: float = 6.5e-10) -> CaseSettings:
    return CaseSettings(
        key="test",
        temperature_c=60.0,
        ambient_moisture=0.0,
        initial_moisture=3.8,
        initial_temperature_c=25.0,
        radius_m=0.00551,
        diffusivity_m2_s=diffusivity,
        mass_transfer_coefficient_m_s=0.008,
        heat_transfer_coefficient_w_m2k=60.0,
        intervals=20,
        end_time_s=3600.0,
    )


def test_case_is_bounded_and_monotone() -> None:
    """预测的平均含水率必须单调不增、落在环境值与初值之间，且表面先干。"""

    trajectory = run_case(_short_settings())
    assert np.all(np.diff(trajectory.moisture_mean) <= 1.0e-9)
    assert np.isclose(trajectory.moisture_mean[0], 3.8)
    assert trajectory.moisture_mean[-1] < trajectory.moisture_mean[0]
    assert np.all(trajectory.moisture_mean <= 3.8 + 1.0e-9)
    assert np.all(trajectory.moisture_mean >= -1.0e-9)
    assert trajectory.moisture_surface[-1] < trajectory.moisture_centre[-1]
    assert np.all(trajectory.moisture_surface <= trajectory.moisture_centre + 1.0e-9)


def test_sampling_is_idempotent() -> None:
    """在自身时间点上插值不应改变轨迹。"""

    trajectory = run_case(_short_settings())
    resampled = sample_at(trajectory, trajectory.time_s)
    assert np.allclose(resampled.moisture_mean, trajectory.moisture_mean)


def test_diffusivity_inversion_recovers_synthetic_value() -> None:
    """用已知 D 生成合成曲线，反演必须回到该值附近（±10%）。"""

    settings = _short_settings(diffusivity=6.5e-10)
    truth = run_case(settings)
    times_min = np.array([0.0, 15.0, 30.0, 45.0, 60.0])
    sampled = sample_at(truth, times_min * 60.0)
    fitted, rmse = fit_diffusivity(times_min, sampled.moisture_mean, settings)
    assert abs(fitted / 6.5e-10 - 1.0) < 0.10
    assert rmse < 0.01


def test_newman_slab_factor_is_bounded_and_decreasing() -> None:
    """端面解析因子必须落在 (0,1] 且随时间下降。"""

    early = newman_slab_mean_factor(6.5e-10, 0.0516, 3600.0)
    late = newman_slab_mean_factor(6.5e-10, 0.0516, 10800.0)
    assert 0.0 < late < early <= 1.0
    assert late > 0.9
