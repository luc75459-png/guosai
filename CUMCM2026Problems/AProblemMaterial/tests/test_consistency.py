"""干物质守恒自洽检验的单元测试。"""

from __future__ import annotations

import numpy as np
import pytest

from aproblem.consistency import (
    check_volume_consistency,
    dry_density_question4,
    radial_volume_weights,
)


def test_radial_volume_weights_are_normalised_and_match_annulus_areas() -> None:
    """环面积权重必须等于控制体横截面积且总和为一。"""

    fraction = np.linspace(0.0, 1.0, 21)
    weights = radial_volume_weights(fraction)
    assert np.isclose(np.sum(weights), 1.0)

    inner = np.empty_like(fraction)
    inner[0] = 0.0
    inner[1:] = 0.5 * (fraction[:-1] + fraction[1:])
    outer = np.empty_like(fraction)
    outer[:-1] = inner[1:]
    outer[-1] = 1.0
    areas = np.pi * (outer**2 - inner**2)
    assert np.allclose(weights, areas / np.sum(areas))


def test_dry_density_is_apparent_density_divided_by_one_plus_moisture() -> None:
    """干物质密度换算公式必须与附录4一致。"""

    moisture = np.array([0.0, 0.15, 2.55])
    assert np.allclose(dry_density_question4(moisture), (760.0 + 90.0 * moisture) / (1.0 + moisture))


def test_uniform_profile_reproduces_the_analytic_volume_ratio() -> None:
    """均匀含水率剖面的检验结果应与解析式完全一致。"""

    initial_moisture = 2.55
    final_moisture = 0.15
    radius_ratio = 0.6
    report = check_volume_consistency(
        np.full(21, final_moisture),
        initial_moisture,
        radius_ratio,
    )

    expected_volume_ratio = float(
        dry_density_question4(np.array([initial_moisture]))[0]
        / dry_density_question4(np.array([final_moisture]))[0]
    )
    assert np.isclose(report.required_volume_ratio, expected_volume_ratio)
    assert np.isclose(report.radial_area_ratio, radius_ratio**2)
    assert np.isclose(
        report.implied_length_ratio,
        expected_volume_ratio / radius_ratio**2,
    )


def test_nonuniform_profile_requires_between_the_extremes() -> None:
    """非均匀剖面要求的体积比应落在最干与最湿两层之间。"""

    initial_moisture = 2.55
    profile = np.linspace(0.15, 0.05, 21)
    report = check_volume_consistency(profile, initial_moisture, 0.6)

    ratio_wet = float(
        dry_density_question4(np.array([initial_moisture]))[0]
        / dry_density_question4(np.array([0.15]))[0]
    )
    ratio_dry = float(
        dry_density_question4(np.array([initial_moisture]))[0]
        / dry_density_question4(np.array([0.05]))[0]
    )
    assert min(ratio_wet, ratio_dry) < report.required_volume_ratio < max(
        ratio_wet,
        ratio_dry,
    )


def test_invalid_inputs_are_rejected() -> None:
    """形状不匹配或非正半径比必须报错。"""

    with pytest.raises(ValueError):
        check_volume_consistency(np.zeros(5), 2.55, 0.6, radial_fraction=np.zeros(4))
    with pytest.raises(ValueError):
        check_volume_consistency(np.zeros(5), 2.55, 0.0)
