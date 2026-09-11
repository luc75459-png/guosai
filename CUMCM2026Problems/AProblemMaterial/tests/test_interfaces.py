"""界面物性取值策略的单元测试。"""

from __future__ import annotations

import numpy as np
import pytest

from aproblem.interfaces import (
    INTERFACE_STRATEGY_NAMES,
    HarmonicInterface,
    MidpointInterface,
    SeriesResistanceInterface,
    resolve_interface_strategy,
)


def test_all_registered_strategies_resolve_by_name() -> None:
    """注册表里的每个名字都应能解析成同名策略实例。"""

    for name in INTERFACE_STRATEGY_NAMES:
        assert resolve_interface_strategy(name).name == name
    # 校准型混合权重不是合法的界面策略，必须被拒绝。
    with pytest.raises(ValueError):
        resolve_interface_strategy("calibrated_blend")


def test_equal_values_return_the_same_value_in_every_strategy() -> None:
    """两侧物性相等时，三套策略都退化为该常数。"""

    radii = np.linspace(0.0, 0.02, 21)
    values = np.full(radii.shape, 3.0e-9)
    temperature = np.linspace(28.0, 50.0, 21)
    moisture = np.linspace(2.55, 0.2, 21)

    for strategy in (
        MidpointInterface(),
        HarmonicInterface(),
        SeriesResistanceInterface(),
    ):
        faces = strategy.face_property(
            node_property=values,
            node_temperature=temperature,
            node_moisture=moisture,
            node_radii_m=radii,
            property_function=lambda _t, _c: np.full(radii.shape, 3.0e-9),
        )
        assert np.allclose(faces, 3.0e-9)


def test_series_reduces_to_harmonic_mean_on_a_flat_slab() -> None:
    """格距远小于半径时（平板极限），圆柱串联阻力退化为调和平均。

    二维模型的轴向方向是平板几何，用远离原点的等距"半径"表示，
    此时对数项趋于常数，公式必须收敛到调和平均。
    """

    intervals = 32
    spacing = 0.25 / intervals
    radii = 0.5 + np.arange(intervals + 1) * spacing
    values = np.geomspace(1.0e-9, 1.0e-5, intervals + 1)
    temperature = np.zeros(intervals + 1)
    moisture = np.zeros(intervals + 1)

    series = SeriesResistanceInterface().face_property(
        node_property=values,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
    )
    harmonic = HarmonicInterface().face_property(
        node_property=values,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
    )
    assert np.allclose(series, harmonic, rtol=1.0e-6)


def test_series_reproduces_exact_annular_series_resistance() -> None:
    """界面等效物性必须与两段圆环串联阻力的解析结果一致。"""

    inner_radius = 0.004
    outer_radius = 0.010
    face_radius = 0.5 * (inner_radius + outer_radius)
    left_value = 2.0e-9
    right_value = 8.0e-7
    radii = np.array([inner_radius, outer_radius])
    values = np.array([left_value, right_value])

    face_value = SeriesResistanceInterface().face_property(
        node_property=values,
        node_temperature=np.zeros(2),
        node_moisture=np.zeros(2),
        node_radii_m=radii,
    )[0]

    # 解析串联阻力：两段圆环的稳态径向导热热阻之和。
    analytic_resistance = np.log(face_radius / inner_radius) / (
        2.0 * np.pi * left_value
    ) + np.log(outer_radius / face_radius) / (2.0 * np.pi * right_value)
    # 有限体积把通量写成 k_face * 2π r_f * ΔT / (r_out - r_in)。
    equivalent_resistance = (outer_radius - inner_radius) / (
        face_value * 2.0 * np.pi * face_radius
    )
    assert np.isclose(equivalent_resistance, analytic_resistance)


def test_small_side_dominates_harmonic_and_series_but_not_midpoint() -> None:
    """两侧物性相差四个数量级时，调和与串联阻力被小值支配。"""

    radii = np.array([0.01, 0.02])
    values = np.array([1.0e-9, 1.0e-5])
    temperature = np.array([30.0, 30.0])
    moisture = np.array([1.0, 1.0])

    harmonic = HarmonicInterface().face_property(
        node_property=values,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
    )[0]
    series = SeriesResistanceInterface().face_property(
        node_property=values,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
    )[0]
    midpoint = MidpointInterface().face_property(
        node_property=values,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
        property_function=lambda _t, c: c,
    )[0]

    assert harmonic < 2.0e-9
    assert series < 2.0e-9
    assert np.isclose(midpoint, 1.0)


def test_series_uses_harmonic_fallback_on_the_centre_face() -> None:
    """半径为零的中心面没有左侧半格，必须退回调和平均而不是发散。"""

    intervals = 20
    radii = np.linspace(0.0, 0.02, intervals + 1)
    values = np.geomspace(1.0e-9, 1.0e-5, intervals + 1)
    temperature = np.zeros(intervals + 1)
    moisture = np.zeros(intervals + 1)

    series = SeriesResistanceInterface().face_property(
        node_property=values,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
    )
    assert np.isfinite(series).all()
    expected_centre = 2.0 * values[0] * values[1] / (values[0] + values[1])
    assert np.isclose(series[0], expected_centre)
