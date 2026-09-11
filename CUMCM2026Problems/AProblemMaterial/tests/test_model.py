"""热湿有限体积模型基础物理性质的单元测试。"""

import numpy as np

from aproblem.grid import RadialGrid
from aproblem.model import DryingModel
from aproblem.physics import Question1Properties


def test_uniform_state_equal_to_environment_is_stationary() -> None:
    """内部状态等于环境边界时，温度和含水率均不应发生变化。"""

    grid = RadialGrid(intervals=20)
    model = DryingModel(
        grid=grid,
        properties=Question1Properties(),
        environment=lambda _time: (28.0, 2.55),
        radius=lambda _time: 0.02,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
    )
    state = np.concatenate((np.full(21, 28.0), np.full(21, 2.55)))
    assert np.allclose(model.rhs(0.0, state), 0.0)


def test_end_face_equivalent_sources_match_the_derived_formula() -> None:
    """开启端面修正前后的 RHS 差值应等于两个端面的等效体积源项。"""

    common_arguments = dict(
        grid=RadialGrid(intervals=20),
        properties=Question1Properties(),
        environment=lambda _time: (40.0, 0.20),
        radius=lambda _time: 0.02,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
        cylinder_length_m=0.25,
    )
    model_with_ends = DryingModel(**common_arguments, include_end_faces=True)
    model_without_ends = DryingModel(**common_arguments, include_end_faces=False)
    temperature = np.full(21, 30.0)
    moisture = np.full(21, 1.00)
    state = np.concatenate((temperature, moisture))

    rhs_difference = model_with_ends.rhs(0.0, state) - model_without_ends.rhs(0.0, state)
    expected_temperature_rate = 2.0 * 25.0 / 0.25 * (40.0 - 30.0) / (820.0 * 2600.0)
    expected_moisture_rate = 2.0 * 8.0e-7 / 0.25 * (0.20 - 1.00)

    assert np.allclose(rhs_difference[:21], expected_temperature_rate)
    assert np.allclose(rhs_difference[21:], expected_moisture_rate)


def test_uniform_average_rate_uses_all_cylinder_surfaces() -> None:
    """均匀初态的体积平均变化率应满足完整圆柱的 A/V 解析关系。"""

    radius_m = 0.02
    length_m = 0.25
    grid = RadialGrid(intervals=20)
    model = DryingModel(
        grid=grid,
        properties=Question1Properties(),
        environment=lambda _time: (40.0, 0.20),
        radius=lambda _time: radius_m,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
        cylinder_length_m=length_m,
        include_end_faces=True,
    )
    state = np.concatenate((np.full(21, 30.0), np.full(21, 1.00)))
    derivative = model.rhs(0.0, state)
    volumes = grid.geometry(radius_m).volumes_per_length_m2
    area_volume_ratio = 2.0 / radius_m + 2.0 / length_m

    average_temperature_rate = np.average(derivative[:21], weights=volumes)
    average_moisture_rate = np.average(derivative[21:], weights=volumes)
    expected_temperature_rate = (
        25.0 * area_volume_ratio * (40.0 - 30.0) / (820.0 * 2600.0)
    )
    expected_moisture_rate = 8.0e-7 * area_volume_ratio * (0.20 - 1.00)

    assert np.isclose(average_temperature_rate, expected_temperature_rate)
    assert np.isclose(average_moisture_rate, expected_moisture_rate)
