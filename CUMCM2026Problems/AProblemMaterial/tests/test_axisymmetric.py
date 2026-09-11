"""二维轴对称材料坐标模型的单元测试。"""

from __future__ import annotations

import numpy as np

from aproblem.axisymmetric import AxisymmetricModel, axisymmetric_geometry
from aproblem.config import SimulationConfig, ProjectPaths
from aproblem.grid import RadialGrid
from aproblem.integrator import integrate_heun
from aproblem.interfaces import resolve_interface_strategy
from aproblem.model import DryingModel
from aproblem.physics import Question1Properties, Question4Properties


def test_control_volumes_and_surfaces_cover_the_cylinder() -> None:
    """控制体体积、侧面面积和端面面积必须与解析值一致。"""

    radius_m = 0.02
    length_m = 0.25
    geometry = axisymmetric_geometry(20, 16, radius_m, length_m)

    assert np.isclose(np.sum(geometry.volumes_m3), np.pi * radius_m**2 * length_m)
    assert np.isclose(np.sum(geometry.side_areas_m2), 2.0 * np.pi * radius_m * length_m)
    assert np.isclose(np.sum(geometry.end_areas_m2), np.pi * radius_m**2)
    assert np.isclose(geometry.cross_section_area_m2, np.pi * radius_m**2)


def test_uniform_state_equal_to_environment_is_stationary() -> None:
    """全场状态等于环境边界时，二维模型的右端项应为零。"""

    shape = (11, 9)
    model = AxisymmetricModel(
        radial_intervals=shape[0] - 1,
        axial_intervals=shape[1] - 1,
        properties=Question1Properties(),
        environment=lambda _time: (28.0, 2.55),
        radius=lambda _time: 0.02,
        length=lambda _time: 0.25,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
    )
    state = model.pack(np.full(shape, 28.0), np.full(shape, 2.55))
    assert np.allclose(model.rhs(0.0, state), 0.0)


def _run_two_dimensional(
    insulated_ends: bool,
    end_time_s: float,
    radius_m: float = 0.02,
    length_m: float = 0.25,
    intervals: int = 20,
    axial_intervals: int = 4,
) -> tuple[AxisymmetricModel, np.ndarray]:
    """在恒定半径和恒定环境下短时积分二维模型，返回末态含水率场。"""

    model = AxisymmetricModel(
        radial_intervals=intervals,
        axial_intervals=axial_intervals,
        properties=Question1Properties(),
        environment=lambda _time: (45.0, 0.05),
        radius=lambda _time: radius_m,
        length=lambda _time: length_m,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
        interface_strategy=resolve_interface_strategy("midpoint"),
        insulated_ends=insulated_ends,
    )
    initial = model.pack(
        np.full(model.shape, 28.0),
        np.full(model.shape, 2.55),
    )
    result = integrate_heun(
        rhs=model.rhs,
        initial_state=initial,
        end_time_s=end_time_s,
        dt_s=0.5,
        save_every_s=end_time_s,
        max_dt_fn=lambda time_s, state: model.stable_time_step(time_s, state),
    )
    states = np.asarray(result.state)
    moisture = states[-1, model.node_count :].reshape(model.shape)
    return model, moisture


def _run_one_dimensional(
    include_end_faces: bool,
    end_time_s: float,
    radius_m: float = 0.02,
    length_m: float = 0.25,
    intervals: int = 20,
) -> np.ndarray:
    """用一维径向模型跑同样的短时算例，返回末态含水率剖面。"""

    model = DryingModel(
        grid=RadialGrid(intervals=intervals),
        properties=Question1Properties(),
        environment=lambda _time: (45.0, 0.05),
        radius=lambda _time: radius_m,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
        cylinder_length_m=length_m,
        include_end_faces=include_end_faces,
        interface_strategy=resolve_interface_strategy("midpoint"),
    )
    state = np.concatenate(
        (np.full(model.node_count, 28.0), np.full(model.node_count, 2.55))
    )
    result = integrate_heun(
        rhs=model.rhs,
        initial_state=state,
        end_time_s=end_time_s,
        dt_s=0.5,
        save_every_s=end_time_s,
        max_dt_fn=lambda time_s, value: model.stable_time_step(time_s, value),
    )
    return np.asarray(result.state)[-1, model.node_count :]


def test_insulated_ends_degenerate_to_the_one_dimensional_model() -> None:
    """端面绝热绝质时，二维解必须与一维"仅侧面"模型逐节点一致。"""

    model, two_dimensional = _run_two_dimensional(
        insulated_ends=True,
        end_time_s=1800.0,
    )
    one_dimensional = _run_one_dimensional(
        include_end_faces=False,
        end_time_s=1800.0,
    )

    # 端面绝热时轴向没有梯度，任意轴向位置都应与一维剖面相同。
    for axial_index in range(model.shape[1]):
        assert np.allclose(two_dimensional[:, axial_index], one_dimensional)
    assert np.allclose(
        two_dimensional.std(axis=1),
        0.0,
        atol=1.0e-12,
    )


def test_open_ends_dry_faster_than_insulated_ends() -> None:
    """同一时刻，计入端面的二维模型失水必须多于端面绝热的模型。"""

    _, insulated = _run_two_dimensional(insulated_ends=True, end_time_s=1800.0)
    _, open_ends = _run_two_dimensional(insulated_ends=False, end_time_s=1800.0)

    assert open_ends.mean() < insulated.mean()
    # 端面只占总表面积的 8%，差距应为一个有限的小量而不是数量级差异。
    loss_with_ends = 2.55 - open_ends.mean()
    loss_without_ends = 2.55 - insulated.mean()
    assert loss_with_ends > loss_without_ends
    assert loss_with_ends / loss_without_ends < 1.5


def test_axisymmetric_model_builds_from_project_paths() -> None:
    """二维模型可以直接从题目附件和配置组装起来。"""

    from aproblem.study import build_axisymmetric_model

    paths = ProjectPaths.discover()
    config = SimulationConfig(radial_intervals=12, interface_strategy="midpoint")
    model = build_axisymmetric_model(4, paths, config, 12, 8)

    assert model.shape == (13, 9)
    assert isinstance(model.properties, Question4Properties)
    assert np.isclose(model.radius(0.0), 0.02)
    assert np.isclose(model.length(0.0), 0.25)
    assert np.isclose(model.radius(259200.0), model.radius(259200.0))
