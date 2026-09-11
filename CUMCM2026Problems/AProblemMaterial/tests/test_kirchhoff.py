"""基尔霍夫通量势与 kirchhoff 界面策略的单元测试。"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from aproblem.interfaces import (
    INTERFACE_STRATEGY_NAMES,
    HarmonicInterface,
    KirchhoffInterface,
    MidpointInterface,
    resolve_interface_strategy,
)
from aproblem.physics import (
    Question1Properties,
    Question23Properties,
    Question4Properties,
)


PROPERTY_MODELS = (
    Question1Properties(),
    Question23Properties(),
    Question4Properties(),
)


def test_kirchhoff_strategy_is_registered_and_resolvable() -> None:
    """新策略必须进入注册表，才能被配置和命令行按名字选中。"""

    assert "kirchhoff" in INTERFACE_STRATEGY_NAMES
    assert resolve_interface_strategy("kirchhoff").name == "kirchhoff"


def test_questions_two_and_three_default_to_kirchhoff() -> None:
    """问题2/3 必须默认走基尔霍夫通量势，显式指定时以显式值为准。"""

    from aproblem.config import (
        ProjectPaths,
        SimulationConfig,
        default_interface_strategy,
    )
    from aproblem.scenarios import build_question

    assert default_interface_strategy(1) == "midpoint"
    assert default_interface_strategy(2) == "kirchhoff"
    assert default_interface_strategy(3) == "kirchhoff"
    assert default_interface_strategy(4) == "midpoint"

    paths = ProjectPaths.discover()
    for question in (2, 3):
        setup = build_question(
            question,
            paths,
            SimulationConfig(radial_intervals=8),
        )
        assert setup.model.interface_strategy.name == "kirchhoff"

    overridden = build_question(
        2,
        paths,
        SimulationConfig(radial_intervals=8, interface_strategy="midpoint"),
    )
    assert overridden.model.interface_strategy.name == "midpoint"


@pytest.mark.parametrize("properties", PROPERTY_MODELS)
def test_flux_potential_derivative_recovers_diffusivity(properties) -> None:
    """解析通量势对含水率的导数必须逐点还原题目给出的 D(C,T)。

    这是整个推导的支点：只有 dΦ/dC = D 成立，ΔΦ/ΔC 才等价于真实通量。
    """

    for temperature_c in (28.0, 50.0):
        for moisture in (0.15, 0.5, 1.0, 2.55):
            temperature = np.array([temperature_c])
            centre = np.array([moisture])
            step = 1.0e-6 * moisture
            derivative = (
                properties.flux_potential(temperature, centre + step)
                - properties.flux_potential(temperature, centre - step)
            ) / (2.0 * step)
            expected = properties.diffusivity(temperature, centre)
            assert np.isclose(derivative[0], expected[0], rtol=1.0e-4)


@pytest.mark.parametrize("properties", PROPERTY_MODELS)
def test_flux_potential_matches_numerical_integration(properties) -> None:
    """解析通量势应与独立数值积分一致。"""

    temperature_c = 50.0
    for moisture in (0.05, 0.15, 0.5, 2.55):
        reference, _ = quad(
            lambda c: float(
                properties.diffusivity(
                    np.array([temperature_c]),
                    np.array([c]),
                )[0]
            ),
            0.0,
            moisture,
            limit=200,
        )
        value = float(
            properties.flux_potential(
                np.array([temperature_c]),
                np.array([moisture]),
            )[0]
        )
        assert np.isclose(value, reference, rtol=1.0e-5, atol=1.0e-30)


@pytest.mark.parametrize("properties", PROPERTY_MODELS)
def test_kirchhoff_face_value_matches_quadrature_of_diffusivity(properties) -> None:
    """界面值必须等于 ∫D dC / ΔC，而不是任何端点组合。"""

    temperature = np.array([50.0, 50.0])
    moisture = np.array([0.10, 0.35])
    faces = KirchhoffInterface().face_property(
        node_property=properties.diffusivity(temperature, moisture),
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=np.array([0.0, 0.001]),
        property_function=properties.diffusivity,
        flux_potential_function=properties.flux_potential,
    )
    reference, _ = quad(
        lambda c: float(
            properties.diffusivity(np.array([50.0]), np.array([c]))[0]
        ),
        0.10,
        0.35,
        limit=200,
    )
    reference /= 0.35 - 0.10
    assert np.isclose(faces[0], reference, rtol=1.0e-5)


def test_kirchhoff_is_far_from_harmonic_in_the_dry_shell() -> None:
    """干壳区必须复现"调和平均只剩真值零头"这一现象。

    这是纯调和平均在粗网格上不收敛的物理根源：它把表面薄壳当成一层真实
    的分层材料，人为制造出巨大的传质阻力。
    """

    properties = Question23Properties()
    temperature = np.array([50.0, 50.0])
    moisture = np.array([0.10, 0.05])
    arguments = dict(
        node_property=properties.diffusivity(temperature, moisture),
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=np.array([0.0, 0.001]),
        property_function=properties.diffusivity,
    )
    kirchhoff = KirchhoffInterface().face_property(
        **arguments,
        flux_potential_function=properties.flux_potential,
    )[0]
    harmonic = HarmonicInterface().face_property(**arguments)[0]
    assert harmonic < 0.2 * kirchhoff


def test_kirchhoff_reduces_to_the_point_value_for_constant_diffusivity() -> None:
    """扩散系数为常数时，通量势斜率必须精确退化为该常数。"""

    constant = 3.0e-9
    face = KirchhoffInterface().face_property(
        node_property=np.full(2, constant),
        node_temperature=np.array([30.0, 30.0]),
        node_moisture=np.array([1.0, 2.0]),
        node_radii_m=np.array([0.0, 0.001]),
        property_function=lambda _t, c: np.full(np.shape(c), constant),
        flux_potential_function=lambda _t, c: constant * np.asarray(c, dtype=float),
    )
    assert np.isclose(face[0], constant)


def test_kirchhoff_falls_back_to_midpoint_without_a_potential() -> None:
    """传热方程没有通量势，此时必须与主方案的中点值逐面一致。"""

    radii = np.linspace(0.0, 0.02, 21)
    temperature = np.linspace(28.0, 50.0, 21)
    moisture = np.linspace(2.55, 0.20, 21)
    conductivity = np.linspace(0.12, 0.32, 21)
    arguments = dict(
        node_property=conductivity,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
        property_function=lambda _t, c: 0.12 + 0.20 * c / (c + 1.0),
    )
    assert np.allclose(
        KirchhoffInterface().face_property(**arguments),
        MidpointInterface().face_property(**arguments),
    )


def test_kirchhoff_handles_equal_moisture_without_division_by_zero() -> None:
    """全场含水率均匀时 ΔC=0，必须退回点值而不是产生 nan。"""

    properties = Question23Properties()
    temperature = np.full(3, 50.0)
    moisture = np.full(3, 1.0)
    faces = KirchhoffInterface().face_property(
        node_property=properties.diffusivity(temperature, moisture),
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=np.array([0.0, 0.001, 0.002]),
        property_function=properties.diffusivity,
        flux_potential_function=properties.flux_potential,
    )
    expected = properties.diffusivity(temperature, moisture)[0]
    assert np.all(np.isfinite(faces))
    assert np.allclose(faces, expected)


def test_kirchhoff_runs_through_the_one_and_two_dimensional_solvers() -> None:
    """基尔霍夫策略必须原样跑通一维和二维两条求解通道。

    只做短时积分：目的是确认通量势接口在两条通道上都接得通、状态保持有限，
    长时程的结果由 interface-scan 的网格对照负责。
    """

    from aproblem.axisymmetric import AxisymmetricModel
    from aproblem.grid import RadialGrid
    from aproblem.integrator import integrate_heun
    from aproblem.model import DryingModel

    strategy = resolve_interface_strategy("kirchhoff")
    properties = Question4Properties()

    def environment(_time_s: float) -> tuple[float, float]:
        return 50.0, 0.05

    one_dimensional = DryingModel(
        grid=RadialGrid(8),
        properties=properties,
        environment=environment,
        radius=lambda _time: 0.02,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
        interface_strategy=strategy,
    )
    one_dimensional_state = np.concatenate(
        (np.full(9, 28.0), np.full(9, 2.55))
    )
    one_dimensional_result = integrate_heun(
        rhs=one_dimensional.rhs,
        initial_state=one_dimensional_state,
        end_time_s=600.0,
        dt_s=0.5,
        save_every_s=600.0,
        max_dt_fn=lambda time_s, value: one_dimensional.stable_time_step(
            time_s, value
        ),
    )
    assert np.all(np.isfinite(one_dimensional_result.state))

    two_dimensional = AxisymmetricModel(
        radial_intervals=6,
        axial_intervals=4,
        properties=properties,
        environment=environment,
        radius=lambda _time: 0.02,
        length=lambda _time: 0.25,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
        interface_strategy=strategy,
    )
    two_dimensional_state = two_dimensional.pack(
        np.full(two_dimensional.shape, 28.0),
        np.full(two_dimensional.shape, 2.55),
    )
    two_dimensional_result = integrate_heun(
        rhs=two_dimensional.rhs,
        initial_state=two_dimensional_state,
        end_time_s=600.0,
        dt_s=0.5,
        save_every_s=600.0,
        max_dt_fn=lambda time_s, value: two_dimensional.stable_time_step(
            time_s, value
        ),
    )
    assert np.all(np.isfinite(two_dimensional_result.state))
