"""L0 解析基准（verification）的单元测试。"""

from __future__ import annotations

import numpy as np
from scipy.special import j0, j1

from aproblem.analytic import (
    BENCHMARK_AMBIENT,
    BENCHMARK_DIFFUSIVITY,
    BENCHMARK_INITIAL,
    BENCHMARK_RADIUS_M,
    CylinderRobinSolution,
    cylinder_eigenvalues,
    observed_order,
    run_steady_case,
    run_time_step_case,
    run_transient_case,
)


def test_eigenvalues_satisfy_robin_condition() -> None:
    """特征值必须满足 βJ₁(β)=Bi·J₀(β)，且严格递增。"""

    biot = 5.0
    eigenvalues = cylinder_eigenvalues(biot, 8)
    residual = eigenvalues * j1(eigenvalues) - biot * j0(eigenvalues)
    assert np.max(np.abs(residual)) < 1.0e-10
    assert np.all(np.diff(eigenvalues) > 0.0)


def test_series_truncation_is_converged() -> None:
    """60 项与 120 项级数解必须一致到远优于空间误差的水平。"""

    common = dict(
        radius_m=BENCHMARK_RADIUS_M,
        diffusivity_m2_s=BENCHMARK_DIFFUSIVITY,
        transfer_coefficient_m_s=5.0 * BENCHMARK_DIFFUSIVITY / BENCHMARK_RADIUS_M,
        initial_value=BENCHMARK_INITIAL,
        ambient_value=BENCHMARK_AMBIENT,
    )
    short = CylinderRobinSolution(**common, terms=60)
    long = CylinderRobinSolution(**common, terms=120)
    time_s = 0.05 * BENCHMARK_RADIUS_M**2 / BENCHMARK_DIFFUSIVITY
    assert np.allclose(short.value(time_s, [0.0, 0.5, 1.0]), long.value(time_s, [0.0, 0.5, 1.0]), rtol=1e-10)


def test_mean_matches_quadrature() -> None:
    """闭式体积平均必须与 2∫θξdξ 的数值积分一致。"""

    solution = CylinderRobinSolution(
        radius_m=BENCHMARK_RADIUS_M,
        diffusivity_m2_s=BENCHMARK_DIFFUSIVITY,
        transfer_coefficient_m_s=5.0 * BENCHMARK_DIFFUSIVITY / BENCHMARK_RADIUS_M,
        initial_value=BENCHMARK_INITIAL,
        ambient_value=BENCHMARK_AMBIENT,
    )
    time_s = 0.1 * BENCHMARK_RADIUS_M**2 / BENCHMARK_DIFFUSIVITY
    xi = np.linspace(0.0, 1.0, 20001)
    values = solution.value(time_s, xi)
    quadrature = 2.0 * np.trapz(values * xi, xi)
    assert np.isclose(solution.mean(time_s), quadrature, rtol=1.0e-8)


def test_transient_converges_at_second_order() -> None:
    """生产 RHS 对解析解的误差必须按二阶收敛（Bi=5，N=40→80，水分与温度）。"""

    coarse = run_transient_case(5.0, 40)
    fine = run_transient_case(5.0, 80)
    order = observed_order(coarse.moisture_surface, fine.moisture_surface, 2.0)
    assert 1.7 <= order <= 2.4
    assert coarse.moisture_surface < 5.0e-3
    assert fine.moisture_surface < 1.5e-3
    temperature_order = observed_order(
        coarse.temperature_surface,
        fine.temperature_surface,
        2.0,
    )
    assert 1.7 <= temperature_order <= 2.4
    assert abs(coarse.temperature_surface - coarse.moisture_surface) < 1.0e-3


def test_heun_time_order() -> None:
    """固定步长 Heun 必须按二阶自收敛（粗网格版，逐档差商比 2）。"""

    intervals = 40
    step = 0.4 * 1.25  # N=40 时 Δt_stable 约 1.25 s，取 40% 与 20% 两档
    coarse = run_time_step_case(step, intervals)
    medium = run_time_step_case(0.5 * step, intervals)
    fine = run_time_step_case(0.25 * step, intervals)
    difference_coarse = float(np.max(np.abs(coarse - medium)))
    difference_fine = float(np.max(np.abs(medium - fine)))
    order = observed_order(difference_coarse, difference_fine, 2.0)
    assert 1.6 <= order <= 2.4


def test_steady_kirchhoff_identity_and_ordering() -> None:
    """基尔霍夫恒等式精确成立，且优于调和与串联策略。"""

    rows = {row["strategy"]: row for row in run_steady_case(40)}
    identity = rows["kirchhoff"]["kirchhoff_identity_residual"]
    assert identity < 1.0e-12
    assert rows["kirchhoff"]["max_relative_flux_error"] < 1.0e-3
    assert rows["harmonic"]["max_relative_flux_error"] > rows["kirchhoff"]["max_relative_flux_error"]
    assert rows["series"]["max_relative_flux_error"] > rows["kirchhoff"]["max_relative_flux_error"]
