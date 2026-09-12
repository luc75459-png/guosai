"""L0 解析基准（verification）：瞬态圆柱 Robin、Heun 时间阶、变系数稳态通量。

本模块只做**代码验证**（verification），不参与任何正式结果生成，也不修改
生产模块。三个算例对应预注册判据（``docs/06`` §9.1–§9.3）：

* L0-A  无限长圆柱、常数物性、侧面 Robin 的解析级数解 vs 生产 ``DryingModel``；
* L0-B  同一半离散右端项下 Heun 时间推进的观测阶；
* L0-C  ``D(C)=A e^{-b/C}`` 下四套界面策略的稳态通量对拍（基尔霍夫恒等式）。

用法：``python -m aproblem.analytic``，产物写入 ``outputs/study/analytic/``。
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import j0, j1

from .config import ProjectPaths
from .grid import RadialGrid
from .integrator import integrate_heun
from .interfaces import (
    INTERFACE_STRATEGY_NAMES,
    InterfaceStrategy,
    resolve_interface_strategy,
)
from .model import DryingModel
from .physics import MaterialProperties, PropertyModel, Question23Properties


# --------------------------------------------------------------------------
# 基础工具：常数物性模型与圆柱 Robin 解析解
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ConstantProperties:
    """常数物性模型，用于让生产求解器退化到可解析的情形。"""

    density_kg_m3: float
    heat_capacity_j_kgk: float
    conductivity_w_mk: float
    diffusivity_m2_s: float

    def _full(self, reference: np.ndarray, value: float) -> np.ndarray:
        return np.full(np.shape(reference), float(value))

    def density(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        return self._full(temperature_c, self.density_kg_m3)

    def heat_capacity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        return self._full(temperature_c, self.heat_capacity_j_kgk)

    def conductivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        return self._full(temperature_c, self.conductivity_w_mk)

    def diffusivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        return self._full(moisture, self.diffusivity_m2_s)

    def flux_potential(
        self,
        temperature_c: np.ndarray,
        moisture: np.ndarray,
    ) -> np.ndarray:
        """常数 D 的通量势 Φ=D·C（解析且精确）。"""

        return self.diffusivity_m2_s * np.asarray(moisture, dtype=float)

    def evaluate(
        self,
        temperature_c: np.ndarray,
        moisture: np.ndarray,
    ) -> MaterialProperties:
        return MaterialProperties(
            density=self.density(temperature_c, moisture),
            heat_capacity=self.heat_capacity(temperature_c, moisture),
            conductivity=self.conductivity(temperature_c, moisture),
            diffusivity=self.diffusivity(temperature_c, moisture),
        )


def cylinder_eigenvalues(biot: float, count: int) -> np.ndarray:
    """求圆柱 Robin 问题的前 ``count`` 个特征值 β（βJ₁(β)=Bi·J₀(β)）。

    特征值两两交错于 J₁ 的零点之间；用密集扫描定位符号变化区间，再用
    Brent 法二分到机器精度。扫描步长 0.01 对前 60 个根已有充分余量。
    """

    if biot <= 0.0:
        raise ValueError("Biot 数必须为正数")
    if count < 1:
        raise ValueError("特征值个数必须为正整数")

    def residual(beta: float) -> float:
        return float(beta * j1(beta) - biot * j0(beta))

    upper = (count + 2) * np.pi
    grid = np.arange(1.0e-6, upper, 0.01)
    values = grid * j1(grid) - biot * j0(grid)
    roots: list[float] = []
    for index in range(grid.size - 1):
        left_value = values[index]
        right_value = values[index + 1]
        if left_value == 0.0:
            root = float(grid[index])
        elif left_value * right_value < 0.0:
            root = brentq(residual, float(grid[index]), float(grid[index + 1]), xtol=1e-14)
        else:
            continue
        if root > 1.0e-6:
            roots.append(root)
        if len(roots) >= count:
            break
    if len(roots) < count:
        raise RuntimeError(f"只找到 {len(roots)} 个特征值，少于要求的 {count} 个")
    return np.asarray(roots, dtype=float)


@dataclass(frozen=True)
class CylinderRobinSolution:
    """无限长圆柱、均匀初值、常数环境、侧面 Robin 的解析级数解。

    θ(ξ,τ)=Σ Aₙ J₀(βₙξ)exp(−βₙ²τ)，τ=Dt/R²，Bi=hR/D，
    Aₙ=2(C₀−C∞)J₁(βₙ)/[βₙ(J₀²(βₙ)+J₁²(βₙ))]。
    """

    radius_m: float
    diffusivity_m2_s: float
    transfer_coefficient_m_s: float
    initial_value: float
    ambient_value: float
    terms: int = 60

    @property
    def biot(self) -> float:
        return self.transfer_coefficient_m_s * self.radius_m / self.diffusivity_m2_s

    @property
    def span(self) -> float:
        return self.initial_value - self.ambient_value

    def __post_init__(self) -> None:
        if self.radius_m <= 0.0 or self.diffusivity_m2_s <= 0.0:
            raise ValueError("半径与扩散系数必须为正数")
        eigenvalues = cylinder_eigenvalues(self.biot, self.terms)
        amplitudes = (
            2.0
            * self.span
            * j1(eigenvalues)
            / (eigenvalues * (j0(eigenvalues) ** 2 + j1(eigenvalues) ** 2))
        )
        object.__setattr__(self, "eigenvalues", eigenvalues)
        object.__setattr__(self, "amplitudes", amplitudes)

    def value(self, time_s: float, radial_fraction) -> np.ndarray:
        """返回给定时刻、给定归一化半径 ξ=r/R 处的解析值。"""

        xi = np.atleast_1d(np.asarray(radial_fraction, dtype=float))
        fourier = self.diffusivity_m2_s * float(time_s) / self.radius_m**2
        decay = np.exp(-(self.eigenvalues**2) * fourier)
        series = np.sum(
            self.amplitudes[:, None]
            * j0(self.eigenvalues[:, None] * xi[None, :])
            * decay[:, None],
            axis=0,
        )
        return self.ambient_value + series

    def centre(self, time_s: float) -> float:
        return float(self.value(time_s, 0.0)[0])

    def surface(self, time_s: float) -> float:
        return float(self.value(time_s, 1.0)[0])

    def mean(self, time_s: float) -> float:
        """体积平均（权重 ξdξ）：⟨θ⟩=Σ Aₙ e^{−βₙ²τ}·2J₁(βₙ)/βₙ。"""

        fourier = self.diffusivity_m2_s * float(time_s) / self.radius_m**2
        decay = np.exp(-(self.eigenvalues**2) * fourier)
        series = np.sum(self.amplitudes * decay * 2.0 * j1(self.eigenvalues) / self.eigenvalues)
        return float(self.ambient_value + series)


# --------------------------------------------------------------------------
# L0-A：瞬态圆柱 Robin 空间收敛阶
# --------------------------------------------------------------------------

# 基准参数（无量纲化后与本题物性无关，只检验离散格式与边界实现）。
BENCHMARK_RADIUS_M = 0.01
BENCHMARK_DIFFUSIVITY = 1.0e-8
BENCHMARK_ALPHA = 1.0e-8
BENCHMARK_DENSITY = 1000.0
BENCHMARK_HEAT_CAPACITY = 1.0
BENCHMARK_INITIAL = 1.0
BENCHMARK_AMBIENT = 0.0
BENCHMARK_INITIAL_TEMPERATURE = 100.0
BENCHMARK_AMBIENT_TEMPERATURE = 20.0

TRANSIT_BIOT = (0.5, 5.0, 50.0)
TRANSIT_FOURIER = (0.01, 0.05, 0.1, 0.2, 0.5)
TRANSIT_GRIDS = (10, 20, 40, 80, 160)
TRANSIT_RADIAL_FRACTION = (0.0, 0.5, 1.0)

BDF_RTOL = 1.0e-11
BDF_ATOL = 1.0e-13


def build_benchmark_model(biot: float, intervals: int) -> DryingModel:
    """构造 L0-A 的基准模型（常数物性、固定半径、仅侧面 Robin）。"""

    radius = BENCHMARK_RADIUS_M
    conductivity = BENCHMARK_DENSITY * BENCHMARK_HEAT_CAPACITY * BENCHMARK_ALPHA
    properties = ConstantProperties(
        density_kg_m3=BENCHMARK_DENSITY,
        heat_capacity_j_kgk=BENCHMARK_HEAT_CAPACITY,
        conductivity_w_mk=conductivity,
        diffusivity_m2_s=BENCHMARK_DIFFUSIVITY,
    )
    return DryingModel(
        grid=RadialGrid(intervals),
        properties=properties,
        environment=lambda _time: (BENCHMARK_AMBIENT_TEMPERATURE, BENCHMARK_AMBIENT),
        radius=lambda _time: radius,
        heat_transfer_coefficient=biot * conductivity / radius,
        mass_transfer_coefficient=biot * BENCHMARK_DIFFUSIVITY / radius,
        include_end_faces=False,
        interface_strategy=resolve_interface_strategy("kirchhoff"),
    )


def _moisture_right_hand_side(
    model: DryingModel,
    temperature_c: float,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """从生产右端项中切出水分方程（常数物性下与温度解耦）。"""

    node_count = model.node_count
    fixed_temperature = np.full(node_count, temperature_c)

    def rhs(time_s: float, moisture: np.ndarray) -> np.ndarray:
        state = np.concatenate((fixed_temperature, np.asarray(moisture, dtype=float)))
        return model.rhs(time_s, state)[node_count:]

    return rhs


def _temperature_right_hand_side(
    model: DryingModel,
    moisture: float,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """从生产右端项中切出温度方程（常数物性下与含水率解耦）。"""

    node_count = model.node_count
    fixed_moisture = np.full(node_count, moisture)

    def rhs(time_s: float, temperature: np.ndarray) -> np.ndarray:
        state = np.concatenate((np.asarray(temperature, dtype=float), fixed_moisture))
        return model.rhs(time_s, state)[:node_count]

    return rhs


def _solve_bdf(
    rhs: Callable[[float, np.ndarray], np.ndarray],
    initial: np.ndarray,
    times: np.ndarray,
) -> np.ndarray:
    """用紧容差 BDF 积分一维场，返回 ``times`` 时刻的解（隔离空间误差用）。"""

    times = np.asarray(times, dtype=float)
    solution = solve_ivp(
        rhs,
        (0.0, float(times[-1])),
        np.asarray(initial, dtype=float),
        method="BDF",
        rtol=BDF_RTOL,
        atol=BDF_ATOL,
        t_eval=times,
    )
    if not solution.success:
        raise RuntimeError(f"L0 基准的 BDF 积分失败：{solution.message}")
    return np.asarray(solution.y, dtype=float).T


def _volume_weights(intervals: int) -> np.ndarray:
    """节点型圆柱控制体的归一化体积权重（正比于环面积）。"""

    nodes = np.linspace(0.0, 1.0, intervals + 1)
    inner = np.empty_like(nodes)
    outer = np.empty_like(nodes)
    inner[0] = 0.0
    inner[1:] = 0.5 * (nodes[:-1] + nodes[1:])
    outer[:-1] = inner[1:]
    outer[-1] = 1.0
    areas = np.pi * (outer**2 - inner**2)
    return areas / float(np.sum(areas))


@dataclass(frozen=True)
class TransientError:
    """L0-A 单个 (Bi, N) 算例的误差指标。"""

    biot: float
    intervals: int
    moisture_centre: float
    moisture_surface: float
    moisture_mean: float
    temperature_centre: float
    temperature_surface: float


def run_transient_case(biot: float, intervals: int) -> TransientError:
    """跑一个 (Bi, N) 算例，返回最坏相对误差（对 Fo 取最大）。"""

    model = build_benchmark_model(biot, intervals)
    radius = BENCHMARK_RADIUS_M
    times = np.asarray(TRANSIT_FOURIER, dtype=float) * radius**2 / BENCHMARK_DIFFUSIVITY
    node_count = model.node_count
    weights = _volume_weights(intervals)

    moisture_solution = _solve_bdf(
        _moisture_right_hand_side(model, BENCHMARK_AMBIENT_TEMPERATURE),
        np.full(node_count, BENCHMARK_INITIAL),
        times,
    )
    temperature_solution = _solve_bdf(
        _temperature_right_hand_side(model, BENCHMARK_AMBIENT),
        np.full(node_count, BENCHMARK_INITIAL_TEMPERATURE),
        times,
    )

    moisture_reference = CylinderRobinSolution(
        radius_m=radius,
        diffusivity_m2_s=BENCHMARK_DIFFUSIVITY,
        transfer_coefficient_m_s=model.mass_transfer_coefficient,
        initial_value=BENCHMARK_INITIAL,
        ambient_value=BENCHMARK_AMBIENT,
    )
    temperature_reference = CylinderRobinSolution(
        radius_m=radius,
        diffusivity_m2_s=BENCHMARK_ALPHA,
        # 温度方程写成 ∂T/∂t=α∇²T 后，Robin 条件的有效系数是 h/(ρc_p)：
        # -k∂T/∂r=h(T_R−T_air) 两边除以 ρc_p 即得。级数解的特征参数
        # 统一为 (有效系数)·R/α，这里因此不能再直接传 h。
        transfer_coefficient_m_s=(
            model.heat_transfer_coefficient
            / (
                BENCHMARK_DENSITY
                * BENCHMARK_HEAT_CAPACITY
            )
        ),
        initial_value=BENCHMARK_INITIAL_TEMPERATURE,
        ambient_value=BENCHMARK_AMBIENT_TEMPERATURE,
    )

    def worst(series: np.ndarray, reference, span: float, kind: str) -> float:
        errors = []
        for index, time_s in enumerate(times):
            if kind == "centre":
                exact = reference.centre(float(time_s))
            elif kind == "surface":
                exact = reference.surface(float(time_s))
            else:
                exact = reference.mean(float(time_s))
            errors.append(abs(float(series[index]) - exact) / span)
        return float(max(errors))

    mean_values = moisture_solution @ weights
    moisture_span = BENCHMARK_INITIAL - BENCHMARK_AMBIENT
    temperature_span = BENCHMARK_INITIAL_TEMPERATURE - BENCHMARK_AMBIENT_TEMPERATURE

    return TransientError(
        biot=biot,
        intervals=intervals,
        moisture_centre=worst(moisture_solution[:, 0], moisture_reference, moisture_span, "centre"),
        moisture_surface=worst(moisture_solution[:, -1], moisture_reference, moisture_span, "surface"),
        moisture_mean=worst(mean_values, moisture_reference, moisture_span, "mean"),
        temperature_centre=worst(
            temperature_solution[:, 0],
            temperature_reference,
            temperature_span,
            "centre",
        ),
        temperature_surface=worst(
            temperature_solution[:, -1],
            temperature_reference,
            temperature_span,
            "surface",
        ),
    )


def observed_order(coarse: float, fine: float, refinement: float) -> float:
    """由两级误差估计观测阶 p=log(e_coarse/e_fine)/log(refinement)。"""

    if coarse <= 0.0 or fine <= 0.0:
        return float("nan")
    return float(np.log(coarse / fine) / np.log(refinement))


# --------------------------------------------------------------------------
# L0-B：Heun 时间推进阶数
# --------------------------------------------------------------------------


def run_time_step_case(step: float, intervals: int = 160) -> np.ndarray:
    """固定 Bi=5、Fo=0.1，返回固定步长 Heun 的末态（默认生产档 N=160）。"""

    biot = 5.0
    model = build_benchmark_model(biot, intervals)
    radius = BENCHMARK_RADIUS_M
    end_time = 0.1 * radius**2 / BENCHMARK_DIFFUSIVITY
    rhs = _moisture_right_hand_side(model, BENCHMARK_AMBIENT_TEMPERATURE)
    result = integrate_heun(
        rhs=rhs,
        initial_state=np.full(model.node_count, BENCHMARK_INITIAL),
        end_time_s=end_time,
        dt_s=step,
        save_every_s=end_time,
        # 传入常数步长函数，既保持固定步长，又跳过"保存间隔必须整除"的检查。
        max_dt_fn=lambda _time_s, _state: step,
    )
    return np.asarray(result.state[-1], dtype=float)


def run_time_step_study() -> list[dict[str, float]]:
    """在 ``c·Δt_stable`` 四档步长上用**自收敛**（逐档差商）测量 Heun 观测阶。

    逐档差商 ``d_k=max|u_k−u_{k+1}|/span`` 中空间误差成对抵消（四档用同一
    空间网格），因此只反映时间积分误差；这避免了"拿另一个时间积分器作参照"
    时被参照解自身容差压住的地板效应。
    """

    model = build_benchmark_model(5.0, 160)
    stable = model.stable_time_step(
        0.0,
        np.concatenate(
            (
                np.full(model.node_count, BENCHMARK_INITIAL_TEMPERATURE),
                np.full(model.node_count, BENCHMARK_INITIAL),
            )
        ),
        safety_factor=1.0,
    )
    states: list[np.ndarray] = []
    rows: list[dict[str, float]] = []
    for fraction in (0.4, 0.2, 0.1, 0.05):
        step = fraction * stable
        states.append(run_time_step_case(step))
        rows.append(
            {
                "fraction_of_stable": fraction,
                "step_s": step,
            }
        )
    span = BENCHMARK_INITIAL - BENCHMARK_AMBIENT
    for index in range(len(rows) - 1):
        rows[index]["successive_difference"] = float(
            np.max(np.abs(states[index] - states[index + 1])) / span
        )
    rows[-1]["successive_difference"] = float("nan")
    for index in range(len(rows) - 2):
        rows[index]["order_to_next"] = observed_order(
            rows[index]["successive_difference"],
            rows[index + 1]["successive_difference"],
            2.0,
        )
    rows[-2]["order_to_next"] = float("nan")
    rows[-1]["order_to_next"] = float("nan")
    return rows


# --------------------------------------------------------------------------
# L0-C：变系数稳态通量恒等式
# --------------------------------------------------------------------------

STEADY_TEMPERATURE_C = 50.0
STEADY_INNER_RADIUS_M = 0.005
STEADY_OUTER_RADIUS_M = 0.020
STEADY_INNER_MOISTURE = 2.5
STEADY_OUTER_MOISTURE = 0.15
STEADY_GRIDS = (20, 40, 80)


def steady_annular_flux_reference(
    properties: PropertyModel,
    temperature_c: float,
    inner_radius_m: float,
    outer_radius_m: float,
    inner_moisture: float,
    outer_moisture: float,
):
    """返回 (Φ₁, Φ₂, 精确总通量 J, 精确剖面函数 C(r))。"""

    if inner_moisture <= outer_moisture:
        raise ValueError("内半径含水率必须高于外半径含水率")
    temperature = np.array([temperature_c], dtype=float)
    phi_inner = float(properties.flux_potential(temperature, np.array([inner_moisture]))[0])
    phi_outer = float(properties.flux_potential(temperature, np.array([outer_moisture]))[0])
    flux = 2.0 * np.pi * (phi_inner - phi_outer) / np.log(outer_radius_m / inner_radius_m)

    def target_potential(radii: np.ndarray) -> np.ndarray:
        return phi_inner + (phi_outer - phi_inner) * np.log(radii / inner_radius_m) / np.log(
            outer_radius_m / inner_radius_m
        )

    def profile(radii: np.ndarray) -> np.ndarray:
        radii = np.atleast_1d(np.asarray(radii, dtype=float))
        targets = target_potential(radii)

        def invert(target: float) -> float:
            return float(
                brentq(
                    lambda c: float(properties.flux_potential(temperature, np.array([c]))[0]) - target,
                    1.0e-8,
                    max(10.0, 2.0 * inner_moisture),
                    xtol=1e-14,
                    rtol=1e-15,
                )
            )

        return np.asarray([invert(float(target)) for target in targets], dtype=float)

    return phi_inner, phi_outer, float(flux), profile


def face_fluxes(
    strategy: InterfaceStrategy,
    properties: PropertyModel,
    temperature_c: float,
    radii: np.ndarray,
    moisture: np.ndarray,
) -> np.ndarray:
    """按指定界面策略计算逐面径向通量（单位长度，正为向外）。"""

    temperature = np.full_like(moisture, temperature_c)
    node_property = properties.diffusivity(temperature, moisture)
    faces = strategy.face_property(
        node_property=node_property,
        node_temperature=temperature,
        node_moisture=moisture,
        node_radii_m=radii,
        property_function=properties.diffusivity,
        flux_potential_function=properties.flux_potential,
    )
    spacing = radii[1] - radii[0]
    face_radius = 0.5 * (radii[:-1] + radii[1:])
    # 向外为正：q = -D dC/dr，故对 (C_{i+1}-C_i) 取负号。
    return -faces * (moisture[1:] - moisture[:-1]) / spacing * 2.0 * np.pi * face_radius


def run_steady_case(intervals: int) -> list[dict[str, float]]:
    """跑一个网格档位，返回各策略的通量误差与基尔霍夫恒等式残差。"""

    properties = Question23Properties()
    radii = np.linspace(STEADY_INNER_RADIUS_M, STEADY_OUTER_RADIUS_M, intervals + 1)
    _, _, exact_flux, profile = steady_annular_flux_reference(
        properties,
        STEADY_TEMPERATURE_C,
        STEADY_INNER_RADIUS_M,
        STEADY_OUTER_RADIUS_M,
        STEADY_INNER_MOISTURE,
        STEADY_OUTER_MOISTURE,
    )
    moisture = profile(radii)
    temperature = np.full_like(moisture, STEADY_TEMPERATURE_C)

    rows: list[dict[str, float]] = []
    for name in INTERFACE_STRATEGY_NAMES:
        strategy = resolve_interface_strategy(name)
        fluxes = face_fluxes(strategy, properties, STEADY_TEMPERATURE_C, radii, moisture)
        error = float(np.max(np.abs(fluxes - exact_flux)) / exact_flux)
        row = {
            "intervals": float(intervals),
            "strategy": name,
            "max_relative_flux_error": error,
        }
        if name == "kirchhoff":
            potentials = properties.flux_potential(temperature, moisture)
            faces = strategy.face_property(
                node_property=properties.diffusivity(temperature, moisture),
                node_temperature=temperature,
                node_moisture=moisture,
                node_radii_m=radii,
                property_function=properties.diffusivity,
                flux_potential_function=properties.flux_potential,
            )
            delta_potential = potentials[1:] - potentials[:-1]
            delta_moisture = moisture[1:] - moisture[:-1]
            residual = np.abs(faces * delta_moisture - delta_potential) / np.abs(delta_potential)
            row["kirchhoff_identity_residual"] = float(np.max(residual))
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# 报告与命令行
# --------------------------------------------------------------------------


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _markdown_table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |"]
    lines.append("|" + "|".join(["---"] * len(header)) + "|")
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def _plot_transient(rows: list[TransientError], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    figure, axis = plt.subplots(figsize=(7.2, 4.6))
    first_errors = None
    for biot in TRANSIT_BIOT:
        subset = [row for row in rows if row.biot == biot]
        grids = np.asarray([row.intervals for row in subset], dtype=float)
        errors = np.asarray([row.moisture_surface for row in subset], dtype=float)
        axis.loglog(grids, errors, marker="o", linewidth=2.0, label=f"Bi={biot:g}（表面）")
        if first_errors is None:
            first_errors = errors
    reference = np.asarray(TRANSIT_GRIDS, dtype=float)
    axis.loglog(
        reference,
        first_errors[0] * (reference / reference[0]) ** -2.0,
        linestyle="--",
        color="#555555",
        label="二阶参考斜率",
    )
    axis.set(
        title="L0-A 圆柱 Robin 基准的空间收敛（表面含水率）",
        xlabel="径向区间数 N",
        ylabel="相对误差（对 Fo 取最大）",
    )
    axis.grid(alpha=0.3, which="both", linestyle="--")
    axis.legend()
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _plot_steady(rows: list[dict[str, float]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    figure, axis = plt.subplots(figsize=(7.2, 4.6))
    for name in INTERFACE_STRATEGY_NAMES:
        subset = [row for row in rows if row["strategy"] == name]
        grids = np.asarray([row["intervals"] for row in subset], dtype=float)
        errors = np.asarray([row["max_relative_flux_error"] for row in subset], dtype=float)
        axis.loglog(grids, errors, marker="o", linewidth=2.0, label=name)
    axis.set(
        title="L0-C 变系数稳态通量的界面策略误差",
        xlabel="径向区间数 N",
        ylabel="最大相对通量误差",
    )
    axis.grid(alpha=0.3, which="both", linestyle="--")
    axis.legend()
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _evaluate_transient_criteria(rows: list[TransientError]) -> list[tuple[str, str, str]]:
    """按冻结判据逐条判定，返回 (判据, 实测, 结论)。"""

    results: list[tuple[str, str, str]] = []
    for biot in TRANSIT_BIOT:
        subset = {row.intervals: row for row in rows if row.biot == biot}
        for quantity in ("moisture_centre", "moisture_surface"):
            errors = {n: getattr(subset[n], quantity) for n in TRANSIT_GRIDS}
            orders = [
                observed_order(errors[TRANSIT_GRIDS[i]], errors[TRANSIT_GRIDS[i + 1]], 2.0)
                for i in range(len(TRANSIT_GRIDS) - 1)
            ]
            finest = orders[-1]
            average = float(np.mean(orders))
            passed = 1.8 <= finest <= 2.2 and 1.8 <= average <= 2.2
            results.append(
                (
                    f"L0-A 观测阶（Bi={biot:g}，{quantity}）",
                    f"最细 {finest:.3f}；平均 {average:.3f}",
                    "通过" if passed else "未通过",
                )
            )
        for grid, limit in ((40, 5.0e-3), (160, 5.0e-4)):
            worst = max(
                getattr(subset[grid], "moisture_centre"),
                getattr(subset[grid], "moisture_surface"),
                getattr(subset[grid], "temperature_centre"),
                getattr(subset[grid], "temperature_surface"),
            )
            results.append(
                (
                    f"L0-A 最大误差（Bi={biot:g}，N={grid}）",
                    f"{worst:.3e}（阈值 {limit:.1e}）",
                    "通过" if worst < limit else "未通过",
                )
            )
    return results


def _evaluate_time_criteria(rows: list[dict[str, float]]) -> list[tuple[str, str, str]]:
    orders = [
        row["order_to_next"]
        for row in rows
        if np.isfinite(row["order_to_next"])
    ]
    average = float(np.mean(orders))
    finest_difference = float(rows[-2]["successive_difference"])
    return [
        (
            "L0-B Heun 观测阶",
            f"逐级 {[f'{value:.3f}' for value in orders]}；平均 {average:.3f}",
            "通过" if 1.8 <= average <= 2.2 else "未通过",
        ),
        (
            "L0-B 最细两档相对差",
            f"{finest_difference:.3e}",
            "通过" if finest_difference < 1.0e-4 else "未通过",
        ),
    ]


def _evaluate_steady_criteria(rows: list[dict[str, float]]) -> list[tuple[str, str, str]]:
    results: list[tuple[str, str, str]] = []
    identity = max(
        row["kirchhoff_identity_residual"]
        for row in rows
        if "kirchhoff_identity_residual" in row
    )
    results.append(
        (
            "L0-C 基尔霍夫恒等式残差",
            f"{identity:.3e}（阈值 1e-12）",
            "通过" if identity < 1.0e-12 else "未通过",
        )
    )
    thresholds = {20: 2.0e-3, 40: 5.0e-4, 80: 1.5e-4}
    for grid, limit in thresholds.items():
        value = next(
            row["max_relative_flux_error"]
            for row in rows
            if row["strategy"] == "kirchhoff" and int(row["intervals"]) == grid
        )
        results.append(
            (
                f"L0-C kirchhoff 通量误差（N={grid}）",
                f"{value:.3e}（阈值 {limit:.1e}）",
                "通过" if value < limit else "未通过",
            )
        )
    for name in ("harmonic", "series"):
        value = next(
            row["max_relative_flux_error"]
            for row in rows
            if row["strategy"] == name and int(row["intervals"]) == 40
        )
        results.append(
            (
                f"L0-C {name} 未收敛证据（N=40）",
                f"{value:.3e}（阈值 >5e-2）",
                "通过" if value > 5.0e-2 else "未通过",
            )
        )
    return results


def run_all(output_dir: Path) -> Path:
    """执行 L0 三个算例，写出 CSV、图与判据报告，返回报告路径。"""

    rows: list[TransientError] = []
    for biot in TRANSIT_BIOT:
        for intervals in TRANSIT_GRIDS:
            rows.append(run_transient_case(biot, intervals))
    time_rows = run_time_step_study()
    steady_rows: list[dict[str, float]] = []
    for intervals in STEADY_GRIDS:
        steady_rows.extend(run_steady_case(intervals))

    _write_csv(
        output_dir / "analytic_transient.csv",
        [row.__dict__ for row in rows],
    )
    _write_csv(output_dir / "analytic_heun_time.csv", time_rows)
    _write_csv(output_dir / "analytic_steady_flux.csv", steady_rows)
    _plot_transient(rows, output_dir / "figures" / "L0A_transient_convergence.png")
    _plot_steady(steady_rows, output_dir / "figures" / "L0C_steady_flux.png")

    criteria = (
        _evaluate_transient_criteria(rows)
        + _evaluate_time_criteria(time_rows)
        + _evaluate_steady_criteria(steady_rows)
    )
    passed = sum(1 for _, _, verdict in criteria if verdict == "通过")

    parts: list[str] = []
    parts.append("# L0 解析基准验证报告")
    parts.append("")
    parts.append(
        "本报告由 `python -m aproblem.analytic` 自动生成，判据见 "
        "[`docs/06 §9`](../../../docs/06_代码与复现指南.md)（先于计算冻结）。"
    )
    parts.append("")
    parts.append(
        f"**判据通过 {passed}/{len(criteria)} 条。** 未通过项按原样记录，不调阈值。"
    )
    parts.append("")
    parts.append("## 1. 判据判定")
    parts.append("")
    parts.append(
        _markdown_table(
            ["判据", "实测", "结论"],
            [[name, value, verdict] for name, value, verdict in criteria],
        )
    )
    parts.append("")
    parts.append("## 2. L0-A 空间收敛（对 Fo 取最大相对误差）")
    parts.append("")
    transient_rows = [
        [
            f"{row.biot:g}",
            str(row.intervals),
            f"{row.moisture_centre:.3e}",
            f"{row.moisture_surface:.3e}",
            f"{row.moisture_mean:.3e}",
            f"{row.temperature_centre:.3e}",
            f"{row.temperature_surface:.3e}",
        ]
        for row in rows
    ]
    parts.append(
        _markdown_table(
            ["Bi", "N", "水分中心", "水分表面", "水分均值", "温度中心", "温度表面"],
            transient_rows,
        )
    )
    parts.append("")
    parts.append("## 3. L0-B Heun 时间阶（N=160，Bi=5，Fo=0.1）")
    parts.append("")
    parts.append(
        _markdown_table(
            ["稳定步长占比", "Δt（s）", "与下一档的最大相对差", "到下一档观测阶"],
            [
                [
                    f"{row['fraction_of_stable']:g}",
                    f"{row['step_s']:.3e}",
                    (
                        "—"
                        if not np.isfinite(row["successive_difference"])
                        else f"{row['successive_difference']:.3e}"
                    ),
                    "—" if not np.isfinite(row["order_to_next"]) else f"{row['order_to_next']:.3f}",
                ]
                for row in time_rows
            ],
        )
    )
    parts.append("")
    parts.append("## 4. L0-C 变系数稳态通量（Φ 沿 ln r 线性）")
    parts.append("")
    parts.append(
        _markdown_table(
            ["N", "界面策略", "最大相对通量误差", "基尔霍夫恒等式残差"],
            [
                [
                    str(int(row["intervals"])),
                    str(row["strategy"]),
                    f"{row['max_relative_flux_error']:.3e}",
                    (
                        f"{row['kirchhoff_identity_residual']:.3e}"
                        if "kirchhoff_identity_residual" in row
                        else "—"
                    ),
                ]
                for row in steady_rows
            ],
        )
    )
    parts.append("")
    parts.append("## 5. 图")
    parts.append("")
    parts.append("- `figures/L0A_transient_convergence.png`")
    parts.append("- `figures/L0C_steady_flux.png`")

    report_path = output_dir / "analytic_report.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return report_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="L0 解析基准验证")
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = ProjectPaths.discover()
    output_dir = (args.output_dir or paths.project_root / "outputs" / "study" / "analytic").resolve()
    report_path = run_all(output_dir)
    print(f"L0 报告：{report_path}")


if __name__ == "__main__":
    main()
