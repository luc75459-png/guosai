"""L1 外部实测对照（validation）：Teleken (2025) 虾体热风干燥。

本模块只做**外部验证**，不参与四问正式结果生成，也不修改生产模块。
设计见预注册判据 ``docs/06 §9.4``：

* L1-A 预测式：直接用论文的 ``D``、``h_m``、几何与实测初值预测 ``X_db(t)``，
  与数字化实测点比较（零拟合参数）；
* L1-B 反演式：固定几何与初值，反演 ``D`` 并检验其是否落在论文报告值附近；
* 敏感性：``h_m×0.1/×10``、环境平衡值、初值、几何映射；
* 端面诊断：用有限圆柱乘积解（Newman）的轴向因子给出解析界。

用法：``python -m aproblem.external``，产物写入
``outputs/study/external_validation/``。
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

from .config import ProjectPaths
from .grid import RadialGrid
from .interfaces import resolve_interface_strategy
from .model import DryingModel
from .physics import MaterialProperties, PropertyModel


DATA_RELATIVE = Path("data") / "external" / "teleken2025"
INTEGRATION_RTOL = 1.0e-8
INTEGRATION_ATOL = 1.0e-10


# --------------------------------------------------------------------------
# 数据读取
# --------------------------------------------------------------------------


def load_conditions(path: Path) -> dict:
    """读取 ``conditions.json``。"""

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_measured(path: Path) -> dict[str, np.ndarray]:
    """读取数字化实测 CSV，返回按列组织的字典（含 ``included`` 掩码）。"""

    times: list[float] = []
    values: list[float] = []
    low: list[float] = []
    high: list[float] = []
    sigma: list[float] = []
    included: list[bool] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            times.append(float(row["time_min"]))
            values.append(float(row["x_db"]))
            low.append(float(row["err_low"]))
            high.append(float(row["err_high"]))
            sigma.append(float(row["digitization_sigma"]))
            included.append(str(row["included"]).strip().lower() == "true")
    return {
        "time_min": np.asarray(times, dtype=float),
        "x_db": np.asarray(values, dtype=float),
        "err_low": np.asarray(low, dtype=float),
        "err_high": np.asarray(high, dtype=float),
        "sigma": np.asarray(sigma, dtype=float),
        "included": np.asarray(included, dtype=bool),
    }


# --------------------------------------------------------------------------
# 物性模型（论文 Table 1 + 式 (7)(8)(9)）
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class TelekenProperties:
    """论文给出的虾体物性；水分扩散系数按工况取常数并固定。"""

    diffusivity_m2_s: float

    @staticmethod
    def _water_density(temperature_c: float) -> float:
        return 997.18 + 3.1439e-3 * temperature_c - 3.7574e-3 * temperature_c**2

    @staticmethod
    def _protein_density(temperature_c: float) -> float:
        return 1329.9 - 0.51840 * temperature_c

    @staticmethod
    def _water_heat_capacity(temperature_c: float) -> float:
        return 4176.2 - 0.0909 * temperature_c + 5.4731e-3 * temperature_c**2

    @staticmethod
    def _protein_heat_capacity(temperature_c: float) -> float:
        return 2008.2 + 1.2089 * temperature_c - 1.3129e-3 * temperature_c**2

    @staticmethod
    def _water_conductivity(temperature_c: float) -> float:
        return 0.57109 + 1.762e-3 * temperature_c - 6.7036e-6 * temperature_c**2

    @staticmethod
    def _protein_conductivity(temperature_c: float) -> float:
        return 0.17881 + 1.958e-3 * temperature_c - 2.7178e-6 * temperature_c**2

    def density(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        temperature = np.asarray(temperature_c, dtype=float)
        c = np.maximum(np.asarray(moisture, dtype=float), 1.0e-8)
        x_w = c / (1.0 + c)
        x_p = 1.0 - x_w
        return 1.0 / (
            x_w / self._water_density(temperature)
            + x_p / self._protein_density(temperature)
        )

    def heat_capacity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        temperature = np.asarray(temperature_c, dtype=float)
        c = np.maximum(np.asarray(moisture, dtype=float), 1.0e-8)
        x_w = c / (1.0 + c)
        return x_w * self._water_heat_capacity(temperature) + (1.0 - x_w) * self._protein_heat_capacity(
            temperature
        )

    def conductivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        temperature = np.asarray(temperature_c, dtype=float)
        c = np.maximum(np.asarray(moisture, dtype=float), 1.0e-8)
        x_w = c / (1.0 + c)
        x_p = 1.0 - x_w
        density = self.density(temperature, c)
        k_w = self._water_conductivity(temperature)
        k_p = self._protein_conductivity(temperature)
        v_w = x_w * density / self._water_density(temperature)
        v_p = x_p * density / self._protein_density(temperature)
        return 0.5 * (v_w * k_w + v_p * k_p) + 0.5 / (v_w / k_w + v_p / k_p)

    def diffusivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        reference = np.asarray(moisture, dtype=float)
        return np.full(np.shape(reference), float(self.diffusivity_m2_s))

    def flux_potential(
        self,
        temperature_c: np.ndarray,
        moisture: np.ndarray,
    ) -> np.ndarray:
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


# --------------------------------------------------------------------------
# 模型组装与积分
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class CaseSettings:
    """一次 L1 算例的全部映射参数（默认取 conditions.json 的主算例值）。"""

    key: str
    temperature_c: float
    ambient_moisture: float
    initial_moisture: float
    initial_temperature_c: float
    radius_m: float
    diffusivity_m2_s: float
    mass_transfer_coefficient_m_s: float
    heat_transfer_coefficient_w_m2k: float
    intervals: int = 40
    end_time_s: float = 10800.0


def build_model(settings: CaseSettings) -> DryingModel:
    """组装一维仅侧面模型（与交付口径一致）。"""

    return DryingModel(
        grid=RadialGrid(settings.intervals),
        properties=TelekenProperties(settings.diffusivity_m2_s),
        environment=lambda _time: (settings.temperature_c, settings.ambient_moisture),
        radius=lambda _time: settings.radius_m,
        heat_transfer_coefficient=settings.heat_transfer_coefficient_w_m2k,
        mass_transfer_coefficient=settings.mass_transfer_coefficient_m_s,
        include_end_faces=False,
        interface_strategy=resolve_interface_strategy("kirchhoff"),
    )


@dataclass(frozen=True)
class CaseTrajectory:
    """采样时刻的预测值（体积平均与中心/表面、温度）。"""

    time_s: np.ndarray
    moisture_mean: np.ndarray
    moisture_centre: np.ndarray
    moisture_surface: np.ndarray
    temperature_centre: np.ndarray
    temperature_mean: np.ndarray


def run_case(settings: CaseSettings) -> CaseTrajectory:
    """运行一个 L1-A 算例，返回采样时刻的预测轨迹。"""

    model = build_model(settings)
    node_count = model.node_count
    initial_state = np.concatenate(
        (
            np.full(node_count, settings.initial_temperature_c),
            np.full(node_count, settings.initial_moisture),
        )
    )
    samples = np.linspace(0.0, settings.end_time_s, 13)
    solution = solve_ivp(
        model.rhs,
        (0.0, settings.end_time_s),
        initial_state,
        method="BDF",
        rtol=INTEGRATION_RTOL,
        atol=INTEGRATION_ATOL,
        t_eval=samples,
    )
    if not solution.success:
        raise RuntimeError(f"L1 积分失败：{solution.message}")

    values = np.asarray(solution.y, dtype=float).T
    temperature = values[:, :node_count]
    moisture = values[:, node_count:]
    weights = model.grid.geometry(settings.radius_m).volumes_per_length_m2
    weights = weights / float(np.sum(weights))
    return CaseTrajectory(
        time_s=np.asarray(solution.t, dtype=float),
        moisture_mean=moisture @ weights,
        moisture_centre=moisture[:, 0],
        moisture_surface=moisture[:, -1],
        temperature_centre=temperature[:, 0],
        temperature_mean=temperature @ weights,
    )


def sample_at(trajectory: CaseTrajectory, times_s: np.ndarray) -> CaseTrajectory:
    """把轨迹线性插值到给定时刻。"""

    times_s = np.asarray(times_s, dtype=float)
    return CaseTrajectory(
        time_s=times_s,
        moisture_mean=np.interp(times_s, trajectory.time_s, trajectory.moisture_mean),
        moisture_centre=np.interp(times_s, trajectory.time_s, trajectory.moisture_centre),
        moisture_surface=np.interp(times_s, trajectory.time_s, trajectory.moisture_surface),
        temperature_centre=np.interp(times_s, trajectory.time_s, trajectory.temperature_centre),
        temperature_mean=np.interp(times_s, trajectory.time_s, trajectory.temperature_mean),
    )


# --------------------------------------------------------------------------
# 指标与判据
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ExternalMetrics:
    rmse: float
    r2: float
    bias: float
    max_abs: float


def compute_metrics(measured: np.ndarray, predicted: np.ndarray) -> ExternalMetrics:
    """计算 RMSE、R²、偏差与最大绝对偏差。"""

    measured = np.asarray(measured, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    residual = predicted - measured
    total = float(np.sum((measured - measured.mean()) ** 2))
    return ExternalMetrics(
        rmse=float(np.sqrt(np.mean(residual**2))),
        r2=float(1.0 - np.sum(residual**2) / total) if total > 0 else float("nan"),
        bias=float(np.mean(residual)),
        max_abs=float(np.max(np.abs(residual))),
    )


def time_to_value(times_min: np.ndarray, values: np.ndarray, target: float) -> float:
    """在单调下降的曲线上线性插值出达到 ``target`` 的时刻（min）。"""

    order = np.argsort(values)
    return float(np.interp(target, values[order], np.asarray(times_min, dtype=float)[order]))


def newman_slab_mean_factor(
    diffusivity_m2_s: float,
    length_m: float,
    time_s: float,
    terms: int = 60,
) -> float:
    """有限平板（两端 Dirichlet）体积平均的乘积解轴向因子。

    圆柱的平均含水率约为「无限长圆柱」与「平板」两个因子的乘积；把两端
    忽略掉会高估平均含水率 ``1 - factor``。
    """

    fourier = diffusivity_m2_s * time_s / length_m**2
    factor = 0.0
    for index in range(terms):
        eigenvalue = (2 * index + 1) * np.pi / 2.0
        factor += 2.0 / eigenvalue**2 * np.exp(-(eigenvalue**2) * fourier)
    return float(factor)


def fit_diffusivity(
    measured_times_min: np.ndarray,
    measured_values: np.ndarray,
    base: CaseSettings,
    lower: float = 1.0e-11,
    upper: float = 1.0e-7,
) -> tuple[float, float]:
    """在固定几何/初值/边界系数下反演水分扩散系数。

    返回 ``(反演 D, 该 D 处的 RMSE)``。用对数坐标的有界标量优化，目标函数是
    体积平均含水率对实测点的 RMSE。
    """

    times_s = np.asarray(measured_times_min, dtype=float) * 60.0
    target = np.asarray(measured_values, dtype=float)

    def objective(log_diffusivity: float) -> float:
        settings = CaseSettings(
            **{
                **base.__dict__,
                "diffusivity_m2_s": float(np.exp(log_diffusivity)),
            }
        )
        trajectory = sample_at(run_case(settings), times_s)
        return float(np.sqrt(np.mean((trajectory.moisture_mean - target) ** 2)))

    result = minimize_scalar(
        objective,
        bounds=(float(np.log(lower)), float(np.log(upper))),
        method="bounded",
        options={"xatol": 1.0e-3},
    )
    return float(np.exp(result.x)), float(result.fun)


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


def _plot_moisture(cases: dict[str, tuple[dict[str, np.ndarray], CaseTrajectory]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    figure, axes = plt.subplots(1, len(cases), figsize=(6.0 * len(cases), 4.6))
    for axis, (key, (measured, trajectory)) in zip(np.atleast_1d(axes), cases.items()):
        axis.errorbar(
            measured["time_min"],
            measured["x_db"],
            yerr=[measured["err_low"], measured["err_high"]],
            fmt="o",
            color="#0072B2",
            capsize=3,
            label="实测（数字化，含误差棒）",
        )
        axis.plot(
            trajectory.time_s / 60.0,
            trajectory.moisture_mean,
            color="#D55E00",
            linewidth=2.0,
            label="本模型预测（零拟合参数）",
        )
        axis.set(
            title=f"Teleken {key}：干基含水率",
            xlabel="时间（min）",
            ylabel="X_db（kg/kg）",
        )
        axis.grid(alpha=0.3, linestyle="--")
        axis.legend(fontsize=9)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _plot_temperature(cases: dict[str, tuple[dict[str, np.ndarray], CaseTrajectory]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    figure, axis = plt.subplots(figsize=(7.6, 4.6))
    for key, (_, trajectory) in cases.items():
        axis.plot(
            trajectory.time_s / 60.0,
            trajectory.temperature_centre,
            linewidth=2.0,
            label=f"{key} 中心温度（本模型）",
        )
    axis.axhline(46.0, color="#009E73", linestyle=":", linewidth=1.5, label="论文 60 ℃ 中心平台 46 ℃")
    axis.axhline(53.0, color="#CC79A7", linestyle=":", linewidth=1.5, label="论文 70 ℃ 中心平台 53 ℃")
    axis.set(
        title="温度诊断（本模型不含蒸发潜热）",
        xlabel="时间（min）",
        ylabel="温度（℃）",
    )
    axis.grid(alpha=0.3, linestyle="--")
    axis.legend(fontsize=9)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def run_all(output_dir: Path, data_dir: Path) -> Path:
    """执行 L1-A / L1-B 与敏感性，写出产物并返回报告路径。"""

    conditions = load_conditions(data_dir / "conditions.json")
    geometry = conditions["geometry"]
    sampling_min = np.asarray(conditions["sampling_times_min"], dtype=float)

    cases: dict[str, tuple[dict[str, np.ndarray], CaseTrajectory]] = {}
    summary_rows: list[dict[str, object]] = []
    curve_rows: list[dict[str, object]] = []

    for key, case in conditions["cases"].items():
        measured = load_measured(data_dir / f"figure4_{key}.csv")
        settings = CaseSettings(
            key=key,
            temperature_c=float(case["temperature_c"]),
            ambient_moisture=0.0,
            initial_moisture=float(case["initial_moisture_measured"]),
            initial_temperature_c=float(conditions["initial_temperature_c"]),
            radius_m=float(geometry["radius_from_diameter_m"]),
            diffusivity_m2_s=float(case["diffusivity_m2_s"]),
            mass_transfer_coefficient_m_s=float(case["mass_transfer_coefficient_m_s"]),
            heat_transfer_coefficient_w_m2k=float(case["heat_transfer_coefficient_w_m2k"]),
            end_time_s=float(conditions["end_time_min"]) * 60.0,
        )
        trajectory = sample_at(run_case(settings), sampling_min * 60.0)
        cases[key] = (measured, trajectory)

        mask = measured["included"]
        metrics = compute_metrics(measured["x_db"][mask], trajectory.moisture_mean[mask])
        falling = measured["time_min"] >= 60.0
        falling_bias = float(
            np.mean(
                (trajectory.moisture_mean - measured["x_db"])[
                    mask & falling
                ]
            )
        )
        measured_time = time_to_value(measured["time_min"], measured["x_db"], 1.0)
        predicted_time = time_to_value(
            trajectory.time_s / 60.0,
            trajectory.moisture_mean,
            1.0,
        )
        plateau = float(np.mean(trajectory.temperature_centre[-3:]))

        summary_rows.append(
            {
                "case": key,
                "rmse_moisture": round(metrics.rmse, 4),
                "r2_moisture": round(metrics.r2, 4),
                "bias_all": round(metrics.bias, 4),
                "bias_falling": round(falling_bias, 4),
                "max_abs": round(metrics.max_abs, 4),
                "time_to_1_measured_min": round(measured_time, 1),
                "time_to_1_predicted_min": round(predicted_time, 1),
                "time_to_1_relative_error": round(
                    (predicted_time - measured_time) / measured_time,
                    4,
                ),
                "centre_plateau_c": round(plateau, 2),
                "reported_plateau_c": float(case["centre_temperature_plateau_c"]),
                "rmse_temperature_c": round(
                    abs(plateau - float(case["centre_temperature_plateau_c"])),
                    2,
                ),
                "reported_rmse_moisture": float(case["reported"]["rmse_moisture_kg_kg"]),
                "reported_r2_moisture": float(case["reported"]["r2_moisture"]),
            }
        )
        for index, time_min in enumerate(measured["time_min"]):
            curve_rows.append(
                {
                    "case": key,
                    "time_min": float(time_min),
                    "measured": float(measured["x_db"][index]),
                    "predicted": round(float(trajectory.moisture_mean[index]), 4),
                    "residual": round(
                        float(trajectory.moisture_mean[index] - measured["x_db"][index]),
                        4,
                    ),
                    "measured_err_low": float(measured["err_low"][index]),
                    "measured_err_high": float(measured["err_high"][index]),
                    "digitization_sigma": float(measured["sigma"][index]),
                    "included": bool(measured["included"][index]),
                }
            )

    # ---- 敏感性 ------------------------------------------------------
    sensitivity_rows: list[dict[str, object]] = []
    for key, case in conditions["cases"].items():
        measured = load_measured(data_dir / f"figure4_{key}.csv")
        variants = {
            "主线（h_m 原值，X_eq=0，实测初值，R=直径换算）": {},
            "h_m×0.1": {"mass_transfer_coefficient_m_s": float(case["mass_transfer_coefficient_m_s"]) * 0.1},
            "h_m×10": {"mass_transfer_coefficient_m_s": float(case["mass_transfer_coefficient_m_s"]) * 10.0},
            "环境值 X_eq=0.02": {"ambient_moisture": 0.02},
            "环境值 X_eq=0.05": {"ambient_moisture": 0.05},
            "预注册初值 X0=4.0": {"initial_moisture": float(case["initial_moisture_plan"])},
            "体积等效半径 R=4.92 mm": {"radius_m": float(geometry["radius_volume_equivalent_m"])},
        }
        for label, override in variants.items():
            parameters = dict(
                key=key,
                temperature_c=float(case["temperature_c"]),
                ambient_moisture=0.0,
                initial_moisture=float(case["initial_moisture_measured"]),
                initial_temperature_c=float(conditions["initial_temperature_c"]),
                radius_m=float(geometry["radius_from_diameter_m"]),
                diffusivity_m2_s=float(case["diffusivity_m2_s"]),
                mass_transfer_coefficient_m_s=float(case["mass_transfer_coefficient_m_s"]),
                heat_transfer_coefficient_w_m2k=float(case["heat_transfer_coefficient_w_m2k"]),
                end_time_s=float(conditions["end_time_min"]) * 60.0,
            )
            parameters.update(override)
            if label.startswith("预注册初值"):
                parameters["ambient_moisture"] = 0.0
            trajectory = sample_at(run_case(CaseSettings(**parameters)), sampling_min * 60.0)
            metrics = compute_metrics(measured["x_db"][measured["included"]], trajectory.moisture_mean[measured["included"]])
            sensitivity_rows.append(
                {
                    "case": key,
                    "variant": label,
                    "rmse_moisture": round(metrics.rmse, 4),
                    "r2_moisture": round(metrics.r2, 4),
                    "bias": round(metrics.bias, 4),
                }
            )

    # ---- L1-B 反演 ----------------------------------------------------
    inversion_rows: list[dict[str, object]] = []
    for key, case in conditions["cases"].items():
        measured = load_measured(data_dir / f"figure4_{key}.csv")
        mask = measured["included"]

        base = CaseSettings(
            key=key,
            temperature_c=float(case["temperature_c"]),
            ambient_moisture=0.0,
            initial_moisture=float(case["initial_moisture_measured"]),
            initial_temperature_c=float(conditions["initial_temperature_c"]),
            radius_m=float(geometry["radius_from_diameter_m"]),
            diffusivity_m2_s=float(case["diffusivity_m2_s"]),
            mass_transfer_coefficient_m_s=float(case["mass_transfer_coefficient_m_s"]),
            heat_transfer_coefficient_w_m2k=float(case["heat_transfer_coefficient_w_m2k"]),
            end_time_s=float(conditions["end_time_min"]) * 60.0,
        )
        fitted, fitted_rmse = fit_diffusivity(
            measured["time_min"][mask],
            measured["x_db"][mask],
            base,
        )
        reported = float(case["diffusivity_m2_s"])
        inversion_rows.append(
            {
                "case": key,
                "reported_diffusivity": reported,
                "fitted_diffusivity": round(fitted, 12),
                "ratio_fitted_over_reported": round(fitted / reported, 3),
                "rmse_at_fit": round(fitted_rmse, 4),
            }
        )

    # ---- 端面解析界 ----------------------------------------------------
    end_face_rows: list[dict[str, object]] = []
    for key, case in conditions["cases"].items():
        factor = newman_slab_mean_factor(
            float(case["diffusivity_m2_s"]),
            float(geometry["length_m"]),
            float(conditions["end_time_min"]) * 60.0,
        )
        end_face_rows.append(
            {
                "case": key,
                "slab_mean_factor_at_180min": round(factor, 4),
                "overestimate_if_ends_ignored": round(1.0 - factor, 4),
            }
        )

    _write_csv(output_dir / "teleken_summary.csv", summary_rows)
    _write_csv(output_dir / "teleken_curves.csv", curve_rows)
    _write_csv(output_dir / "teleken_sensitivity.csv", sensitivity_rows)
    _write_csv(output_dir / "teleken_inversion.csv", inversion_rows)
    _write_csv(output_dir / "teleken_end_faces.csv", end_face_rows)
    _plot_moisture(cases, output_dir / "figures" / "teleken_moisture.png")
    _plot_temperature(cases, output_dir / "figures" / "teleken_temperature.png")

    # ---- 判据判定 ------------------------------------------------------
    criteria: list[tuple[str, str, str]] = []
    thresholds = {"60C": 0.35, "70C": 0.30}
    for row in summary_rows:
        key = str(row["case"])
        criteria.append(
            (
                f"L1 水分 R²（{key}）",
                f"{row['r2_moisture']:.4f}（阈值 ≥0.90）",
                "通过" if float(row["r2_moisture"]) >= 0.90 else "未通过",
            )
        )
        limit = thresholds[key]
        criteria.append(
            (
                f"L1 水分 RMSE（{key}）",
                f"{row['rmse_moisture']:.4f}（阈值 ≤{limit}）",
                "通过" if float(row["rmse_moisture"]) <= limit else "未通过",
            )
        )
        criteria.append(
            (
                f"L1 降速段偏差（{key}）",
                f"{row['bias_falling']:+.4f}（阈值 |·|≤0.10）",
                "通过" if abs(float(row["bias_falling"])) <= 0.10 else "未通过",
            )
        )
        criteria.append(
            (
                f"L1 到达 X=1.0 时间（{key}）",
                f"{row['time_to_1_relative_error']:+.2%}（阈值 |·|≤20%）",
                "通过" if abs(float(row["time_to_1_relative_error"])) <= 0.20 else "未通过",
            )
        )
        deviation = abs(float(row["centre_plateau_c"]) - float(row["reported_plateau_c"]))
        criteria.append(
            (
                f"L1 中心温度平台（{key}，诊断项）",
                f"{row['centre_plateau_c']:.2f} ℃ vs {row['reported_plateau_c']:.1f} ℃（偏差 {deviation:.2f}，阈值 ≤5）",
                "通过" if deviation <= 5.0 else "未通过",
            )
        )
    for row in inversion_rows:
        ratio = float(row["ratio_fitted_over_reported"])
        criteria.append(
            (
                f"L1-B 反演 D 比值（{row['case']}）",
                f"{ratio:.3f}（阈值 0.5–2.0）",
                "通过" if 0.5 <= ratio <= 2.0 else "未通过",
            )
        )

    passed = sum(1 for _, _, verdict in criteria if verdict == "通过")
    parts: list[str] = []
    parts.append("# L1 外部实测对照报告（Teleken 2025）")
    parts.append("")
    parts.append(
        "本报告由 `python -m aproblem.external` 自动生成；映射规则与判据见 "
        "[`docs/06 §9.4`](../../../docs/06_代码与复现指南.md)（先于计算冻结，含执行前偏差记录）。"
    )
    parts.append("")
    parts.append(f"**判据通过 {passed}/{len(criteria)} 条。** 未通过项按原样记录，不调参。")
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
    parts.append("## 2. 零拟合参数预测 vs 数字化实测")
    parts.append("")
    parts.append(
        _markdown_table(
            [
                "工况",
                "RMSE_X",
                "R²_X",
                "降速段偏差",
                "最大偏差",
                "到达 X=1 实测/预测(min)",
                "中心平台(℃)",
            ],
            [
                [
                    str(row["case"]),
                    f"{row['rmse_moisture']:.4f}",
                    f"{row['r2_moisture']:.4f}",
                    f"{row['bias_falling']:+.4f}",
                    f"{row['max_abs']:.4f}",
                    f"{row['time_to_1_measured_min']:.0f}/{row['time_to_1_predicted_min']:.0f}",
                    f"{row['centre_plateau_c']:.2f}",
                ]
                for row in summary_rows
            ],
        )
    )
    parts.append("")
    parts.append(
        f"论文自报的 RMSE_X 为 0.22（60 ℃）/0.16（70 ℃）；本次数字化不确定度约 "
        f"{float(curve_rows[0]['digitization_sigma']):.3f} kg/kg，是残差的主导下限。"
    )
    parts.append("")
    parts.append("## 3. 敏感性（RMSE_X，kg/kg）")
    parts.append("")
    parts.append(
        _markdown_table(
            ["工况", "变体", "RMSE_X", "R²_X", "偏差"],
            [
                [
                    str(row["case"]),
                    str(row["variant"]),
                    f"{row['rmse_moisture']:.4f}",
                    f"{row['r2_moisture']:.4f}",
                    f"{row['bias']:+.4f}",
                ]
                for row in sensitivity_rows
            ],
        )
    )
    parts.append("")
    parts.append("## 4. L1-B 反演 D")
    parts.append("")
    parts.append(
        _markdown_table(
            ["工况", "论文 D（m²/s）", "反演 D（m²/s）", "比值", "拟合 RMSE"],
            [
                [
                    str(row["case"]),
                    f"{row['reported_diffusivity']:.3e}",
                    f"{row['fitted_diffusivity']:.3e}",
                    f"{row['ratio_fitted_over_reported']:.3f}",
                    f"{row['rmse_at_fit']:.4f}",
                ]
                for row in inversion_rows
            ],
        )
    )
    parts.append("")
    parts.append("## 5. 端面影响的解析界（Newman 乘积解）")
    parts.append("")
    parts.append(
        _markdown_table(
            ["工况", "180 min 平板平均因子", "忽略两端的高估幅度"],
            [
                [
                    str(row["case"]),
                    f"{row['slab_mean_factor_at_180min']:.4f}",
                    f"{row['overestimate_if_ends_ignored']:.2%}",
                ]
                for row in end_face_rows
            ],
        )
    )
    parts.append("")
    parts.append(
        "有限圆柱的平均含水率约为「无限长圆柱」与「平板」两个因子的乘积；"
        "上表给出仅算侧面时的最大高估幅度，用于替代预注册的二维轴对称诊断"
        "（理由见 §9.4 执行前偏差记录）。"
    )
    parts.append("")
    parts.append("## 6. 图")
    parts.append("")
    parts.append("- `figures/teleken_moisture.png`")
    parts.append("- `figures/teleken_temperature.png`")

    report_path = output_dir / "external_validation.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return report_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="L1 外部实测对照（Teleken 2025）")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--data-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = ProjectPaths.discover()
    data_dir = (args.data_dir or paths.project_root / DATA_RELATIVE).resolve()
    output_dir = (
        args.output_dir or paths.project_root / "outputs" / "study" / "external_validation"
    ).resolve()
    report_path = run_all(output_dir, data_dir)
    print(f"L1 报告：{report_path}")


if __name__ == "__main__":
    main()
