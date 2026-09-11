"""生成温度场与含水率场的中文径向剖面图。"""

from __future__ import annotations

from pathlib import Path

import matplotlib

# 使用无界面的 Agg 后端，保证比赛电脑或终端没有 Tk 图形环境时仍能保存图片。
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .integrator import SimulationResult
from .model import DryingModel
from .outputs import split_result


def _configure_chinese_font() -> None:
    """设置常见中文字体回退顺序，并保证坐标轴负号正常显示。"""

    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False


def plot_final_profiles(
    path: Path,
    model: DryingModel,
    result: SimulationResult,
) -> None:
    """绘制最终时刻的径向温度、含水率剖面并保存为 PNG。

    两条数据线使用不同颜色：橙红色表示温度，蓝色表示干基含水率。
    图题、坐标轴、图例以及文字注释均使用中文。
    """

    _configure_chinese_font()
    temperature, moisture = split_result(model, result)
    radius_cm = model.grid.geometry(model.radius(result.time_s[-1])).nodes_m * 100.0
    final_time_h = result.time_s[-1] / 3600.0
    end_face_text = "含端面等效修正" if model.include_end_faces else "仅考虑圆柱侧面"

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(
        radius_cm,
        temperature[-1],
        marker="o",
        color="#D55E00",
        linewidth=2.0,
        markersize=4,
        label="温度",
    )
    axes[0].set(
        title="最终径向温度分布",
        xlabel="距药材中心的距离（cm）",
        ylabel="温度（℃）",
    )
    axes[0].legend()
    axes[0].grid(alpha=0.3, linestyle="--")

    axes[1].plot(
        radius_cm,
        moisture[-1],
        marker="s",
        color="#0072B2",
        linewidth=2.0,
        markersize=4,
        label="干基含水率",
    )
    axes[1].set(
        title="最终径向含水率分布",
        xlabel="距药材中心的距离（cm）",
        ylabel="干基含水率（kg/kg）",
    )
    axes[1].legend()
    axes[1].grid(alpha=0.3, linestyle="--")

    figure.suptitle(
        f"药材最终径向温湿分布（t = {final_time_h:.3f} h，{end_face_text}）"
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_centre_surface_history(
    path: Path,
    time_hours: np.ndarray,
    centre_moisture: np.ndarray,
    surface_moisture: np.ndarray,
    radius_cm: np.ndarray,
    event_hours: float | None = None,
    threshold: float = 0.15,
) -> None:
    """绘制中心与表面含水率随时间变化，并附上半径收缩曲线。"""

    _configure_chinese_font()
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    axes[0].plot(
        time_hours,
        centre_moisture,
        color="#0072B2",
        linewidth=2.0,
        label="中心含水率",
    )
    axes[0].plot(
        time_hours,
        surface_moisture,
        color="#D55E00",
        linewidth=2.0,
        label="表面含水率",
    )
    axes[0].axhline(
        threshold,
        color="#555555",
        linestyle="--",
        linewidth=1.2,
        label=f"达标阈值 {threshold} kg/kg",
    )
    if event_hours is not None:
        axes[0].axvline(
            event_hours,
            color="#009E73",
            linestyle=":",
            linewidth=1.5,
            label=f"达标时刻 {event_hours:.3f} h",
        )
    axes[0].set(
        title="中心与表面含水率演化",
        xlabel="时间（h）",
        ylabel="干基含水率（kg/kg）",
    )
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.3, linestyle="--")

    axes[1].plot(
        time_hours,
        radius_cm,
        color="#CC79A7",
        linewidth=2.0,
        label="药材外半径",
    )
    axes[1].set(
        title="药材半径收缩过程",
        xlabel="时间（h）",
        ylabel="半径（cm）",
    )
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3, linestyle="--")

    figure.suptitle("问题4 收缩条件下的水分演化与几何变化")
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_axial_moisture_map(
    path: Path,
    moisture_final: np.ndarray,
    radius_cm: float,
    length_cm: float,
    event_hours: float,
    axial_spread: float,
) -> None:
    """绘制二维轴对称模型末态的轴向—径向含水率分布热图。"""

    _configure_chinese_font()
    figure, axis = plt.subplots(figsize=(7.6, 4.6))
    radial_cm = np.linspace(0.0, radius_cm, moisture_final.shape[0])
    axial_cm = np.linspace(0.0, length_cm, moisture_final.shape[1])
    mesh = axis.pcolormesh(
        axial_cm,
        radial_cm,
        moisture_final,
        shading="nearest",
        cmap="viridis",
    )
    bar = figure.colorbar(mesh, ax=axis)
    bar.set_label("干基含水率（kg/kg）")
    axis.set(
        title=(
            "二维轴对称末态含水率分布"
            f"（达标时刻 {event_hours:.2f} h，轴向极差 {axial_spread:.4f}）"
        ),
        xlabel="轴向位置 z（cm）",
        ylabel="距中心距离 r（cm）",
    )
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_interface_convergence(
    path: Path,
    grids: tuple[int, ...],
    series: dict[str, list[float]],
    reference_hours: float | None = None,
) -> None:
    """绘制不同界面策略的网格收敛曲线。"""

    _configure_chinese_font()
    colors = {
        "midpoint": "#0072B2",
        "harmonic": "#D55E00",
        "series": "#009E73",
        "kirchhoff": "#CC79A7",
    }
    markers = {"midpoint": "o", "harmonic": "s", "series": "^", "kirchhoff": "D"}

    figure, axis = plt.subplots(figsize=(7.6, 4.6))
    for name, values in series.items():
        axis.plot(
            grids,
            values,
            marker=markers.get(name, "o"),
            color=colors.get(name, "#333333"),
            linewidth=2.0,
            label=name,
        )
    if reference_hours is not None:
        axis.axhline(
            reference_hours,
            color="#555555",
            linestyle="--",
            linewidth=1.2,
            label="中点点值 N=80 参考值",
        )
    axis.set(
        title="界面扩散系数策略的网格收敛性",
        xlabel="径向区间数 N",
        ylabel="达标时刻（h）",
    )
    axis.set_xticks(list(grids))
    axis.legend()
    axis.grid(alpha=0.3, linestyle="--")
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)
