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
