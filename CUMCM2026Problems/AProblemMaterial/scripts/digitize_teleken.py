"""从 Teleken et al. (2025) Figure 4 数字化实测干燥曲线。

两种**互相独立**的算法各跑一遍，用于满足预注册的数字化质量门
（`docs/06 §9.4`）：

* 算法 A：列窗口法——在预期采样时刻的像素列邻域内，按"连续多行达到标记宽度"
  定位标记中心；
* 算法 B：连通域法——对颜色掩膜做连通域标记，按宽高比与"宽行"结构筛掉文字
  （Stage 1/2/3 标注与蓝色/红色文字），再用最近邻匹配到采样时刻。

两次结果逐点比较，`|Δ| ≤ max(1% 读数, 0.03 kg/kg)` 才判为通过；并输出叠加
校验图供目视确认。用法：

    python scripts/digitize_teleken.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURE = (
    PROJECT_ROOT
    / "data"
    / "external"
    / "teleken2025"
    / "figures"
    / "figure4_drying_curves.jpg"
)
DATA_DIR = FIGURE.parents[1]

SAMPLE_TIMES_MIN = (0, 15, 30, 45, 60, 90, 120, 150, 180)

# 由坐标轴刻度标签中心线性标定（见 docs/06 §9.4 与 README）。
X_AT_ZERO_MIN = 61.0
X_AT_180_MIN = 337.0


@dataclass(frozen=True)
class Panel:
    """一个面板的标定与颜色定义。"""

    key: str
    temperature_c: float
    row_start: int
    row_end: int
    y_at_zero: float
    y_at_max: float
    y_max_value: float
    colour: str

    def x_pixel(self, time_min: float) -> float:
        return X_AT_ZERO_MIN + (X_AT_180_MIN - X_AT_ZERO_MIN) * time_min / 180.0

    def value_at_row(self, row: float) -> float:
        return (self.y_at_zero - row) / (self.y_at_zero - self.y_at_max) * self.y_max_value

    def row_at_value(self, value: float) -> float:
        return self.y_at_zero - value / self.y_max_value * (self.y_at_zero - self.y_at_max)


PANELS = (
    Panel("60C", 60.0, 0, 210, 174.0, 30.5, 4.0, "blue"),
    Panel("70C", 70.0, 250, 428, 419.5, 274.0, 4.0, "red"),
)

# 质量门阈值（预注册）：|Δ| ≤ max(1% 读数, 0.03 kg/kg)。
RELATIVE_GATE = 0.01
ABSOLUTE_GATE = 0.03


def colour_mask(image: np.ndarray, colour: str) -> np.ndarray:
    """返回布尔颜色掩膜（阈值由图像实测中位色确定）。"""

    red = image[:, :, 0].astype(int)
    green = image[:, :, 1].astype(int)
    blue = image[:, :, 2].astype(int)
    if colour == "blue":
        return (blue > 120) & (blue - red > 50) & (green < 190) & (green > 60)
    if colour == "red":
        return (red > 140) & (red - blue > 60) & (green < 110)
    raise ValueError(f"未知颜色：{colour}")


def _contiguous_rows(rows: np.ndarray, gap: int = 2) -> list[tuple[int, int]]:
    """把行号数组切成允许 ``gap`` 像素间断的连续段。"""

    if rows.size == 0:
        return []
    blocks: list[tuple[int, int]] = []
    start = previous = int(rows[0])
    for value in rows[1:]:
        value = int(value)
        if value - previous > gap:
            blocks.append((start, previous))
            start = value
        previous = value
    blocks.append((start, previous))
    return blocks


def _block_fill(window: np.ndarray, block: tuple[int, int]) -> tuple[float, int, int]:
    """返回候选块的（填充率、宽度、高度）。

    实心标记在同宽高比下的填充率约 0.7–0.8；文字笔画只有 0.3–0.5，
    因此填充率是把"蓝色文字"与"蓝色标记"分开的关键判据。
    """

    rows = window[block[0] : block[1] + 1, :]
    columns = np.where(rows.any(axis=0))[0]
    if columns.size == 0:
        return 0.0, 0, 0
    width = int(columns.max() - columns.min() + 1)
    height = int(block[1] - block[0] + 1)
    area = int(rows[:, columns.min() : columns.max() + 1].sum())
    return area / float(width * height), width, height


def digitize_window(
    mask: np.ndarray,
    panel: Panel,
) -> tuple[dict[int, tuple[float, float, float]], list[int]]:
    """算法 A：列窗口法。返回 ({时间: (值, 下误差, 上误差)}, 标记体宽度列表)。"""

    region = mask[panel.row_start : panel.row_end, :]
    results: dict[int, tuple[float, float, float]] = {}
    body_widths: list[int] = []
    for time_min in SAMPLE_TIMES_MIN:
        centre_column = int(round(panel.x_pixel(time_min)))
        window = region[:, centre_column - 5 : centre_column + 6]
        widths = window.sum(axis=1)
        wide_rows = np.where(widths >= 5)[0]
        blocks = _contiguous_rows(wide_rows, gap=1)
        if not blocks:
            raise RuntimeError(f"{panel.key} t={time_min} 未找到标记")
        # 标记体是持续多行、填充率高的实心块；误差棒端帽只有 1–2 行，
        # 文字笔画宽度够但填充率低。
        scored = []
        for block in blocks:
            fill, width, height = _block_fill(window, block)
            if width >= 5 and 3 <= height <= 12 and fill >= 0.55:
                scored.append((fill, block))
        if not scored:
            raise RuntimeError(f"{panel.key} t={time_min} 没有通过形状判据的标记")
        chosen = max(scored, key=lambda item: item[0])[1]
        body_widths.append(_block_fill(window, chosen)[1])
        centre_row = 0.5 * (chosen[0] + chosen[1])

        coloured_rows = np.where(widths >= 1)[0]
        containing = [
            block
            for block in _contiguous_rows(coloured_rows, gap=2)
            if block[0] <= centre_row <= block[1]
        ]
        top_row, bottom_row = containing[0] if containing else chosen

        offset = panel.row_start
        value = panel.value_at_row(centre_row + offset)
        high = panel.value_at_row(top_row + offset)
        low = panel.value_at_row(bottom_row + offset)
        results[time_min] = (value, value - low, high - value)
    return results, body_widths


def digitize_components(
    mask: np.ndarray,
    panel: Panel,
) -> dict[int, tuple[float, float, float]]:
    """算法 B：连通域法。返回 {时间: (值, 下误差, 上误差)}。"""

    region = mask[panel.row_start : panel.row_end, :]
    labels, count = ndimage.label(region)
    candidates: list[tuple[float, float, float, float]] = []
    for index in range(1, count + 1):
        rows, columns = np.where(labels == index)
        if rows.size < 15:
            continue
        width = int(columns.max() - columns.min() + 1)
        height = int(rows.max() - rows.min() + 1)
        if width > 16 or height > 70:
            continue  # 文字词块或标注线
        row_min, row_max = int(rows.min()), int(rows.max())
        component = (
            labels[row_min : row_max + 1, columns.min() : columns.max() + 1] == index
        )
        profile = np.bincount(rows - row_min, minlength=component.shape[0])
        wide = np.where(profile >= max(4, 0.7 * profile.max()))[0]
        blocks = _contiguous_rows(wide, gap=1)
        scored = []
        for block in blocks:
            fill, block_width, block_height = _block_fill(component, block)
            if block_width >= 5 and 3 <= block_height <= 12 and fill >= 0.55:
                scored.append((fill, block))
        if not scored:
            continue
        body = max(scored, key=lambda item: item[0])[1]
        centre_row = row_min + 0.5 * (body[0] + body[1])
        centre_column = 0.5 * (columns.min() + columns.max())
        offset = panel.row_start
        value = panel.value_at_row(centre_row + offset)
        high = panel.value_at_row(rows.min() + offset)
        low = panel.value_at_row(rows.max() + offset)
        candidates.append((centre_column, value, value - low, high - value))

    results: dict[int, tuple[float, float, float]] = {}
    for time_min in SAMPLE_TIMES_MIN:
        target = panel.x_pixel(time_min)
        near = [item for item in candidates if abs(item[0] - target) <= 6.0]
        if not near:
            raise RuntimeError(f"{panel.key} t={time_min} 连通域匹配失败")
        best = min(near, key=lambda item: abs(item[0] - target))
        results[time_min] = (best[1], best[2], best[3])
    return results


def overlay(path: Path, image: np.ndarray, panel: Panel, digits: dict[int, tuple[float, float, float]]) -> None:
    """把数字化结果画回原图，供目视 QA。"""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    figure, axis = plt.subplots(figsize=(9.0, 6.2))
    axis.imshow(image)
    for time_min, (value, low_error, high_error) in digits.items():
        x = panel.x_pixel(time_min)
        y = panel.row_at_value(value)
        axis.plot([x], [y], marker="+", color="#00A000", markersize=10, markeredgewidth=2.0)
        axis.plot(
            [x, x],
            [panel.row_at_value(value + high_error), panel.row_at_value(value - low_error)],
            color="#FF00FF",
            linewidth=1.0,
        )
    axis.set_title(f"{panel.key} 数字化 QA：绿色十字=标记中心，洋红线=误差棒")
    axis.set_xlim(0, image.shape[1])
    axis.set_ylim(image.shape[0], 0)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    image = np.asarray(Image.open(FIGURE).convert("RGB")).astype(int)
    rows: list[dict[str, object]] = []
    for panel in PANELS:
        mask = colour_mask(image, panel.colour)
        pass_a, body_widths = digitize_window(mask, panel)
        pass_b = digitize_components(mask, panel)
        overlay(
            DATA_DIR / f"digitization_overlay_{panel.key}.png",
            image,
            panel,
            pass_a,
        )
        # 数字化不确定度 = 两次差值的一半 + 标记半径换算值（预注册定义）。
        pixels_per_unit = (panel.y_at_zero - panel.y_at_max) / panel.y_max_value
        marker_sigma = 0.5 * float(np.median(body_widths)) / pixels_per_unit
        for time_min in SAMPLE_TIMES_MIN:
            value_a, low_a, high_a = pass_a[time_min]
            value_b, low_b, high_b = pass_b[time_min]
            difference = abs(value_a - value_b)
            gate = max(RELATIVE_GATE * abs(value_a), ABSOLUTE_GATE)
            rows.append(
                {
                    "temperature_c": panel.temperature_c,
                    "time_min": time_min,
                    "x_db": round(0.5 * (value_a + value_b), 4),
                    "err_low": round(0.5 * (low_a + low_b), 4),
                    "err_high": round(0.5 * (high_a + high_b), 4),
                    "pass_a": round(value_a, 4),
                    "pass_b": round(value_b, 4),
                    "difference": round(difference, 4),
                    "gate": round(gate, 4),
                    "digitization_sigma": round(0.5 * difference + marker_sigma, 4),
                    "included": bool(difference <= gate),
                }
            )

    for temperature_c, name in ((60.0, "figure4_60C.csv"), (70.0, "figure4_70C.csv")):
        subset = [row for row in rows if row["temperature_c"] == temperature_c]
        path = DATA_DIR / name
        with path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "temperature_c",
                    "time_min",
                    "x_db",
                    "err_low",
                    "err_high",
                    "pass_a",
                    "pass_b",
                    "difference",
                    "gate",
                    "digitization_sigma",
                    "included",
                ],
            )
            writer.writeheader()
            writer.writerows(subset)
        failed = [row for row in subset if not row["included"]]
        print(
            f"{name}: {len(subset)} 点，未过质量门 {len(failed)} 点"
            + ("" if not failed else f" → {[row['time_min'] for row in failed]}")
        )


if __name__ == "__main__":
    main()
