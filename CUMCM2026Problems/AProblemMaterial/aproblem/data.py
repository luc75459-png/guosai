#附件读取、环境插值、半径 PCHIP
"""读取竞赛附件，并构造随时间变化的环境边界与收缩半径函数。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator


@dataclass(frozen=True)
class EnvironmentBoundary:
    """烘房温度与含水率边界条件的分段线性插值器。"""

    time_s: np.ndarray
    temperature_c: np.ndarray
    moisture: np.ndarray
    plateau_temperature_c: float
    plateau_moisture: float

    def __call__(self, time_s: float) -> tuple[float, float]:
        """返回给定时刻的烘房温度和含水率。

        附件时间范围外：左侧保持首值，4 小时以后保持配置的恒温段值。
        """

        temperature = float(
            np.interp(
                time_s,
                self.time_s,
                self.temperature_c,
                left=self.temperature_c[0],
                right=self.plateau_temperature_c,
            )
        )
        moisture = float(
            np.interp(
                time_s,
                self.time_s,
                self.moisture,
                left=self.moisture[0],
                right=self.plateau_moisture,
            )
        )
        return temperature, moisture


@dataclass(frozen=True)
class RadiusBoundary:
    """问题4所用的药材半径 PCHIP 保形插值器。"""

    time_s: np.ndarray
    radius_m: np.ndarray
    _interpolator: PchipInterpolator

    @classmethod
    def from_arrays(cls, time_s: np.ndarray, radius_m: np.ndarray) -> "RadiusBoundary":
        """由已校验的时间和半径数组建立保形插值。"""

        interpolator = PchipInterpolator(time_s, radius_m, extrapolate=False)
        return cls(time_s=time_s, radius_m=radius_m, _interpolator=interpolator)

    def __call__(self, time_s: float) -> float:
        """返回给定时刻的半径；超出附件范围时保持端点半径。"""

        clipped_time = float(np.clip(time_s, self.time_s[0], self.time_s[-1]))
        return float(self._interpolator(clipped_time))


def _read_numeric_columns(path: Path, expected_columns: int) -> np.ndarray:
    """读取 Excel 的前若干列，并检查数值类型和时间单调性。"""

    if not path.exists():
        raise FileNotFoundError(f"未找到输入工作簿：{path}")
    frame = pd.read_excel(path)
    if frame.shape[1] < expected_columns:
        raise ValueError(f"{path.name} 至少需要 {expected_columns} 列数据")
    values = frame.iloc[:, :expected_columns].apply(pd.to_numeric, errors="raise").to_numpy(float)
    if np.any(np.diff(values[:, 0]) <= 0):
        raise ValueError(f"{path.name} 的时间列必须严格递增")
    return values


def load_environment(
    path: Path,
    plateau_temperature_c: float = 50.0,
    plateau_moisture: float = 0.05,
) -> EnvironmentBoundary:
    """读取附件1并生成烘房环境边界函数。"""

    values = _read_numeric_columns(path, expected_columns=3)
    return EnvironmentBoundary(
        time_s=values[:, 0],
        temperature_c=values[:, 1],
        moisture=values[:, 2],
        plateau_temperature_c=plateau_temperature_c,
        plateau_moisture=plateau_moisture,
    )


def load_radius(path: Path) -> RadiusBoundary:
    """读取附件2，将 cm 转为 m，并生成问题4的半径函数。"""

    values = _read_numeric_columns(path, expected_columns=2)
    radius_m = values[:, 1] * 0.01
    if np.any(np.diff(radius_m) > 1.0e-12):
        raise ValueError("附件2中的药材半径必须保持不增")
    return RadiusBoundary.from_arrays(values[:, 0], radius_m)
