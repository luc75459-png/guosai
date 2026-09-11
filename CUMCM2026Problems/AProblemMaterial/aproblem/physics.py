#问题 1、2/3、4 的物性接口
"""实现题目附录给出的三套热物性和水分扩散系数公式。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


@dataclass(frozen=True)
class MaterialProperties:
    """同一组径向节点上的密度、比热、导热系数和扩散系数。"""

    density: np.ndarray
    heat_capacity: np.ndarray
    conductivity: np.ndarray
    diffusivity: np.ndarray


class PropertyModel(Protocol):
    """所有问题物性模型必须遵循的统一调用接口。"""

    def evaluate(self, temperature_c: np.ndarray, moisture: np.ndarray) -> MaterialProperties:
        """根据当前温度和含水率返回各节点物性。"""

        ...


def _safe_moisture(moisture: np.ndarray) -> np.ndarray:
    """仅在物性公式分母中设置极小正下限，避免数值除零。"""

    return np.maximum(np.asarray(moisture, dtype=float), 1.0e-8)


@dataclass(frozen=True)
class Question1Properties:
    """问题1：常密度、常比热、常导热系数，扩散系数只依赖含水率。"""

    def evaluate(self, temperature_c: np.ndarray, moisture: np.ndarray) -> MaterialProperties:
        """计算问题1各径向节点的物性数组。"""

        c = _safe_moisture(moisture)
        shape = np.asarray(temperature_c, dtype=float).shape
        return MaterialProperties(
            density=np.full(shape, 820.0),
            heat_capacity=np.full(shape, 2600.0),
            conductivity=np.full(shape, 0.36),
            diffusivity=7.0e-9 * np.exp(-0.89 / c),
        )


@dataclass(frozen=True)
class Question23Properties:
    """问题2和问题3共用的含水率相关热物性与温湿耦合扩散系数。"""

    def evaluate(self, temperature_c: np.ndarray, moisture: np.ndarray) -> MaterialProperties:
        """计算问题2/3物性；Arrhenius 指数中的温度使用开尔文。"""

        c = _safe_moisture(moisture)
        temperature_k = np.asarray(temperature_c, dtype=float) + 273.15
        return MaterialProperties(
            density=650.0 + 128.0 * c,
            heat_capacity=1450.0 + 2736.0 * c / (c + 1.0),
            conductivity=0.21 + 0.38 * c / (c + 1.0),
            diffusivity=2.4e-3 * np.exp(-0.45 / c) * np.exp(-3850.0 / temperature_k),
        )


@dataclass(frozen=True)
class Question4Properties:
    """问题4收缩条件下使用的热物性和水分扩散系数。"""

    def evaluate(self, temperature_c: np.ndarray, moisture: np.ndarray) -> MaterialProperties:
        """计算问题4物性；Arrhenius 指数中的温度使用开尔文。"""

        c = _safe_moisture(moisture)
        temperature_k = np.asarray(temperature_c, dtype=float) + 273.15
        return MaterialProperties(
            density=760.0 + 90.0 * c,
            heat_capacity=1850.0 + 2150.0 * c / (c + 1.0),
            conductivity=0.12 + 0.20 * c / (c + 1.0),
            diffusivity=4.2e-4 * np.exp(-0.30 / c) * np.exp(-3850.0 / temperature_k),
        )
