#问题 1、2/3、4 的物性接口
"""实现题目附录给出的三套热物性和水分扩散系数公式。

每个物性模型同时提供"四个字段单独计算"和"一次性返回全部字段"两种
调用方式。前者供界面中点策略在面两侧的中点状态上求点值，后者供有限
体积右端项按节点取值，两者共用同一组公式，不会出现两处写法不一致。
"""

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

    def conductivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """单独返回导热系数，供界面中点策略调用。"""

        ...

    def diffusivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """单独返回水分扩散系数，供界面中点策略调用。"""

        ...


def _safe_moisture(moisture: np.ndarray) -> np.ndarray:
    """仅在物性公式分母中设置极小正下限，避免数值除零。"""

    return np.maximum(np.asarray(moisture, dtype=float), 1.0e-8)


class _BasePropertyModel:
    """把四个物性字段组合成统一返回结构，子类只需实现各自公式。"""

    def density(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """返回密度；由子类实现。"""

        raise NotImplementedError

    def heat_capacity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """返回比热容；由子类实现。"""

        raise NotImplementedError

    def conductivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """返回导热系数；由子类实现。"""

        raise NotImplementedError

    def diffusivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """返回水分扩散系数；由子类实现。"""

        raise NotImplementedError

    def evaluate(
        self,
        temperature_c: np.ndarray,
        moisture: np.ndarray,
    ) -> MaterialProperties:
        """按节点状态一次性返回全部物性。"""

        return MaterialProperties(
            density=np.asarray(self.density(temperature_c, moisture), dtype=float),
            heat_capacity=np.asarray(
                self.heat_capacity(temperature_c, moisture),
                dtype=float,
            ),
            conductivity=np.asarray(
                self.conductivity(temperature_c, moisture),
                dtype=float,
            ),
            diffusivity=np.asarray(
                self.diffusivity(temperature_c, moisture),
                dtype=float,
            ),
        )


@dataclass(frozen=True)
class Question1Properties(_BasePropertyModel):
    """问题1：常密度、常比热、常导热系数，扩散系数只依赖含水率。"""

    def density(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """问题1的密度为常数 820 kg/m3。"""

        return np.full(np.shape(temperature_c), 820.0)

    def heat_capacity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """问题1的比热容为常数 2600 J/(kg·K)。"""

        return np.full(np.shape(temperature_c), 2600.0)

    def conductivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """问题1的导热系数为常数 0.36 W/(m·K)。"""

        return np.full(np.shape(temperature_c), 0.36)

    def diffusivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """问题1的扩散系数只依赖含水率。"""

        c = _safe_moisture(moisture)
        return 7.0e-9 * np.exp(-0.89 / c)


@dataclass(frozen=True)
class Question23Properties(_BasePropertyModel):
    """问题2和问题3共用的含水率相关热物性与温湿耦合扩散系数。"""

    def density(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录3的密度公式。"""

        c = _safe_moisture(moisture)
        return 650.0 + 128.0 * c

    def heat_capacity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录3的比热容公式。"""

        c = _safe_moisture(moisture)
        return 1450.0 + 2736.0 * c / (c + 1.0)

    def conductivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录3的导热系数公式。"""

        c = _safe_moisture(moisture)
        return 0.21 + 0.38 * c / (c + 1.0)

    def diffusivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录3的扩散系数；Arrhenius 指数中的温度使用开尔文。"""

        c = _safe_moisture(moisture)
        temperature_k = np.asarray(temperature_c, dtype=float) + 273.15
        return 2.4e-3 * np.exp(-0.45 / c) * np.exp(-3850.0 / temperature_k)


@dataclass(frozen=True)
class Question4Properties(_BasePropertyModel):
    """问题4收缩条件下使用的热物性和水分扩散系数。"""

    def density(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录4的密度公式。"""

        c = _safe_moisture(moisture)
        return 760.0 + 90.0 * c

    def heat_capacity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录4的比热容公式。"""

        c = _safe_moisture(moisture)
        return 1850.0 + 2150.0 * c / (c + 1.0)

    def conductivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录4的导热系数公式。"""

        c = _safe_moisture(moisture)
        return 0.12 + 0.20 * c / (c + 1.0)

    def diffusivity(self, temperature_c: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """附录4的扩散系数；Arrhenius 指数中的温度使用开尔文。"""

        c = _safe_moisture(moisture)
        temperature_k = np.asarray(temperature_c, dtype=float) + 273.15
        return 4.2e-4 * np.exp(-0.30 / c) * np.exp(-3850.0 / temperature_k)
