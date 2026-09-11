#问题 1、2/3、4 的物性接口
"""实现题目附录给出的三套热物性和水分扩散系数公式。

每个物性模型同时提供"四个字段单独计算"和"一次性返回全部字段"两种
调用方式。前者供界面中点策略在面两侧的中点状态上求点值，后者供有限
体积右端项按节点取值，两者共用同一组公式，不会出现两处写法不一致。

此外每个模型还提供水分扩散的通量势 ``flux_potential``：

    Φ(C, T) = ∫₀^C D(c, T) dc

三套附录的 D 都是 ``A(T)·exp(−b/C)`` 形式，因此 Φ 有解析原函数
``A(T)·[C·exp(−b/C) − b·E₁(b/C)]``（E₁ 为指数积分），且 dΦ/dC = D。
基尔霍夫界面策略直接用 ΔΦ/ΔC 作为界面扩散系数，不需要对 D 做任何
加权平均。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from scipy.special import exp1


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

    def flux_potential(
        self,
        temperature_c: np.ndarray,
        moisture: np.ndarray,
    ) -> np.ndarray:
        """返回水分扩散通量势 Φ(C)=∫₀^C D dc，供基尔霍夫界面策略调用。"""

        ...


def _safe_moisture(moisture: np.ndarray) -> np.ndarray:
    """仅在物性公式分母中设置极小正下限，避免数值除零。"""

    return np.maximum(np.asarray(moisture, dtype=float), 1.0e-8)


# 通量势积分下限对应的极小含水率，避免 C=0 处除零；Φ 在该量级已为 0。
_POTENTIAL_MIN_MOISTURE = 1.0e-8


class _BasePropertyModel:
    """把四个物性字段组合成统一返回结构，子类只需实现各自公式。

    子类通过三个类属性声明自己的 D 公式形状 ``A·exp(−b/C)``：

    ``diffusivity_prefactor``  A 的前因子（不含温度因子）
    ``diffusivity_decay``      b
    ``diffusivity_uses_temperature``  是否还要乘 exp(−3850/T)，T 用开尔文
    """

    # 默认值仅用于兜底，三个子类都会显式覆盖。
    diffusivity_prefactor = 0.0
    diffusivity_decay = 1.0
    diffusivity_uses_temperature = False

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

    def flux_potential(
        self,
        temperature_c: np.ndarray,
        moisture: np.ndarray,
    ) -> np.ndarray:
        """返回解析通量势 Φ(C)=A(T)·[C·exp(−b/C) − b·E₁(b/C)]。

        推导：令 u=b/C，则

            ∫exp(−b/C)dC = C·exp(−b/C) + b·∫exp(−u)/u du
                         = C·exp(−b/C) − b·E₁(b/C) + const

        下限 C→0 时两项都趋于 0，所以该式就是 ∫₀^C D dc 本身。
        对 C 求导恰好还原 D(C)，因此 ΔΦ/ΔC 与真实通量严格一致。
        """

        c = np.maximum(
            np.asarray(moisture, dtype=float),
            _POTENTIAL_MIN_MOISTURE,
        )
        prefactor = np.full(np.shape(c), float(self.diffusivity_prefactor))
        if self.diffusivity_uses_temperature:
            temperature_k = np.asarray(temperature_c, dtype=float) + 273.15
            prefactor = prefactor * np.exp(-3850.0 / temperature_k)
        decay = float(self.diffusivity_decay)
        ratio = decay / c
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            integral = c * np.exp(-ratio) - decay * exp1(ratio)
        return prefactor * integral

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

    diffusivity_prefactor = 7.0e-9
    diffusivity_decay = 0.89
    diffusivity_uses_temperature = False

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

    diffusivity_prefactor = 2.4e-3
    diffusivity_decay = 0.45
    diffusivity_uses_temperature = True

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

    diffusivity_prefactor = 4.2e-4
    diffusivity_decay = 0.30
    diffusivity_uses_temperature = True

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
