# 界面物性取值策略
"""定义有限体积界面上的导热系数和扩散系数取值方式。

本题的水分扩散系数 D 在干燥过程中跨越约四个数量级，表面先干、内部
后干，界面上的 D 取多大直接决定干燥时长。因此界面取值不能靠"结果像
不像"来选，必须做成可切换的策略，再用网格收敛性判断。

四套策略：

midpoint
    用面两侧节点状态的中点 (T_f, C_f) 代入物性公式求点值。连续、同质
    材料中这是最自然的解释，不引入任何分层假设。
harmonic
    调和平均 2xy/(x+y)，经典的串联阻力结果。
series
    圆柱径向串联阻力，含 ln(r_f/r_i) 对数项：
    x_f = (r2-r1) / ( r_f [ ln(r_f/r1)/x1 + ln(r2/r_f)/x2 ] )。
    等距平板时该式退化为调和平均，因此同一策略也能用于二维模型的
    轴向界面。
kirchhoff
    通量势（基尔霍夫变换）取法：D_face = ΔΦ/ΔC，Φ(C)=∫₀^C D dc。
    它是唯一能让两点有限体积通量格式与真实通量严格相等的取值，
    不含任何自由参数。传热侧没有通量势输入时退回中点值。

series 在圆心那一格需要特殊说明：中心控制体是一块实心小圆盘，它的
内侧没有界面，半径从 0 出发，对数项在 r=0 处发散。薄层极限下该半格
的串联阻力就是调和平均，因此中心面退回调和平均，误差随网格加密按
O(dr^2) 衰减。

kirchhoff 的推导：无源稳态一维通量满足

    D_face·(C_R − C_L) = Φ(C_R) − Φ(C_L)

所以 ΔΦ/ΔC 与真实通量逐点相等。调和平均与算术平均分别对应"把 D 当成
阶梯跳变"和"把 D 当成线性变化"两种重构假设，都只是它的近似；当 D 在
一层网格内跨越若干数量级时，调和平均会被小值支配，形成人为的分层阻力。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol

import numpy as np


# 物性函数：输入温度（摄氏度）和干基含水率（kg/kg），返回同形状的物性数组。
PropertyFunction = Callable[[np.ndarray, np.ndarray], np.ndarray]

# 防止除零的下限；真实物性远大于该值，不影响结果。
TINY_VALUE = 1.0e-30

# 通量势斜率 ΔΦ/ΔC 在两侧含水率几乎相等时会退化为 0/0，
# 低于该差值时改用中点状态点值兜底（此时两者在浮点意义下本就相等）。
MIN_MOISTURE_DIFFERENCE = 1.0e-9


class InterfaceStrategy(Protocol):
    """所有界面取值策略必须遵循的调用接口。"""

    name: str

    def face_property(
        self,
        node_property: np.ndarray,
        node_temperature: np.ndarray,
        node_moisture: np.ndarray,
        node_radii_m: np.ndarray,
        property_function: PropertyFunction | None = None,
        flux_potential_function: PropertyFunction | None = None,
    ) -> np.ndarray:
        """返回相邻节点之间所有界面的物性值。

        所有输入数组沿**最后一维**排列节点；因此二维模型只要把要取界面的
        方向放到最后一维即可复用同一套策略。``node_radii_m`` 是一维数组，
        长度等于节点个数。

        ``property_function(T, C)`` 返回点物性，``flux_potential_function
        (T, C)`` 返回通量势 Φ。后者只在水分方程上有定义；传热方程的调用
        不传入它，此时需要通量势的策略会退回各自的安全取法。
        """

        ...


def harmonic_pair(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """成对的调和平均；分母为零时返回零。"""

    denominator = left + right
    return np.divide(
        2.0 * left * right,
        denominator,
        out=np.zeros_like(denominator),
        where=denominator > 0.0,
    )


@dataclass(frozen=True)
class MidpointInterface:
    """中点状态点值策略。"""

    name: str = "midpoint"

    def face_property(
        self,
        node_property: np.ndarray,
        node_temperature: np.ndarray,
        node_moisture: np.ndarray,
        node_radii_m: np.ndarray,
        property_function: PropertyFunction | None = None,
        flux_potential_function: PropertyFunction | None = None,
    ) -> np.ndarray:
        """在相邻节点的中点状态上直接调用物性公式。"""

        if property_function is None:
            raise ValueError("midpoint 策略需要物性函数来计算界面中点状态")
        temperature = np.asarray(node_temperature, dtype=float)
        moisture = np.asarray(node_moisture, dtype=float)
        face_temperature = 0.5 * (temperature[..., :-1] + temperature[..., 1:])
        face_moisture = 0.5 * (moisture[..., :-1] + moisture[..., 1:])
        return np.asarray(
            property_function(face_temperature, face_moisture),
            dtype=float,
        )


@dataclass(frozen=True)
class HarmonicInterface:
    """调和平均策略。"""

    name: str = "harmonic"

    def face_property(
        self,
        node_property: np.ndarray,
        node_temperature: np.ndarray,
        node_moisture: np.ndarray,
        node_radii_m: np.ndarray,
        property_function: PropertyFunction | None = None,
        flux_potential_function: PropertyFunction | None = None,
    ) -> np.ndarray:
        """对相邻节点的物性取调和平均。"""

        values = np.asarray(node_property, dtype=float)
        return harmonic_pair(values[..., :-1], values[..., 1:])


@dataclass(frozen=True)
class SeriesResistanceInterface:
    """圆柱径向串联阻力策略。"""

    name: str = "series"

    def face_property(
        self,
        node_property: np.ndarray,
        node_temperature: np.ndarray,
        node_moisture: np.ndarray,
        node_radii_m: np.ndarray,
        property_function: PropertyFunction | None = None,
        flux_potential_function: PropertyFunction | None = None,
    ) -> np.ndarray:
        """按对数半径的串联阻力求等效界面物性。

        对每个界面同时计算左、右半段的阻力。半径为零的中心面没有左侧
        半格，退回调和平均；等距平板时全部界面都退化为调和平均。
        """

        values = np.maximum(np.asarray(node_property, dtype=float), TINY_VALUE)
        radii = np.asarray(node_radii_m, dtype=float)
        if radii.ndim != 1 or radii.shape[0] != values.shape[-1]:
            raise ValueError("节点半径必须是一维数组，且长度等于节点个数")

        left_value = values[..., :-1]
        right_value = values[..., 1:]
        left_radius = radii[:-1]
        right_radius = radii[1:]
        face_radius = np.maximum(0.5 * (left_radius + right_radius), TINY_VALUE)

        # 半径为零表示中心实心格的退化界面，改用调和平均。
        annular = left_radius > 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            left_term = np.log(face_radius / np.maximum(left_radius, TINY_VALUE))
            right_term = np.log(right_radius / face_radius)
        resistance = np.where(
            annular,
            left_term / left_value + right_term / right_value,
            0.0,
        )

        thickness = right_radius - left_radius
        denominator = face_radius * resistance
        series_value = np.divide(
            thickness,
            denominator,
            out=np.zeros(np.broadcast_shapes(thickness.shape, left_value.shape)),
            where=denominator > 0.0,
        )
        return np.where(
            annular,
            series_value,
            harmonic_pair(left_value, right_value),
        )


@dataclass(frozen=True)
class KirchhoffInterface:
    """通量势（基尔霍夫变换）策略：D_face = ΔΦ/ΔC。

    无源稳态一维通量在界面两侧相等，对 D·dC/dr 积分得

        D_face·(C_R − C_L)/Δr = (Φ(C_R) − Φ(C_L))/Δr

    因此 ΔΦ/ΔC 是唯一能让有限体积两点通量格式与真实通量严格相等的界面
    取值。它不含任何可调参数，也不需要假设 D 在网格内是阶梯状或线性。

    传热方程没有对应的通量势（热流是 k·∇T，不是 k·∇C），调用方不会传入
    ``flux_potential_function``，此时本策略退回中点状态点值，与主方案保持
    一致。
    """

    name: str = "kirchhoff"

    def face_property(
        self,
        node_property: np.ndarray,
        node_temperature: np.ndarray,
        node_moisture: np.ndarray,
        node_radii_m: np.ndarray,
        property_function: PropertyFunction | None = None,
        flux_potential_function: PropertyFunction | None = None,
    ) -> np.ndarray:
        """用通量势斜率作为界面物性；无通量势时退回中点值。"""

        if flux_potential_function is None:
            return MidpointInterface().face_property(
                node_property=node_property,
                node_temperature=node_temperature,
                node_moisture=node_moisture,
                node_radii_m=node_radii_m,
                property_function=property_function,
            )
        if property_function is None:
            raise ValueError("kirchhoff 策略需要物性函数作为小差值的兜底取值")

        temperature = np.asarray(node_temperature, dtype=float)
        moisture = np.asarray(node_moisture, dtype=float)
        face_temperature = 0.5 * (temperature[..., :-1] + temperature[..., 1:])
        left_moisture = moisture[..., :-1]
        right_moisture = moisture[..., 1:]

        # 通量势在同一个面温度上求值，使 ΔΦ/ΔC 与面通量严格配对。
        potential_left = np.asarray(
            flux_potential_function(face_temperature, left_moisture),
            dtype=float,
        )
        potential_right = np.asarray(
            flux_potential_function(face_temperature, right_moisture),
            dtype=float,
        )
        delta_moisture = right_moisture - left_moisture
        delta_potential = potential_right - potential_left

        # 两侧含水率几乎相等时 0/0，用中点状态点值兜底（两者极限相同）。
        distinguishable = np.abs(delta_moisture) > MIN_MOISTURE_DIFFERENCE
        fallback = np.asarray(
            property_function(
                face_temperature,
                0.5 * (left_moisture + right_moisture),
            ),
            dtype=float,
        )
        return np.where(
            distinguishable,
            np.divide(
                delta_potential,
                delta_moisture,
                out=np.zeros_like(delta_potential),
                where=distinguishable,
            ),
            fallback,
        )


# 策略名称到实例的映射，供配置和命令行按名字选择。
INTERFACE_STRATEGIES: dict[str, InterfaceStrategy] = {
    "midpoint": MidpointInterface(),
    "harmonic": HarmonicInterface(),
    "series": SeriesResistanceInterface(),
    "kirchhoff": KirchhoffInterface(),
}

# 策略的名字集合，供参数校验使用。
INTERFACE_STRATEGY_NAMES = tuple(INTERFACE_STRATEGIES)


def resolve_interface_strategy(name: str) -> InterfaceStrategy:
    """按名称返回界面策略实例。"""

    try:
        return INTERFACE_STRATEGIES[name]
    except KeyError as error:
        raise ValueError(
            f"未知的界面策略 {name!r}；可选：{', '.join(INTERFACE_STRATEGY_NAMES)}"
        ) from error
