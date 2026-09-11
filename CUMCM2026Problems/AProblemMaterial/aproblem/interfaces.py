# 界面物性取值策略
"""定义有限体积界面上的导热系数和扩散系数取值方式。

本题的水分扩散系数 D 在干燥过程中跨越约四个数量级，表面先干、内部
后干，界面上的 D 取多大直接决定干燥时长。因此界面取值不能靠"结果像
不像"来选，必须做成可切换的策略，再用网格收敛性判断。

三套策略：

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

series 在圆心那一格需要特殊说明：中心控制体是一块实心小圆盘，它的
内侧没有界面，半径从 0 出发，对数项在 r=0 处发散。薄层极限下该半格
的串联阻力就是调和平均，因此中心面退回调和平均，误差随网格加密按
O(dr^2) 衰减。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol

import numpy as np


# 物性函数：输入温度（摄氏度）和干基含水率（kg/kg），返回同形状的物性数组。
PropertyFunction = Callable[[np.ndarray, np.ndarray], np.ndarray]

# 防止除零的下限；真实物性远大于该值，不影响结果。
TINY_VALUE = 1.0e-30


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
    ) -> np.ndarray:
        """返回相邻节点之间所有界面的物性值。

        所有输入数组沿**最后一维**排列节点；因此二维模型只要把要取界面的
        方向放到最后一维即可复用同一套策略。``node_radii_m`` 是一维数组，
        长度等于节点个数。
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


# 策略名称到实例的映射，供配置和命令行按名字选择。
INTERFACE_STRATEGIES: dict[str, InterfaceStrategy] = {
    "midpoint": MidpointInterface(),
    "harmonic": HarmonicInterface(),
    "series": SeriesResistanceInterface(),
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
