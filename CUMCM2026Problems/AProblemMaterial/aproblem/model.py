#热湿有限体积 RHS
"""建立温度—含水率耦合的一维圆柱有限体积半离散模型。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .grid import RadialGrid
from .interfaces import HarmonicInterface, InterfaceStrategy
from .physics import PropertyModel


EnvironmentFunction = Callable[[float], tuple[float, float]]
RadiusFunction = Callable[[float], float]


@dataclass(frozen=True)
class DryingModel:
    """封装网格、物性、侧面边界、端面等效源项和热湿有限体积右端项。

    状态向量前半段为温度，后半段为干基含水率。问题4通过固定节点编号
    配合随时间变化的物理半径实现材料坐标网格，不额外添加收缩拖曳项。

    界面上的导热系数和扩散系数由 ``interface_strategy`` 决定；右端项和
    稳定步长估计共用同一个策略对象，避免两处取值不一致。
    """

    grid: RadialGrid
    properties: PropertyModel
    environment: EnvironmentFunction
    radius: RadiusFunction
    heat_transfer_coefficient: float
    mass_transfer_coefficient: float
    cylinder_length_m: float = 0.25
    include_end_faces: bool = True
    interface_strategy: InterfaceStrategy = HarmonicInterface()
    length_function: Callable[[float], float] | None = None

    def __post_init__(self) -> None:
        """检查端面等效源项所需的圆柱长度。"""

        if self.cylinder_length_m <= 0:
            raise ValueError("药材长度 cylinder_length_m 必须为正数")

    def length_at(self, time_s: float) -> float:
        """返回当前时刻的药材长度。

        ``length_function`` 为空表示长度保持 ``cylinder_length_m`` 不变；
        传入函数则支持各向同性同步收缩等敏感性工况。两种情形下长度都
        必须为正，否则无法定义端面等效源项。
        """

        if self.length_function is None:
            return self.cylinder_length_m
        length_m = float(self.length_function(time_s))
        if length_m <= 0.0:
            raise ValueError("药材长度必须为正数")
        return length_m

    @property
    def node_count(self) -> int:
        """返回包含圆心和表面的径向节点总数。"""

        return self.grid.intervals + 1

    def split_state(self, state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """将一维状态向量拆分为温度场和含水率场。"""

        if state.shape != (2 * self.node_count,):
            raise ValueError(
                f"状态向量形状应为 {(2 * self.node_count,)}，实际为 {state.shape}"
            )
        return state[: self.node_count], state[self.node_count :]

    def rhs(self, time_s: float, state: np.ndarray) -> np.ndarray:
        """计算当前时刻温度和含水率在所有节点上的变化率。"""

        temperature, moisture = self.split_state(np.asarray(state, dtype=float))
        geometry = self.grid.geometry(self.radius(time_s))
        material = self.properties.evaluate(temperature, moisture)
        environment_temperature, environment_moisture = self.environment(time_s)

        # 先计算全部内部界面热通量，再以相反符号计入相邻控制体，保证守恒。
        heat_rate = np.zeros(self.node_count)
        conductivity_faces = self.interface_strategy.face_property(
            node_property=material.conductivity,
            node_temperature=temperature,
            node_moisture=moisture,
            node_radii_m=geometry.nodes_m,
            property_function=self.properties.conductivity,
        )
        heat_conductance = (
            conductivity_faces
            * geometry.interface_areas_per_length_m
            / geometry.spacing_m
        )
        delta_temperature = temperature[1:] - temperature[:-1]
        heat_rate[:-1] += heat_conductance * delta_temperature
        heat_rate[1:] -= heat_conductance * delta_temperature
        # 表面采用第三类边界条件：烘房温度高于表面时，热量流入药材。
        heat_rate[-1] += (
            self.heat_transfer_coefficient
            * geometry.surface_area_per_length_m
            * (environment_temperature - temperature[-1])
        )

        if self.include_end_faces:
            # 将两个端面的轴向对流通量除以长度，折算为一维径向方程中的体积热源。
            # 该闭合近似假设端面状态可由同一半径处的轴向平均状态表示。
            length_m = self.length_at(time_s)
            end_heat_source = (
                2.0
                * self.heat_transfer_coefficient
                / length_m
                * (environment_temperature - temperature)
            )
            heat_rate += end_heat_source * geometry.volumes_per_length_m2

        d_temperature = heat_rate / (
            material.density * material.heat_capacity * geometry.volumes_per_length_m2
        )

        # 水分通量与热通量使用同一套控制体几何和符号约定。
        moisture_rate = np.zeros(self.node_count)
        diffusivity_faces = self.interface_strategy.face_property(
            node_property=material.diffusivity,
            node_temperature=temperature,
            node_moisture=moisture,
            node_radii_m=geometry.nodes_m,
            property_function=self.properties.diffusivity,
        )
        moisture_conductance = (
            diffusivity_faces
            * geometry.interface_areas_per_length_m
            / geometry.spacing_m
        )
        delta_moisture = moisture[1:] - moisture[:-1]
        moisture_rate[:-1] += moisture_conductance * delta_moisture
        moisture_rate[1:] -= moisture_conductance * delta_moisture
        # 烘房含水率更低时，该项为负，表示药材从表面失水。
        moisture_rate[-1] += (
            self.mass_transfer_coefficient
            * geometry.surface_area_per_length_m
            * (environment_moisture - moisture[-1])
        )

        if self.include_end_faces:
            # 两个端面的失水通量等效为体积水分汇；环境更干时该项为负。
            length_m = self.length_at(time_s)
            end_moisture_source = (
                2.0
                * self.mass_transfer_coefficient
                / length_m
                * (environment_moisture - moisture)
            )
            moisture_rate += end_moisture_source * geometry.volumes_per_length_m2

        d_moisture = moisture_rate / geometry.volumes_per_length_m2

        return np.concatenate((d_temperature, d_moisture))

    def stable_time_step(
        self,
        time_s: float,
        state: np.ndarray,
        safety_factor: float = 0.8,
    ) -> float:
        """返回当前离散系统允许的显式稳定时间步上限。

        该方法采用与控制体离散一致的总流出系数估计。界面导热系数和扩散
        系数继续使用调和平均；侧面 Robin 通量和端面等效源项均纳入稳定
        性限制。问题 4 收缩后网格间距变小，因此该方法会在每个成功时间步
        重新计算。
        """

        if not 0.0 < safety_factor <= 1.0:
            raise ValueError("时间步安全系数必须位于 (0, 1] 区间")

        temperature, moisture = self.split_state(np.asarray(state, dtype=float))
        geometry = self.grid.geometry(self.radius(time_s))
        material = self.properties.evaluate(temperature, moisture)
        node_count = self.node_count

        conductivity_faces = self.interface_strategy.face_property(
            node_property=material.conductivity,
            node_temperature=temperature,
            node_moisture=moisture,
            node_radii_m=geometry.nodes_m,
            property_function=self.properties.conductivity,
        )
        heat_faces = (
            conductivity_faces
            * geometry.interface_areas_per_length_m
            / geometry.spacing_m
        )
        heat_outflow = np.zeros(node_count)
        heat_outflow[:-1] += heat_faces
        heat_outflow[1:] += heat_faces
        heat_outflow[-1] += (
            self.heat_transfer_coefficient
            * geometry.surface_area_per_length_m
        )

        heat_capacity_per_length = (
            material.density
            * material.heat_capacity
            * geometry.volumes_per_length_m2
        )
        if self.include_end_faces:
            length_m = self.length_at(time_s)
            heat_outflow += (
                2.0 * self.heat_transfer_coefficient / length_m
            ) * geometry.volumes_per_length_m2
        temperature_rates = heat_outflow / heat_capacity_per_length

        diffusivity_faces = self.interface_strategy.face_property(
            node_property=material.diffusivity,
            node_temperature=temperature,
            node_moisture=moisture,
            node_radii_m=geometry.nodes_m,
            property_function=self.properties.diffusivity,
        )
        moisture_faces = (
            diffusivity_faces
            * geometry.interface_areas_per_length_m
            / geometry.spacing_m
        )
        moisture_outflow = np.zeros(node_count)
        moisture_outflow[:-1] += moisture_faces
        moisture_outflow[1:] += moisture_faces
        moisture_outflow[-1] += (
            self.mass_transfer_coefficient
            * geometry.surface_area_per_length_m
        )
        if self.include_end_faces:
            length_m = self.length_at(time_s)
            moisture_outflow += (
                2.0 * self.mass_transfer_coefficient / length_m
            ) * geometry.volumes_per_length_m2
        moisture_rates = moisture_outflow / geometry.volumes_per_length_m2

        maximum_rate = float(np.max(np.concatenate((temperature_rates, moisture_rates))))
        if not np.isfinite(maximum_rate):
            raise FloatingPointError("稳定步长估计中出现非有限值")
        if maximum_rate <= 0.0:
            return float("inf")
        return safety_factor / maximum_rate
