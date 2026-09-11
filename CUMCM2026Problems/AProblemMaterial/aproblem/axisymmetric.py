# 二维轴对称材料坐标模型
"""圆柱（含两个端面）的二维轴对称材料坐标有限体积模型。

本模块只用于校验：一维径向模型把两个端面折算成体积源项 ``2h/L`` 时，
内含"端面状态约等于同半径处轴向平均状态"的闭合假设。二维模型把端面
当作真实边界解出来，因此能给出轴向梯度究竟有多大、降维近似误差究竟
是多少。

坐标与方程
----------
令 ``xi = r/R(t)``、``eta = z/L(t)``，两个方向都跟随材料，因此
``partial/partial t`` 就是材料导数，不需要再补拖曳项。拉普拉斯算子
在两个方向各缩一个尺度：

    dT/dt = (1/(xi R^2)) d/dxi(xi k/(rho cp) dT/dxi) + (1/L^2) d/deta(k/(rho cp) dT/deta)
    dC/dt = (1/(xi R^2)) d/dxi(xi D dC/dxi) + (1/L^2) d/deta(D dC/deta)

边界条件：``xi=0`` 轴对称；``xi=1`` 是圆柱侧面；``eta=0, 1`` 是两个
端面，后两者均为第三类（Robin）换热与传质边界。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .interfaces import HarmonicInterface, InterfaceStrategy
from .physics import PropertyModel


EnvironmentFunction = Callable[[float], tuple[float, float]]
RadiusFunction = Callable[[float], float]
LengthFunction = Callable[[float], float]


@dataclass(frozen=True)
class AxisymmetricGeometry:
    """某一时刻的二维轴对称控制体几何量（真实物理面积与体积）。"""

    radial_nodes_m: np.ndarray
    axial_nodes_m: np.ndarray
    radial_spacing_m: float
    axial_spacing_m: float
    volumes_m3: np.ndarray
    radial_face_areas_m2: np.ndarray
    axial_face_areas_m2: np.ndarray
    side_areas_m2: np.ndarray
    end_areas_m2: np.ndarray
    cross_section_area_m2: float


def axisymmetric_geometry(
    radial_intervals: int,
    axial_intervals: int,
    radius_m: float,
    length_m: float,
) -> AxisymmetricGeometry:
    """构造节点中心型二维控制体几何。

    径向网格在 ``xi`` 上等距，所以物理半径 ``r = xi R`` 随半径缩放；
    轴向同理。圆心与圆柱表面、两个端面都是半控制体。
    """

    if radial_intervals < 2 or axial_intervals < 2:
        raise ValueError("二维网格的径向和轴向区间数都至少为 2")
    if radius_m <= 0.0 or length_m <= 0.0:
        raise ValueError("半径和长度必须为正数")

    radial_nodes = np.linspace(0.0, radius_m, radial_intervals + 1)
    axial_nodes = np.linspace(0.0, length_m, axial_intervals + 1)
    radial_spacing = radius_m / radial_intervals
    axial_spacing = length_m / axial_intervals

    # 径向控制体的内外边界：首个固定在圆心，末个固定在圆柱侧面。
    inner_faces = np.empty_like(radial_nodes)
    outer_faces = np.empty_like(radial_nodes)
    inner_faces[0] = 0.0
    inner_faces[1:] = 0.5 * (radial_nodes[:-1] + radial_nodes[1:])
    outer_faces[:-1] = inner_faces[1:]
    outer_faces[-1] = radius_m

    # 每个径向环带的横截面积；它同时是轴向界面的面积。
    annulus_areas = np.pi * (outer_faces**2 - inner_faces**2)
    radial_face_radii = 0.5 * (radial_nodes[:-1] + radial_nodes[1:])

    # 两个端面节点只控制半格，其余节点控制整格；这样轴向体积自然加满。
    axial_extents = np.full(axial_intervals + 1, axial_spacing)
    axial_extents[0] = 0.5 * axial_spacing
    axial_extents[-1] = 0.5 * axial_spacing

    return AxisymmetricGeometry(
        radial_nodes_m=radial_nodes,
        axial_nodes_m=axial_nodes,
        radial_spacing_m=radial_spacing,
        axial_spacing_m=axial_spacing,
        volumes_m3=annulus_areas[:, None] * axial_extents[None, :],
        radial_face_areas_m2=(
            2.0 * np.pi * radial_face_radii[:, None] * axial_extents[None, :]
        ),
        axial_face_areas_m2=annulus_areas,
        side_areas_m2=2.0 * np.pi * radius_m * axial_extents,
        end_areas_m2=annulus_areas,
        cross_section_area_m2=float(np.sum(annulus_areas)),
    )


@dataclass(frozen=True)
class AxisymmetricModel:
    """二维轴对称热湿耦合有限体积模型。

    ``insulated_ends=True`` 把两个端面的 Robin 边界换成齐次 Neumann
    边界，用于退化测试：此时二维解应与一维"仅侧面"模型完全一致。
    """

    radial_intervals: int
    axial_intervals: int
    properties: PropertyModel
    environment: EnvironmentFunction
    radius: RadiusFunction
    length: LengthFunction
    heat_transfer_coefficient: float
    mass_transfer_coefficient: float
    interface_strategy: InterfaceStrategy = HarmonicInterface()
    insulated_ends: bool = False

    @property
    def shape(self) -> tuple[int, int]:
        """温度场和含水率场的二维形状。"""

        return (self.radial_intervals + 1, self.axial_intervals + 1)

    @property
    def node_count(self) -> int:
        """单场节点总数。"""

        return self.shape[0] * self.shape[1]

    def geometry(self, radius_m: float, length_m: float) -> AxisymmetricGeometry:
        """返回给定半径和长度下的控制体几何。"""

        return axisymmetric_geometry(
            self.radial_intervals,
            self.axial_intervals,
            radius_m,
            length_m,
        )

    def split(self, state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """把扁平状态向量还原成二维温度场和含水率场。"""

        values = np.asarray(state, dtype=float)
        if values.shape != (2 * self.node_count,):
            raise ValueError(
                f"状态向量形状应为 {(2 * self.node_count,)}，实际为 {values.shape}"
            )
        temperature = values[: self.node_count].reshape(self.shape)
        moisture = values[self.node_count :].reshape(self.shape)
        return temperature, moisture

    def pack(self, temperature: np.ndarray, moisture: np.ndarray) -> np.ndarray:
        """把二维温度场和含水率场拼成扁平状态向量。"""

        return np.concatenate((temperature.ravel(), moisture.ravel()))

    def _faces(
        self,
        material,
        temperature: np.ndarray,
        moisture: np.ndarray,
        geometry: AxisymmetricGeometry,
        field: str,
    ) -> tuple[np.ndarray, np.ndarray]:
        """返回径向和轴向界面的物性值。

        径向界面把数组转置，使径向成为最后一维；轴向界面本身就是
        最后一维。两个方向共用同一个界面策略对象。
        """

        node_values = (
            material.conductivity if field == "temperature" else material.diffusivity
        )
        property_function = (
            self.properties.conductivity
            if field == "temperature"
            else self.properties.diffusivity
        )
        # 只有水分方程存在通量势；传热侧传 None，基尔霍夫策略会退回中点值。
        flux_potential_function = (
            None if field == "temperature" else self.properties.flux_potential
        )
        radial_faces = self.interface_strategy.face_property(
            node_property=node_values.T,
            node_temperature=temperature.T,
            node_moisture=moisture.T,
            node_radii_m=geometry.radial_nodes_m,
            property_function=property_function,
            flux_potential_function=flux_potential_function,
        )
        axial_faces = self.interface_strategy.face_property(
            node_property=node_values,
            node_temperature=temperature,
            node_moisture=moisture,
            node_radii_m=geometry.axial_nodes_m,
            property_function=property_function,
            flux_potential_function=flux_potential_function,
        )
        return radial_faces, axial_faces

    def _diffusion_flux(
        self,
        field: np.ndarray,
        radial_faces: np.ndarray,
        axial_faces: np.ndarray,
        geometry: AxisymmetricGeometry,
    ) -> np.ndarray:
        """把两个方向的界面通量累加成一个净流入率。"""

        net = np.zeros_like(field)
        # 径向界面值按 (轴向, 径向) 排列，转置成 (径向, 轴向) 后与面积相乘。
        radial_conductance = (
            radial_faces.T
            * geometry.radial_face_areas_m2
            / geometry.radial_spacing_m
        )
        # 与一维模型同一约定：流入本控制体为正。
        radial_flux = radial_conductance * (field[1:, :] - field[:-1, :])
        net[:-1, :] += radial_flux
        net[1:, :] -= radial_flux

        axial_conductance = (
            axial_faces
            * geometry.axial_face_areas_m2[:, None]
            / geometry.axial_spacing_m
        )
        axial_flux = axial_conductance * (field[:, 1:] - field[:, :-1])
        net[:, :-1] += axial_flux
        net[:, 1:] -= axial_flux
        return net

    def rhs(self, time_s: float, state: np.ndarray) -> np.ndarray:
        """计算二维温度场和含水率场的时间导数。"""

        temperature, moisture = self.split(state)
        radius_m = float(self.radius(time_s))
        length_m = float(self.length(time_s))
        geometry = self.geometry(radius_m, length_m)
        material = self.properties.evaluate(temperature, moisture)
        environment_temperature, environment_moisture = self.environment(time_s)

        heat = np.zeros_like(temperature)
        radial_faces, axial_faces = self._faces(
            material,
            temperature,
            moisture,
            geometry,
            "temperature",
        )
        heat += self._diffusion_flux(temperature, radial_faces, axial_faces, geometry)
        # 圆柱侧面的第三类边界，作用在最外侧径向节点的每一个轴向位置上。
        heat[-1, :] += (
            self.heat_transfer_coefficient
            * geometry.side_areas_m2
            * (environment_temperature - temperature[-1, :])
        )
        if not self.insulated_ends:
            # 两个端面是真实边界，按各自环带的横截面积加权。
            for index in (0, -1):
                heat[:, index] += (
                    self.heat_transfer_coefficient
                    * geometry.end_areas_m2
                    * (environment_temperature - temperature[:, index])
                )
        d_temperature = heat / (
            material.density * material.heat_capacity * geometry.volumes_m3
        )

        moisture_rate = np.zeros_like(moisture)
        radial_faces, axial_faces = self._faces(
            material,
            temperature,
            moisture,
            geometry,
            "moisture",
        )
        moisture_rate += self._diffusion_flux(
            moisture,
            radial_faces,
            axial_faces,
            geometry,
        )
        moisture_rate[-1, :] += (
            self.mass_transfer_coefficient
            * geometry.side_areas_m2
            * (environment_moisture - moisture[-1, :])
        )
        if not self.insulated_ends:
            for index in (0, -1):
                moisture_rate[:, index] += (
                    self.mass_transfer_coefficient
                    * geometry.end_areas_m2
                    * (environment_moisture - moisture[:, index])
                )
        d_moisture = moisture_rate / geometry.volumes_m3

        return self.pack(d_temperature, d_moisture)

    def stable_time_step(
        self,
        time_s: float,
        state: np.ndarray,
        safety_factor: float = 0.8,
    ) -> float:
        """按节点局部总流出系数估计显式稳定步长上限。"""

        if not 0.0 < safety_factor <= 1.0:
            raise ValueError("时间步安全系数必须位于 (0, 1] 区间")

        temperature, moisture = self.split(state)
        geometry = self.geometry(
            float(self.radius(time_s)),
            float(self.length(time_s)),
        )
        material = self.properties.evaluate(temperature, moisture)

        radial_faces, axial_faces = self._faces(
            material,
            temperature,
            moisture,
            geometry,
            "temperature",
        )
        heat_outflow = np.zeros_like(temperature)
        radial_conductance = (
            radial_faces.T
            * geometry.radial_face_areas_m2
            / geometry.radial_spacing_m
        )
        heat_outflow[:-1, :] += radial_conductance
        heat_outflow[1:, :] += radial_conductance
        axial_conductance = (
            axial_faces
            * geometry.axial_face_areas_m2[:, None]
            / geometry.axial_spacing_m
        )
        heat_outflow[:, :-1] += axial_conductance
        heat_outflow[:, 1:] += axial_conductance
        heat_outflow[-1, :] += (
            self.heat_transfer_coefficient * geometry.side_areas_m2
        )
        if not self.insulated_ends:
            for index in (0, -1):
                heat_outflow[:, index] += (
                    self.heat_transfer_coefficient * geometry.end_areas_m2
                )
        temperature_rates = heat_outflow / (
            material.density * material.heat_capacity * geometry.volumes_m3
        )

        radial_faces, axial_faces = self._faces(
            material,
            temperature,
            moisture,
            geometry,
            "moisture",
        )
        moisture_outflow = np.zeros_like(moisture)
        radial_conductance = (
            radial_faces.T
            * geometry.radial_face_areas_m2
            / geometry.radial_spacing_m
        )
        moisture_outflow[:-1, :] += radial_conductance
        moisture_outflow[1:, :] += radial_conductance
        axial_conductance = (
            axial_faces
            * geometry.axial_face_areas_m2[:, None]
            / geometry.axial_spacing_m
        )
        moisture_outflow[:, :-1] += axial_conductance
        moisture_outflow[:, 1:] += axial_conductance
        moisture_outflow[-1, :] += (
            self.mass_transfer_coefficient * geometry.side_areas_m2
        )
        if not self.insulated_ends:
            for index in (0, -1):
                moisture_outflow[:, index] += (
                    self.mass_transfer_coefficient * geometry.end_areas_m2
                )
        moisture_rates = moisture_outflow / geometry.volumes_m3

        maximum_rate = float(
            np.max(np.concatenate((temperature_rates.ravel(), moisture_rates.ravel())))
        )
        if not np.isfinite(maximum_rate):
            raise FloatingPointError("稳定步长估计中出现非有限值")
        if maximum_rate <= 0.0:
            return float("inf")
        return safety_factor / maximum_rate
