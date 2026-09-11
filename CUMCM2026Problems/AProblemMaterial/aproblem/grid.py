#圆柱节点型控制体几何
"""构造一维圆柱径向节点网格及有限体积几何量。"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RadialGeometry:
    """某一给定半径下的径向节点、控制体积和界面面积。"""

    radius_m: float
    nodes_m: np.ndarray
    inner_faces_m: np.ndarray
    outer_faces_m: np.ndarray
    volumes_per_length_m2: np.ndarray
    interface_areas_per_length_m: np.ndarray
    surface_area_per_length_m: float
    spacing_m: float


@dataclass(frozen=True)
class RadialGrid:
    """用固定区间数生成圆柱径向节点型有限体积网格。"""

    intervals: int

    def geometry(self, radius_m: float) -> RadialGeometry:
        """计算单位圆柱长度下的完整有限体积几何信息。

        节点包含圆心和药材表面；圆心与表面节点对应半控制体。
        问题4调用本方法时传入随时间变化的半径，即可得到收缩网格。
        """

        if radius_m <= 0:
            raise ValueError("药材半径 radius_m 必须为正数")
        if self.intervals < 2:
            raise ValueError("径向区间数 intervals 至少为 2")

        nodes = np.linspace(0.0, radius_m, self.intervals + 1)
        spacing = radius_m / self.intervals

        # 控制体边界取相邻节点中点；首尾边界分别固定在圆心和药材表面。
        inner_faces = np.empty_like(nodes)
        outer_faces = np.empty_like(nodes)
        inner_faces[0] = 0.0
        inner_faces[1:] = 0.5 * (nodes[:-1] + nodes[1:])
        outer_faces[:-1] = inner_faces[1:]
        outer_faces[-1] = radius_m

        # 以下量均按单位圆柱长度计算，轴向长度会在通量与容量中约去。
        volumes = np.pi * (outer_faces**2 - inner_faces**2)
        interface_radii = 0.5 * (nodes[:-1] + nodes[1:])
        interface_areas = 2.0 * np.pi * interface_radii
        surface_area = 2.0 * np.pi * radius_m

        return RadialGeometry(
            radius_m=radius_m,
            nodes_m=nodes,
            inner_faces_m=inner_faces,
            outer_faces_m=outer_faces,
            volumes_per_length_m2=volumes,
            interface_areas_per_length_m=interface_areas,
            surface_area_per_length_m=surface_area,
            spacing_m=spacing,
        )
