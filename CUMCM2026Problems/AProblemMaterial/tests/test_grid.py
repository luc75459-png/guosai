"""径向有限体积网格的几何单元测试。"""

import numpy as np

from aproblem.grid import RadialGrid


def test_annular_control_volumes_cover_the_cross_section() -> None:
    """所有控制体面积之和应恰好覆盖完整圆形横截面。"""

    radius = 0.02
    geometry = RadialGrid(intervals=20).geometry(radius)
    assert len(geometry.nodes_m) == 21
    assert np.isclose(geometry.nodes_m[1] - geometry.nodes_m[0], 0.001)
    assert np.isclose(np.sum(geometry.volumes_per_length_m2), np.pi * radius**2)
    assert geometry.inner_faces_m[0] == 0.0
    assert geometry.outer_faces_m[-1] == radius
