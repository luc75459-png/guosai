#基础范围检查
"""对仿真结果执行有限性、非负性和物理范围等基础检查。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .integrator import SimulationResult
from .model import DryingModel
from .outputs import split_result


@dataclass(frozen=True)
class ValidationReport:
    """汇总一次仿真的基础数值与物理检查结果。"""

    all_finite: bool
    nonnegative_moisture: bool
    material_mesh_consistent: bool
    maximum_radius_error_m: float
    minimum_node_spacing_m: float
    minimum_step_s: float
    maximum_step_s: float
    temperature_bounds_c: tuple[float, float]
    moisture_bounds: tuple[float, float]
    final_max_moisture: float

    def __str__(self) -> str:
        """以中文输出验证摘要，便于终端阅读和记录。"""

        finite_text = "是" if self.all_finite else "否"
        nonnegative_text = "是" if self.nonnegative_moisture else "否"
        mesh_text = "是" if self.material_mesh_consistent else "否"
        return (
            f"数值全部有限：{finite_text}；含水率均非负：{nonnegative_text}；"
            f"材料网格一致：{mesh_text}；"
            f"最大外半径误差：{self.maximum_radius_error_m:.3e} m；"
            f"最小节点间距：{self.minimum_node_spacing_m:.6g} m；"
            f"实际步长范围：{self.minimum_step_s:.6g}～"
            f"{self.maximum_step_s:.6g} s；"
            f"温度范围：{self.temperature_bounds_c[0]:.4f}～"
            f"{self.temperature_bounds_c[1]:.4f} ℃；"
            f"含水率范围：{self.moisture_bounds[0]:.6f}～"
            f"{self.moisture_bounds[1]:.6f} kg/kg；"
            f"最终全场最大含水率：{self.final_max_moisture:.6f} kg/kg"
        )


def validate_result(model: DryingModel, result: SimulationResult) -> ValidationReport:
    """根据温度和含水率矩阵生成基础验证报告。"""

    temperature, moisture = split_result(model, result)
    maximum_radius_error = 0.0
    minimum_spacing = float("inf")
    material_mesh_consistent = True
    expected_fraction = np.linspace(0.0, 1.0, model.node_count)
    for time_s in result.time_s:
        radius_m = model.radius(float(time_s))
        geometry = model.grid.geometry(radius_m)
        maximum_radius_error = max(
            maximum_radius_error,
            abs(float(geometry.nodes_m[-1] - radius_m)),
        )
        spacing = np.diff(geometry.nodes_m)
        minimum_spacing = min(minimum_spacing, float(np.min(spacing)))
        if (
            not np.all(np.isfinite(geometry.nodes_m))
            or np.any(spacing <= 0.0)
            or not np.allclose(
                geometry.nodes_m / radius_m,
                expected_fraction,
                rtol=0.0,
                atol=1.0e-12,
            )
        ):
            material_mesh_consistent = False

    if result.step_s is None or result.step_s.size == 0:
        minimum_step = 0.0
        maximum_step = 0.0
    else:
        minimum_step = float(np.min(result.step_s))
        maximum_step = float(np.max(result.step_s))

    return ValidationReport(
        all_finite=bool(np.all(np.isfinite(result.state))),
        nonnegative_moisture=bool(np.all(moisture >= 0.0)),
        material_mesh_consistent=material_mesh_consistent,
        maximum_radius_error_m=maximum_radius_error,
        minimum_node_spacing_m=minimum_spacing,
        minimum_step_s=minimum_step,
        maximum_step_s=maximum_step,
        temperature_bounds_c=(float(np.min(temperature)), float(np.max(temperature))),
        moisture_bounds=(float(np.min(moisture)), float(np.max(moisture))),
        final_max_moisture=float(np.max(moisture[-1])),
    )
