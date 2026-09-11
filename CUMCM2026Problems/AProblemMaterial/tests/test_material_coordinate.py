"""问题 4 材料坐标、动态步长和正式输出测试。"""

from __future__ import annotations

import numpy as np
import openpyxl
from openpyxl import Workbook

from aproblem.crosscheck import integrate_bdf
from aproblem.grid import RadialGrid
from aproblem.integrator import SimulationResult, integrate_heun
from aproblem.model import DryingModel
from aproblem.outputs import map_material_profile_to_physical, write_result4_workbook
from aproblem.physics import Question4Properties


def test_uniform_material_geometry_scales_with_radius() -> None:
    """统一半径缩放应使长度、面积和体积分别按 R、R² 和 R² 缩放。"""

    grid = RadialGrid(intervals=20)
    full = grid.geometry(0.02)
    half = grid.geometry(0.01)

    assert np.allclose(half.nodes_m, 0.5 * full.nodes_m)
    assert np.allclose(
        half.interface_areas_per_length_m,
        0.5 * full.interface_areas_per_length_m,
    )
    assert np.allclose(
        half.volumes_per_length_m2,
        0.25 * full.volumes_per_length_m2,
    )
    assert np.isclose(sum(half.volumes_per_length_m2), np.pi * 0.01**2)


def test_question4_properties_match_appendix_formulas() -> None:
    """问题 4 物性应逐项匹配附录 4 的经验公式。"""

    moisture = np.array([2.55])
    temperature_c = np.array([28.0])
    material = Question4Properties().evaluate(temperature_c, moisture)
    c = moisture[0]
    temperature_k = temperature_c[0] + 273.15

    assert np.isclose(material.density[0], 760.0 + 90.0 * c)
    assert np.isclose(material.heat_capacity[0], 1850.0 + 2150.0 * c / (c + 1.0))
    assert np.isclose(material.conductivity[0], 0.12 + 0.20 * c / (c + 1.0))
    assert np.isclose(
        material.diffusivity[0],
        4.2e-4 * np.exp(-0.30 / c) * np.exp(-3850.0 / temperature_k),
    )


def test_stable_time_step_decreases_on_finer_material_grid() -> None:
    """更细的径向网格应给出更严格的显式稳定步长。"""

    def build_model(intervals: int) -> DryingModel:
        return DryingModel(
            grid=RadialGrid(intervals),
            properties=Question4Properties(),
            environment=lambda _time: (50.0, 0.05),
            radius=lambda _time: 0.02,
            heat_transfer_coefficient=25.0,
            mass_transfer_coefficient=8.0e-7,
        )

    state_n20 = np.concatenate((np.full(21, 50.0), np.full(21, 0.15)))
    state_n40 = np.concatenate((np.full(41, 50.0), np.full(41, 0.15)))
    stable_n20 = build_model(20).stable_time_step(0.0, state_n20)
    stable_n40 = build_model(40).stable_time_step(0.0, state_n40)

    assert stable_n20 > 0.0
    assert stable_n40 > 0.0
    assert stable_n40 < stable_n20


def test_adaptive_heun_respects_dynamic_step_limit_and_save_times() -> None:
    """动态步长必须服从稳定性上限，同时精确保留指定输出时刻。"""

    result = integrate_heun(
        rhs=lambda _time, state: -state,
        initial_state=np.array([1.0]),
        end_time_s=1.0,
        dt_s=0.4,
        save_every_s=0.2,
        max_dt_fn=lambda _time, _state: 0.05,
    )

    assert result.step_s is not None
    assert np.all(result.step_s <= 0.05 + 1.0e-12)
    assert np.allclose(result.time_s, np.linspace(0.0, 1.0, 6))
    assert np.isclose(result.state[-1, 0], np.exp(-1.0), rtol=2.0e-2)


def test_bdf_crosscheck_locates_same_exponential_event() -> None:
    """BDF 独立对照应能找到与主求解器相同的连续事件时刻。"""

    result = integrate_bdf(
        rhs=lambda _time, state: -state,
        initial_state=np.array([1.0]),
        end_time_s=2.0,
        event=lambda _time, state: float(state[0] - 0.5),
    )

    assert result.event_time_s is not None
    assert np.isclose(result.event_time_s, np.log(2.0), rtol=1.0e-6)


def test_material_profile_mapping_leaves_region_outside_radius_empty() -> None:
    """超过当前药材半径的固定物理位置应映射为 nan。"""

    field = np.array([1.0, 0.8, 0.5])
    mapped, surface = map_material_profile_to_physical(field, radius_m=0.01)

    assert mapped[0] == 1.0
    assert np.isclose(mapped[10], 0.5)
    assert np.all(np.isnan(mapped[11:]))
    assert np.isclose(surface, 0.5)


def test_result4_workbook_uses_template_and_writes_rounded_values(tmp_path) -> None:
    """正式工作簿应保留模板样式、留空药材外区域并保留四位小数。"""

    template = tmp_path / "result4_template.xlsx"
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "Sheet1"
    worksheet.cell(row=1, column=1, value="时间\\到药材中心的距离")
    worksheet.cell(row=1, column=2, value=0.0)
    worksheet.cell(row=2, column=1, value=60)
    worksheet.cell(row=2, column=2, value=0.0)
    workbook.save(template)
    workbook.close()

    model = DryingModel(
        grid=RadialGrid(2),
        properties=Question4Properties(),
        environment=lambda _time: (50.0, 0.05),
        radius=lambda _time: 0.01,
        heat_transfer_coefficient=25.0,
        mass_transfer_coefficient=8.0e-7,
    )
    result = SimulationResult(
        time_s=np.array([0.0, 60.0, 120.0]),
        state=np.array(
            [
                [50.0, 50.0, 50.0, 1.23456, 1.23456, 1.23456],
                [50.0, 50.0, 50.0, 0.30001, 0.30001, 0.30001],
                [50.0, 50.0, 50.0, 0.15004, 0.15004, 0.15004],
            ]
        ),
        event_time_s=120.75,
    )

    output = write_result4_workbook(
        template,
        tmp_path / "result4.xlsx",
        model,
        result,
    )
    written = openpyxl.load_workbook(output)
    worksheet = written.active

    assert worksheet.cell(row=1, column=22).value == "药材表面"
    assert worksheet.cell(row=1, column=8).value == 0.6
    assert worksheet.cell(row=2, column=1).value == 60
    assert worksheet.cell(row=2, column=2).value == 0.3
    assert worksheet.cell(row=2, column=12).value == 0.3
    assert worksheet.cell(row=2, column=13).value is None
    assert worksheet.cell(row=3, column=1).value == 120
    assert worksheet.cell(row=4, column=1).value == 120.75
    assert worksheet.cell(row=4, column=22).value == 0.15
    assert worksheet.cell(row=2, column=2).number_format == "0.0000"
    written.close()
