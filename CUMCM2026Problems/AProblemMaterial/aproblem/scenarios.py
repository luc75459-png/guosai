#四问组装
"""定义四个问题的差异，并把数据、物性、模型和积分器组装起来。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ProjectPaths, SimulationConfig
from .data import load_environment, load_radius
from .grid import RadialGrid
from .integrator import SimulationResult, integrate_heun
from .model import DryingModel
from .physics import Question1Properties, Question23Properties, Question4Properties


@dataclass(frozen=True)
class QuestionSpec:
    """单个问题的终止时间、保存间隔、收缩和阈值事件配置。"""

    end_time_s: float
    save_every_s: float
    shrinking: bool
    stop_at_threshold: bool


QUESTION_SPECS = {
    1: QuestionSpec(1800.0, 1.0, False, False),
    2: QuestionSpec(10800.0, 1.0, False, False),
    3: QuestionSpec(259200.0, 60.0, False, True),
    4: QuestionSpec(259200.0, 60.0, True, True),
}


def _property_model(question: int):
    """根据问题编号选择题目规定的物性公式。"""

    if question == 1:
        return Question1Properties()
    if question in (2, 3):
        return Question23Properties()
    if question == 4:
        return Question4Properties()
    raise ValueError(f"不支持的问题编号：{question}")


def run_question(
    question: int,
    paths: ProjectPaths,
    config: SimulationConfig | None = None,
) -> tuple[DryingModel, SimulationResult]:
    """组装并运行指定问题，返回模型对象和完整仿真结果。"""

    if question not in QUESTION_SPECS:
        raise ValueError("问题编号 question 必须是 1、2、3、4 之一")
    cfg = config or SimulationConfig()
    cfg.validate()
    spec = QUESTION_SPECS[question]

    environment = load_environment(
        paths.environment_file,
        plateau_temperature_c=cfg.plateau_temperature_c,
        plateau_moisture=cfg.plateau_moisture,
    )
    # 只有问题4读取收缩半径；其余问题始终使用初始固定半径。
    radius_data = load_radius(paths.radius_file) if spec.shrinking else None
    radius_function = radius_data if radius_data is not None else lambda _time: cfg.radius_m

    grid = RadialGrid(cfg.radial_intervals)
    model = DryingModel(
        grid=grid,
        properties=_property_model(question),
        environment=environment,
        radius=radius_function,
        heat_transfer_coefficient=cfg.heat_transfer_coefficient,
        mass_transfer_coefficient=cfg.mass_transfer_coefficient,
        cylinder_length_m=cfg.cylinder_length_m,
        include_end_faces=cfg.include_end_faces,
    )
    node_count = model.node_count
    initial_state = np.concatenate(
        (
            np.full(node_count, cfg.initial_temperature_c),
            np.full(node_count, cfg.initial_moisture),
        )
    )

    event = None
    if spec.stop_at_threshold:
        # 使用全场最大含水率判断“药材各处均低于阈值”，不预设中心必为最大值。
        event = lambda _time, state: float(np.max(state[node_count:]) - cfg.moisture_threshold)

    result = integrate_heun(
        rhs=model.rhs,
        initial_state=initial_state,
        end_time_s=spec.end_time_s,
        dt_s=cfg.dt_s,
        save_every_s=spec.save_every_s,
        event=event,
        max_dt_fn=(
            lambda time_s, state: model.stable_time_step(
                time_s,
                state,
                safety_factor=cfg.dt_safety_factor,
            )
            if cfg.adaptive_dt
            else None
        ),
    )
    return model, result
