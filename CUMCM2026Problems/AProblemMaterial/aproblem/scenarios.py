#四问组装
"""定义四个问题的差异，并把数据、物性、模型和积分器组装起来。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import (
    ProjectPaths,
    SimulationConfig,
    default_interface_strategy,
)
from .data import load_environment, load_radius
from .grid import RadialGrid
from .integrator import EventFunction, SimulationResult, integrate_heun
from .interfaces import resolve_interface_strategy
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


@dataclass(frozen=True)
class QuestionSetup:
    """一个问题在积分之前的完整组装结果。"""

    spec: QuestionSpec
    model: DryingModel
    initial_state: np.ndarray
    threshold_event: EventFunction | None

    def max_time_step(self, config: SimulationConfig):
        """返回可直接传给积分器的动态步长上限函数；关闭自适应时为 None。"""

        if not config.adaptive_dt:
            return None
        return lambda time_s, state: self.model.stable_time_step(
            time_s,
            state,
            safety_factor=config.dt_safety_factor,
        )


def build_question(
    question: int,
    paths: ProjectPaths,
    config: SimulationConfig | None = None,
) -> QuestionSetup:
    """组装指定问题的模型、初值和终止事件，但不进行时间积分。

    把组装与积分分开，是为了让收敛性研究可以对同一套模型自由选择
    Heun 或 BDF 两种时间推进方式。
    """

    if question not in QUESTION_SPECS:
        raise ValueError("问题编号 question 必须是 1、2、3、4 之一")
    cfg = config or SimulationConfig()
    cfg.validate()
    spec = QUESTION_SPECS[question]

    environment = load_environment(
        paths.environment_file,
        plateau_temperature_c=cfg.plateau_temperature_c,
        plateau_moisture=cfg.plateau_moisture,
        hold_last_value=(cfg.plateau_mode == "last_value"),
    )
    # 只有问题4读取收缩半径；其余问题始终使用初始固定半径。
    radius_data = load_radius(paths.radius_file) if spec.shrinking else None
    radius_function = radius_data if radius_data is not None else lambda _time: cfg.radius_m

    # 长度模式只对收缩问题有意义：各向同性假设下长度按同一比例缩短。
    length_function = None
    if spec.shrinking and cfg.length_scaling == "isotropic":
        initial_radius_m = radius_function(0.0)
        length_function = (
            lambda time_s: cfg.cylinder_length_m
            * radius_function(time_s)
            / initial_radius_m
        )

    grid = RadialGrid(cfg.radial_intervals)
    # 未显式指定界面策略时按题号取默认（问题2/3 默认基尔霍夫通量势）。
    interface_strategy_name = (
        cfg.interface_strategy
        if cfg.interface_strategy is not None
        else default_interface_strategy(question)
    )
    model = DryingModel(
        grid=grid,
        properties=_property_model(question),
        environment=environment,
        radius=radius_function,
        heat_transfer_coefficient=cfg.heat_transfer_coefficient,
        mass_transfer_coefficient=cfg.mass_transfer_coefficient,
        cylinder_length_m=cfg.cylinder_length_m,
        include_end_faces=cfg.include_end_faces,
        interface_strategy=resolve_interface_strategy(interface_strategy_name),
        length_function=length_function,
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

    return QuestionSetup(
        spec=spec,
        model=model,
        initial_state=initial_state,
        threshold_event=event,
    )


def run_question(
    question: int,
    paths: ProjectPaths,
    config: SimulationConfig | None = None,
) -> tuple[DryingModel, SimulationResult]:
    """组装并用 Heun 方法运行指定问题，返回模型对象和完整仿真结果。"""

    cfg = config or SimulationConfig()
    setup = build_question(question, paths, cfg)
    result = integrate_heun(
        rhs=setup.model.rhs,
        initial_state=setup.initial_state,
        end_time_s=setup.spec.end_time_s,
        dt_s=cfg.dt_s,
        save_every_s=setup.spec.save_every_s,
        event=setup.threshold_event,
        max_dt_fn=setup.max_time_step(cfg),
    )
    return setup.model, result
