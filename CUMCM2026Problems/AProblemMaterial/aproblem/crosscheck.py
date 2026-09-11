# BDF 独立对拍
"""使用 SciPy 的 BDF 求解器对同一半离散方程进行独立时间推进。"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .integrator import EventFunction, RightHandSide, SimulationResult


def integrate_bdf(
    rhs: RightHandSide,
    initial_state: np.ndarray,
    end_time_s: float,
    event: EventFunction | None = None,
    relative_tolerance: float = 1.0e-7,
    absolute_tolerance: float = 1.0e-9,
) -> SimulationResult:
    """使用 BDF 方法积分同一个有限体积右端项。

    BDF 仅作为问题 4 的独立数值对照，不替代二阶显式 Heun 主求解器。
    若设置事件函数，采用其中由正穿负的根作为连续事件时刻。
    """

    if end_time_s <= 0.0:
        raise ValueError("终止时间必须为正数")
    if relative_tolerance <= 0.0 or absolute_tolerance <= 0.0:
        raise ValueError("BDF 相对和绝对容差必须为正数")

    event_function = None
    if event is not None:
        def terminal_event(time_s: float, state: np.ndarray) -> float:
            return float(event(time_s, state))

        terminal_event.terminal = True
        terminal_event.direction = -1.0
        event_function = terminal_event

    solution = solve_ivp(
        fun=rhs,
        t_span=(0.0, end_time_s),
        y0=np.asarray(initial_state, dtype=float),
        method="BDF",
        rtol=relative_tolerance,
        atol=absolute_tolerance,
        events=event_function,
    )
    if not solution.success:
        raise RuntimeError(f"BDF 对拍失败：{solution.message}")

    event_time_s = None
    event_state = None
    if solution.t_events and solution.t_events[0].size:
        event_time_s = float(solution.t_events[0][0])
        event_state = np.asarray(solution.y_events[0][0], dtype=float)

    return SimulationResult(
        time_s=np.asarray(solution.t, dtype=float),
        state=np.asarray(solution.y.T, dtype=float),
        event_time_s=event_time_s,
        event_state=event_state,
        step_s=np.diff(np.asarray(solution.t, dtype=float)),
    )
