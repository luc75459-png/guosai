# Heun 推进和事件定位
"""提供 Heun 时间推进器、动态稳定步长和含水率阈值事件定位。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


RightHandSide = Callable[[float, np.ndarray], np.ndarray]
EventFunction = Callable[[float, np.ndarray], float]
MaxTimeStepFunction = Callable[[float, np.ndarray], float]


@dataclass(frozen=True)
class SimulationResult:
    """保存指定时刻的完整状态，以及可选的连续事件时刻和事件状态。"""

    time_s: np.ndarray
    state: np.ndarray
    event_time_s: float | None = None
    event_state: np.ndarray | None = None
    step_s: np.ndarray | None = None


def integrate_heun(
    rhs: RightHandSide,
    initial_state: np.ndarray,
    end_time_s: float,
    dt_s: float,
    save_every_s: float,
    event: EventFunction | None = None,
    max_dt_fn: MaxTimeStepFunction | None = None,
) -> SimulationResult:
    """用二阶显式 Heun 方法积分半离散方程。

    ``dt_s`` 表示最大内部时间步。传入 ``max_dt_fn`` 后，每个时间步还会
    查询当前状态允许的稳定步长，并取两者中的较小值。输出时刻仍严格落在
    ``save_every_s`` 的整数倍上。

    若事件函数由正变为非正，使用跨步线性插值估计首次越过阈值的连续时刻，
    并立即停止积分。
    """

    if end_time_s <= 0 or dt_s <= 0 or save_every_s <= 0:
        raise ValueError("终止时间、内部时间步和保存间隔都必须为正数")
    if max_dt_fn is None:
        save_ratio = save_every_s / dt_s
        if not np.isclose(save_ratio, round(save_ratio), atol=1.0e-10):
            raise ValueError("固定步长模式下，保存间隔必须是 dt_s 的整数倍")

    state = np.asarray(initial_state, dtype=float).copy()
    saved_times = [0.0]
    saved_states = [state.copy()]
    accepted_steps: list[float] = []

    previous_event_value = event(0.0, state) if event is not None else None
    event_time = None
    event_state = None
    if previous_event_value is not None and previous_event_value <= 0.0:
        return SimulationResult(
            time_s=np.asarray(saved_times),
            state=np.vstack(saved_states),
            event_time_s=0.0,
            event_state=state.copy(),
            step_s=np.asarray(accepted_steps),
        )

    time_s = 0.0
    next_save_time = save_every_s
    time_tolerance = max(1.0e-12, 1.0e-10 * max(1.0, end_time_s))

    def append_unique(time_value: float, state_value: np.ndarray) -> None:
        """追加保存点，并避免事件与输出时刻完全重合时重复。"""

        if saved_times and np.isclose(
            saved_times[-1],
            time_value,
            atol=time_tolerance,
            rtol=0.0,
        ):
            saved_times[-1] = float(time_value)
            saved_states[-1] = np.asarray(state_value, dtype=float).copy()
            return
        saved_times.append(float(time_value))
        saved_states.append(np.asarray(state_value, dtype=float).copy())

    while time_s < end_time_s - time_tolerance:
        requested_dt = dt_s
        if max_dt_fn is not None:
            stable_dt = float(max_dt_fn(time_s, state))
            if np.isnan(stable_dt) or stable_dt <= 0.0:
                raise ValueError("稳定步长必须为正数")
            requested_dt = min(requested_dt, stable_dt)

        next_target = (
            next_save_time
            if next_save_time < end_time_s - time_tolerance
            else end_time_s
        )
        actual_dt = min(requested_dt, next_target - time_s)
        if actual_dt <= time_tolerance:
            if next_save_time <= time_s + time_tolerance:
                append_unique(time_s, state)
                while next_save_time <= time_s + time_tolerance:
                    next_save_time += save_every_s
                continue
            raise RuntimeError("时间步推进失败：实际步长过小")

        # Heun 方法：先用欧拉预测，再用起点和预测点斜率的平均值校正。
        k1 = rhs(time_s, state)
        predictor = state + actual_dt * k1
        k2 = rhs(time_s + actual_dt, predictor)
        next_state = state + 0.5 * actual_dt * (k1 + k2)
        next_time = time_s + actual_dt

        if not np.all(np.isfinite(next_state)):
            raise FloatingPointError(f"在 t={next_time:.6g} s 出现非有限状态")

        if event is not None:
            next_event_value = event(next_time, next_state)
            if previous_event_value is not None and previous_event_value > 0 >= next_event_value:
                # 在线性近似下定位真正的阈值穿越时刻，避免只报告离散步终点。
                fraction = previous_event_value / (previous_event_value - next_event_value)
                event_time = time_s + fraction * actual_dt
                event_state = state + fraction * (next_state - state)
                accepted_steps.append(actual_dt)
                append_unique(event_time, event_state)
                break
            previous_event_value = next_event_value

        state = next_state
        time_s = next_time
        accepted_steps.append(actual_dt)

        if time_s >= next_save_time - time_tolerance:
            append_unique(time_s, state)
            while next_save_time <= time_s + time_tolerance:
                next_save_time += save_every_s
        elif time_s >= end_time_s - time_tolerance:
            append_unique(time_s, state)

    return SimulationResult(
        time_s=np.asarray(saved_times),
        state=np.vstack(saved_states),
        event_time_s=event_time,
        event_state=event_state,
        step_s=np.asarray(accepted_steps),
    )
