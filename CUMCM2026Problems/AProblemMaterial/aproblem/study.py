# 收敛性与敏感性批量研究驱动
"""批量运行问题 4 的网格、时间步、界面策略和模型假设对比。

本模块只负责"跑一批算例并把结果整理成表格"，不参与正式结果导出。
所有结论都来自同一套半离散模型，只是替换网格、界面策略、积分器或
模型假设，便于把数值误差和模型假设的不确定度分开观察。
"""

from __future__ import annotations

import argparse
import csv
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .axisymmetric import AxisymmetricModel
from .config import (
    ProjectPaths,
    SimulationConfig,
    default_interface_strategy,
)
from .crosscheck import integrate_bdf
from .data import load_environment, load_radius
from .integrator import integrate_heun
from .interfaces import INTERFACE_STRATEGY_NAMES, resolve_interface_strategy
from .outputs import split_result
from .physics import Question1Properties, Question23Properties, Question4Properties
from .scenarios import QUESTION_SPECS, build_question


@dataclass(frozen=True)
class CaseOutcome:
    """单个算例的关键结果，用于生成收敛表和敏感性表。"""

    label: str
    question: int
    intervals: int
    interface_strategy: str
    include_end_faces: bool
    length_scaling: str
    plateau_mode: str
    integrator: str
    end_time_s: float
    event_time_s: float
    event_hours: float
    # 达标事件是否真的在终止时刻前触发；BDF 触发时末时刻正好等于事件时刻，
    # 因此不能用 event_time_s 与 end_time_s 是否相等来判断。
    event_reached: bool
    elapsed_s: float
    step_count: int
    min_step_s: float
    max_step_s: float
    final_max_moisture: float
    moisture_min: float
    temperature_max_c: float
    moisture_monotone: bool

    def as_row(self) -> dict[str, object]:
        """转换成可直接写入 CSV 的键值对。"""

        return asdict(self)


def _monotone_decreasing(values: np.ndarray, tolerance: float = 1.0e-9) -> bool:
    """判断序列是否单调不增（允许给定容差内的数值抖动）。"""

    return bool(np.all(np.diff(values) <= tolerance))


def run_case(
    *,
    question: int,
    intervals: int,
    interface_strategy: str,
    paths: ProjectPaths,
    include_end_faces: bool = False,
    length_scaling: str = "constant",
    plateau_mode: str = "fixed",
    integrator: str = "heun",
    dt_cap_s: float = 0.5,
    label: str = "",
) -> CaseOutcome:
    """运行单个算例并返回关键指标。"""

    config = SimulationConfig(
        radial_intervals=intervals,
        dt_s=dt_cap_s,
        adaptive_dt=(integrator == "heun"),
        include_end_faces=include_end_faces,
        interface_strategy=interface_strategy,
        length_scaling=length_scaling,
        plateau_mode=plateau_mode,
    )
    setup = build_question(question, paths, config)
    spec = QUESTION_SPECS[question]

    started = time.perf_counter()
    if integrator == "heun":
        result = integrate_heun(
            rhs=setup.model.rhs,
            initial_state=setup.initial_state,
            end_time_s=spec.end_time_s,
            dt_s=config.dt_s,
            save_every_s=spec.save_every_s,
            event=setup.threshold_event,
            max_dt_fn=setup.max_time_step(config),
        )
    elif integrator == "bdf":
        result = integrate_bdf(
            setup.model.rhs,
            setup.initial_state,
            end_time_s=spec.end_time_s,
            event=setup.threshold_event,
        )
    else:
        raise ValueError(f"不支持的积分器：{integrator}")
    elapsed_s = time.perf_counter() - started

    _, moisture = split_result(setup.model, result)
    temperature, _ = split_result(setup.model, result)
    event_time_s = (
        float(result.event_time_s)
        if result.event_time_s is not None
        else float(result.time_s[-1])
    )
    event_reached = result.event_time_s is not None
    steps = np.asarray(result.step_s if result.step_s is not None else [])

    return CaseOutcome(
        label=label or f"N{intervals}-{interface_strategy}-{integrator}",
        question=question,
        intervals=intervals,
        interface_strategy=interface_strategy,
        include_end_faces=include_end_faces,
        length_scaling=length_scaling,
        plateau_mode=plateau_mode,
        integrator=integrator,
        end_time_s=float(result.time_s[-1]),
        event_time_s=event_time_s,
        event_hours=event_time_s / 3600.0,
        event_reached=event_reached,
        elapsed_s=elapsed_s,
        step_count=int(steps.size),
        min_step_s=float(np.min(steps)) if steps.size else 0.0,
        max_step_s=float(np.max(steps)) if steps.size else 0.0,
        final_max_moisture=float(np.max(moisture[-1])),
        moisture_min=float(np.min(moisture)),
        temperature_max_c=float(np.max(temperature)),
        moisture_monotone=_monotone_decreasing(moisture.max(axis=1)),
    )


def _relative_change(previous: float, current: float) -> float:
    """返回相邻网格结果相对前一档的变化率。"""

    if previous == 0.0:
        return float("nan")
    return abs(current - previous) / abs(previous)


def _markdown_table(header: list[str], rows: list[list[str]]) -> str:
    """把二维字符串拼成 Markdown 表格。"""

    lines = ["| " + " | ".join(header) + " |"]
    lines.append("|" + "|".join(["---"] * len(header)) + "|")
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def _write_outcomes(path: Path, outcomes: list[object]) -> None:
    """把算例结果写成 CSV，便于后续作图。"""

    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [asdict(outcome) for outcome in outcomes]
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def interface_scan(
    paths: ProjectPaths,
    grids: tuple[int, ...] = (20, 40, 80),
    integrator: str = "bdf",
    strategies: tuple[str, ...] = INTERFACE_STRATEGY_NAMES,
    question: int = 4,
    include_end_faces: bool = True,
) -> list[CaseOutcome]:
    """遍历"界面策略 × 网格"，输出达标时刻和相邻网格相对变化。"""

    outcomes: list[CaseOutcome] = []
    for strategy in strategies:
        for intervals in grids:
            outcome = run_case(
                question=question,
                intervals=intervals,
                interface_strategy=strategy,
                paths=paths,
                integrator=integrator,
                include_end_faces=include_end_faces,
                label=f"{strategy}-N{intervals}",
            )
            print(
                f"[interface-scan] {outcome.label}: "
                f"{outcome.event_hours:.4f} h，耗时 {outcome.elapsed_s:.1f} s",
                flush=True,
            )
            outcomes.append(outcome)
    return outcomes


SENSITIVITY_CASES: tuple[tuple[str, dict[str, str]], ...] = (
    ("主线：长度恒定 + 50 ℃/0.05 保持", {"length_scaling": "constant", "plateau_mode": "fixed"}),
    ("敏感性①：长度与半径同步收缩", {"length_scaling": "isotropic", "plateau_mode": "fixed"}),
    ("敏感性②：4 h 后保持附件末值", {"length_scaling": "constant", "plateau_mode": "last_value"}),
    ("敏感性③：长度收缩且保持附件末值", {"length_scaling": "isotropic", "plateau_mode": "last_value"}),
)


def run_sensitivity(
    paths: ProjectPaths,
    *,
    intervals: int = 40,
    interface_strategy: str = "midpoint",
    integrator: str = "bdf",
) -> list[CaseOutcome]:
    """运行长度假设和 4 h 后边界延拓两项敏感性分析。"""

    outcomes: list[CaseOutcome] = []
    for label, options in SENSITIVITY_CASES:
        outcome = run_case(
            question=4,
            intervals=intervals,
            interface_strategy=interface_strategy,
            paths=paths,
            integrator=integrator,
            include_end_faces=True,
            label=label,
            **options,
        )
        print(
            f"[sensitivity] {label}: {outcome.event_hours:.4f} h，"
            f"耗时 {outcome.elapsed_s:.1f} s",
            flush=True,
        )
        outcomes.append(outcome)
    return outcomes


def render_sensitivity(outcomes: list[CaseOutcome]) -> str:
    """把敏感性结果渲染成相对主线的偏差表。"""

    baseline = outcomes[0].event_time_s
    header = ["工况", "达标时刻", "相对主线偏差", "长度模式", "边界延拓"]
    rows = []
    for outcome in outcomes:
        rows.append(
            [
                outcome.label,
                f"{outcome.event_hours:.4f} h",
                f"{(outcome.event_time_s - baseline) / baseline:+.3%}",
                outcome.length_scaling,
                outcome.plateau_mode,
            ]
        )
    return _markdown_table(header, rows)


def render_time_convergence(outcomes: list[CaseOutcome]) -> str:
    """把时间步收敛结果渲染成表，最后一列是与最细步长的相对偏差。"""

    finest = outcomes[-1].event_time_s
    header = ["内部步长上限", "积分器", "达标时刻", "相对最细步长偏差", "实际步数", "耗时"]
    rows = []
    for outcome in outcomes:
        rows.append(
            [
                f"{outcome.label}",
                outcome.integrator,
                f"{outcome.event_hours:.4f} h",
                f"{(outcome.event_time_s - finest) / finest:+.4%}",
                f"{outcome.step_count}",
                f"{outcome.elapsed_s:.1f} s",
            ]
        )
    return _markdown_table(header, rows)


def render_interface_scan(outcomes: list[CaseOutcome], grids: tuple[int, ...]) -> str:
    """把界面策略扫描结果渲染成 Markdown 收敛表。"""

    strategies = []
    for outcome in outcomes:
        if outcome.interface_strategy not in strategies:
            strategies.append(outcome.interface_strategy)

    header = ["界面策略"] + [f"N={n}" for n in grids] + [
        f"N={grids[-2]}→{grids[-1]} 相对变化" if len(grids) > 1 else "相对变化"
    ]
    rows: list[list[str]] = []

    def event_cell(outcome: CaseOutcome) -> str:
        """达标事件在终止时刻前触发才给出定值，否则标注未达标。"""

        if outcome.event_reached:
            return f"{outcome.event_hours:.4f} h"
        return f">{outcome.event_hours:.4f} h（未达标）"

    for strategy in strategies:
        by_grid = {
            outcome.intervals: outcome
            for outcome in outcomes
            if outcome.interface_strategy == strategy
        }
        cells = [
            event_cell(by_grid[n]) if n in by_grid else "—"
            for n in grids
        ]
        both_reached = (
            len(grids) > 1
            and grids[-1] in by_grid
            and grids[-2] in by_grid
            and by_grid[grids[-1]].event_reached
            and by_grid[grids[-2]].event_reached
        )
        if both_reached:
            change = _relative_change(
                by_grid[grids[-2]].event_time_s,
                by_grid[grids[-1]].event_time_s,
            )
            cells.append(f"{change:.4%}")
        else:
            cells.append("—")
        rows.append([strategy, *cells])
    return _markdown_table(header, rows)


def property_model_for(question: int):
    """返回指定问题的物性模型实例。"""

    if question == 1:
        return Question1Properties()
    if question in (2, 3):
        return Question23Properties()
    if question == 4:
        return Question4Properties()
    raise ValueError(f"不支持的问题编号：{question}")


def build_axisymmetric_model(
    question: int,
    paths: ProjectPaths,
    config: SimulationConfig,
    radial_intervals: int,
    axial_intervals: int,
    insulated_ends: bool = False,
) -> AxisymmetricModel:
    """按配置组装二维轴对称模型。

    只有收缩问题（问题4）读附件2 的半径演化；问题1~3 半径固定为
    ``config.radius_m``。长度按 ``length_scaling`` 取恒定值或与半径同步
    缩短。二维模型里端面是真实边界，不再使用 ``2h/L`` 体积源项。
    """

    environment = load_environment(
        paths.environment_file,
        plateau_temperature_c=config.plateau_temperature_c,
        plateau_moisture=config.plateau_moisture,
        hold_last_value=(config.plateau_mode == "last_value"),
    )
    if QUESTION_SPECS[question].shrinking:
        radius_data = load_radius(paths.radius_file)
    else:
        radius_data = None
    initial_radius_m = (
        float(radius_data(0.0)) if radius_data is not None else config.radius_m
    )
    if radius_data is not None and config.length_scaling == "isotropic":
        length_function = (
            lambda time_s: config.cylinder_length_m
            * float(radius_data(time_s))
            / initial_radius_m
        )
    else:
        length_function = lambda _time: config.cylinder_length_m
    radius_function = (
        radius_data
        if radius_data is not None
        else (lambda _time: config.radius_m)
    )

    return AxisymmetricModel(
        radial_intervals=radial_intervals,
        axial_intervals=axial_intervals,
        properties=property_model_for(question),
        environment=environment,
        radius=radius_function,
        length=length_function,
        heat_transfer_coefficient=config.heat_transfer_coefficient,
        mass_transfer_coefficient=config.mass_transfer_coefficient,
        interface_strategy=resolve_interface_strategy(
            config.interface_strategy
            if config.interface_strategy is not None
            else default_interface_strategy(question)
        ),
        insulated_ends=insulated_ends,
    )


@dataclass(frozen=True)
class AxisymmetricOutcome:
    """二维轴对称算例的关键结果。"""

    label: str
    radial_intervals: int
    axial_intervals: int
    dt_cap_s: float
    insulated_ends: bool
    event_time_s: float
    event_hours: float
    elapsed_s: float
    step_count: int
    min_step_s: float
    max_step_s: float
    final_max_moisture: float
    axial_spread_final: float
    centre_moisture_final: float
    surface_moisture_final: float


def run_axisymmetric_case(
    *,
    question: int,
    radial_intervals: int,
    axial_intervals: int,
    paths: ProjectPaths,
    interface_strategy: str = "midpoint",
    length_scaling: str = "constant",
    plateau_mode: str = "fixed",
    dt_cap_s: float = 0.5,
    insulated_ends: bool = False,
    label: str = "",
    output_path: Path | None = None,
) -> tuple[AxisymmetricOutcome, AxisymmetricModel, np.ndarray, np.ndarray, np.ndarray]:
    """运行单个二维轴对称算例。

    返回算例指标、模型、时间序列和温度、含水率三维数组，便于后续
    计算轴向平均剖面和绘制空间分布。
    """

    config = SimulationConfig(
        radial_intervals=radial_intervals,
        dt_s=dt_cap_s,
        adaptive_dt=True,
        interface_strategy=interface_strategy,
        length_scaling=length_scaling,
        plateau_mode=plateau_mode,
        include_end_faces=True,
    )
    model = build_axisymmetric_model(
        question,
        paths,
        config,
        radial_intervals,
        axial_intervals,
        insulated_ends=insulated_ends,
    )
    initial_temperature = np.full(model.shape, config.initial_temperature_c)
    initial_moisture = np.full(model.shape, config.initial_moisture)
    initial_state = model.pack(initial_temperature, initial_moisture)
    spec = QUESTION_SPECS[question]

    def event(_time_s: float, state: np.ndarray) -> float:
        _, moisture = model.split(state)
        return float(np.max(moisture) - config.moisture_threshold)

    started = time.perf_counter()
    result = integrate_heun(
        rhs=model.rhs,
        initial_state=initial_state,
        end_time_s=spec.end_time_s,
        dt_s=dt_cap_s,
        save_every_s=spec.save_every_s,
        event=event,
        max_dt_fn=lambda time_s, state: model.stable_time_step(
            time_s,
            state,
            safety_factor=config.dt_safety_factor,
        ),
    )
    elapsed_s = time.perf_counter() - started

    states = np.asarray(result.state, dtype=float)
    node_count = model.node_count
    temperature = states[:, :node_count].reshape(
        (states.shape[0], *model.shape)
    )
    moisture = states[:, node_count:].reshape((states.shape[0], *model.shape))
    radii_m = np.asarray([model.radius(float(value)) for value in result.time_s])
    lengths_m = np.asarray([model.length(float(value)) for value in result.time_s])
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            output_path,
            time_s=result.time_s,
            temperature_c=temperature,
            moisture=moisture,
            radius_m=radii_m,
            length_m=lengths_m,
            event_time_s=np.nan if result.event_time_s is None else result.event_time_s,
        )
    steps = np.asarray(result.step_s if result.step_s is not None else [])
    event_time_s = (
        float(result.event_time_s)
        if result.event_time_s is not None
        else float(result.time_s[-1])
    )
    final_moisture = moisture[-1]
    axial_mean = final_moisture.mean(axis=1, keepdims=True)

    outcome = AxisymmetricOutcome(
        label=label or f"N{radial_intervals}x{axial_intervals}",
        radial_intervals=radial_intervals,
        axial_intervals=axial_intervals,
        dt_cap_s=dt_cap_s,
        insulated_ends=insulated_ends,
        event_time_s=event_time_s,
        event_hours=event_time_s / 3600.0,
        elapsed_s=elapsed_s,
        step_count=int(steps.size),
        min_step_s=float(np.min(steps)) if steps.size else 0.0,
        max_step_s=float(np.max(steps)) if steps.size else 0.0,
        final_max_moisture=float(np.max(final_moisture)),
        axial_spread_final=float(np.max(np.abs(final_moisture - axial_mean))),
        centre_moisture_final=float(final_moisture[0].mean()),
        surface_moisture_final=float(final_moisture[-1].mean()),
    )
    return outcome, model, result.time_s, temperature, moisture


def parse_args() -> argparse.Namespace:
    """解析研究脚本的命令行参数。"""

    parser = argparse.ArgumentParser(description="问题4收敛性与敏感性批量研究")
    parser.add_argument(
        "study",
        choices=(
            "interface-scan",
            "grid-convergence",
            "time-convergence",
            "sensitivity",
            "axisymmetric-scan",
        ),
        help="要执行的研究类型",
    )
    parser.add_argument("--attachments", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--intervals",
        type=int,
        nargs="+",
        default=None,
        help="要测试的径向区间数列表",
    )
    parser.add_argument("--radial-intervals", type=int, default=40)
    parser.add_argument(
        "--axial-intervals",
        type=int,
        nargs="+",
        default=None,
        help="二维模型要测试的轴向区间数列表",
    )
    parser.add_argument(
        "--interface-strategy",
        default="all",
        help="界面策略名，或 all 表示遍历全部三套策略",
    )
    parser.add_argument("--integrator", choices=("heun", "bdf"), default="bdf")
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument(
        "--dt-list",
        type=float,
        nargs="+",
        default=None,
        help="时间收敛研究要测试的内部步长上限列表",
    )
    parser.add_argument(
        "--length-scaling",
        choices=("constant", "isotropic"),
        default="constant",
    )
    parser.add_argument(
        "--plateau-mode",
        choices=("fixed", "last_value"),
        default="fixed",
    )
    parser.add_argument(
        "--insulated-ends",
        action="store_true",
        help="二维模型把端面改成齐次 Neumann 边界，用于退化测试",
    )
    parser.add_argument(
        "--include-end-faces",
        action="store_true",
        help="一维 interface-scan 计入 2h/L 端面等效源；默认只算圆柱侧面",
    )
    parser.add_argument(
        "--ignore-end-faces",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--question",
        type=int,
        choices=(1, 2, 3, 4),
        default=4,
        help="研究针对的问题编号，默认为问题4",
    )
    return parser.parse_args()


def main() -> None:
    """执行指定的研究并把结果写入输出目录。"""

    args = parse_args()
    paths = ProjectPaths.discover(attachments=args.attachments)
    output_dir = (args.output_dir or paths.project_root / "outputs" / "study").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.study == "interface-scan":
        grids = tuple(args.intervals or (20, 40, 80))
        if args.interface_strategy == "all":
            strategies = INTERFACE_STRATEGY_NAMES
        elif args.interface_strategy in INTERFACE_STRATEGY_NAMES:
            strategies = (args.interface_strategy,)
        else:
            raise SystemExit(f"未知的界面策略：{args.interface_strategy}")
        outcomes = interface_scan(
            paths,
            grids=grids,
            integrator=args.integrator,
            strategies=strategies,
            question=args.question,
            include_end_faces=args.include_end_faces,
        )
        table = render_interface_scan(outcomes, grids)
        (output_dir / "interface_scan.md").write_text(
            "# 界面扩散系数策略 × 网格收敛表\n\n" + table + "\n",
            encoding="utf-8",
        )
        _write_outcomes(output_dir / "interface_scan.csv", outcomes)
        print(table)
    elif args.study == "axisymmetric-scan":
        axial_intervals = tuple(args.axial_intervals or (20, 40, 80))
        rows: list[list[str]] = []
        outcomes: list[AxisymmetricOutcome] = []
        for nz in axial_intervals:
            outcome, _, _, _, _ = run_axisymmetric_case(
                question=args.question,
                radial_intervals=args.radial_intervals,
                axial_intervals=nz,
                paths=paths,
                interface_strategy=args.interface_strategy,
                length_scaling=args.length_scaling,
                plateau_mode=args.plateau_mode,
                dt_cap_s=args.dt,
                insulated_ends=args.insulated_ends,
                label=f"N{args.radial_intervals}x{nz}",
                output_path=output_dir / f"axisym_N{args.radial_intervals}x{nz}.npz",
            )
            print(
                f"[axisymmetric-scan] {outcome.label}: {outcome.event_hours:.4f} h，"
                f"轴向极差 {outcome.axial_spread_final:.4f}，耗时 {outcome.elapsed_s:.1f} s",
                flush=True,
            )
            rows.append(
                [
                    outcome.label,
                    f"{outcome.event_hours:.4f} h",
                    f"{outcome.axial_spread_final:.5f}",
                    f"{outcome.centre_moisture_final:.4f}",
                    f"{outcome.surface_moisture_final:.4f}",
                    f"{outcome.step_count}",
                    f"{outcome.elapsed_s:.1f} s",
                ]
            )
            outcomes.append(outcome)
        table = _markdown_table(
            ["网格", "达标时刻", "轴向极差", "中心含水率", "表面含水率", "步数", "耗时"],
            rows,
        )
        (output_dir / "axisymmetric_scan.md").write_text(
            "# 二维轴对称端面校验\n\n" + table + "\n",
            encoding="utf-8",
        )
        _write_outcomes(output_dir / "axisymmetric_scan.csv", outcomes)
        print(table)
    elif args.study == "time-convergence":
        dt_list = tuple(args.dt_list or (0.5, 0.25, 0.125))
        # 时间步收敛必须用固定步长上限的主求解器 Heun；BDF 是自适应容差，
        # 没有"步长上限"这一档，所以不作为步长收敛的对照。
        integration = "heun"
        outcomes = []
        for dt_cap in dt_list:
            outcome = run_case(
                question=args.question,
                intervals=args.intervals[0] if args.intervals else 40,
                interface_strategy=(
                    "midpoint"
                    if args.interface_strategy == "all"
                    else args.interface_strategy
                ),
                paths=paths,
                integrator=integration,
                dt_cap_s=dt_cap,
                include_end_faces=not args.insulated_ends,
                label=f"dt≤{dt_cap:g} s",
            )
            print(
                f"[time-convergence] dt≤{dt_cap:g}: {outcome.event_hours:.4f} h，"
                f"步数 {outcome.step_count}，耗时 {outcome.elapsed_s:.1f} s",
                flush=True,
            )
            outcomes.append(outcome)
        table = render_time_convergence(outcomes)
        (output_dir / "time_convergence.md").write_text(
            "# 时间步收敛表\n\n" + table + "\n",
            encoding="utf-8",
        )
        _write_outcomes(output_dir / "time_convergence.csv", outcomes)
        print(table)
    elif args.study == "sensitivity":
        outcomes = run_sensitivity(
            paths,
            intervals=args.intervals[0] if args.intervals else 40,
            interface_strategy=(
                "midpoint"
                if args.interface_strategy == "all"
                else args.interface_strategy
            ),
            integrator=args.integrator,
        )
        table = render_sensitivity(outcomes)
        (output_dir / "sensitivity.md").write_text(
            "# 长度假设与边界延拓敏感性分析\n\n" + table + "\n",
            encoding="utf-8",
        )
        _write_outcomes(output_dir / "sensitivity.csv", outcomes)
        print(table)
    else:
        raise SystemExit(f"研究类型 {args.study!r} 尚未实现")


if __name__ == "__main__":
    main()
