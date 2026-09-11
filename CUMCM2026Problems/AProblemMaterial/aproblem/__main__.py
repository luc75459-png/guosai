"""命令行入口：组装参数、运行指定问题并生成预览结果与图像。"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .config import ProjectPaths, SimulationConfig
from .crosscheck import integrate_bdf
from .interfaces import INTERFACE_STRATEGY_NAMES
from .outputs import (
    write_preview_files,
    write_result4_workbook,
    write_result_workbook,
)
from .plotting import plot_final_profiles
from .scenarios import QUESTION_SPECS, run_question
from .validation import (
    validate_result,
    verify_result4_workbook,
    verify_result_workbook,
)


def parse_args() -> argparse.Namespace:
    """解析命令行参数，并限制题号只能为 1～4。"""

    parser = argparse.ArgumentParser(description="CUMCM 2026 A题药材烘干数值模拟")
    parser.add_argument(
        "--question",
        type=int,
        choices=(1, 2, 3, 4),
        default=1,
        help="要计算的问题编号，默认为问题1",
    )
    parser.add_argument(
        "--attachments",
        type=Path,
        default=None,
        help="附件目录；不填写时使用项目约定的默认目录",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="结果输出目录；不填写时使用项目内的 outputs 目录",
    )
    parser.add_argument("--dt", type=float, default=0.5, help="内部时间步长（秒）")
    parser.add_argument(
        "--intervals",
        type=int,
        default=20,
        help="径向有限体积区间数；生产计算建议使用 40",
    )
    parser.add_argument(
        "--fixed-dt",
        action="store_true",
        help="关闭动态稳定步长，严格使用 --dt 指定的固定步长",
    )
    parser.add_argument(
        "--no-result",
        "--no-result4",
        dest="no_result",
        action="store_true",
        help="只生成预览文件，不写正式 resultN.xlsx",
    )
    parser.add_argument(
        "--bdf-check",
        action="store_true",
        help="使用 BDF 对同一模型做独立时间积分并比较事件时刻",
    )
    parser.add_argument(
        "--include-end-faces",
        action="store_true",
        help="计入两个端面的 2h/L 等效源；主模型默认不计入，只算圆柱侧面",
    )
    parser.add_argument(
        "--ignore-end-faces",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--interface-strategy",
        choices=INTERFACE_STRATEGY_NAMES,
        default=None,
        help=(
            "水界面上扩散系数的取值方式；不填写时按题号取默认"
            "（问题1/4 为 midpoint，问题2/3 为 kirchhoff）"
        ),
    )
    parser.add_argument(
        "--length-scaling",
        choices=("constant", "isotropic"),
        default="constant",
        help="药材长度模式：constant 保持 0.25 m，isotropic 与半径同步收缩",
    )
    parser.add_argument(
        "--plateau-mode",
        choices=("fixed", "last_value"),
        default="fixed",
        help="附件时间范围以外的烘房条件：fixed 用 50 ℃/0.05，last_value 保持附件末值",
    )
    return parser.parse_args()


def main() -> None:
    """执行一次完整仿真，并输出结果路径和基础物理检查报告。"""

    args = parse_args()
    paths = ProjectPaths.discover(attachments=args.attachments)
    output_dir = (args.output_dir or paths.outputs).resolve()
    config = SimulationConfig(
        dt_s=args.dt,
        radial_intervals=args.intervals,
        adaptive_dt=not args.fixed_dt,
        include_end_faces=args.include_end_faces and not args.ignore_end_faces,
        interface_strategy=args.interface_strategy,
        length_scaling=args.length_scaling,
        plateau_mode=args.plateau_mode,
    )
    config.validate()

    # 四个问题共用相同入口，只在场景组装阶段切换物性、时长和半径函数。
    model, result = run_question(args.question, paths, config)
    files = list(write_preview_files(output_dir, args.question, model, result))
    if not args.no_result:
        if args.question == 4:
            files.append(
                write_result4_workbook(
                    paths.result_template(4),
                    output_dir / "result4.xlsx",
                    model,
                    result,
                )
            )
        else:
            files.append(
                write_result_workbook(
                    paths.result_template(args.question),
                    output_dir / f"result{args.question}.xlsx",
                    args.question,
                    model,
                    result,
                )
            )
    bdf_result = None
    if args.bdf_check:
        event = None
        # 问题3/4 都以"全场最大含水率降到阈值"为终止事件；
        # 问题1/2 是定时长输出，没有事件。
        if args.question in (3, 4):
            event = lambda _time, state: float(
                np.max(state[model.node_count :]) - config.moisture_threshold
            )
        bdf_result = integrate_bdf(
            model.rhs,
            np.asarray(result.state[0], dtype=float),
            # 用该题的最大时间窗，而不是 Heun 的事件时刻：否则 BDF 的事件
            # 只要比 Heun 晚一点就落不到窗口内，会被误报成"没有事件"。
            end_time_s=QUESTION_SPECS[args.question].end_time_s,
            event=event,
        )
    figure_path = output_dir / f"question{args.question}_final_profiles.png"
    plot_final_profiles(figure_path, model, result)
    report = validate_result(model, result)

    print(f"问题 {args.question} 计算完成，最终时刻：{result.time_s[-1]:.3f} s")
    print(
        f"界面策略：{model.interface_strategy.name}；"
        f"长度模式：{args.length_scaling}；"
        f"边界延拓：{args.plateau_mode}；"
        f"端面：{'计入' if model.include_end_faces else '忽略'}"
    )
    if result.step_s is not None and result.step_s.size:
        print(
            "实际内部时间步范围："
            f"{np.min(result.step_s):.6g}～{np.max(result.step_s):.6g} s"
        )
    if bdf_result is not None:
        print(
            "BDF 对照："
            f"最终时刻={bdf_result.time_s[-1]:.3f} s；"
            f"事件时刻={bdf_result.event_time_s}"
        )
        if result.event_time_s is not None and bdf_result.event_time_s is not None:
            relative_difference = abs(
                result.event_time_s - bdf_result.event_time_s
            ) / result.event_time_s
            print(f"BDF 与 Heun 事件时刻相对差异：{relative_difference:.6%}")
    if result.event_time_s is not None:
        print(
            "达到含水率阈值的时刻："
            f"{result.event_time_s:.3f} s（{result.event_time_s / 3600:.6f} h）"
        )
    print(f"基础验证报告：{report}")
    result_workbook = output_dir / f"result{args.question}.xlsx"
    if result_workbook.exists():
        if args.question == 4:
            print(verify_result4_workbook(result_workbook))
        else:
            print(verify_result_workbook(result_workbook, args.question))
    print("已生成文件：")
    for path in [*files, figure_path]:
        print(f"- {path}")


if __name__ == "__main__":
    main()
