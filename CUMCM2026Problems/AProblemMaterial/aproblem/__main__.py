"""命令行入口：组装参数、运行指定问题并生成预览结果与图像。"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .config import ProjectPaths, SimulationConfig
from .crosscheck import integrate_bdf
from .outputs import write_preview_files, write_result4_workbook
from .plotting import plot_final_profiles
from .scenarios import run_question
from .validation import validate_result


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
        "--no-result4",
        action="store_true",
        help="问题 4 仅生成预览文件，不写正式 result4.xlsx",
    )
    parser.add_argument(
        "--bdf-check",
        action="store_true",
        help="使用 BDF 对同一模型做独立时间积分并比较事件时刻",
    )
    parser.add_argument(
        "--ignore-end-faces",
        action="store_true",
        help="关闭两个端面的轴向平均等效源项，仅计算圆柱侧面换热和传质",
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
        include_end_faces=not args.ignore_end_faces,
    )

    # 四个问题共用相同入口，只在场景组装阶段切换物性、时长和半径函数。
    model, result = run_question(args.question, paths, config)
    files = list(write_preview_files(output_dir, args.question, model, result))
    if args.question == 4 and not args.no_result4:
        files.append(
            write_result4_workbook(
                paths.result_template(4),
                output_dir / "result4.xlsx",
                model,
                result,
            )
        )
    bdf_result = None
    if args.bdf_check:
        event = None
        if args.question == 4:
            event = lambda _time, state: float(
                np.max(state[model.node_count :]) - config.moisture_threshold
            )
        bdf_result = integrate_bdf(
            model.rhs,
            np.asarray(result.state[0], dtype=float),
            end_time_s=result.time_s[-1],
            event=event,
        )
    figure_path = output_dir / f"question{args.question}_final_profiles.png"
    plot_final_profiles(figure_path, model, result)
    report = validate_result(model, result)

    print(f"问题 {args.question} 计算完成，最终时刻：{result.time_s[-1]:.3f} s")
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
    print("已生成文件：")
    for path in [*files, figure_path]:
        print(f"- {path}")


if __name__ == "__main__":
    main()
