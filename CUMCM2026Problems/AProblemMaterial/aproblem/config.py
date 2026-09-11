#路径、半径、步长和阈值配置
"""项目路径与仿真参数配置。

路径配置负责定位竞赛附件和输出目录；仿真配置集中保存四问共享的
几何、时间步、初值、边界换热/传质系数及终止阈值。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .interfaces import INTERFACE_STRATEGY_NAMES


@dataclass(frozen=True)
class ProjectPaths:
    """项目中输入、输出目录的统一路径集合。"""

    project_root: Path
    attachments: Path
    outputs: Path

    @classmethod
    def discover(
        cls,
        project_root: Path | None = None,
        attachments: Path | None = None,
    ) -> "ProjectPaths":
        """根据项目根目录和可选附件目录生成绝对路径。"""

        root = (project_root or Path(__file__).resolve().parents[1]).resolve()
        default_attachments = root.parent / "CUMCM2026Problems" / "A题" / "附件"
        return cls(
            project_root=root,
            attachments=(attachments or default_attachments).resolve(),
            outputs=root / "outputs",
        )

    @property
    def environment_file(self) -> Path:
        """返回烘房温度和含水率附件的路径。"""

        return self.attachments / "附件1.xlsx"

    @property
    def radius_file(self) -> Path:
        """返回药材半径随时间变化附件的路径。"""

        return self.attachments / "附件2.xlsx"

    def result_template(self, question: int) -> Path:
        """返回指定问题的官方 Excel 结果模板。"""

        if question not in (1, 2, 3, 4):
            raise ValueError("问题编号必须为 1、2、3、4 之一")
        return self.attachments / "附件3" / f"result{question}.xlsx"


# 各题默认的界面扩散系数取法。
#
# 问题2/3 的干燥时长由 D(C,T) 跨越约四个数量级主导，界面取值直接决定答案，
# 因此默认用通量势（基尔霍夫）取法 D_face = ΔΦ/ΔC：它是唯一能让两点有限体积
# 通量格式与真实通量严格相等的取值，不含自由参数，网格收敛也最快。
#
# 问题1 只有 30 min、水分几乎没有移动（渗透深度约 0.4 cm），界面取法对结果
# 影响可忽略，保持既有的中点点值；问题4 的正式结果已按中点值定稿，暂不改动。
DEFAULT_INTERFACE_STRATEGY_BY_QUESTION: dict[int, str] = {
    1: "midpoint",
    2: "kirchhoff",
    3: "kirchhoff",
    4: "midpoint",
}


def default_interface_strategy(question: int) -> str:
    """返回指定问题的默认界面策略名称。"""

    try:
        return DEFAULT_INTERFACE_STRATEGY_BY_QUESTION[question]
    except KeyError as error:
        raise ValueError(f"不支持的问题编号：{question}") from error


@dataclass(frozen=True)
class SimulationConfig:
    """四问共享的数值模拟配置，内部单位统一采用 SI 制。

    ``dt_s`` 始终表示允许的最大内部时间步。启用 ``adaptive_dt`` 后，
    积分器还会根据当前几何和物性计算稳定步长上限，并取两者中的较小值。
    """

    radius_m: float = 0.02
    cylinder_length_m: float = 0.25
    radial_intervals: int = 20
    dt_s: float = 0.5
    adaptive_dt: bool = True
    dt_safety_factor: float = 0.8
    initial_temperature_c: float = 28.0
    initial_moisture: float = 2.55
    heat_transfer_coefficient: float = 25.0
    mass_transfer_coefficient: float = 8.0e-7
    moisture_threshold: float = 0.15
    plateau_temperature_c: float = 50.0
    plateau_moisture: float = 0.05
    # 主模型一律只算圆柱侧面，不计入两个端面的 2h/L 等效源。
    # 问题4 的二维轴对称校验已经证明该闭合把端面影响高估约 45%
    # （中截面按对称性没有轴向通量，不该获得端面失水的好处）；
    # 问题1~3 用同一口径，保证四问的几何假设一致。
    include_end_faces: bool = False
    # 留空表示"按题号取默认"，见 DEFAULT_INTERFACE_STRATEGY_BY_QUESTION；
    # 显式赋值（含命令行 --interface-strategy）时以显式值为准。
    interface_strategy: str | None = None
    length_scaling: str = "constant"
    plateau_mode: str = "fixed"

    def validate(self) -> None:
        """在运行前检查会导致求解失败的基础参数。"""

        if self.radius_m <= 0:
            raise ValueError("药材半径 radius_m 必须为正数")
        if self.cylinder_length_m <= 0:
            raise ValueError("药材长度 cylinder_length_m 必须为正数")
        if self.radial_intervals < 2:
            raise ValueError("径向区间数 radial_intervals 至少为 2")
        if self.dt_s <= 0:
            raise ValueError("内部时间步长 dt_s 必须为正数")
        if not 0.0 < self.dt_safety_factor <= 1.0:
            raise ValueError("时间步安全系数 dt_safety_factor 必须位于 (0, 1] 区间")
        if self.initial_moisture <= 0:
            raise ValueError("初始干基含水率 initial_moisture 必须为正数")
        if (
            self.interface_strategy is not None
            and self.interface_strategy not in INTERFACE_STRATEGY_NAMES
        ):
            raise ValueError(
                "界面策略 interface_strategy 必须为 "
                + "、".join(INTERFACE_STRATEGY_NAMES)
            )
        if self.length_scaling not in ("constant", "isotropic"):
            raise ValueError("长度模式 length_scaling 必须为 constant 或 isotropic")
        if self.plateau_mode not in ("fixed", "last_value"):
            raise ValueError("边界延拓模式 plateau_mode 必须为 fixed 或 last_value")
