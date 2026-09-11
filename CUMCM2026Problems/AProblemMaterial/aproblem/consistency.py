# 干物质守恒自洽检验
"""用附录4的密度公式和附件2的半径数据反推长度演化是否自洽。

思路
----
附录4给出的是药材的**表观湿密度** ``rho(C)``（总质量除以当前体积）。
干基含水率 ``C`` 的定义是水/干物质，因此该处的干物质密度为

    rho_d(C) = rho(C) / (1 + C)

干物质在烘干过程中不增不减，所以同一材料层在初末两态必须满足

    rho_d(C_末) * V_末 = rho_d(C_初) * V_初

即每一层的体积收缩倍数是唯一确定的。均匀收缩下 ``V ∝ R^2 L``，
于是可以从附件2的半径数据和本次算出的含水率剖面反推出长度比
``L_末/L_初``，再与"长度不变"和"各向同性同步收缩"两个假设对照。

这是一个**数据一致性检验**，不是新模型：它用的是题目已经给出的
密度公式和半径数据，不含任何自由参数。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def apparent_density_question4(moisture: np.ndarray) -> np.ndarray:
    """附录4给出的表观湿密度 rho(C) = 760 + 90C。"""

    return 760.0 + 90.0 * np.asarray(moisture, dtype=float)


def dry_density_question4(moisture: np.ndarray) -> np.ndarray:
    """由表观湿密度换算的干物质密度 rho(C)/(1+C)。"""

    values = np.asarray(moisture, dtype=float)
    return apparent_density_question4(values) / (1.0 + values)


@dataclass(frozen=True)
class VolumeConsistencyReport:
    """干物质守恒自洽检验的结果。"""

    initial_dry_density_kg_m3: float
    final_dry_density_range_kg_m3: tuple[float, float]
    required_volume_ratio: float
    radial_area_ratio: float
    implied_length_ratio: float
    isotropic_volume_ratio: float
    mismatch_constant_length: float
    mismatch_isotropic_length: float

    def __str__(self) -> str:
        """输出中文摘要，便于写进论文和验证报告。"""

        return (
            f"初态干物质密度 {self.initial_dry_density_kg_m3:.2f} kg/m³，"
            f"末态 {self.final_dry_density_range_kg_m3[0]:.2f}～"
            f"{self.final_dry_density_range_kg_m3[1]:.2f} kg/m³；"
            f"干物质守恒要求的体积比 {self.required_volume_ratio:.4f}，"
            f"附件2给出的径向面积比 {self.radial_area_ratio:.4f}，"
            f"反推长度比 L末/L初 = {self.implied_length_ratio:.4f}；"
            f"长度不变的相对偏差 {self.mismatch_constant_length:+.2%}，"
            f"各向同性收缩的相对偏差 {self.mismatch_isotropic_length:+.2%}"
        )


def radial_volume_weights(radial_fraction: np.ndarray) -> np.ndarray:
    """返回圆柱节点型控制体的归一化体积权重（正比于环面积）。"""

    fraction = np.asarray(radial_fraction, dtype=float)
    inner = np.empty_like(fraction)
    outer = np.empty_like(fraction)
    inner[0] = 0.0
    inner[1:] = 0.5 * (fraction[:-1] + fraction[1:])
    outer[:-1] = inner[1:]
    outer[-1] = 1.0
    areas = np.pi * (outer**2 - inner**2)
    return areas / float(np.sum(areas))


def check_volume_consistency(
    final_moisture: np.ndarray,
    initial_moisture: float,
    radius_ratio: float,
    radial_fraction: np.ndarray | None = None,
) -> VolumeConsistencyReport:
    """执行干物质守恒自洽检验。

    ``final_moisture`` 是末态各材料层的干基含水率，``radius_ratio``
    是末态外半径与初始外半径之比。体积收缩倍数按干物质守恒逐层计算，
    再用体积权重取平均，最后反推长度比。
    """

    values = np.asarray(final_moisture, dtype=float)
    fraction = (
        np.linspace(0.0, 1.0, values.size)
        if radial_fraction is None
        else np.asarray(radial_fraction, dtype=float)
    )
    if fraction.shape != values.shape:
        raise ValueError("径向坐标数组必须与含水率剖面形状一致")
    if radius_ratio <= 0.0:
        raise ValueError("半径比必须为正数")

    weights = radial_volume_weights(fraction)
    initial_dry_density = float(dry_density_question4(np.array([initial_moisture]))[0])
    final_dry_density = dry_density_question4(values)

    # 每一层要求的体积收缩倍数；干物质守恒下它是唯一确定的。
    shrink_factor = final_dry_density / initial_dry_density
    required_volume_ratio = 1.0 / float(np.sum(weights * shrink_factor))
    radial_area_ratio = radius_ratio**2
    implied_length_ratio = required_volume_ratio / radial_area_ratio
    isotropic_volume_ratio = radial_area_ratio * radius_ratio

    return VolumeConsistencyReport(
        initial_dry_density_kg_m3=initial_dry_density,
        final_dry_density_range_kg_m3=(
            float(np.min(final_dry_density)),
            float(np.max(final_dry_density)),
        ),
        required_volume_ratio=required_volume_ratio,
        radial_area_ratio=radial_area_ratio,
        implied_length_ratio=implied_length_ratio,
        isotropic_volume_ratio=isotropic_volume_ratio,
        mismatch_constant_length=(required_volume_ratio - radial_area_ratio)
        / radial_area_ratio,
        mismatch_isotropic_length=(required_volume_ratio - isotropic_volume_ratio)
        / isotropic_volume_ratio,
    )
