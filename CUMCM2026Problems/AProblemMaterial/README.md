# CUMCM 2026 A 题：材料坐标数值求解版

本目录是 `AProblem` 的独立扩展版本，不修改队友原有代码。核心求解器、
物性模型、附件读取、绘图和测试结构尽量保持原样，主要补全问题 4 所需的
材料坐标收缩、动态稳定步长、BDF 对拍和正式 `result4.xlsx` 导出。

完整推导、模型边界和验收方案见：

```text
../问题4_材料坐标工作流.md
```

## 问题 4 主模型

令：

```text
ξ = r/R(t)
r(ξ,t) = ξR(t)
```

其中 `ξ` 是固定的材料层坐标，`R(t)` 来自附件 2 的 PCHIP 保形插值。
圆心始终为 0，表面始终为 `R(t)`。在该材料坐标下，干基含水率的材料
时间导数直接等于固定 `ξ` 的时间导数，因此不额外添加朗道伪对流项。

均匀收缩下：

```text
节点间距 dr ∝ R
界面面积 A ∝ R
控制体体积 V ∝ R²
径向扩散通量密度 ∝ 1/R²
```

控制体几何仍在每个时间步根据 `R(t)` 重建，侧边 Robin 边界和两个端面
的轴向平均等效源项沿用原模型。

## 与原版的区别

- 原版代码及文档保持不动，本目录可以独立运行和提交。
- Heun 主求解器支持动态稳定步长；`--dt` 表示最大步长而不是强行固定步长。
- 问题 4 自动基于官方模板生成 `result4.xlsx`。
- 增加材料网格一致性、实际步长范围和 BDF 独立对拍。
- `result4.xlsx` 只写固定物理距离 `0.0～1.9 cm` 和单独的药材表面列；
  超过当前半径的位置为空白。

## 目录结构

```text
AProblemMaterial/
├── aproblem/
│   ├── __main__.py       # 命令行入口与 result4 导出
│   ├── config.py         # 路径、物性和时间步配置
│   ├── crosscheck.py     # SciPy BDF 独立时间积分
│   ├── data.py           # 附件读取、环境插值和 PCHIP 半径
│   ├── grid.py           # 材料坐标对应的圆柱控制体几何
│   ├── integrator.py     # Heun、动态步长和事件定位
│   ├── model.py          # 热湿有限体积 RHS 和稳定步长估计
│   ├── outputs.py        # NPZ/CSV 预览与 result4.xlsx
│   ├── physics.py        # 四问经验物性公式
│   ├── scenarios.py      # 四问组装逻辑
│   └── validation.py     # 材料网格与数值范围检查
├── tests/
└── outputs/
```

## 运行方式

快速验证问题 4：

```powershell
python -m aproblem --question 4 --intervals 20 --output-dir "outputs\问题4\快速验证"
```

正式计算建议使用 `N=40`：

```powershell
python -m aproblem --question 4 --intervals 40 --output-dir "outputs\问题4\考虑端面"
python -m aproblem --question 4 --intervals 40 --ignore-end-faces --output-dir "outputs\问题4\不考虑端面"
```

加入 BDF 事件时刻对拍：

```powershell
python -m aproblem --question 4 --intervals 40 --bdf-check --output-dir "outputs\问题4\BDF对拍"
```

关闭动态步长并严格使用 `--dt`：

```powershell
python -m aproblem --question 4 --intervals 40 --fixed-dt --dt 0.25 --output-dir "outputs\问题4\固定步长"
```

运行全部测试：

```powershell
python -m pytest -q
```

## 已完成的验证

问题 4 使用 `N=20`、Heun 主求解器的结果：

```text
达到 max(C)=0.15 的时间：约 39.892696 h
```

同一 RHS 使用 SciPy BDF 的结果：

```text
达到 max(C)=0.15 的时间：约 39.892690 h
相对差异：约 1.6×10⁻⁵%
```

问题 4 使用 `N=40`、Heun 主求解器的候选结果为：

```text
达到 max(C)=0.15 的时间：约 36.100211 h
同网格 BDF：约 36.100220 h
```

`N=20` 与 `N=40` 的结束时间相差约 3.79 h，说明论文正式结果不能直接
使用快速版。当前 `N=40` 结果只是候选值，仍需完成 `N=80` 网格收敛和
时间步收敛后才能定稿。

## 模型边界

本版本只把“均匀、各向同性的自由径向收缩”作为问题 4 主模型。非均匀
收缩、孔隙率演化、硬壳力学和完整热-湿-力耦合不在主模型中，原因是附件
没有提供内部位移或力学参数。后续如需增加非均匀收缩，只能作为敏感性
分析，不应混入 `result4.xlsx` 的正式结果。
