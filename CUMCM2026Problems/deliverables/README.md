# 官方结果文件

本目录是四个问题的 `resultN.xlsx` 汇总，与 `附件3` 的官方模板格式一致。
源文件由求解器直接生成，这里是便于交付的副本。

| 文件 | 题目 | 工作表 | 数据行 | 生成配置 |
|---|---|---|---|---|
| `result1.xlsx` | 问题1 | 温度、水分浓度 | 1800（1–1800 s，1 s 步长） | 附录2 物性；固定半径；**不计端面**；`midpoint`；N=40 |
| `result2.xlsx` | 问题2 | 温度、水分浓度 | 10800（1–10800 s，1 s 步长） | 附录3 物性；固定半径；**不计端面**；**`kirchhoff`**；N=40 |
| `result3.xlsx` | 问题3 | Sheet1（水分浓度） | 3449（60 s 步长 + 末行连续事件时刻） | 附录3 物性；固定半径；**不计端面**；**`kirchhoff`**；N=40 |
| `result4.xlsx` | 问题4 | Sheet1（水分浓度） | 3066（60 s 步长 + 末行连续事件时刻） | 附录4 物性；**收缩**（附件2）；**不计端面**、报告**中截面口径**；**`kirchhoff`**；一维 N=40 |
| `result4_对照_拟一维端面模型.xlsx` | 问题4 对照 | Sheet1 | 2113 | 一维拟一维端面模型（`2h/L` 等效源），**仅作对照** |

## 关键结论

- **问题3 烘干时长 = 57.48298 h**（Heun）/ 57.48304 h（BDF），相对差 1.1×10⁻⁶。
  最后一个达标点是圆心，符合径向干燥的物理。
- 问题2 在 3 h 末的中心含水率 1.7662、表面 1.0080；温度已接近环境（49.85–49.97 ℃）。
- 问题4：达标 **51.0967 h**（一维仅侧面 + `kirchhoff`，N=40，中截面口径），与队友最新版
  **逐格最大差 0.0001**（仅 52/64400 个格点在末位差 1）；端面影响由二维轴对称模型复验。
- 界面扩散系数统一采用**基尔霍夫通量势** `D_face = ΔΦ/ΔC`，推导见
  [`../docs/01_论文推导稿_通用篇.md`](../docs/01_论文推导稿_通用篇.md) §6，证据与数字索引见
  [`../docs/03_证据与验证总表.md`](../docs/03_证据与验证总表.md)。

## 统一口径

- **端面**：主模型一律只算圆柱侧面，报告**中截面口径**的径向分布；问题4 的端面影响由
  二维轴对称模型（真实 Robin 边界）复验，与一维"仅侧面"解一致到 10⁻⁷ 量级。
- **边界延拓**：附件1 只给到 4 h，之后取恒温干燥段的 50 ℃、0.05 kg/kg。
- **网格**：径向 N=40；问题3 的 N=40→80 相对变化 0.0137%。
- **精度**：四位小数；超出药材边界的位置留空。

## 复现

在 `../AProblemMaterial/` 下执行：

```powershell
python -m aproblem --question 1 --intervals 40 --output-dir "outputs\问题1\主模型_midpoint_N40"
python -m aproblem --question 2 --intervals 40 --output-dir "outputs\问题2\主模型_kirchhoff_N40"
python -m aproblem --question 3 --intervals 40 --output-dir "outputs\问题3\主模型_kirchhoff_N40" --bdf-check
python -m aproblem --question 4 --intervals 40 --interface-strategy kirchhoff --output-dir "outputs\问题4\生产_kirchhoff_N40"

# 问题4 端面复验（二维轴对称，kirchhoff）
python -m aproblem.study axisymmetric-scan --radial-intervals 40 --axial-intervals 20 40 80 --interface-strategy kirchhoff
```

每个文件生成时都会自动回读校验（工作表名、表头、时间间隔、非空、四位小数格式），
校验不通过会直接报错。

## 待办

- 问题4 的二维 `kirchhoff` 端面复验（三个网格）正在运行，跑完回填
  [`../docs/03_证据与验证总表.md`](../docs/03_证据与验证总表.md) 的 Q4-4/Q4-8。
- 论文表 1～表 6 已可自动抽取：`python -m aproblem.paper_tables --question N`
  （用法见 [`../docs/06_代码与复现指南.md`](../docs/06_代码与复现指南.md) §3）。
