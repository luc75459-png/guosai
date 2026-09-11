# CUMCM 2026 A 题 · 中药材烘干 · 项目总入口

本目录是整个项目的根。所有路径都相对于本文件。

---

## 1. 目录结构

```text
CUMCM2026Problems/
├── README.md              ← 本文件
├── docs/                  ← 全部工作文档（索引见 docs/00_文档索引.md）
│   └── 历史/               ← 被取代的早期文档，不作为实现依据
├── AProblemMaterial/      ← 【正式工程】问题1~4 的求解与结果导出
├── AProblem/              ← 队友原始工程（只读参考，不改动）
├── CUMCM2026Problems/     ← 题目与附件（A~E 题、官方模板）
└── 11.pdf                 ← 一份外部建模推导参考
```

> **注意**：`CUMCM2026Problems/CUMCM2026Problems/` 这个重复的目录名是压缩包
> 解压留下的，它里面才是真正的题目与附件。该路径同时被本工程和队友的
> `AProblem` 写死，改名会同时破坏两边，因此保留不动。

---

## 2. 当前主模型配置

四个问题共用同一套有限体积求解内核，只在物性公式、是否收缩、终止条件上不同。

| 问题 | 物性 | 半径 | 端面 | 界面扩散系数取法 | 网格 | 终止条件 |
|---|---|---|---|---|---|---|
| 1 | 附录2 | 固定 | 不计 | `midpoint` | N=40 | 1800 s 定时 |
| 2 | 附录3 | 固定 | 不计 | **`kirchhoff`** | N=40 | 10800 s 定时 |
| 3 | 附录3 | 固定 | 不计 | **`kirchhoff`** | N=40 | 全场 `max C = 0.15` |
| 4 | 附录4 | 收缩（附件2） | 二维模型为真实边界 | `midpoint` | 径向40×轴向40 | 全场 `max C = 0.15` |

统一口径：

- **端面**：主模型一律只算圆柱侧面。问题4 的二维轴对称校验已证明 `2h/L`
  等效源把端面影响高估约 45%（中截面按对称性没有轴向通量）。
- **界面取法**：`D_face = ΔΦ/ΔC`（基尔霍夫通量势），推导与证据见
  [docs/03](docs/03_界面扩散系数取法与物理依据.md)。传热方程不适用通量势，
  该策略在传热侧退回中点点值。
- 问题1 的界面取法对结果影响低于输出精度（温度严格不变，水分差 <5×10⁻⁵），
  因此沿用中点值；问题4 已定稿，暂未改动。

---

## 3. 交付物

官方结果文件由求解器直接生成，位置如下：

| 文件 | 说明 | 数据行 |
|---|---|---|
| `AProblemMaterial/outputs/问题1/主模型_midpoint_N40/result1.xlsx` | 温度 + 水分浓度两个工作表 | 1800 |
| `AProblemMaterial/outputs/问题2/主模型_kirchhoff_N40/result2.xlsx` | 温度 + 水分浓度两个工作表 | 10800 |
| `AProblemMaterial/outputs/问题3/主模型_kirchhoff_N40/result3.xlsx` | 水分浓度单表 | 3449 |
| `AProblemMaterial/outputs/问题4/result4.xlsx` | 问题4 正式结果（二维轴对称模型） | 3072 |
| `AProblemMaterial/outputs/问题4/result4_对照_拟一维端面模型.xlsx` | 问题4 一维对照 | 2113 |

全部文件在生成时都会做**回读校验**（工作表名、表头、时间列间隔、非空、
四位小数格式），校验不通过会直接报错。

---

## 4. 复现命令

在 `AProblemMaterial/` 目录下执行，Python 用 Anaconda 环境：

```powershell
# 环境与测试
python -m pytest -q

# 问题1 / 2 / 3 正式计算（默认配置，直接产出 resultN.xlsx）
python -m aproblem --question 1 --intervals 40 --output-dir "outputs\问题1\主模型_midpoint_N40"
python -m aproblem --question 2 --intervals 40 --output-dir "outputs\问题2\主模型_kirchhoff_N40"
python -m aproblem --question 3 --intervals 40 --output-dir "outputs\问题3\主模型_kirchhoff_N40" --bdf-check

# 问题4 正式结果与验证报告
python -m aproblem.report

# 界面策略 × 网格收敛对照
python -m aproblem.study interface-scan --question 3 --intervals 20 40 80 --integrator bdf
```

命令行开关：

- `--interface-strategy {midpoint,harmonic,series,kirchhoff}`：覆盖按题号默认
- `--include-end-faces`：计入 `2h/L` 端面等效源（主模型默认不计入）
- `--bdf-check`：用 SciPy BDF 独立对拍事件时刻
- `--no-result`：只出预览文件，不写 `resultN.xlsx`

---

## 5. 当前状态

| 问题 | 状态 |
|---|---|
| 1 | 结果文件已生成并校验 |
| 2 | 结果文件已生成并校验；界面已切基尔霍夫 |
| 3 | 结果文件已生成并校验；达标 **57.48298 h**（Heun）/ 57.48304 h（BDF），相对差 1.1×10⁻⁶ |
| 4 | 结果文件已生成并校验；正式结果取二维轴对称模型，**尚未按新默认（不计端面）重跑一维对照** |

未完成事项：

1. 论文表 1～表 5 的自动抽取（目前需从 `resultN.xlsx` 手工取点）。
2. 问题4 的一维对照与文档尚未同步到"主模型不计端面"的新口径。
