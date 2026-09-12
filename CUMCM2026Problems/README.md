# CUMCM 2026 A 题 · 中药材烘干 · 项目总入口

> **一句话**：正式版本只有一个——`AProblemMaterial/`（`cumcm-a-drying-material` v0.2.0）。
> 队友的 `AProblem/`、早期独立实现 `代码骨架/`、以及各种备选口径**都只是参考或对照**，
> 不进交付文件。完整配置表见 [`docs/06 §2.1`](docs/06_代码与复现指南.md)，全部文档索引见
> [`docs/00_文档索引.md`](docs/00_文档索引.md)。

---

## 1. 采用的版本（一眼速览）

### 1.1 四问配置与交付文件

| 问题 | 物性 | 半径 | 端面 | 界面扩散系数 | 网格 | 终止条件 | 交付文件 |
|---|---|---|---|---|---|---|---|
| 1 | 附录2 | 固定 `R=0.02 m` | 不计 | `midpoint` | N=40 | `t=1800 s` 定时 | `result1.xlsx` |
| 2 | 附录3 | 固定 | 不计 | **`kirchhoff`** | N=40 | `t=10800 s` 定时 | `result2.xlsx` |
| 3 | 附录3 | 固定 | 不计 | **`kirchhoff`** | N=40 | 全场 `max C=0.15` | `result3.xlsx` |
| 4 | 附录4 | **收缩**（附件2） | 不计（二维模型复验端面口径） | **`kirchhoff`** | N=40 | 全场 `max C=0.15` | `result4.xlsx` |

### 1.2 统一口径三条

1. **端面**：一维主模型一律只算圆柱侧面（`include_end_faces=False`），报告口径为**中截面**的径向分布。
   问题4 的端面影响由**二维轴对称模型**（真实 Robin 边界）复验，与一维"仅侧面"解一致到 10⁻⁷ 量级；
   `2h/L` 等效源已证明会把端面作用高估约 45%，**只作对照、不作结果**。
2. **界面扩散系数**：按题号取用。问题2/3 用基尔霍夫通量势 `D_face=ΔΦ/ΔC`（由通量守恒唯一确定、
   无自由参数）；问题1 因 `D` 全场仅变化 1.27 倍而沿用中点值；问题4 自 2026-09-12 起
   **全队统一到 `kirchhoff`**（与问题2/3 同口径、与队友一致）。传热侧的界面导热系数
   统一取中点状态点值。
3. **时间与输出**：二阶显式 Heun + 按节点局部总流出系数的动态稳定步长（安全系数 0.8）；
   终止事件用跨步线性插值定位；结果文件四位小数、由求解器直接生成并自动回读校验。

### 1.3 版本谱系（存在哪些版本、用不用）

| 版本 | 位置 | 状态 | 说明 |
|---|---|---|---|
| **`AProblemMaterial/` v0.2.0** | 本目录 | ✅ **采用** | 四问正式求解与交付：一维有限体积内核 + 问题4 二维轴对称模型 |
| `AProblem/` | 本目录 | ❌ 只读参考 | 队友原始工程：`0.7247/0.8705` 校准权重、`2h/L` 端面源、固定 `dt`；不修改、不引用其结果 |
| `代码骨架/` | 工作区根目录（未入库） | ⚪ 仅作交叉验证 | 早期独立实现，用于"两套代码互证"；不作交付 |
| 二维轴对称材料坐标模型 | `AProblemMaterial/aproblem/axisymmetric.py` | ✅ 问题4 正式结果 | 端面为真实边界；同时用于端面口径验证 |
| 拟一维 `2h/L` 端面模型 | 同工程，`--include-end-faces` | ❌ 仅对照 | 高估端面作用约 45%，结果文件后缀标注"对照" |
| 备选口径：问题1 N=160 | 本地 `AProblemMaterial/outputs/问题1/主模型_midpoint_N160/`（**未入库**） | ⚪ 已备查、本轮不采用 | 表面含水率末位更稳（见 `docs/03 §4.2`） |
| 旧口径：问题4 `midpoint` + 轴向平均 | 本地 `AProblemMaterial/outputs/问题4/生产_midpoint_N40*`（**未入库**）；交付件的旧版备份在 `outputs/_归档/` | ❌ 已弃用（2026-09-12 统一口径） | 51.1837 h，作为对照留档；与采用口径差 0.17% |

### 1.4 权威文档在哪

| 需要什么 | 去哪 |
|---|---|
| 公式与推导（唯一来源） | [`docs/01 通用篇`](docs/01_论文推导稿_通用篇.md)、[`docs/02 四问篇`](docs/02_论文推导稿_四问篇.md) |
| 哪个数字该写多少、有什么证据 | [`docs/03 证据与验证总表`](docs/03_证据与验证总表.md)、[`docs/05 问题4 数据档案`](docs/05_问题4_数据与技术档案.md) |
| 文献条目与等级 | [`docs/04 文献与外部材料`](docs/04_文献与外部材料.md) |
| 跑代码、复现任何数字 | [`docs/06 代码与复现指南`](docs/06_代码与复现指南.md) |
| **正式写论文** | [`docs/07 论文写作手册`](docs/07_论文写作手册.md)（大纲、图表清单、逐章素材模板、分工与时间表） |
| **新队友 / 新 AI 接手** | [`docs/08 工作进度与交接`](docs/08_工作进度与交接.md)（现状、决策、待办、仓库信息） |
| **验证实验（L0/L1）** | [`docs/09 验证实验报告`](docs/09_验证实验报告.md)（完整实验报告；判据原件在 `docs/06 §9`） |

---

## 2. 目录结构

```text
CUMCM2026Problems/
├── README.md              ← 本文件（版本速览）
├── docs/                  ← 全部工作文档（索引见 docs/00_文档索引.md）
│   ├── 00～06             ← 7 份活跃文档（推导稿 / 证据 / 文献 / 数据档案 / 复现指南）
│   ├── 文献/              ← 检索审计原件（只读）
│   ├── 源材料/            ← 题面原文提取与队友侧原件（只读）
│   └── 历史/              ← 早期文档与本次合并归档（已合并/，不作为实现依据）
├── AProblemMaterial/      ← 【正式工程】问题1~4 的求解与结果导出
├── AProblem/              ← 队友原始工程（只读参考，不改动）
├── deliverables/          ← 官方 result1~4.xlsx 交付副本
├── CUMCM2026Problems/     ← 题目与附件（A~E 题、官方模板）
└── 11.pdf                 ← 一份外部建模推导参考
```

> **注意**：`CUMCM2026Problems/CUMCM2026Problems/` 这个重复的目录名是压缩包解压留下的，
> 里面才是真正的题目与附件。该路径同时被本工程和队友的 `AProblem` 写死，改名会同时破坏两边，
> 因此保留不动。

---

## 3. 交付物

| 文件 | 说明 | 数据行 |
|---|---|---|
| `deliverables/result1.xlsx`（源自 `AProblemMaterial/outputs/问题1/主模型_midpoint_N40/`） | 温度 + 水分浓度两个工作表 | 1800 |
| `deliverables/result2.xlsx`（源自 `outputs/问题2/主模型_kirchhoff_N40/`） | 温度 + 水分浓度两个工作表 | 10800 |
| `deliverables/result3.xlsx`（源自 `outputs/问题3/主模型_kirchhoff_N40/`） | 水分浓度单表 | 3449 |
| `deliverables/result4.xlsx`（源自 `AProblemMaterial/outputs/问题4/生产_kirchhoff_N40/`，一维仅侧面 + `kirchhoff` + 中截面口径） | 水分浓度单表 | 3066 |
| `AProblemMaterial/outputs/问题4/result4_对照_拟一维端面模型.xlsx` | 问题4 对照（**非交付**） | 2113 |

全部文件在生成时都会做**回读校验**（工作表名、表头、时间列间隔、非空、四位小数格式），
校验不通过会直接报错。论文表 1～表 6 可用
`python -m aproblem.paper_tables --question N` 从这些文件自动抽取，不必手工取点。

---

## 4. 复现命令

在 `AProblemMaterial/` 目录下执行，Python 用 Anaconda 环境（`C:\Users\luzho\anaconda3\python.exe`）：

```powershell
# 环境与测试
python -m pytest -q

# 问题1 / 2 / 3 正式计算
python -m aproblem --question 1 --intervals 40 --output-dir "outputs\问题1\主模型_midpoint_N40"
python -m aproblem --question 2 --intervals 40 --output-dir "outputs\问题2\主模型_kirchhoff_N40"
python -m aproblem --question 3 --intervals 40 --output-dir "outputs\问题3\主模型_kirchhoff_N40" --bdf-check

# 问题4 正式结果（二维轴对称）与验证报告
python -m aproblem.report

# 批量研究：interface-scan / grid-convergence / time-convergence / sensitivity / axisymmetric-scan
python -m aproblem.study interface-scan --intervals 20 40 80 --integrator bdf

# 论文表格
python -m aproblem.paper_tables --question 3
```

命令行开关：`--interface-strategy {midpoint,harmonic,series,kirchhoff}`（覆盖按题号默认）、
`--include-end-faces`（计入 `2h/L` 端面等效源，仅作对照）、`--bdf-check`（BDF 独立对拍）、
`--no-result`（只出预览文件）。详细说明见 [`docs/06`](docs/06_代码与复现指南.md)。

---

## 5. 当前状态

| 问题 | 状态 |
|---|---|
| 1 | 结果已生成并校验（N=40）。表面含水率末位不确定度约 `±2×10⁻⁴`，论文不确定度一节可直接引用 `docs/03 §4.4` |
| 2 | 结果已生成并校验；界面取基尔霍夫 |
| 3 | 结果已生成并校验；达标时刻（N=40）见 `docs/03 §2`；端面影响定量正在进行（二维对照算例） |
| 4 | 结果已生成并校验；**2026-09-12 起统一为 `kirchhoff` + 中截面口径**（一维仅侧面 N=40，51.0967 h），与队友最新版逐点一致；二维 `kirchhoff` 三网格端面复验进行中 |

**验证实验（2026-09-12 新增）**：L0 解析基准判据 **20/20 通过**（圆柱 Robin 观测阶
1.99–2.18、Heun 时间阶 2.000–2.002、基尔霍夫恒等式残差 1.5×10⁻¹⁶）；L1 外部实测对照判据
**10/12 通过**（零拟合参数 `RMSE_X=0.104/0.090`、`R²=0.988/0.991`，反演 `D` 为论文值的
0.935/0.972 倍）；未通过的两条是温度平台，已定量归因于"忽略蒸发潜热"（+14/+17 ℃）。
完整报告见 [`docs/09 验证实验报告`](docs/09_验证实验报告.md)。

**已决定（2026-09-12）**：① **问题4 统一口径**——全队改用 `kirchhoff` + 中截面口径，
`result4.xlsx` 由一维仅侧面 N=40 重出（与队友最新版逐点一致），二维 `kirchhoff` 三网格做端面复验；
旧 `midpoint` + 轴向平均版本留档对照。② **问题1 维持 N=40**，不重出交付文件，
末位不确定度按 `docs/03 §4.2` 的表述写进论文。详见 `docs/03` 缺口 T1/T2/T8/T9/T10 与 `docs/05 §8`。

**未完成**：① 问题3 的端面影响定量（T3，二维对照算例已暂停，需要时按 `docs/06 §3` 重跑）；
② 问题4 二维 `kirchhoff` 端面复验（T8，同样暂停）；③ 外部文献实测数据对照（首选 Teleken 2025，方案见 `docs/06 §6`）。

---

## 6. 仓库与协作

| 项 | 值 |
|---|---|
| **本仓库（上传目标）** | `https://github.com/luc75459-png/guosai.git`（remote 名 `mine`） |
| 队友仓库 | `https://github.com/Shyyyyyy-0920/26-.git`（remote 名 `origin`，只读参考） |
| 本地克隆 | `E:\大三\2-26国赛\26-`（`main` 分支跟踪 `mine/main`） |
| 推送须知 | **推送会直接改队友仓库的 `main`，推送前必须先确认**；或用分支 + PR |

文档结构与阅读思路见 [`docs/00_文档索引.md`](docs/00_文档索引.md)；
新接手的人或 AI 请先读 [`docs/08 工作进度与交接`](docs/08_工作进度与交接.md)。
