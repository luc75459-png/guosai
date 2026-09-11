# CUMCM 2026 A 题：问题 4 材料坐标数值求解版

本目录是队友原始工程 `AProblem` 的独立扩展版本，**不修改**原目录。
这里只做问题 4：收缩药材的烘干时长与 `result4.xlsx`。

所有工作文档集中在 `../docs/`，索引见
[`../docs/00_文档索引.md`](../docs/00_文档索引.md)：

```text
../docs/03_界面扩散系数取法与物理依据.md   ← 界面取法的推导、文献与数值证据
../docs/04_问题4_完整技术文档.md           ← 问题4 模型、代码逐模块解说与正式结果
../docs/05_问题4_对照分析.md               ← 对照算例、误差拆解与敏感性分析
../docs/06_问题4_材料坐标工作流.md          ← 问题4 工作流与收敛试验记录
```

## 主方案（当前定稿）

| 项目 | 取值 | 依据 |
|---|---|---|
| 坐标 | 材料坐标 ξ = r/R(t) | 网格跟随材料收缩，表面永远是 ξ=1 |
| 界面扩散系数 | `midpoint` | 网格收敛扫描，N=40→80 只差 0.046% |
| **端面** | **二维轴对称模型中的真实 Robin 边界** | 唯一不使用 `2h/L` 闭合假设的做法 |
| 长度 | L = 0.25 m 恒定 | 与附录4密度公式的干物质守恒检验相容 |
| 网格 | 径向 40 × 轴向 40 | 40×40 与 40×80 相对差 1.6×10⁻⁸ |
| 时间推进 | 二阶 Heun + 动态稳定步长 | 与 BDF 对拍相对差 5×10⁻⁷ |
| 4 h 后边界 | 50 ℃、0.05 kg/kg 保持 | 敏感性分析显示只影响 0.36% |

**正式结果：达标（全场 max(C) = 0.15 kg/kg）时刻 = 51.1837 h**
（184261.270 s，二维轴对称模型）。

一维拟一维端面模型（端面用 `2h/L` 体积源）给出 35.2133 h。它数值上同样
收敛，但其端面闭合假设已被二维模型证伪，**只作为对照保留**，理由见
[`05_问题4_对照分析.md`](../docs/05_问题4_对照分析.md) 第 5 节。

## 关键结论：界面扩散系数不能用调和平均

题目给的 `D(C,T)` 在干燥过程中跨越约四个数量级，界面上取什么值是本题
最敏感的数值选择。队友原工程用一个 AI 估算出来的混合权重（`0.8705`）
在校准结果，**本目录完全没有这个参数**，改用三套有物理依据的策略并
用网格收敛性定夺：

| 界面策略 | N=20 | N=40 | N=80 | N=40→80 相对变化 |
|---|---|---|---|---|
| **midpoint（中点状态点值）** | 35.2503 h | 35.2132 h | 35.1971 h | **0.046%** |
| harmonic（调和平均） | 39.8927 h | 36.1002 h | 35.3700 h | 2.02% |
| series（圆柱径向串联阻力） | 39.7687 h | 36.0886 h | 35.3691 h | 1.99% |

调和平均与圆柱串联阻力式几乎完全重合，说明两者同源：它们都把连续变化
的 `D` 当成两层不同材料串联，等于人为插进一层不存在的分层阻力，因此
粗网格上严重高估传质阻力、系统偏保守。中点状态点值不引入任何分层
假设，粗网格上就已收敛。

## 目录结构

```text
AProblemMaterial/
├── aproblem/
│   ├── __main__.py        # 命令行入口与 result4 导出
│   ├── config.py          # 路径、物性、界面策略、长度和边界延拓配置
│   ├── data.py            # 附件读取、环境插值和 PCHIP 半径
│   ├── grid.py            # 圆柱径向控制体几何
│   ├── interfaces.py      # 界面物性策略：midpoint / harmonic / series
│   ├── physics.py         # 四问经验物性公式（含逐字段接口）
│   ├── model.py           # 一维热湿有限体积 RHS 与端面等效源
│   ├── axisymmetric.py    # 二维轴对称材料坐标模型（端面校验用）
│   ├── integrator.py      # Heun、动态步长和事件定位
│   ├── crosscheck.py      # SciPy BDF 独立时间积分
│   ├── outputs.py         # NPZ/CSV 预览与 result4.xlsx
│   ├── validation.py      # 物理范围检查与 result4 回读校验
│   ├── consistency.py     # 干物质守恒自洽检验（长度假设）
│   ├── study.py           # 收敛性与敏感性批量研究驱动
│   ├── report.py          # 汇总图表与验证报告
│   └── plotting.py        # 中文图表
├── tests/
└── outputs/
```

## 运行方式

以下命令都在本目录执行，Python 使用 Anaconda 环境。

### 正式计算

```powershell
# 第一步：一维对照（端面用 2h/L 等效源）
python -m aproblem --question 4 --intervals 40 --interface-strategy midpoint `
    --output-dir "outputs\问题4\生产_midpoint_N40" --bdf-check

# 第二步：消融对照（忽略两个端面）
python -m aproblem --question 4 --intervals 40 --interface-strategy midpoint `
    --ignore-end-faces --output-dir "outputs\问题4\生产_midpoint_N40_不含端面"

# 第三步：二维轴对称（正式结果来源），三个轴向网格建议并行执行
python -m aproblem.study axisymmetric-scan --radial-intervals 40 `
    --axial-intervals 40 --interface-strategy midpoint `
    --output-dir "outputs\study\axisym_40"

# 第四步：汇总，生成正式 result4.xlsx、对照文件、图表与验证报告
python -m aproblem.report
```

### 命令行参数

| 参数 | 说明 |
|---|---|
| `--question` | 问题编号，问题 4 用 `4` |
| `--intervals` | 径向区间数，生产用 40 |
| `--interface-strategy` | `midpoint`（默认）/ `harmonic` / `series` |
| `--length-scaling` | `constant`（默认，L=0.25 m）/ `isotropic`（L∝R） |
| `--plateau-mode` | `fixed`（默认，50 ℃/0.05）/ `last_value`（保持附件末值） |
| `--ignore-end-faces` | 关闭端面等效源，只算圆柱侧面 |
| `--fixed-dt` / `--dt` | 关闭动态步长 / 指定最大内部步长 |
| `--bdf-check` | 用 BDF 对同一条半离散方程做独立对拍 |
| `--no-result4` | 只出预览文件，不写 result4.xlsx |

### 收敛性与敏感性研究

```powershell
# 界面策略 × 网格收敛表（BDF 快，适合扫描）
python -m aproblem.study interface-scan --intervals 20 40 80 --integrator bdf

# 时间步收敛（Heun，主求解器）
python -m aproblem.study time-convergence --intervals 40 --dt-list 0.5 0.25 0.125

# 长度假设与边界延拓敏感性
python -m aproblem.study sensitivity --intervals 40 --integrator bdf

# 二维轴对称端面校验
python -m aproblem.study axisymmetric-scan --radial-intervals 40 `
    --axial-intervals 20 40 80 --interface-strategy midpoint
```

### 汇总交付物

```powershell
python -m aproblem.report
```

生成 `outputs/问题4/result4.xlsx`、`outputs/问题4/figures/*.png`
和 `outputs/问题4/问题4_验证报告.md`。

### 测试

```powershell
python -m pytest -q
```

## 问题 4 的主模型

### 材料坐标

```text
ξ = r/R(t),   r(ξ,t) = ξR(t)
```

`ξ` 标记固定材料层，圆心 `ξ=0` 与表面 `ξ=1` 永远不动。干基含水率
跟随固体骨架，因此从材料导数出发变换后**不出现额外的拖曳项**；
只有从固定空间坐标方程出发才需要保留迁移项，两者不能混用。

均匀收缩下的几何缩放：

```text
节点间距 dr ∝ R    界面面积 A ∝ R
控制体体积 V ∝ R²  径向扩散项 ∝ 1/R²
```

### 端面

正式模型用**二维轴对称材料坐标模型**（`aproblem/axisymmetric.py`）把两个
端面当作真实边界求解：

```text
ξ = r/R(t)，η = z/L(t)
∂T/∂t = (1/(ξR²))∂_ξ(ξ·α·∂_ξT) + (1/L²)∂_η(α·∂_ηT)
∂C/∂t = (1/(ξR²))∂_ξ(ξ·D·∂_ξC) + (1/L²)∂_η(D·∂_ηC)
边界：ξ=0 对称；ξ=1 侧面 Robin；η=0,1 两个端面 Robin
```

一维代码里保留了另一种做法——把端面对流通量沿轴向积分平均、折算成体积
源项：

```text
S_T = (2 h_T / L)(T_air − T)
S_C = (2 h_m / L)(C_air − C)
```

该式的闭合假设是"端面状态约等于同半径处的轴向平均状态"。二维校验表明
**这个假设在本题几何下不成立**：端面到轴向中截面有 `L/2 = 12.5 cm`，
而轴向扩散在整个烘干过程中的穿透深度只有 `0.8～2.1 cm`，因此端面只能
烘干两端各约 1～2 cm 的区域，够不到决定达标时刻的中截面。把端面失水
摊到每一个径向控制体上会把达标时刻从 51.18 h 压到 35.21 h。

所以一维 `2h/L` 结果只作为对照，正式结果取二维模型。
完整对照见 [`05_问题4_对照分析.md`](../docs/05_问题4_对照分析.md) 第 5 节。

### 长度

附件 2 只给半径，不给长度。主方案取 `L = 0.25 m` 恒定，并用附录 4 的
密度公式做干物质守恒自洽检验：把 `rho(C)` 换算成干物质密度
`rho_d = rho(C)/(1+C)`，干物质守恒要求每层的体积收缩倍数唯一确定，
反推得到的长度比与"长度不变"相容，而各向同性同步收缩要求药材多压缩
一倍以上，与密度公式冲突。同步收缩仍作为敏感性上界报告
（达标时刻 30.31 h，比主线短约 14%）。

## 输出文件

正式结果与对照结果都由 `python -m aproblem.report` 生成，格式完全一致：

```text
outputs/问题4/
├── result4.xlsx                    正式：二维轴对称模型
├── result4_对照_拟一维端面模型.xlsx   对照：一维 2h/L 等效源
└── 问题4_验证报告.md
```

两份文件均按官方模板生成：

```text
列：0.0, 0.1, ..., 1.9 cm + 药材表面（共 21 个数据列）
行：60 s 的整数倍，末行是连续事件时刻
精度：四位小数
超出当前 R(t) 的固定距离位置留空，表面列始终有值
```

`aproblem.validation` 里的 `verify_result4_workbook()` 会自动检查以上
全部条目，`python -m aproblem.report` 每次运行时都会执行。

## 模型边界

- 只把"均匀、各向同性的自由径向收缩"作为问题 4 主模型。非均匀收缩、
  孔隙率演化和硬壳力学不在主模型中，因为附件没有提供内部位移或力学参数。
- 题目给出的 `D` 视为干基含水率方程的有效扩散系数；潜热与收缩功
  因缺少闭合参数不进入主模型。
- 二维轴对称模型只用于量化端面降维误差，不参与 `result4.xlsx` 生成。
- 不允许为了匹配某个预期时间而校准物理参数。
