# CUMCM 2026 A题：药材烘干数值模型

本项目采用一维圆柱径向有限体积法，将温度场和干基含水率场统一离散，并使用固定步长 Heun（显式二阶 Runge-Kutta）方法推进。根据最新物理模型，两个端面的对流换热和传质通过轴向积分平均折算为一维体积源/汇，默认参与四问计算，同时保留仅考虑圆柱侧面的消融开关。代码中的模块说明、关键类、函数和主要算法步骤均使用中文注释，便于队员阅读、复核和论文推导。

完整建模、算法选择和验收流程见：[02_建模与代码工作流.md](../docs/02_建模与代码工作流.md)。

## 目录

```text
AProblem/
├── aproblem/
│   ├── __main__.py       # 命令行入口
│   ├── config.py         # 路径、圆柱尺寸、端面开关和仿真配置
│   ├── data.py           # 附件读取与边界插值
│   ├── grid.py           # 圆柱径向控制体几何
│   ├── physics.py        # 三组物性公式
│   ├── model.py          # 热湿有限体积半离散方程与端面等效源项
│   ├── integrator.py     # Heun推进与终止事件
│   ├── scenarios.py      # 问题1-4的组装逻辑
│   ├── outputs.py        # NPZ/CSV预览输出
│   ├── plotting.py       # 剖面图与热力图
│   └── validation.py     # 数值和物理检查
├── tests/                # 不依赖竞赛附件的单元测试
└── outputs/              # 本地运行输出，不覆盖官方模板
```

## 环境要求与依赖安装

建议使用 **Python 3.11 或更高版本**，不要直接使用系统全局环境。项目所需包已经写在 `pyproject.toml` 中：

| 环境包 | 最低版本 | 用途 |
|---|---:|---|
| `numpy` | 1.26 | 数组运算、有限体积状态和结果保存 |
| `pandas` | 2.1 | 读取附件 Excel、生成 CSV 预览 |
| `scipy` | 1.11 | PCHIP 半径插值及后续 BDF 对照 |
| `matplotlib` | 3.8 | 生成中文温度与含水率图像 |
| `openpyxl` | 3.1 | 供 pandas 读取 `.xlsx`，并用于正式 Excel 导出 |
| `pytest` | 8.0 | 运行项目测试，属于开发依赖 |

### 方案一：使用 Conda 环境

队员已经有 Anaconda/Miniconda 时，推荐为比赛单独创建环境：

```powershell
conda create -n cumcm-a python=3.11 -y  # 创建名为 cumcm-a 的独立 Python 3.11 环境
conda activate cumcm-a                  # 激活该环境，后续命令都在该环境中执行
cd "E:\数学建模\CUMCM2026Problems\AProblem"  # 进入项目根目录，路径不同的队员应替换成自己的路径
python -m pip install --upgrade pip setuptools  # 更新安装工具，减少可编辑安装失败的概率
python -m pip install -e ".[dev]"       # 安装项目及 numpy、pandas、scipy、matplotlib、openpyxl、pytest
python -m pytest                        # 使用当前环境的 Python 运行全部测试
```

如果想继续使用已有环境，例如 `si100`，无需重新创建：

```powershell
conda activate si100                    # 激活已有的 si100 环境
cd "E:\数学建模\CUMCM2026Problems\AProblem"  # 进入项目根目录
python -m pip install -e ".[dev]"       # 补齐该环境缺失的全部项目依赖，包括 openpyxl
python -m pytest                        # 确认安装后的环境能够通过测试
```

### 方案二：使用 Python 自带虚拟环境

没有 Conda、但已安装 Python 3.11 以上版本时使用：

```powershell
cd "E:\数学建模\CUMCM2026Problems\AProblem"  # 进入项目根目录
python -m venv .venv                    # 在项目内创建名为 .venv 的独立虚拟环境
.\.venv\Scripts\Activate.ps1           # 激活虚拟环境
python -m pip install --upgrade pip setuptools  # 更新 pip 和项目安装工具
python -m pip install -e ".[dev]"       # 一次安装程序运行和测试所需的全部环境包
python -m pytest                        # 运行 tests/ 中的全部测试
```

安装完成后可以一次检查所有关键包：

```powershell
python -c "import numpy, pandas, scipy, matplotlib, openpyxl, pytest; print('全部依赖安装成功')"  # 任一包缺失时会直接显示包名
```

## 快速开始

完成上面的任意一种环境安装方案后，在 `AProblem` 目录执行。每个问题都必须分别计算“考虑端面”和“不考虑端面”两种工况；以下以问题1为例：

```powershell
python -m pytest  # 再次确认当前代码和环境通过全部单元测试
python -m aproblem --question 1 --output-dir "outputs\问题1\考虑端面"  # 计算问题1的复杂工况：计入两个端面，结果单独保存
python -m aproblem --question 1 --ignore-end-faces --output-dir "outputs\问题1\不考虑端面"  # 计算问题1的简单工况：只计圆柱侧面，结果单独保存
```

其中，不添加 `--ignore-end-faces` 时，程序默认**考虑两个端面**；添加 `--ignore-end-faces` 时，程序**不考虑两个端面**。两次计算必须指定不同的 `--output-dir`，否则同名结果文件可能被后一次运行覆盖。

如果 PowerShell 因执行策略阻止激活脚本，可以先在当前终端临时执行 `Set-ExecutionPolicy -Scope Process Bypass`，关闭该终端后设置会自动失效。

### 常见环境报错

如果出现 `No module named 'openpyxl'`，说明当前正在运行程序的 Python 环境缺少 Excel 引擎。在项目目录执行：

```powershell
python -m pip install -e ".[dev]"       # 推荐：按项目配置补齐所有依赖，而不只安装单个缺失包
python -c "import openpyxl; print(openpyxl.__version__)"  # 验证 openpyxl 已安装到当前环境
```

如果旧版本项目安装时出现 `Multiple top-level packages discovered`，说明 `setuptools` 把 `outputs/` 误判成了 Python 包。当前 `pyproject.toml` 已明确只打包 `aproblem`；拉取或复制最新文件后重新执行：

```powershell
python -m pip install -e ".[dev]"       # 使用修正后的打包配置安装项目和全部依赖
```

比赛现场如果只需要立刻补装 Excel 读取包，也可以先绕过项目安装：

```powershell
python -m pip install "openpyxl>=3.1"   # 仅安装缺失的 openpyxl；版本表达式要放在引号内
```

如果安装后仍提示缺包，通常是 `pip` 和 `python` 指向了不同环境。使用下面命令检查：

```powershell
where.exe python                        # 查看终端实际找到的所有 python.exe 路径
python -c "import sys; print(sys.executable)"  # 显示当前运行程序使用的 Python 路径
python -m pip --version                 # 显示 pip 所属的 Python 环境，路径应与上一行一致
```

始终使用 `python -m pip ...` 和 `python -m pytest`，可以最大限度避免 Conda、系统 Python 和用户目录中的包互相混用。

默认附件目录为：

```text
../CUMCM2026Problems/A题/附件
```

也可以显式指定：

```powershell
python -m aproblem --question 1 --attachments "完整的附件目录" --output-dir "outputs\问题1\考虑端面"  # 从指定目录读取附件，并计算考虑端面的工况
python -m aproblem --question 1 --attachments "完整的附件目录" --ignore-end-faces --output-dir "outputs\问题1\不考虑端面"  # 从同一附件目录读取数据，并计算不考虑端面的工况
```

### 四个问题的两种工况

建议始终使用下面统一的目录命名。每个问题运行两条命令，共运行八次：

```powershell
# 问题1：常物性、固定半径、计算3小时
python -m aproblem --question 1 --output-dir "outputs\问题1\考虑端面"  # 问题1复杂工况：计入两个端面
python -m aproblem --question 1 --ignore-end-faces --output-dir "outputs\问题1\不考虑端面"  # 问题1简单工况：忽略两个端面

# 问题2：变物性、固定半径、计算3小时
python -m aproblem --question 2 --output-dir "outputs\问题2\考虑端面"  # 问题2复杂工况：计入两个端面
python -m aproblem --question 2 --ignore-end-faces --output-dir "outputs\问题2\不考虑端面"  # 问题2简单工况：忽略两个端面

# 问题3：变物性、达到全场含水率阈值时停止
python -m aproblem --question 3 --output-dir "outputs\问题3\考虑端面"  # 问题3复杂工况：计入两个端面
python -m aproblem --question 3 --ignore-end-faces --output-dir "outputs\问题3\不考虑端面"  # 问题3简单工况：忽略两个端面

# 问题4：变物性、半径随附件2收缩、达到阈值时停止
python -m aproblem --question 4 --output-dir "outputs\问题4\考虑端面"  # 问题4复杂工况：计入两个端面
python -m aproblem --question 4 --ignore-end-faces --output-dir "outputs\问题4\不考虑端面"  # 问题4简单工况：忽略两个端面
```

运行完成后的目录结构如下：

```text
outputs/
├── 问题1/
│   ├── 考虑端面/
│   └── 不考虑端面/
├── 问题2/
│   ├── 考虑端面/
│   └── 不考虑端面/
├── 问题3/
│   ├── 考虑端面/
│   └── 不考虑端面/
└── 问题4/
    ├── 考虑端面/
    └── 不考虑端面/
```

需要试验其他内部时间步时，两种工况应使用相同的 `--dt`。例如，问题1使用 0.25 秒内部步长时：

```powershell
python -m aproblem --question 1 --dt 0.25 --output-dir "outputs\问题1\时间步0.25秒\考虑端面"  # 使用0.25秒步长重算考虑端面的工况
python -m aproblem --question 1 --dt 0.25 --ignore-end-faces --output-dir "outputs\问题1\时间步0.25秒\不考虑端面"  # 使用0.25秒步长重算不考虑端面的工况
```

## 两个端面的处理

设圆柱长度为 `L=0.25 m`。二维轴向方程沿长度积分后，用同一半径处的轴向平均状态近似端面状态，得到

```text
温度体积源：S_T = (2 h_T / L) · (T_air - T)
水分体积汇：S_C = (2 h_m / L) · (C_air - C)
```

两个源项作用于每个径向控制体；原有最外节点 Robin 边界仍只表示圆柱侧面，因此不会重复计算。该方法是轴向平均的拟一维闭合，不等价于完整二维轴对称模型。正式论文应同时给出默认模型和 `--ignore-end-faces` 的对比结果。

## 当前约定

- 所有内部计算使用 SI 单位：m、s、K/Celsius 差值、kg/kg。
- 圆柱长度默认为 0.25 m；四问默认计入两个端面的等效换热和失水贡献。
- 固定半径问题使用 20 个径向区间，即 21 个节点，恰好对应 0.0-2.0 cm、间隔 0.1 cm。
- 默认内部步长 0.5 s。问题1、2分别每1 s保存，问题3、4每60 s保存。
- 问题2/3的界面扩散系数采用调和平均与算术平均的校准混合，算术平均权重默认为 0.7247；问题1和问题4仍使用调和平均。
- 问题3、4在附件1结束后将烘房条件保持为 50°C 和 0.05 kg/kg。
- `outputs` 目前输出便于检查的 NPZ 和 CSV；官方 result*.xlsx 的模板填充单独实现，避免早期调试覆盖模板。
- 问题4的预览 CSV 使用归一化半径 `xi=r/R(t)` 作为列；NPZ 同时保存每个输出时刻的实际半径，正式导出时再插值到题目要求的物理距离。

## 图表约定

- 所有图题、坐标轴和图例均使用中文。
- 温度曲线使用橙红色圆点线，干基含水率曲线使用蓝色方点线，避免黑白之外难以区分。
- 横坐标统一为“距药材中心的距离（cm）”；纵坐标分别为“温度（℃）”和“干基含水率（kg/kg）”。
- 绘图模块会依次尝试微软雅黑、黑体、思源黑体等中文字体；若本机仍出现方框，请安装其中任意一种字体后重新运行。

## 建议开发顺序

1. 完成并核验问题1。
2. 对问题1运行含端面/仅侧面的消融对比，核验等效源项。
3. 在同一固定半径内核上切换问题2/3物性。
4. 加入问题3的含水率阈值事件。
5. 使用归一化径向坐标处理问题4的收缩半径。
6. 最后接入官方 Excel 模板和论文图表。
