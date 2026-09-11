"""CUMCM 2026 A 题药材烘干数值模拟包。

该包对外暴露路径配置、仿真配置、仿真结果和四问统一运行入口，
便于命令行脚本、测试代码和后续结果导出模块复用同一套核心逻辑。
"""

from .config import ProjectPaths, SimulationConfig
from .integrator import SimulationResult
from .scenarios import run_question

__all__ = ["ProjectPaths", "SimulationConfig", "SimulationResult", "run_question"]
