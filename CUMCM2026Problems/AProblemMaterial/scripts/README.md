# scripts 目录

存放**一次性/辅助**脚本，不参与正式求解，也不进入 `aproblem` 包。

| 脚本 | 用途 | 相关文档 |
|---|---|---|
| `digitize_teleken.py` | 从 Teleken et al. (2025) Figure 4 数字化实测干燥曲线：两遍独立算法 + 质量门 + QA 叠加图 | `docs/06 §9.4`、`docs/09 §4`、`data/external/teleken2025/README.md` |

运行方式（在 `AProblemMaterial/` 下）：

```powershell
python scripts/digitize_teleken.py
```

输出：

- `data/external/teleken2025/figure4_60C.csv`、`figure4_70C.csv`
- `data/external/teleken2025/digitization_overlay_60C.png`、`_70C.png`

脚本依赖 `Pillow`、`numpy`、`scipy`（连通域）与 `matplotlib`（仅用于画 QA 图）。
数字化结果**不要手工修改**：任何改动都会使 `docs/09` 里的不确定度与质量门结论失效。
