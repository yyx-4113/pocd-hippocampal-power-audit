# Step 5e — pseudobulk 检验选择敏感性分析（审稿意见 P2 回应）

> 目的：回应审稿人"pseudobulk 检验选择未论证"。主分析用 Welch t（log1p 均值 pseudobulk），
> 本步补 limma-voom（edgeR TMM → voom → lmFit → eBayes）作为标准替代路线。
> 环境：R 4.4.3 + limma 3.62.2 + edgeR（Bioconductor 二进制，本沙箱安装成功）。

## 1. 数据

- 输入：QC 后 18,328 细胞，9 个细胞类型，来自**原始 10x 整数 counts**（`GSE267933.h5ad`，20,684 × 27,998）
  与 QC 后 celltype 注释对齐（overlap 18,328/18,328）。
- 逐 (细胞类型 × 样本) counts 求和 → 19,269–17,849 非零基因 × 6 样本（3 Control / 3 Surgery）。
- 产物：`sc_out/pbcounts/pbcounts_<celltype>.csv`、`sample_meta.csv`。

## 2. 结果：limma-voom vs Welch（主分析）

低表达过滤（至少 3 样本 counts≥10）后共 72,204 个基因级检验（主分析未过滤为 173,421）。

| 细胞类型 | voom 检验数 | voom p<0.05 | 期望 | **voom adjP<0.05** | voom min adjP | Welch p<0.05 |
|---|---|---|---|---|---|---|
| Astrocyte | 9,214 | 592 | 460.7 | **2** | 0.0021 | 521 |
| Endothelial | 7,167 | 422 | 358.4 | 0 | 0.354 | 339 |
| Ependymal | 10,974 | 2,030 | 548.7 | 0 | 0.206 | 1,167 |
| Microglia | 11,309 | 1,018 | 565.5 | 0 | 0.071 | 1,056 |
| Neuron | 6,128 | 553 | 306.4 | **2** | 0.040 | 377 |
| OPC | 9,706 | 866 | 485.3 | 0 | 0.061 | 543 |
| Oligodendrocyte | 11,102 | 741 | 555.1 | **1** | 0.046 | 558 |
| Pericyte | 2,750 | 208 | 137.5 | **1** | 0.030 | 414 |
| VLMC | 3,854 | 306 | 192.7 | **5** | 0.0031 | 398 |
| **合计** | **72,204** | | | **11** | | **0 / 173,421** |

两种归一化（lib.size vs edgeR TMM）给出**相同的 11 个显著基因数**，稳健。

## 3. 11 个显著基因的身份（决定解释）

全部 `logFC < 0`（11/11 下调）。**Ttr 占 5 个**。

| 细胞类型 | 基因 | logFC | AveExpr (log2CPM) | P | adjP |
|---|---|---|---|---|---|
| Astrocyte | **Ttr** | −4.64 | 7.0 | 2.3e-07 | 0.0021 |
| Astrocyte | Utp14b | −1.13 | 7.0 | 9.8e-07 | 0.0045 |
| VLMC | **Ttr** | −4.59 | 7.0 | 1.6e-06 | 0.0031 |
| VLMC | Il1b | −1.71 | 11.5 | 1.2e-06 | 0.0031 |
| VLMC | Rgs1 | −3.86 | 7.2 | 9.1e-06 | 0.0090 |
| VLMC | Bcl2a1d | −1.93 | 8.1 | 9.3e-06 | 0.0090 |
| VLMC | Fosb | −1.18 | 9.7 | 5.8e-05 | 0.044 |
| Neuron | **Ttr** | −5.70 | 6.8 | 1.1e-05 | 0.040 |
| Neuron | Nrgn | −3.09 | 8.7 | 1.3e-05 | 0.040 |
| Oligodendrocyte | **Ttr** | −4.93 | 6.5 | 4.2e-06 | 0.046 |
| Pericyte | **Ttr** | −4.39 | 7.0 | 1.1e-05 | 0.030 |

四机制线命中：仅 **Fosb**（Set A）与 **Il1b**（Set B），各 1 个，均在最小的 VLMC 类型中。

## 4. 关键解释：Ttr 是 ambient RNA（soup）指纹，不是细胞自主信号

Ttr（transthyretin）由**室管膜/脉络丛**高表达（本数据 Ependymal 的 mean log2CPM = 16.28），
其余细胞类型仅 6.5–8.5（背景水平）。观察到的模式：

| 细胞类型 | Ttr 表达水平 (log2CPM) | Ttr logFC | Ttr P |
|---|---|---|---|
| Ependymal（真实来源） | **16.28** | **−3.06** | 1.3e-03 |
| Astrocyte | 6.97 | −4.64 | 2.3e-07 |
| VLMC | 7.02 | −4.59 | 1.6e-06 |
| Neuron | 6.84 | **−5.70** | 1.1e-05 |
| Oligodendrocyte | 6.52 | −4.93 | 4.2e-06 |
| Pericyte | 6.99 | −4.39 | 1.1e-05 |
| Microglia | 6.79 | −4.69 | 6.3e-06 (adjP 0.071) |

**在几乎不表达 Ttr 的细胞类型中，Ttr 的变化幅度（−4.4 ~ −5.7）与显著性都超过其真实来源细胞
（Ependymal，−3.06）**。这一"非表达者变化更大"的模式是游离 mRNA 污染（ambient/soup RNA）的
教科书指纹，无法由细胞自主转录调控解释。它与本数据中室管膜比例收缩
（10.7% → 3.4%，动物级 p=0.038）自洽：室管膜细胞减少 → 悬液中游离 Ttr 减少 →
所有液滴的 Ttr 背景下降 → 每一类细胞的 pseudobulk Ttr 都下降。

## 5. 结论（写入稿件）

1. 两种主流 pseudobulk 路线（Welch/lognorm 与 limma-voom/TMM）在**定性结论上一致**：
   3v3 下不存在跨细胞类型、跨机制线可复现的细胞状态位移。
2. 定量上有差异，且必须诚实报告：BH 后显著基因数 0/173,421（Welch）vs 11/72,204（voom，0.015%）。
3. 差异可完全归因于三点技术因素：
   (a) 低表达过滤把检验负担从 173,421 降到 72,204（BH 阈值放宽约 2.4 倍）；
   (b) voom 对高丰度基因赋更高权重，放大了 Ttr 这类 soup 主导转录本的信号（5/11）；
   (c) 11 个命中全部为负、无跨类型一致性、仅 2 个落在机制线上，且集中在最小的 VLMC。
4. raw p 分布也需诚实处理：Ependymal 在 voom 下 raw p<0.05 = 2,030（期望 549，3.7× 富集），
   与 Welch 下 1,167（期望 964）相比偏斜更大；但两条路线在 BH 层面无基因存活。
   → 因此"检验保守"的表述仅在 Welch 路线成立，voom 路线应表述为"p 值分布偏斜、BH 后无存活"。

## 6. 产物

- `sc_out/pbcounts/pbcounts_<celltype>.csv`、`sample_meta.csv`
- `sc_out/step5e_limma_voom_stats.csv`（逐类型汇总）
- `sc_out/step5e_voom_significant_genes.csv`（11 个基因）
- `sc_out/step5e_method_comparison.csv`（两法对比）
- `sc_out/pbcounts/voom_<celltype>.csv`（全基因 voom 结果）
- 脚本：`step5e_export_pseudobulk_counts.py`、`step5e_limma_voom.R`、`step5e_report_numbers.py`
