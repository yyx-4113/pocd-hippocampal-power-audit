# Step 5–6 报告 — GSE267933 术后海马单细胞定位与深入分析

> 数据：GSE267933（小鼠海马，Control 3 vs Surgery 3，手术诱导 POCD 模型）
> 原始 20,684 细胞 × 27,998 基因 → QC 后 **18,328 细胞**（保留率 88.6%）
> 流程脚本：`step5_sc_analysis.py`（Step 5/6）、`step5b_sc_pseudobulk_geneset.py`（5b）、`step5b2_sensitivity.py`（5b-2）
> 产物目录：`sc_out/`

---

## 1. 技术路线与质量控制

| 步骤 | 方法与参数 | 结果 |
|---|---|---|
| QC | min_genes=200、min_cells=3、pct_mt<20% | 20,684 → 18,328 细胞（中位线粒体比 7.9%，质量良好） |
| 归一化 | normalize_total(1e4) + log1p | — |
| 高变基因 | 2,000 HVG（batch_key=sample） | 19269 非 HVG 基因保留于 `raw` 供定位用 |
| 降维 | PCA 50 PC（arpack） | — |
| **去批次** | **Harmony（by sample），14 次迭代收敛** | `X_pca_harmony` (18328×50) |
| 聚类 | Leiden res=0.6，igraph flavor | **21 簇** |
| 注释 | 12 类 marker 基因集打分取最大（`sc.tl.score_genes`） | **9 种细胞类型** |

> 环境要点：`scanpy.external.pp.harmony_integrate` 与 `harmonypy 0.2.0`（PyTorch 后端）不兼容，本流程**直接调用 `harmonypy.run_harmony` 并显式 `.cpu().detach().numpy()`**，已在日志中确认 `[Harmony] rep=X_pca_harmony`（非回退）。

## 2. 细胞类型构成与术后比例变化（核心结果）

| 细胞类型 | Control % | Surgery % | Fisher p | BH p | 判读 |
|---|---|---|---|---|---|
| **Microglia 小胶质** | 39.53 | **47.56** | 6.1e-28 | **2.7e-27** | **显著扩张 ↑** |
| Oligodendrocyte 少突 | 29.82 | 29.71 | 0.87 | 0.91 | 无变化 |
| Astrocyte 星形胶质 | 9.72 | 7.86 | 9.8e-6 | **2.2e-5** | 下降 ↓ |
| Ependymal 室管膜 | 10.72 | **3.36** | 1.1e-88 | **1.0e-87** | **大幅丧失 ↓↓** |
| Endothelial 内皮 | 2.17 | 3.62 | 4.5e-9 | **1.3e-8** | 上升 ↑ |
| OPC | 4.03 | 3.98 | 0.91 | 0.91 | 无变化 |
| Neuron 神经元 | 1.11 | 1.61 | 4.1e-3 | **7.5e-3** | 名义上升（n 极小，见注） |
| VLMC / Pericyte | 1.70 / 1.20 | 1.35 / 0.94 | 0.053 / 0.098 | n.s. | — |

**读法**：术后海马最稳健的单细胞层变化是**细胞组成重塑**——小胶质占比扩张 8 个百分点、室管膜几乎崩解、星形胶质轻度减少。

## 3. 六特征基因的单细胞定位

| 基因 | 最高表达细胞类型（平均 lognorm） | 术后 pseudobulk log2FC（该类型） |
|---|---|---|
| **Fos** | Microglia **3.03** | Microglia +0.22；Ependymal +0.65 |
| **Junb** | Microglia **3.69** | Microglia +0.18；Ependymal +0.56 |
| **Sgk1** | Oligodendrocyte **2.44** | **Oligodendrocyte −0.65**（唯一明显下调） |
| **Gltp** | Oligodendrocyte **1.57** | Oligodendrocyte −0.06 |
| **Runx2** | Pericyte 0.20（其余≈0） | ≈0 |
| **Sla2** | **全类型≈0（未检出）** | ≈0 |

- **即刻早期基因（Fos/Junb）信号集中在小胶质**（远超神经元与少突）—— 与 bulk 层"神经活动-IEG"条目落在小胶质倾向一致，但需注意 IEG 对取材/解离应激极敏感。
- **Sgk1 与 Gltp 定位于少突胶质谱系**，与 Set D 髓鞘轴呼应。
- **⚠️ Sla2 在单细胞层完全 dropout（检出率≈0）**：Step 3 的 "稳健核心 3 基因"中，**Sla2 无法在单细胞层定位**，只保留 Junb、Sgk1 可定位。这是必须如实写入限制的。

## 4. 四条机制线的单细胞层验证（Step 5b）

方法：把细胞级降到样本级——逐 (细胞类型 × 样本) 聚合为 pseudobulk，集合内基因跨样本 z 后取均值作为集合得分，3v3 上 Welch t + BH（**细胞不是独立重复单元，故不做细胞级检验**）。

**结果：四条机制线在全部 9 种细胞类型上，BH 后无一显著**（最小 BH_p = 0.457）。

| 机制线 | Microglia ΔZ | Oligodendrocyte ΔZ | Ependymal ΔZ | Astrocyte ΔZ | 方向与 bulk 一致? |
|---|---|---|---|---|---|
| A 补体/突触修剪 | +0.22 | +0.31 | +0.37 | +0.08 | ✗（bulk 为整体下调） |
| B DAM | −0.03 | +0.55 | +0.84 | +0.36 | ✗ |
| C 线粒体 | +0.28 | −0.17 | −0.62 | −0.28 | 部分（少突/星形为负） |
| D 髓鞘/少突 | +0.51 | **−0.20** | +0.70 | +0.10 | **弱一致**（少突为负） |

定向 panel 补充：
- **Oligodendrocyte Myelin_struct ΔZ = −0.58（p=0.038，BH 0.46）** —— 全表中**唯一方向与 bulk 髓鞘轴下调一致**的条目；
- Microglia：DAM_core ΔZ=+0.07、Homeostatic −0.36、Complement_core −0.15 → **未观察到经典稳态→DAM 极化**；
- 名义显著项集中在 **Ependymal**（Myelin +1.36 p=0.0055；Complement +1.35 p=0.045；DAM +0.84 p=0.028）——但室管膜细胞数已崩解（10.7%→3.4%），残余细胞构成高度偏倚，**不作机制结论**。

## 5. 构建方式敏感性验证（Step 5b-2）

担心 log-norm 细胞均值聚合会压掉信号，故用更标准的 **counts 求和 → CPM → log1p** 重做：

| 指标 | 结果 |
|---|---|
| 两构建法 delta_Z 相关性 | **Pearson r = 0.857（p=4e-10）；Spearman ρ = 0.835**，n=32 |
| 名义 p<0.05 计数 | lognorm 5 项 vs counts-sum 3 项 |
| BH<0.05 计数 | **两法均为 0** |
| 方向反转条目（\|Δ\|>0.4） | **无** |

→ **结论稳健**：单细胞层机制线位移极弱并非聚合方式造成的假阴性。两种构建下均唯一稳定为负的是**少突胶质髓鞘集合**（−0.20 / −0.34）。

## 6. 配体-受体轴（探索性，均值乘积打分，非 CellChat）

| 配体→受体 | Sender → Receiver | score |
|---|---|---|
| Apoe → Trem2 | VLMC / Astrocyte / Ependymal → **Microglia** | 17.5 / 15.3 / 9.5 |
| Apoe → Trem2 | Astrocyte → VLMC | 9.4 |
| Tyrobp → Trem2 | VLMC → Microglia | 9.0 |
| **C1qa → Lrp1** | **Microglia → OPC** | 8.4 |
| C1qa → Lrp1 | VLMC → OPC | 8.1 |

- **Apoe→Trem2 是主导轴**，指向小胶质为接收端，与 bulk 层 DAM/吞噬受体轴（Tyrobp↓、Trem2）形成呼应；
- **C1qa→Lrp1（小胶质→OPC）**提供了"补体-吞噬信号作用于少突谱系前体"的候选桥接，与 Set D 髓鞘轴可连成一条机制链。
- 局限：该方法为表达均值乘积，未做通路显著性检验，仅作方向提示；建议稿件中明确标注为 exploratory，或后续以 CellChat/LIANA 复算。

## 7. 少突胶质谱系再聚类与拟时序

- 9 个 OL 亚簇（Control/Surgery 构成见表 `sc_ol_subcluster_composition.csv`），未见明显亚群比例的组间分离；
- **DPT 拟时序根簇 = cluster 7**（OPC_score 最高），`Myelin_score` 沿拟时序上升，轨迹结构成立（图 `fig_sc_ol_dpt.png`）；
- 但术后沿拟时序的成熟度位移微弱，与第 4 节少突髓鞘集合弱下调一致。

---

## 8. 必须写进稿件的四条限制

1. **单细胞层无 BH 显著位移**：所有机制线/panel 在 9 种细胞类型上 BH_p ≥ 0.457；**经两种 pseudobulk 构建法验证**。因此单细胞层的贡献是**定位与方向提示**，不是确证。
2. **细胞比例/状态解耦**：小胶质占比显著扩张（47.6%）而其转录组状态几乎不变（DAM_core ΔZ≈+0.07）→ 提示**募集/增殖**而非状态极化。这一"解耦"本身是可写的故事，但需在 Discuss 中与 bulk 层 DAM 集合下调 (-0.91Z) 谨慎对接。
3. **Sla2 单细胞层 dropout**，Step 3 特征基因签名在单细胞层只能定位 Junb、Sgk1（+Fos/Gltp）。
4. **Neuron 比例差异不可过度解读**：Control 仅 99 个神经元 vs Surgery 152 个，且小鼠海马 scRNA 中神经元占比极低（解离偏倚），p=7.5e-3 不构成生物学结论。

## 9. 与 bulk 层的对接（同轴闭环）

| 层 | 结论 | 强度 |
|---|---|---|
| bulk（三数据集 meta） | 四条机制线**协调整体下调**（D −1.48Z / B −0.91 / C −0.50 / A −0.29），个体基因弱 | 强（Stouffer p 1e-4 ~ 1e-37） |
| sc 组成 | 小胶质↑、室管膜↓↓、星形胶质↓ | 强（BH 1e-5 ~ 1e-87） |
| sc 表达状态 | 极弱；唯一方向一致者为**少突髓鞘下调** | 弱（BH n.s.，两法一致） |

**叙事落点**：bulk 层给出"方向 + 协调性"，单细胞层给出"**发生在哪些细胞（小胶质扩张、少突髓鞘、室管膜崩解）**"，两层共同支持"神经胶质-突触轴"框架，而精确基因签名仍需 Step 7–8 外部/临床样本验证。

---

## 附：产物清单（`geo_meta/sc_out/`）

**图**：`fig_sc_umap.png`（Leiden/细胞类型/分组三联）、`fig_sc_marker_dotplot.png`、`fig_sc_composition.png`、`fig_sc_feature_dotplot.png`（6 特征基因 × 9 类型气泡）、`fig_sc_feature_umap.png`、`fig_sc_pseudobulk_feature.png`、`fig_sc_ol_subclusters.png`、`fig_sc_ol_dpt.png`、`fig_sc_ol_featuregenes.png`、`fig_sc_microglia_DAM.png`、`fig_sc5b_geneset_heatmap.png`、`fig_sc5b_persample_scores.png`、`fig_sc5b2_sensitivity.png`

**表**：`sc_celltype_composition.csv`、`sc_cluster_markers.csv`、`sc_cluster_typescore.csv`、`sc_featuregene_meancount_by_celltype.csv`、`sc_featuregene_pct_by_celltype.csv`、`sc_pseudobulk_featuregenes.csv`、`sc_pseudobulkDEG_<celltype>.csv`（9 个）、`sc_ol_subcluster_composition.csv`、`sc_microglia_scores.csv`、`sc_lr_axis.csv`、`step5b_geneset_stats.csv`、`step5b_pseudobulk_<celltype>.csv`、`step5b2_sensitivity.csv`

**数据**：`GSE267933_processed.h5ad`（18328 细胞，含 raw 与 Harmony 嵌入）
