# Step 2 摘要 — bulk 整合 + 通路基因集接入

## 数据可用性修正
- **GSE178995 不可用**: 其 4 样本 `Sample_type=SRA`, GEO 仅有注释文件 `GSE178995_S_C.anno.txt.gz`, 无处理表达矩阵(原始 fastq 在 SRA SRP325799)。方案原计划三集合并, 实际整合 **GSE276942 + GSE215410 + GSE174412**。
- 另修正: GSE215410 样本为 C1,C2(对照) vs P1,P2(POCD) = **2v2**; GSE174412 = CON1-3 vs PND1-3 = 3v3。

## 单数据集 DEG
- GSE276942: n_post vs n_ctrl = 13v4; BH<0.05 = 7, p<0.05 = 672
- GSE215410: n_post vs n_ctrl = 2v2; BH<0.05 = 0, p<0.05 = 972
- GSE174412: n_post vs n_ctrl = 3v3; BH<0.05 = 0, p<0.05 = 529

## 跨数据集元分析 (Stouffer 加权 Z)
- 可比基因(≥2 数据集): 20348
- meta_p<0.05: 717; meta_padj(BH)<0.05: 0

## ComBat 合并
- 方法: 内置参数化 ComBat (位置 EB, 保留 group)
- 共同基因: 15998; 样本: 27
- 合并矩阵 DEG: BH<0.05=0, p<0.05=843

## 四套通路基因集 (命中本步表达矩阵)
- A_complement_synapse_pruning: 163
- B_DAM_microglia: 118
- C_mitochondrial_OXPHOS: 191
- D_myelin_oligodendrocyte: 75
- 并集: 523

## 候选基因 (通路 ∩ 元分析)
- 共 50 个
- T1 (n=32): Opalin, Sla2, Myo1d, Mog, Gltp, Tmem63a, Pmp22, Bcas1, Prr5l, Enpp6, Mag, Cebpb, Plp1, Cd83, Gpr37, Csf1, Runx2, Sox10, Gjc2, Prr18, Ctss, Ccl6, Siglech, Ndrg1, Lcn2, Ccl12, C5ar1, Cd80, Tyrobp, Cx3cl1, Megf10, Fabp5
- T2 (n=18): Junb, Egr2, Arc, Tlr2, Irf8, Myrf, Nkx2-2, Cnp, Fos, Aif1, Lpar1, Egr1, Irf5, Kctd17, Gpr34, Cx3cr1, Sgk1, Cd9
- T3 (n=0): 

## 基因集水平统计 (Step 2b)
| 集合 | 成员 | 均值 meta_Z | 下调比例 | 单样本 t p | Wilcoxon p | Stouffer p | 单体 meta_p<0.05 |
|---|---|---|---|---|---|---|---|
| A_complement_synapse_pruning | 149 | -0.285 | 0.711 | 9.55e-04 | 4.75e-06 | 5.04e-04 | 11 |
| B_DAM_microglia | 112 | -0.906 | 0.839 | 3.23e-16 | 3.31e-13 | 9.40e-22 | 19 |
| C_mitochondrial_OXPHOS | 187 | -0.504 | 0.797 | 1.79e-20 | 3.34e-17 | 5.26e-12 | 0 |
| D_myelin_oligodendrocyte | 73 | -1.484 | 0.904 | 5.85e-15 | 3.61e-12 | 7.87e-37 | 23 |

**关键结论**: 线粒体 OXPHOS(C)与补体(A)在**个体基因层面**多为弱改变, 但**集合层面**方向高度一致 -> 支持'需基因集/ML 才能捕获协调性弱信号'的论点。

**Step 3 候选池**: 60 基因 (基因水平 + 集合代表), 见 step2_candidate_pool.csv。
