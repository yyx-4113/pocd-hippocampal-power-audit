# Step 5-6 摘要 — GSE267933 单细胞定位与深入分析

- QC 后: 18328 细胞 × 2000 HVG (原始 20684 细胞)
- 整合: X_pca_harmony (Harmony by sample)
- Leiden 21 簇 -> 9 种细胞类型

## 细胞类型分布
celltype
Microglia          8003
Oligodendrocyte    5455
Astrocyte          1606
Ependymal          1270
OPC                 734
Endothelial         535
VLMC                278
Neuron              251
Pericyte            196

## 细胞比例 (Control vs Surgery)
       celltype  n_ctrl  n_surg  pct_ctrl  pct_surg     fisher_p         BH_p
      Microglia    3514    4489     39.53     47.56 6.104862e-28 2.747188e-27
Oligodendrocyte    2651    2804     29.82     29.71 8.716135e-01 9.099853e-01
      Astrocyte     864     742      9.72      7.86 9.848530e-06 2.215919e-05
            OPC     358     376      4.03      3.98 9.099853e-01 9.099853e-01
    Endothelial     193     342      2.17      3.62 4.498986e-09 1.349696e-08
      Ependymal     953     317     10.72      3.36 1.134615e-88 1.021153e-87
         Neuron      99     152      1.11      1.61 4.146163e-03 7.463093e-03
           VLMC     151     127      1.70      1.35 5.303227e-02 7.954841e-02
       Pericyte     107      89      1.20      0.94 9.821560e-02 1.262772e-01

## 特征基因定位 (pseudobulk log2FC, Surgery vs Control)
celltype  Astrocyte  Endothelial  Ependymal  Microglia  Neuron    OPC  Oligodendrocyte  Pericyte   VLMC
gene                                                                                                   
Fos           0.154        0.263      0.653      0.219   0.565  0.145            0.211     0.489 -0.117
Gltp          0.007        0.023      0.090      0.051   0.060  0.047           -0.064    -0.074  0.208
Junb          0.206        0.037      0.561      0.180   0.637  0.402            0.153     0.306 -0.156
Runx2         0.004        0.001      0.062     -0.000  -0.051 -0.001            0.000     0.101 -0.007
Sgk1          0.038       -0.101     -0.166      0.035  -0.042  0.052           -0.650     0.032 -0.031
Sla2          0.000        0.000     -0.001      0.000   0.000 -0.002           -0.001     0.006  0.000

注: 3v3 设计, 单细胞层仅作定位与方向验证; 统计推断依赖 bulk 层与 pseudobulk(逐样本聚合)。
