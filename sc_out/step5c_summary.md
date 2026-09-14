# Step 5c — 单细胞层"无显著"的校准与功效分析

- n=3v3, alpha=0.05, power=0.80 的标准化最小可检测效应 d = **3.07**（非中心 t 分布, df=2n−2，与稿件 Methods 同法）
  > 2026-09-14 更正：此处原记为 3.27，与本文件所述方法不符（正态近似亦仅为 2.29）；已按非中心 t 重算为 3.07，与 `step5c_mde.csv`、`step5f_required_n.csv` 及稿件一致。
- 判定逻辑: 若全基因组 BH 显著基因≈0 且 p 值分布均匀 -> 是**功效极限**;
  若全基因组有大量显著基因而机制线全无 -> 是**集合确实没变**(强结论)。

## 全基因组校准
       celltype  n_genes  padj_lt_0.05  padj_lt_0.10  p_lt_0.05  p_lt_0.01  exp_p_lt_0.05   KS_D          KS_p    min_p top_gene  top_log2FC
      Astrocyte    19269             0             0        521         87          963.5 0.1661  0.000000e+00 0.000354    Nr4a3   -0.033380
    Endothelial    19269             0             0        339         47          963.5 0.1928  0.000000e+00 0.000092  Gm17354   -0.033322
      Ependymal    19269             0             1       1167        199          963.5 0.1669  0.000000e+00 0.000005     Pgk1   -0.341482
      Microglia    19269             0             0       1056        176          963.5 0.0919 3.833926e-142 0.000057    Bola1    0.031709
         Neuron    19269             0             0        377         56          963.5 0.1926  0.000000e+00 0.000050    Camk4   -0.113969
            OPC    19269             0             0        543         91          963.5 0.1707  0.000000e+00 0.000119      Npb    0.032562
Oligodendrocyte    19269             0             0        558         87          963.5 0.1043 8.057155e-183 0.000176      Dbi   -0.229620
       Pericyte    19269             0             0        414         59          963.5 0.3371  0.000000e+00 0.000771    Rpl28   -0.255309
           VLMC    19269             0             0        398         61          963.5 0.3166  0.000000e+00 0.000377    Gpaa1    0.144802

## MDE vs 观测
       celltype                      geneset  sd_score  observed_delta_Z  MDE_delta_Z(80%)  observed_over_MDE     p
      Astrocyte A_complement_synapse_pruning     0.114             0.083             0.373              0.222 0.442
      Astrocyte              B_DAM_microglia     0.379             0.365             1.240              0.294 0.330
      Astrocyte       C_mitochondrial_OXPHOS     0.306            -0.281             1.002              0.281 0.357
      Astrocyte     D_myelin_oligodendrocyte     0.392             0.103             1.283              0.080 0.786
    Endothelial A_complement_synapse_pruning     0.243             0.050             0.794              0.063 0.841
    Endothelial              B_DAM_microglia     0.282             0.112             0.923              0.121 0.682
    Endothelial       C_mitochondrial_OXPHOS     0.459             0.083             1.502              0.055 0.854
    Endothelial     D_myelin_oligodendrocyte     0.313             0.036             1.025              0.035 0.909
      Ependymal A_complement_synapse_pruning     0.265             0.367             0.866              0.424 0.116
      Ependymal              B_DAM_microglia     0.534             0.843             1.749              0.482 0.028
      Ependymal       C_mitochondrial_OXPHOS     0.494            -0.616             1.617              0.381 0.192
      Ependymal     D_myelin_oligodendrocyte     0.443             0.696             1.450              0.480 0.046
      Microglia A_complement_synapse_pruning     0.225             0.224             0.738              0.304 0.302
      Microglia              B_DAM_microglia     0.219            -0.032             0.718              0.045 0.881
      Microglia       C_mitochondrial_OXPHOS     0.211             0.280             0.690              0.406 0.158
      Microglia     D_myelin_oligodendrocyte     0.590             0.507             1.932              0.262 0.391
         Neuron A_complement_synapse_pruning     0.296            -0.196             0.967              0.202 0.501
         Neuron              B_DAM_microglia     0.443             0.541             1.451              0.373 0.199
         Neuron       C_mitochondrial_OXPHOS     0.360            -0.255             1.178              0.216 0.450
         Neuron     D_myelin_oligodendrocyte     0.421            -0.143             1.378              0.104 0.729
            OPC A_complement_synapse_pruning     0.150             0.044             0.492              0.090 0.768
            OPC              B_DAM_microglia     0.310             0.485             1.015              0.478 0.037
            OPC       C_mitochondrial_OXPHOS     0.513             0.113             1.679              0.067 0.823
            OPC     D_myelin_oligodendrocyte     0.224             0.201             0.734              0.273 0.327
Oligodendrocyte A_complement_synapse_pruning     0.403             0.305             1.320              0.231 0.454
Oligodendrocyte              B_DAM_microglia     0.544             0.549             1.782              0.308 0.278
Oligodendrocyte       C_mitochondrial_OXPHOS     0.376            -0.171             1.229              0.139 0.635
Oligodendrocyte     D_myelin_oligodendrocyte     0.576            -0.195             1.885              0.104 0.725
       Pericyte A_complement_synapse_pruning     0.265             0.376             0.868              0.433 0.076
       Pericyte              B_DAM_microglia     0.393             0.590             1.286              0.459 0.066
       Pericyte       C_mitochondrial_OXPHOS     0.184             0.317             0.601              0.527 0.017
       Pericyte     D_myelin_oligodendrocyte     0.292             0.138             0.957              0.144 0.626
           VLMC A_complement_synapse_pruning     0.184             0.299             0.603              0.496 0.034
           VLMC              B_DAM_microglia     0.318             0.046             1.041              0.045 0.881
           VLMC       C_mitochondrial_OXPHOS     0.415             0.627             1.359              0.461 0.061
           VLMC     D_myelin_oligodendrocyte     0.489             0.593             1.602              0.370 0.200
