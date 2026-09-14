# Step 5b — 单细胞逐样本 pseudobulk 基因集统计

方法: 逐 (细胞类型 x 样本) 聚合 lognorm 均值为 pseudobulk; 集合内基因跨样本 z 后取均值
      作为集合得分; 3v3 上 Welch t 检验 + BH。细胞不是重复单元, 故不做细胞级检验。

## 四条机制线 x 细胞类型 (delta Z = Surgery - Control)
       celltype                      geneset  n_genes  mean_Z_ctrl  mean_Z_surg  delta_Z      p   BH_p
      Astrocyte A_complement_synapse_pruning      152      -0.0415       0.0415   0.0830 0.4421 0.7245
      Astrocyte              B_DAM_microglia      116      -0.1824       0.1824   0.3647 0.3302 0.7026
      Astrocyte       C_mitochondrial_OXPHOS      189       0.1407      -0.1407  -0.2815 0.3568 0.7065
      Astrocyte     D_myelin_oligodendrocyte       74      -0.0515       0.0515   0.1031 0.7861 0.8744
    Endothelial A_complement_synapse_pruning      152      -0.0250       0.0250   0.0500 0.8407 0.8854
    Endothelial              B_DAM_microglia      116      -0.0559       0.0559   0.1118 0.6818 0.8250
    Endothelial       C_mitochondrial_OXPHOS      189      -0.0414       0.0414   0.0829 0.8538 0.8898
    Endothelial     D_myelin_oligodendrocyte       74      -0.0178       0.0178   0.0356 0.9088 0.9181
      Ependymal A_complement_synapse_pruning      152      -0.1835       0.1835   0.3671 0.1163 0.5293
      Ependymal              B_DAM_microglia      116      -0.4214       0.4214   0.8427 0.0280 0.4574
      Ependymal       C_mitochondrial_OXPHOS      189       0.3078      -0.3079  -0.6157 0.1918 0.5499
      Ependymal     D_myelin_oligodendrocyte       74      -0.3480       0.3480   0.6960 0.0462 0.4574
      Microglia A_complement_synapse_pruning      152      -0.1121       0.1121   0.2243 0.3020 0.7026
      Microglia              B_DAM_microglia      116       0.0161      -0.0161  -0.0323 0.8814 0.8996
      Microglia       C_mitochondrial_OXPHOS      189      -0.1400       0.1400   0.2801 0.1581 0.5417
      Microglia     D_myelin_oligodendrocyte       74      -0.2533       0.2533   0.5067 0.3910 0.7137
         Neuron A_complement_synapse_pruning      152       0.0978      -0.0978  -0.1956 0.5012 0.7633
         Neuron              B_DAM_microglia      116      -0.2704       0.2704   0.5409 0.1994 0.5499
         Neuron       C_mitochondrial_OXPHOS      189       0.1273      -0.1273  -0.2546 0.4500 0.7245
         Neuron     D_myelin_oligodendrocyte       74       0.0714      -0.0714  -0.1428 0.7286 0.8486
            OPC A_complement_synapse_pruning      152      -0.0220       0.0220   0.0440 0.7680 0.8640
            OPC              B_DAM_microglia      116      -0.2425       0.2425   0.4850 0.0370 0.4574
            OPC       C_mitochondrial_OXPHOS      189      -0.0564       0.0564   0.1128 0.8226 0.8854
            OPC     D_myelin_oligodendrocyte       74      -0.1003       0.1003   0.2006 0.3266 0.7026
Oligodendrocyte A_complement_synapse_pruning      152      -0.1526       0.1526   0.3051 0.4537 0.7245
Oligodendrocyte              B_DAM_microglia      116      -0.2745       0.2745   0.5490 0.2775 0.7026
Oligodendrocyte       C_mitochondrial_OXPHOS      189       0.0855      -0.0855  -0.1710 0.6345 0.8053
Oligodendrocyte     D_myelin_oligodendrocyte       74       0.0976      -0.0976  -0.1951 0.7249 0.8486
       Pericyte A_complement_synapse_pruning      152      -0.1878       0.1878   0.3757 0.0761 0.4707
       Pericyte              B_DAM_microglia      116      -0.2950       0.2950   0.5900 0.0664 0.4696
       Pericyte       C_mitochondrial_OXPHOS      189      -0.1583       0.1583   0.3167 0.0175 0.4574
       Pericyte     D_myelin_oligodendrocyte       74      -0.0689       0.0689   0.1377 0.6262 0.8053
           VLMC A_complement_synapse_pruning      152      -0.1496       0.1496   0.2993 0.0342 0.4574
           VLMC              B_DAM_microglia      116      -0.0232       0.0232   0.0464 0.8806 0.8996
           VLMC       C_mitochondrial_OXPHOS      189      -0.3135       0.3135   0.6270 0.0606 0.4613
           VLMC     D_myelin_oligodendrocyte       74      -0.2965       0.2965   0.5930 0.2000 0.5499

## 定向 panel
       celltype         geneset  n_genes  delta_Z      p   BH_p
      Astrocyte Complement_core        9   0.8158 0.0747 0.4707
      Astrocyte        DAM_core        9  -0.0599 0.8260 0.8854
      Astrocyte     Homeostatic        7   0.8835 0.3691 0.7137
      Astrocyte    IEG_activity        8   0.6085 0.3965 0.7137
      Astrocyte     Mito_OXPHOS        8  -0.9478 0.0101 0.4574
      Astrocyte   Myelin_struct       10   0.5795 0.4959 0.7633
      Astrocyte  Phagocytic_rec        9   0.2745 0.4517 0.7245
    Endothelial Complement_core        9   0.3196 0.3780 0.7137
    Endothelial        DAM_core        9  -0.1537 0.6975 0.8320
    Endothelial     Homeostatic        7   0.3736 0.6834 0.8250
    Endothelial    IEG_activity        8   0.3325 0.6715 0.8250
    Endothelial     Mito_OXPHOS        8  -0.3559 0.4521 0.7245
    Endothelial   Myelin_struct       10   0.3383 0.5666 0.8053
    Endothelial  Phagocytic_rec        9   0.1119 0.8222 0.8854
      Ependymal Complement_core        9   1.3499 0.0445 0.4574
      Ependymal        DAM_core        9   0.5682 0.2290 0.6126
      Ependymal     Homeostatic        7   1.3191 0.1176 0.5293
      Ependymal    IEG_activity        8   1.1985 0.1536 0.5417
      Ependymal     Mito_OXPHOS        8  -0.8136 0.1474 0.5417
      Ependymal   Myelin_struct       10   1.3626 0.0055 0.4574
      Ependymal  Phagocytic_rec        9   0.7094 0.1792 0.5499
      Microglia Complement_core        9  -0.1476 0.7635 0.8640
      Microglia        DAM_core        9   0.0732 0.8389 0.8854
      Microglia     Homeostatic        7  -0.3550 0.6125 0.8053
      Microglia    IEG_activity        8   0.4702 0.2759 0.7026
      Microglia     Mito_OXPHOS        8  -0.1470 0.6777 0.8250
      Microglia   Myelin_struct       10   0.6120 0.5277 0.7893
      Microglia  Phagocytic_rec        9  -0.0925 0.6277 0.8053
         Neuron Complement_core        9   0.6317 0.1645 0.5417
         Neuron        DAM_core        9   0.2891 0.4497 0.7245
         Neuron     Homeostatic        7   1.2830 0.1410 0.5417
         Neuron    IEG_activity        8   1.0290 0.0308 0.4574
         Neuron     Mito_OXPHOS        8  -0.5155 0.1696 0.5417
         Neuron   Myelin_struct       10   0.4639 0.5342 0.7893
         Neuron  Phagocytic_rec        9   0.3112 0.5765 0.8053
            OPC Complement_core        9   0.6375 0.0999 0.5283
            OPC        DAM_core        9   0.1708 0.5931 0.8053
            OPC     Homeostatic        7   1.0671 0.1420 0.5417
            OPC    IEG_activity        8   0.7969 0.1014 0.5283
            OPC     Mito_OXPHOS        8  -0.0453 0.9417 0.9417
            OPC   Myelin_struct       10   1.1144 0.1550 0.5417
            OPC  Phagocytic_rec        9   0.3949 0.3029 0.7026
Oligodendrocyte Complement_core        9   0.6927 0.1662 0.5417
Oligodendrocyte        DAM_core        9   0.4158 0.3336 0.7026
Oligodendrocyte     Homeostatic        7   0.7694 0.4277 0.7245
Oligodendrocyte    IEG_activity        8   0.5204 0.4720 0.7418
Oligodendrocyte     Mito_OXPHOS        8  -0.4093 0.3292 0.7026
Oligodendrocyte   Myelin_struct       10  -0.5827 0.0383 0.4574
Oligodendrocyte  Phagocytic_rec        9   0.5811 0.3491 0.7065
       Pericyte Complement_core        9   0.9019 0.0883 0.5141
       Pericyte        DAM_core        9   0.5406 0.3874 0.7137
       Pericyte     Homeostatic        7   1.4095 0.0587 0.4613
       Pericyte    IEG_activity        8   1.0703 0.1086 0.5293
       Pericyte     Mito_OXPHOS        8   0.3501 0.5440 0.7919
       Pericyte   Myelin_struct       10   0.2116 0.7377 0.8492
       Pericyte  Phagocytic_rec        9   1.0161 0.0578 0.4613
           VLMC Complement_core        9   0.5052 0.2926 0.7026
           VLMC        DAM_core        9   0.2226 0.6315 0.8053
           VLMC     Homeostatic        7   0.1695 0.6345 0.8053
           VLMC    IEG_activity        8  -0.3143 0.6182 0.8053
           VLMC     Mito_OXPHOS        8   0.4606 0.3286 0.7026
           VLMC   Myelin_struct       10   0.8155 0.3538 0.7065
           VLMC  Phagocytic_rec        9   0.5355 0.1833 0.5499

**局限**: n=3v3, 功效极低; 本步为方向一致性验证, 非确证性统计。
