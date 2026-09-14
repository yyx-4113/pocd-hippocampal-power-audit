# Step 1a 摘要 — GSE276942 时序 DEG + 轨迹聚类

- 基因数: 16885; 样本数: 17 (8 时间点)
- BH-padj<0.05 显著基因(并集): **0** (经 16885 基因多重校正后极弱)
- 探索性候选基因(p<0.05 & |log2FC|>0.5, 未校正, 并集): **1306**
  - 24H: up=66, down=51
  - 3D: up=26, down=910
  - 7D: up=302, down=23
  - 1M: up=0, down=0
  - 2M: up=0, down=0
  - 3M: up=0, down=0
  - 4M: up=0, down=0
- 轨迹聚类 K=5, 各簇: C1=778, C2=57, C3=323, C4=88, C5=60

## 各簇代表基因 (top 10 by |z| peak)

### Cluster 1 (n=778, 峰值时间=CON)
Irs2, Fosl2, Dusp5, Ttyh2, Homer1, Kctd17, Egr1, Fa2h, Junb, Fosb

### Cluster 2 (n=57, 峰值时间=CON)
9530082P21Rik, Phtf1os, Gm43536, Rab26os, Hlx, Gm43690, Iigp1, Usp27x, Ccdc67, Nod2

### Cluster 3 (n=323, 峰值时间=CON)
Zfp644, Kiz, Mcm8, Vamp7, Lysmd3, Slfn5, Zfp763, Tmem161b, Matr3, Ankrd49

### Cluster 4 (n=88, 峰值时间=CON)
Tlr2, Spint1, Rgag4, Gm42851, Lrrc9, Cd44, Dnmt3b, Gm15246, Pqlc3, Zfp558

### Cluster 5 (n=60, 峰值时间=CON)
5730522E02Rik, Rpp40, A330069E16Rik, Gbp7, Asf1b, Xlr3b, Nat1, Tnfrsf22, Gm11128, Gm43328
