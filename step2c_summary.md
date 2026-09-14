# Step 2c 摘要 — 扩展 bulk 元分析 (并入 GSE199318 星形胶质)

## 候选筛选结果 (6 个 GSE 实测)
- GSE95426 海马 POCD 组织 6v6 (定制 Agilent 阵列): 平台无基因符号列, GPL22782.annot.gz 沙箱 FTP 全路径 404 -> **排除**
- GSE165798 PND 海马 (circRNA 芯片): 平台表仅 circRNA ID 无 mRNA 符号 -> **排除**
- GSE303920 海马 (标称 bulk): RAW 实为 10x scRNA -> **排除**
- GSE316433 PBMC / GSE234493 snRNA: 非海马 bulk -> **排除**
- **GSE199318 分选海马星形胶质 RNA-seq (Gene Symbol+计数): 唯一干净可用** -> 并入

## 纳入数据集
- 原 3 套组织层: GSE276942, GSE215410, GSE174412
- 新增: GSE199318 星形胶质 (SEV=麻醉对照 n=3, LA=剖腹手术 PND n=3)

## 四机制线基因集统计 (扩展前后)
                         set  n_members  mean_meta_Z_3ds  frac_down_3ds      t_p_3ds  stouffer_p_3ds     BH_p_3ds  mean_meta_Z_4ds  frac_down_4ds      t_p_4ds  stouffer_p_4ds     BH_p_4ds  delta_mean_Z
A_complement_synapse_pruning        149           -0.285          0.711 9.547469e-04    5.038303e-04 9.547469e-04           -0.277          0.703 1.931432e-03    8.458200e-04 1.931432e-03         0.008
             B_DAM_microglia        112           -0.906          0.839 3.228731e-16    9.398802e-22 6.457461e-16           -0.893          0.829 4.718491e-16    4.963638e-21 9.436983e-16         0.013
      C_mitochondrial_OXPHOS        187           -0.504          0.797 1.787426e-20    5.256096e-12 7.149705e-20           -0.512          0.804 9.163965e-21    3.924861e-12 3.665586e-20        -0.008
    D_myelin_oligodendrocyte         73           -1.484          0.904 5.853693e-15    7.869118e-37 7.804924e-15           -1.471          0.903 1.011692e-14    9.350312e-36 1.348922e-14         0.013

## 关键判读
- 见 step2c_geneset_compare.csv 与正文。
