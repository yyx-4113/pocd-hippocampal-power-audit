# Step 3 摘要 — 双机器学习锁靶 + ROC

候选池: 60 基因 (Step 2: 四机制线 A/B/C/D)

## 双模型结果
### discovery_GSE276942 (LASSO C*=0.229, λ*=4.38; Boruta=BorutaPy)
- LASSO 选中 (2): [np.str_('Fos'), np.str_('Junb')]
- Boruta 确认 (9): [np.str_('Arc'), np.str_('Cx3cr1'), np.str_('Egr1'), np.str_('Egr2'), np.str_('Fos'), np.str_('Junb'), np.str_('Myo1d'), np.str_('Opalin'), np.str_('Prr5l')]
- 交集 (LASSO∩Boruta) = [np.str_('Fos'), np.str_('Junb')]
- LASSO bootstrap 稳定性>=0.6 (0): []

### integrated_ComBat (LASSO C*=1.914, λ*=0.52; Boruta=BorutaPy)
- LASSO 选中 (11): [np.str_('Cx3cl1'), np.str_('Gltp'), np.str_('Irf8'), np.str_('Junb'), np.str_('Kctd17'), np.str_('Mog'), np.str_('Prr18'), np.str_('Runx2'), np.str_('Sgk1'), np.str_('Sla2'), np.str_('Tyrobp')]
- Boruta 确认 (13): [np.str_('Bcas1'), np.str_('C5ar1'), np.str_('Enpp6'), np.str_('Gltp'), np.str_('Myo1d'), np.str_('Opalin'), np.str_('Pmp22'), np.str_('Prr5l'), np.str_('Runx2'), np.str_('Sgk1'), np.str_('Sla2'), np.str_('Sox10'), np.str_('Tmem63a')]
- 交集 (LASSO∩Boruta) = [np.str_('Gltp'), np.str_('Runx2'), np.str_('Sgk1'), np.str_('Sla2')]
- LASSO bootstrap 稳定性>=0.6 (3): [np.str_('Junb'), np.str_('Sgk1'), np.str_('Sla2')]

## 特征基因汇总 (两训练集并集, n=6)
[np.str_('Fos'), np.str_('Gltp'), np.str_('Junb'), np.str_('Runx2'), np.str_('Sgk1'), np.str_('Sla2')]
- **稳健核心** (LASSO bootstrap 频率>=0.6): [np.str_('Junb'), np.str_('Sgk1'), np.str_('Sla2')]
- 两训练集特征交集 = 0 个 -> **不稳定**(小样本, 特征选择须谨慎)

## 特征基因稳定性/重要性 (step3_feature_genes.csv)
 gene  is_final  robust_core  discovery_feature  integrated_feature  discovery_lasso_stab  integrated_lasso_stab  discovery_rf_rank  integrated_rf_rank  discovery_boruta  integrated_boruta
 Sgk1      True         True              False                True                 0.175                  0.875                 22                   8             False               True
 Sla2      True         True              False                True                 0.000                  0.780                 20                   1             False               True
 Junb      True         True               True               False                 0.360                  0.730                  2                  24              True              False
Runx2      True        False              False                True                 0.000                  0.455                 58                   7             False               True
 Gltp      True        False              False                True                 0.000                  0.400                 10                  10             False               True
  Fos      True        False               True               False                 0.265                  0.030                  6                  54              True              False

## 重复交叉验证 (5 折 x 20 次)
            dataset  n  n_feat  repeatedCV_AUC_mean  repeatedCV_AUC_sd
discovery_GSE276942 17       6                  1.0                0.0
  integrated_ComBat 27       6                  1.0                0.0

## 跨数据集 train -> test AUC
              train                test  n_test  n_feat   AUC
discovery_GSE276942           GSE215410       4       6 0.000
discovery_GSE276942           GSE174412       6       6 0.889
discovery_GSE276942 discovery_GSE276942      17       6 1.000
discovery_GSE276942   integrated_ComBat      27       6 1.000
  integrated_ComBat           GSE215410       4       6 0.250
  integrated_ComBat           GSE174412       6       6 1.000
  integrated_ComBat discovery_GSE276942      17       6 1.000
  integrated_ComBat   integrated_ComBat      27       6 1.000

## 单基因 ROC AUC
dataset  discovery_GSE276942  integrated_ComBat  validation_GSE174412  validation_GSE215410
gene                                                                                       
Fos                    1.000              0.574                 1.000                  1.00
Gltp                   0.865              0.864                 0.889                  1.00
Junb                   1.000              0.753                 0.556                  1.00
Runx2                  0.712              0.877                 1.000                  1.00
Sgk1                   0.865              0.846                 1.000                  0.75
Sla2                   0.885              0.932                 1.000                  1.00

## ⚠️ 关键局限 (必须写入稿件 Limitation)
1. **样本量极小**: 发现集 n=17, 整合矩阵 n=27, 验证集 n=4/6 -> 所有 ML 结果均为**探索性**, AUC 的高值部分来自过拟合。
2. **选择偏倚/双重浸用(double dipping)**: 候选基因本身由这些数据集的 DEG/元分析选出, 再用同一批数据训练+评估, AUC 必然偏高; 真实泛化须依赖**独立外部队列**(本方案 Step 8 临床血样或未来数据集)。
3. **特征选择不稳定**: 两个训练集选出的特征集交叠很小, 说明在当前样本量下特征选择不可复现; 报告应以"机制方向"而非"精确基因列表"为主。
