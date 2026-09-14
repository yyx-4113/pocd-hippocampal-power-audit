# -*- coding: utf-8 -*-
"""step5h — 计算稿件修订所需的、可复现的不确定性数字:
   (1) 组成效应(动物级)的 Welch 95% CI + Cohen d 95% CI;
   (2) GSE199318 星形胶质 C3 的 log2FC 95% CI;
   (3) Welch-Satterthwaite 最坏情形(方差比极端, df=n-1)下的 MDE 带上界,
       与等方差(df=2n-2)最优情形对照 -> 图 2 阴影带 + 正文数字。
所有数字均从既有产物文件读取, 不硬编码。
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
SCO = os.path.join(HERE, "sc_out")
ALPHA, POWER = 0.05, 0.80


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    v1, v2 = a.var(ddof=1), b.var(ddof=1)
    m1, m2 = a.mean(), b.mean()
    se = np.sqrt(v1 / n1 + v2 / n2)
    df = (v1 / n1 + v2 / n2) ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    t = (m1 - m2) / se
    return m1 - m2, se, df, t


def mde_noncentral(n, df):
    ds = np.linspace(0.01, 14, 8000)
    tcrit = stats.nct.ppf(1 - ALPHA / 2, df, 0)
    nc = ds * np.sqrt(n / 2.0)
    pw = stats.nct.sf(tcrit, df, nc) + stats.nct.cdf(-tcrit, df, nc)
    ok = np.where(pw >= POWER)[0]
    return float(ds[ok[0]]) if len(ok) else np.nan


print("===== (1) 组成效应 动物级 95% CI =====")
d5 = pd.read_csv(os.path.join(SCO, "step5d_composition_animallevel.csv")).set_index("celltype")
for ct in ["Microglia", "Ependymal"]:
    row = d5.loc[ct]
    ctrl = [float(x) for x in str(row.per_animal_ctrl).replace("[", "").replace("]", "").split(";")]
    surg = [float(x) for x in str(row.per_animal_surg).replace("[", "").replace("]", "").split(";")]
    diff, se, df, t = welch(ctrl, surg)
    tc = stats.t.ppf(1 - ALPHA / 2, df)
    lo, hi = diff - tc * se, diff + tc * se
    sp = np.sqrt((np.var(ctrl, ddof=1) + np.var(surg, ddof=1)) / 2.0)  # pooled SD
    d = diff / sp
    dlo, dhi = lo / sp, hi / sp
    print(f"{ct}: diff={diff:+.2f} pp, 95%CI [{lo:.1f}, {hi:.1f}], Welch df={df:.2f}, "
          f"t={t:.2f}, p={row.Welch_p_animal:.3f}")
    print(f"   Cohen d={d:.2f}, 95%CI [{dlo:.2f}, {dhi:.2f}] (from diff CI / pooled SD={sp:.2f})")

print("\n===== (2) GSE199318 C3 log2FC 95% CI =====")
deg = pd.read_csv(os.path.join(HERE, "step2c_per_dataset_DEG.csv"))
c3 = deg[(deg.dataset == "GSE199318_astro") & (deg.gene == "C3")].iloc[0]
fc, p = float(c3.log2FC), float(c3.p)
df_c = 3 + 3 - 2  # Welch approx equal n=3 -> df in [2,4]; use 4 (conservative t)
# recover t from two-sided p
t_abs = stats.t.ppf(1 - p / 2, df_c)
se = abs(fc) / t_abs if t_abs > 0 else np.nan
tc = stats.t.ppf(1 - ALPHA / 2, df_c)
lo, hi = fc - tc * se, fc + tc * se
print(f"C3 log2FC={fc:+.2f}, p={p:.2f}, recovered t={t_abs:.2f}, SE={se:.2f}, "
      f"95%CI [{lo:.1f}, {hi:.1f}] (df={df_c})")

print("\n===== (3) Welch 最坏情形 MDE 带 (df=n-1) vs 等方差 (df=2n-2) =====")
mde_csv = pd.read_csv(os.path.join(SCO, "step5c_mde_by_samplesize.csv"))
for _, r in mde_csv.iterrows():
    n = int(r["每组n"])
    best = float(r["可检出标准化效应 d (80% power)"])
    worst = mde_noncentral(n, max(n - 1, 1))  # 方差比极端 -> df 下界 n-1
    print(f"n={n:2d}: equal-var(df={2*n-2}) d={best:.2f} | Welch-worst(df={max(n-1,1)}) d={worst:.2f}")
