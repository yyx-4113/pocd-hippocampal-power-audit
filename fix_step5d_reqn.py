"""[已作废] 用精确非中心 t 重算 step5d CSV 的 req_n_per_group (原为两样本近似公式, 与 step5f 的精确值不一致)。

本脚本是在 step5d 尚无独立脚本时打的一次性补丁。现 step5d_composition_animallevel.py
已原生输出 req_n_per_group_exact 与 req_n_per_group (正态近似, 保留供追溯),
因此本文件不再需要, 保留仅作历史记录; 请勿再运行。
"""
import numpy as np
import pandas as pd
from scipy.stats import nct

F = "sc_out/step5d_composition_animallevel.csv"
df = pd.read_csv(F)
ALPHA, POWER = 0.05, 0.80


def power_at(ns, d):
    ns = np.asarray(ns, dtype=float)
    dfs = 2 * ns - 2
    tcrit = nct.ppf(1 - ALPHA / 2, dfs, 0)
    nc = d * np.sqrt(ns / 2.0)
    return nct.sf(tcrit, dfs, nc) + nct.cdf(-tcrit, dfs, nc)


def req_n(d):
    if not np.isfinite(d) or d == 0:
        return np.nan
    ns = np.arange(2, 601)
    ok = np.where(power_at(ns, abs(d)) >= POWER)[0]
    return int(ns[ok[0]]) if len(ok) else np.nan


df["req_n_per_group_exact"] = [req_n(d) for d in df["Cohen_d_animal"]]
df["req_n_method"] = "exact non-central t, alpha=0.05, power=0.80"
df.to_csv(F, index=False)
print(df[["celltype", "Cohen_d_animal", "req_n_per_group", "req_n_per_group_exact"]].to_string(index=False))
