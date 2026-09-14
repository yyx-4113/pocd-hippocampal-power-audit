"""
Step 5d — 细胞组成的动物级重算（本稿中心分析）
================================================================
背景
----
step5_sc_analysis.py 产出的 sc_celltype_composition.csv 用 Fisher 精确检验比较
Control vs Surgery 的细胞类型占比，但合计的是**池化细胞**：推断单元是细胞
(n ≈ 18,328)，而不是动物 (n = 3 vs 3)。这是单细胞文献里最具体的一类伪重复——
同一只动物身上的细胞不是独立重复，它们共享技术噪声与个体背景。

本步把推断单元降到动物：
  1. 先算每只动物各类细胞占该动物全部细胞的比例；
  2. 对逐动物比例做 Welch t 检验（不等方差）与 Mann-Whitney U，BH 校正；
  3. 同表并列保留池化细胞级 Fisher 结果作对照，量化"推断单元"本身造成的 p 值落差；
  4. 按非中心 t 分布精确计算 80% 功效所需每组动物数（与 Step 5f / 稿件 Fig 2 同一函数）。

公式说明
--------
- 效应量 Cohen's d = (mean_surg − mean_ctrl) / sqrt((var_surg + var_ctrl)/2)，
  组内方差用 ddof=1，两组 n 相等时该式即经典合并方差形式。输入用未取整的逐动物比例。
- req_n_per_group_exact：非中心 t，nc = d·sqrt(n/2)，df = 2n−2，α=0.05 双侧，power=0.80。
- req_n_per_group：**已被上者取代**的早期正态近似，保留仅为与历史产物对齐，勿引用。

输入
----
sc_out/GSE267933_processed.h5ad      QC + 注释后的单细胞对象（obs 含 celltype/sample/group）
sc_out/sc_celltype_composition.csv   池化细胞级 Fisher 结果（作对照列）

输出
----
sc_out/step5d_composition_animallevel.csv

用法
----
python step5d_composition_animallevel.py [h5ad路径] [输出目录]
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import nct
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
H5AD = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "sc_out", "GSE267933_processed.h5ad")
SCO = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "sc_out")
ALPHA, POWER = 0.05, 0.80


# ---------------------------------------------------------------- 非中心 t 功效
def power_at(ns, d):
    """向量化非中心 t 功效；ns 为每组 n（标量或数组），d 为效应量。"""
    ns = np.asarray(ns, dtype=float)
    dfs = 2 * ns - 2
    tcrit = nct.ppf(1 - ALPHA / 2, dfs, 0)
    nc = d * np.sqrt(ns / 2.0)
    return nct.sf(tcrit, dfs, nc) + nct.cdf(-tcrit, dfs, nc)


def req_n(d):
    """达到 80% 功效所需的最小每组 n（搜索 2..600）。"""
    d = abs(float(d)) if np.isfinite(d) else np.nan
    if not np.isfinite(d) or d == 0:
        return np.nan
    ns = np.arange(2, 601)
    ok = np.where(power_at(ns, d) >= POWER)[0]
    return int(ns[ok[0]]) if len(ok) else np.nan


def req_n_normalapprox(d):
    """早期正态近似（superseded）：2(z_{1-a/2}+z_{1-b})^2 / d^2，保留供追溯。"""
    d = abs(float(d)) if np.isfinite(d) else np.nan
    if not np.isfinite(d) or d == 0:
        return np.nan
    z = stats.norm.ppf(1 - ALPHA / 2) + stats.norm.ppf(POWER)
    return round(2 * z ** 2 / d ** 2, 1)


def cohen_d(a, b):
    """a = surgery，b = control；等 n 平均方差形式，组内 ddof=1。"""
    va = np.var(a, ddof=1) if len(a) > 1 else 0.0
    vb = np.var(b, ddof=1) if len(b) > 1 else 0.0
    sd = np.sqrt((va + vb) / 2.0)
    return (np.mean(a) - np.mean(b)) / sd if sd > 0 else np.nan


def main():
    print(f"[load] {H5AD}")
    import anndata as adata
    A = adata.read_h5ad(H5AD, backed="r")          # 只读 obs，不载入表达矩阵
    obs = A.obs[["celltype", "sample", "group"]].copy()
    A.file.close()
    obs["sample"] = obs["sample"].astype(str)
    obs["group"] = obs["group"].astype(str)

    # ---------------- 逐动物比例 ----------------
    counts = pd.crosstab(obs["celltype"], obs["sample"])
    counts = counts[sorted(counts.columns, key=lambda s: int(s) if s.isdigit() else s)]
    prop = 100.0 * counts / counts.sum(0)          # 细胞类型 x 动物，单位为 %

    grp = obs.drop_duplicates("sample").set_index("sample")["group"]
    ctrl_s = [s for s in prop.columns if grp[s] == "Control"]
    surg_s = [s for s in prop.columns if grp[s] == "Surgery"]
    print(f"[design] Control = {ctrl_s} | Surgery = {surg_s}")

    # 逐样本 QC 后细胞数（稿件 Appendix 的 "per sample 2,384-3,680" 即此表）
    per_sample = counts.sum(0).rename("n_cells_qc").to_frame()
    per_sample["group"] = [grp[s] for s in per_sample.index]
    per_sample.to_csv(os.path.join(SCO, "step5d_per_sample_cells.csv"))
    print(f"[per-sample QC cells] {per_sample.n_cells_qc.min()}-{per_sample.n_cells_qc.max()}"
          f" (total {per_sample.n_cells_qc.sum()})\n"
          + per_sample.to_string())

    # 池化细胞级 Fisher（对照列，来自 step5_sc_analysis.py）
    comp_file = os.path.join(SCO, "sc_celltype_composition.csv")
    cell_level = pd.read_csv(comp_file).set_index("celltype") if os.path.exists(comp_file) else None

    rows = []
    for ct in prop.index:
        c = prop.loc[ct, ctrl_s].values.astype(float)
        s = prop.loc[ct, surg_s].values.astype(float)
        welch_p = float(stats.ttest_ind(s, c, equal_var=False).pvalue)
        mwu_p = float(stats.mannwhitneyu(s, c, alternative="two-sided").pvalue)
        d = float(cohen_d(s, c))
        rows.append({
            "celltype": ct,
            "pct_ctrl_animalmean": round(float(c.mean()), 2),
            "pct_surg_animalmean": round(float(s.mean()), 2),
            "per_animal_ctrl": ";".join(f"{v:.2f}" for v in c),
            "per_animal_surg": ";".join(f"{v:.2f}" for v in s),
            "Welch_p_animal": welch_p,
            "MWU_p_animal": mwu_p,
            "Cohen_d_animal": round(d, 3),
            "req_n_per_group": req_n_normalapprox(d),
            "cellLevel_fisher_p": float(cell_level.loc[ct, "fisher_p"]) if cell_level is not None else np.nan,
            "cellLevel_BH_p": float(cell_level.loc[ct, "BH_p"]) if cell_level is not None else np.nan,
            "req_n_per_group_exact": req_n(d),
            "req_n_method": "exact non-central t, alpha=0.05, power=0.80",
        })

    df = pd.DataFrame(rows)
    df["BH_Welch_animal"] = multipletests(df["Welch_p_animal"], method="fdr_bh")[1]
    df["BH_MWU_animal"] = multipletests(df["MWU_p_animal"], method="fdr_bh")[1]

    cols = ["celltype", "pct_ctrl_animalmean", "pct_surg_animalmean",
            "per_animal_ctrl", "per_animal_surg", "Welch_p_animal", "MWU_p_animal",
            "Cohen_d_animal", "req_n_per_group", "cellLevel_fisher_p", "cellLevel_BH_p",
            "BH_Welch_animal", "BH_MWU_animal", "req_n_per_group_exact", "req_n_method"]
    df = df[cols].sort_values("Welch_p_animal").reset_index(drop=True)

    out = os.path.join(SCO, "step5d_composition_animallevel.csv")
    df.to_csv(out, index=False)
    pd.set_option("display.width", 220)
    print(df[["celltype", "pct_ctrl_animalmean", "pct_surg_animalmean",
              "Welch_p_animal", "BH_Welch_animal", "Cohen_d_animal",
              "req_n_per_group_exact", "cellLevel_BH_p"]].to_string(index=False))
    n_sig_animal = int((df.BH_Welch_animal < 0.05).sum())
    n_sig_cell = int((df.cellLevel_BH_p < 0.05).sum())
    print(f"\n[结果] 动物级 BH<0.05：{n_sig_animal}/9 类；池化细胞级 BH<0.05：{n_sig_cell}/9 类")
    print(f"[结果] 最大动物级效应：{df.iloc[0].celltype} "
          f"d={df.iloc[0].Cohen_d_animal:.2f}（需 n≈{df.iloc[0].req_n_per_group_exact}/组）")
    print("[Done]", out)


if __name__ == "__main__":
    main()
