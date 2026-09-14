"""
Step 5b-2 — pseudobulk 构建方式敏感性验证
========================================
动机:
  step5b 用 (log1p-normalized 均值) 构建 pseudobulk, 结论是"单细胞层机制线位移极弱"。
  该结论很强, 必须排除构建方式造成的假阴性。本步用更标准的
  counts 求和 -> CPM -> log1p 路线重做, 比较两种构建下 delta_Z 是否一致。

输入: GSE267933.h5ad (原始 counts) + sc_out/GSE267933_processed.h5ad (细胞类型标签)
输出: sc_out/step5b2_sensitivity.csv, fig_sc5b2_sensitivity.png
"""
import os, json, warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "sc_out")


def log(*a):
    print(*a, flush=True)


def score_tab(pb, genes):
    """pb: 样本 x 基因 (DataFrame); 返回 集合内基因跨样本 z 的样本均值"""
    sub = pb[genes]
    z = (sub - sub.mean(0)) / sub.std(0, ddof=1).replace(0, np.nan)
    return z.fillna(0.0).mean(1)


def main():
    log("[load] 原始 counts + 处理后标签")
    a0 = sc.read_h5ad(os.path.join(HERE, "GSE267933.h5ad"))
    a1 = sc.read_h5ad(os.path.join(WORK, "GSE267933_processed.h5ad"))
    a0.layers["counts"] = a0.X.copy()

    # 复现 step5 的 QC, 保证细胞集合一致
    sc.pp.filter_cells(a0, min_genes=200)
    sc.pp.filter_genes(a0, min_cells=3)
    a0 = a0[a0.obs.pct_counts_mt < 20].copy()
    common = a1.obs_names.intersection(a0.obs_names)
    log(f"[align] 处理后 {a1.n_obs} 细胞, 复现 QC 后 {a0.n_obs}, 交集 {len(common)}")
    a0 = a0[common].copy()
    a0.obs["celltype"] = a1.obs.loc[common, "celltype"].values
    a0.obs["sample"] = a0.obs["sample"].astype(str)
    a0.obs["group"] = a0.obs["group"].astype(str)

    C = a0.layers["counts"]                       # 细胞 x 全基因 原始 counts (稀疏)
    genes = np.asarray(a0.var_names)
    gidx = {g: i for i, g in enumerate(genes)}

    pw = json.load(open(os.path.join(HERE, "step2_pathway_sets.json"), encoding="utf-8"))
    sets = {k: [g for g in v if g in gidx] for k, v in pw.items()}

    # ---- counts 求和 -> CPM -> log1p 的 pseudobulk ----
    celltypes = [c for c in a0.obs.celltype.value_counts().index
                 if (a0.obs.celltype.values == c).sum() >= 100]
    rows, pbs = [], {}
    for ct in celltypes:
        tot, samples, groups = [], [], []
        for s in sorted(pd.unique(a0.obs["sample"].values)):
            m = (a0.obs.celltype.values == ct) & (a0.obs["sample"].values == s)
            if m.sum() < 20:
                continue
            tot.append(np.asarray(C[m].sum(0)).ravel())
            samples.append(s); groups.append(a0.obs["group"].values[m][0])
        if len(samples) < 6 or len(set(groups)) < 2:
            continue
        T = pd.DataFrame(tot, index=samples, columns=genes)
        cpm = T.div(T.sum(1), axis=0) * 1e6
        pb = np.log1p(cpm)
        pbs[ct] = pb
        grp = pd.Series(groups, index=samples)
        for sname, glist in sets.items():
            gp = [g for g in glist if g in pb.columns]
            if len(gp) < 3:
                continue
            sc_ = score_tab(pb, gp)
            surg = sc_[(grp == "Surgery").values].values
            ctrl = sc_[(grp == "Control").values].values
            t, p = stats.ttest_ind(surg, ctrl, equal_var=False)
            rows.append({"celltype": ct, "geneset": sname, "method": "counts_sum_CPM",
                         "n_genes": len(gp), "delta_Z": float(surg.mean() - ctrl.mean()),
                         "t": float(t), "p": float(p) if np.isfinite(p) else 1.0})
    b2 = pd.DataFrame(rows)

    ref = pd.read_csv(os.path.join(WORK, "step5b_geneset_stats.csv"))
    ref = ref[ref.geneset.str.startswith(("A_", "B_", "C_", "D_"))]
    ref = ref.assign(method="lognorm_mean")[["celltype", "geneset", "method", "n_genes", "delta_Z", "p"]]
    both = pd.concat([ref, b2], ignore_index=True)
    both["BH_p"] = multipletests(both.p.fillna(1), method="fdr_bh")[1]
    both.round(5).to_csv(os.path.join(WORK, "step5b2_sensitivity.csv"), index=False)

    piv = both.pivot_table(index=["geneset", "celltype"], columns="method", values="delta_Z").dropna()
    r = stats.pearsonr(piv["lognorm_mean"], piv["counts_sum_CPM"])
    rho = stats.spearmanr(piv["lognorm_mean"], piv["counts_sum_CPM"]).statistic
    log(f"\n[一致性] 两构建法 delta_Z 相关: Pearson r={r[0]:.3f} (p={r[1]:.2g}), Spearman rho={rho:.3f}, n={len(piv)}")
    log("[各方法名义显著 (p<0.05) 计数] " + str(both[both.p < 0.05].groupby("method").size().to_dict()))
    log("[各方法 BH<0.05 计数] " + str(both[both.BH_p < 0.05].groupby("method").size().to_dict()))

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2))
    axes[0].scatter(piv["lognorm_mean"], piv["counts_sum_CPM"], s=26, c="#4c72b0")
    lim = np.nanmax(np.abs(piv.values)) * 1.1
    axes[0].plot([-lim, lim], [-lim, lim], ls="--", c="grey", lw=0.8)
    axes[0].axhline(0, ls=":", c="grey", lw=0.7); axes[0].axvline(0, ls=":", c="grey", lw=0.7)
    axes[0].set_xlabel("delta_Z (lognorm cell-mean)")
    axes[0].set_ylabel("delta_Z (counts sum -> CPM)")
    axes[0].set_title(f"Pseudobulk construction sensitivity (r={r[0]:.2f})")
    order = ["A_complement_synapse_pruning", "B_DAM_microglia",
             "C_mitochondrial_OXPHOS", "D_myelin_oligodendrocyte"]
    b2o = b2[b2.geneset.isin(order)].copy()
    hp = b2o.pivot_table(index="geneset", columns="celltype", values="delta_Z").reindex(order)
    sns.heatmap(hp, cmap="RdBu_r", center=0, annot=True, fmt=".2f", annot_kws={"size": 6.5},
                ax=axes[1], cbar_kws={"label": "delta_Z (counts sum)"})
    axes[1].set_title("Mechanism lines (counts-sum pseudobulk)")
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc5b2_sensitivity.png"), dpi=150)
    plt.close(fig)
    log("[Done] Step 5b-2 输出写入", WORK)


if __name__ == "__main__":
    main()
