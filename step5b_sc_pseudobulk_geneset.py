"""
Step 5b — 单细胞逐样本 pseudobulk 基因集统计 (与 bulk 层同轴闭环)
=================================================================
动机:
  Step 5-6 的细胞级打分 (DAM_score/Homeo_score 均值) 不能做统计推断 —— 细胞不是独立重复。
  本步把"细胞级"降到"样本级": 每个 (细胞类型 x 样本) 聚合成一条 pseudobulk 表达向量,
  再算基因集水平统计量 (集合内基因 z 均值), 在 3v3 上做 Welch t 检验 + BH。
  这样四条机制线 (A/B/C/D) 在 bulk 与单细胞两层得到同轴、可比较的效应量。

输入: sc_out/GSE267933_processed.h5ad (含 raw, 全基因 lognorm), step2_pathway_sets.json
输出: sc_out/step5b_*.csv, fig_sc5b_*.png, step5b_summary.md
注意: n=3v3, 检验功效极低; 本步定位为"方向一致性验证", 不是确证性统计。
"""
import os, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "sc_out")
SEED = 42

# 微胶质/少突 定向 panel (与 Step 2 的机制线互补, 用于精细定位)
PANELS = {
    "DAM_core":        ["Trem2", "Apoe", "Cst7", "Lpl", "Clec7a", "Itgax", "Spp1", "Gpnmb", "Cd9", "Timp2"],
    "Homeostatic":     ["P2ry12", "Cx3cr1", "Tmem119", "Siglech", "Hexb", "Gpr34", "Mafb"],
    "Phagocytic_rec":  ["C5ar1", "Cx3cr1", "Megf10", "Tyrobp", "Cx3cl1", "Ctss", "Cd68", "Mertk", "Gas6"],
    "Complement_core": ["C1qa", "C1qb", "C1qc", "C3", "C3ar1", "C4b", "Serping1", "Cfb", "Cfh", "Cfi"],
    "Myelin_struct":   ["Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mal", "Mobp", "Opalin", "Cldn11", "Bcas1"],
    "IEG_activity":    ["Fos", "Fosb", "Junb", "Egr1", "Egr2", "Arc", "Jun", "Nr4a1"],
    "Mito_OXPHOS":     ["Ndufa1", "Ndufb2", "Cox7b", "Cox10", "Atp5f1", "Sdha", "Uqcrb", "Atp5a1"],
}


def log(*a):
    print(*a, flush=True)


def pseudobulk_matrix(adx, ct_name, gene_idx):
    """返回 (样本 x 基因) DataFrame, 值 = 该样本该细胞类型内 lognorm 均值"""
    mask = (adx.obs.celltype.values == ct_name)
    if mask.sum() < 20:
        return None, None
    Xs = adx.raw.X[mask][:, gene_idx]
    samples = adx.obs["sample"].values[mask]
    groups = adx.obs["group"].values[mask]
    rows, meta = [], []
    for s in sorted(pd.unique(samples)):
        m = (samples == s)
        rows.append(np.asarray(Xs[m].mean(0)).ravel())
        meta.append({"sample": s, "group": groups[m][0], "n_cells": int(m.sum())})
    mix = pd.DataFrame(meta)
    # 只保留两组均有样本的 (3v3)
    if mix.group.nunique() < 2:
        return None, None
    return pd.DataFrame(rows, index=mix["sample"].values), mix.set_index("sample")


def set_score(pb, genes_present):
    """z 逐基因(跨样本) -> 集合内均值 = 集合水平得分; 同时给 下调比例"""
    sub = pb[genes_present]
    mu = sub.mean(0)
    sd = sub.std(0, ddof=1).replace(0, np.nan)
    z = (sub - mu) / sd
    z = z.fillna(0.0)
    score = z.mean(1)
    frac_down = (z < 0).mean(1)
    return score, frac_down, z


def main():
    log("[load] sc_out/GSE267933_processed.h5ad")
    a = ad.read_h5ad(os.path.join(WORK, "GSE267933_processed.h5ad"))
    RV = list(a.raw.var_names)
    gidx = {g: i for i, g in enumerate(RV)}
    log(f"[info] {a.n_obs} 细胞 x {a.n_vars} HVG; raw {len(RV)} 基因; 类型={sorted(a.obs.celltype.unique())}")

    pw = json.load(open(os.path.join(HERE, "step2_pathway_sets.json"), encoding="utf-8"))
    sets = {k: [g for g in v if g in gidx] for k, v in pw.items()}
    sets.update({k: [g for g in v if g in gidx] for k, v in PANELS.items()})

    all_genes = sorted({g for v in sets.values() for g in v})
    gi = [gidx[g] for g in all_genes]
    gene_pos = {g: i for i, g in enumerate(all_genes)}

    celltypes = [c for c in a.obs.celltype.value_counts().index if (a.obs.celltype.values == c).sum() >= 100]
    log(f"[pseudobulk] 纳入细胞类型 (>=100 细胞): {celltypes}")

    # ---- 逐 (细胞类型, 基因集) 统计 ----
    rows = []
    store = {}          # (ct, set) -> (score series, z table)
    for ct in celltypes:
        pb, mix = pseudobulk_matrix(a, ct, gi)
        if pb is None:
            continue
        pb.columns = all_genes
        log(f"[pseudobulk] {ct}: {pb.shape[0]} 样本 x {pb.shape[1]} 基因 (n_cells={dict(mix.n_cells)})")
        for sname, glist in sets.items():
            gp = [g for g in glist if g in pb.columns]
            if len(gp) < 3:
                continue
            score, frac_down, ztab = set_score(pb, gp)
            surg = score[(mix.group == "Surgery").values].values
            ctrl = score[(mix.group == "Control").values].values
            t, p = stats.ttest_ind(surg, ctrl, equal_var=False)
            rows.append({
                "celltype": ct, "geneset": sname, "n_genes": len(gp),
                "mean_Z_ctrl": float(ctrl.mean()), "mean_Z_surg": float(surg.mean()),
                "delta_Z": float(surg.mean() - ctrl.mean()),
                "frac_down_surg": float(frac_down[(mix.group == "Surgery").values].mean()),
                "t": float(t), "p": float(p) if np.isfinite(p) else 1.0,
            })
            store[(ct, sname)] = (score, ztab, mix)
        # 存逐样本矩阵 (便于复核)
        pb.round(4).to_csv(os.path.join(WORK, f"step5b_pseudobulk_{ct}.csv"))

    res = pd.DataFrame(rows)
    res["BH_p"] = multipletests(res.p.fillna(1), method="fdr_bh")[1]
    res = res.sort_values(["celltype", "geneset"])
    res.round(5).to_csv(os.path.join(WORK, "step5b_geneset_stats.csv"), index=False)
    log("\n[基因集 x 细胞类型 统计]\n" + res[["celltype", "geneset", "n_genes", "delta_Z", "p", "BH_p"]]
        .round(4).to_string(index=False))

    # ---- 图 1: 四机制线 x 细胞类型 热图 (delta_Z) ----
    plan = ["A_complement_synapse_pruning", "B_DAM_microglia",
            "C_mitochondrial_OXPHOS", "D_myelin_oligodendrocyte"]
    sub = res[res.geneset.isin(plan)]
    piv = sub.pivot_table(index="geneset", columns="celltype", values="delta_Z")
    piv = piv.reindex([p for p in plan if p in piv.index])
    pv = sub.pivot_table(index="geneset", columns="celltype", values="BH_p").reindex(piv.index)
    fig, ax = plt.subplots(figsize=(1.25 * piv.shape[1] + 4, 3.0))
    sns.heatmap(piv, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                annot_kws={"size": 7}, ax=ax, cbar_kws={"label": "delta pseudobulk Z (Surg - Ctrl)"})
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if pv.iloc[i, j] < 0.05:
                ax.text(j + 0.5, i + 0.78, "*", ha="center", va="center", fontsize=12, color="black")
    ax.set_title("Step 5b: mechanism lines in scRNA pseudobulk (* BH p<0.05)")
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc5b_geneset_heatmap.png"), dpi=150)
    plt.close(fig)

    # ---- 图 2: 微胶质/少突 定向 panel 逐样本得分 ----
    focus_ct = [c for c in ["Microglia", "Oligodendrocyte", "OPC", "Astrocyte", "Neuron"] if c in celltypes]
    focus_set = ["DAM_core", "Homeostatic", "Phagocytic_rec", "Complement_core",
                 "Myelin_struct", "IEG_activity", "Mito_OXPHOS"]
    focus_set = [s for s in focus_set if s in res.geneset.unique()]
    n = len(focus_ct)
    fig, axes = plt.subplots(1, n, figsize=(3.1 * n, 3.6), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, ct in zip(axes, focus_ct):
        for k, sname in enumerate(focus_set):
            if (ct, sname) not in store:
                continue
            score, _, mix = store[(ct, sname)]
            for grp, color, off in [("Control", "#4c72b0", -0.16), ("Surgery", "#c44e52", 0.16)]:
                m = (mix.group == grp).values
                y = score.values[m]
                ax.scatter(np.full(len(y), k + off) + np.random.uniform(-0.05, 0.05, len(y)),
                           y, s=26, c=color, label=grp if k == 0 else None, zorder=3)
                ax.plot([k + off - 0.08, k + off + 0.08], [y.mean()] * 2, c="black", lw=1.6, zorder=4)
        ax.axhline(0, ls=":", c="grey", lw=0.8)
        ax.set_xticks(range(len(focus_set)))
        ax.set_xticklabels([s.replace("_", "\n") for s in focus_set], fontsize=6.5)
        ax.set_title(ct, fontsize=9)
    axes[0].set_ylabel("pseudobulk gene-set Z")
    axes[0].legend(fontsize=7, loc="lower left")
    fig.suptitle("Step 5b: per-sample pseudobulk gene-set scores (3v3)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc5b_persample_scores.png"), dpi=150)
    plt.close(fig)

    # ---- 摘要 ----
    with open(os.path.join(WORK, "step5b_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Step 5b — 单细胞逐样本 pseudobulk 基因集统计\n\n")
        f.write("方法: 逐 (细胞类型 x 样本) 聚合 lognorm 均值为 pseudobulk; 集合内基因跨样本 z 后取均值\n")
        f.write("      作为集合得分; 3v3 上 Welch t 检验 + BH。细胞不是重复单元, 故不做细胞级检验。\n\n")
        f.write("## 四条机制线 x 细胞类型 (delta Z = Surgery - Control)\n")
        f.write(sub[["celltype", "geneset", "n_genes", "mean_Z_ctrl", "mean_Z_surg",
                     "delta_Z", "p", "BH_p"]].round(4).to_string(index=False) + "\n\n")
        f.write("## 定向 panel\n")
        f.write(res[~res.geneset.isin(plan)][["celltype", "geneset", "n_genes", "delta_Z", "p", "BH_p"]]
                .round(4).to_string(index=False) + "\n\n")
        f.write("**局限**: n=3v3, 功效极低; 本步为方向一致性验证, 非确证性统计。\n")
    log("[Done] Step 5b 输出写入", WORK)


if __name__ == "__main__":
    main()
