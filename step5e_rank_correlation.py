"""
Step 5e (c) -- P2-A: limma-voom vs Welch 全基因排序相关系数
============================================================
动机 (审稿 m6 / P2-A):
  第 65 行构造敏感性只报了 Welch ΔZ vs count-sum 的 r=0.857/ρ=0.835，
  limma-voom 路线只写了 "the same"。需补 limma-voom(计数求和) 与
  Welch(log-mean 细胞均值) 两条路线在**基因排序**上的直接相关系数，
  作为"两路线定性一致"的构造对照。

输入:
  sc_out/sc_pseudobulkDEG_<ct>.csv  (Welch 主分析全基因: gene, log2FC, p)
  sc_out/pbcounts/voom_<ct>.csv     (limma-voom 全基因: gene, logFC, P.Value, ...)

计算:
  对每个细胞类型取两文件基因交集, 在 -log10(P) 上算
    - Spearman ρ (秩相关, 稳健)
    - Pearson  r (数值相关)
  另报 |log2FC| 的 Spearman 作稳健性。汇总逐细胞类型 + 全基因合并 + 逐类型均值。

输出:
  sc_out/step5e_rank_correlation.csv
  sc_out/fig_sc5e_rank_corr.png  (Microglia 散点 + 全类型 ρ 条形)
"""
import os, glob, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "sc_out")


def log(*a):
    print(*a, flush=True)


def main():
    welch_files = glob.glob(os.path.join(WORK, "sc_pseudobulkDEG_*.csv"))
    voom_files = glob.glob(os.path.join(WORK, "pbcounts", "voom_*.csv"))
    welch = {os.path.basename(f).replace("sc_pseudobulkDEG_", "").replace(".csv", ""): f
             for f in welch_files}
    voom = {os.path.basename(f).replace("voom_", "").replace(".csv", ""): f
            for f in voom_files}
    common = sorted(set(welch) & set(voom))
    log(f"[cell types] welch={len(welch)} voom={len(voom)} common={len(common)} -> {common}")

    rows = []
    pooled = []  # 全基因合并 (跨 cell type)
    per_ct_logp = {}
    for ct in common:
        w = pd.read_csv(welch[ct]).rename(columns={"gene": "gene", "log2FC": "w_log2FC", "p": "w_p"})
        v = pd.read_csv(voom[ct]).rename(columns={"Unnamed: 0": "gene", "logFC": "v_log2FC",
                                                  "P.Value": "v_p", "adj.P.Val": "v_adjp"})
        m = w.merge(v, on="gene")
        m = m[np.isfinite(m.w_p) & np.isfinite(m.v_p)]
        m = m[(m.w_p > 0) & (m.v_p > 0)]           # log10(0) 安全
        wl = -np.log10(m.w_p.values)
        vl = -np.log10(m.v_p.values)
        rho = float(stats.spearmanr(wl, vl).statistic)
        r = float(stats.pearsonr(wl, vl)[0])
        rho_fc = float(stats.spearmanr(m.w_log2FC.abs().values, m.v_log2FC.abs().values).statistic)
        rows.append({
            "celltype": ct, "n_intersect": len(m),
            "spearman_rho": round(rho, 4), "pearson_r": round(r, 4),
            "spearman_abslog2FC": round(rho_fc, 4),
            "welch_p_lt_0.05": int((m.w_p < 0.05).sum()),
            "voom_p_lt_0.05": int((m.v_p < 0.05).sum()),
            "welch_top_gene": m.loc[m.w_p.idxmin(), "gene"],
            "voom_top_gene": m.loc[m.v_p.idxmin(), "gene"],
        })
        pooled.append(m[["w_p", "v_p", "w_log2FC", "v_log2FC"]].copy())
        per_ct_logp[ct] = (wl, vl)
        log(f"  {ct:14s} n={len(m):6d}  rho={rho:.3f}  r={r:.3f}  rho|FC|={rho_fc:.3f}  "
            f"welch_p05={int((m.w_p<0.05).sum())} voom_p05={int((m.v_p<0.05).sum())}")

    rc = pd.DataFrame(rows)
    # 全基因合并统一相关
    allm = pd.concat(pooled, ignore_index=True)
    awl = -np.log10(allm.w_p.values)
    avl = -np.log10(allm.v_p.values)
    rho_all = float(stats.spearmanr(awl, avl).statistic)
    r_all = float(stats.pearsonr(awl, avl)[0])
    rho_all_fc = float(stats.spearmanr(allm.w_log2FC.abs().values, allm.v_log2FC.abs().values).statistic)

    mean_rho = float(rc.spearman_rho.mean())
    mean_r = float(rc.pearson_r.mean())
    mean_rho_fc = float(rc.spearman_abslog2FC.mean())

    summary = {
        "n_celltypes": len(rc),
        "mean_spearman_rho": round(mean_rho, 4),
        "mean_pearson_r": round(mean_r, 4),
        "mean_spearman_abslog2FC": round(mean_rho_fc, 4),
        "pooled_n_genes": len(allm),
        "pooled_spearman_rho": round(rho_all, 4),
        "pooled_pearson_r": round(r_all, 4),
        "pooled_spearman_abslog2FC": round(rho_all_fc, 4),
        "min_celltype_rho": round(float(rc.spearman_rho.min()), 4),
        "min_celltype_r": round(float(rc.pearson_r.min()), 4),
    }
    log("\n=== SUMMARY ===")
    for k, val in summary.items():
        log(f"  {k}: {val}")

    # 落盘
    rc.round(4).to_csv(os.path.join(WORK, "step5e_rank_correlation.csv"), index=False)
    with open(os.path.join(WORK, "step5e_rank_correlation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    log("\n[Done] -> sc_out/step5e_rank_correlation.csv + summary.json")

    # 散点图 (Microglia 代表) + 逐类型 ρ 条形
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.2))
    mref = pooled[common.index("Microglia")] if "Microglia" in common else allm
    ax = axes[0]
    ax.scatter(-np.log10(mref.v_p.values), -np.log10(mref.w_p.values), s=10,
               alpha=0.35, c="#4c72b0", edgecolors="none")
    ax.set_xlabel("limma-voom  −log10(P)")
    ax.set_ylabel("Welch  −log10(P)")
    ax.set_title(f"Gene-level −log10(P) agreement (Microglia, n={len(mref)})\n"
                 f"Spearman ρ={stats.spearmanr(-np.log10(mref.w_p.values),-np.log10(mref.v_p.values)).statistic:.2f}")
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], ls="--", c="grey", lw=0.8)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)

    ax = axes[1]
    order = rc.sort_values("spearman_rho")
    ax.barh(order.celltype, order.spearman_rho, color="#55a868")
    ax.axvline(1.0, ls=":", c="grey")
    ax.set_xlabel("Spearman ρ (limma-voom vs Welch, −log10 P)")
    ax.set_title(f"Rank correlation per cell type\nmean ρ={mean_rho:.2f} (min {rc.spearman_rho.min():.2f})")
    ax.set_xlim(0, 1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(WORK, "fig_sc5e_rank_corr.png"), dpi=150)
    plt.close(fig)
    log("[Done] fig -> sc_out/fig_sc5e_rank_corr.png")


if __name__ == "__main__":
    main()
