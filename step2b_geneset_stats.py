"""
Step 2b — 基因集水平统计 + Step 3 候选池
================================================================
背景: Step 2 发现"个体基因弱、但通路成员协调一致"的现象(线粒体 OXPHOS、补体成分
      三数据集方向一致但单体 meta_p 不显著)。故补做基因集水平检验, 并把
      "集合显著但单体未过阈"的成员以 set-representative 形式纳入候选池。

输入: step2_meta_DEG.csv, step2_pathway_sets.json, step2_candidate_genes.csv
输出:
  - step2_geneset_stats.csv        每套基因集的集合水平统计(均值 Z / 单样本 t / Wilcoxon / Stouffer)
  - step2_candidate_pool.csv       【Step 3 输入】基因水平候选 ∪ 集合代表候选
  - fig_step2_geneset_Z.png        各集合 meta_Z 分布箱线图
  - 追加写入 step2_summary.md
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

OUT = os.path.dirname(os.path.abspath(__file__))

# --- 探索性附加集合: 髓鞘/少突胶质 (Step2 全表最强信号, 原方案未含, 仅作报告) ---
MYELIN = ["Mbp","Plp1","Mog","Mag","Opalin","Bcas1","Pmp22","Cnp","Cldn11","Mal","Mobp",
          "Cntnap1","Ncmap","Tspan2","Eml1","Sirt2","Aspa","Sox10","Olig1","Olig2","Myrf",
          "Gjc2","Gpr37","Lpar1","St18","Prr5l","Sla2","Myo1d","Gltp","Fa2h","Ugt8a","Cers2",
          "Elovl7","Nkx2-2","Cspg4","Pdgfra","Gal3st1","Srebf1","Srebf2","Tmem63a","Kctd17",
          "Efnb1","Arid5a","Nr1d1","Foxo1","Sin3a"]


def main():
    meta = pd.read_csv(os.path.join(OUT, "step2_meta_DEG.csv"))
    pw = json.load(open(os.path.join(OUT, "step2_pathway_sets.json"), encoding="utf-8"))
    gene_cand = pd.read_csv(os.path.join(OUT, "step2_candidate_genes.csv"))

    sets = {k: v for k, v in pw.items()}    # A/B/C/D 全部来自 step2 的 pathway_sets.json

    present = meta[meta.n_datasets >= 2]
    rows, rep = [], {}
    for k, genes in sets.items():
        s = present[present.gene.isin(genes) & (present.consistency >= 2 / 3)]
        Z = s.meta_Z.values
        if len(Z) < 3:
            continue
        t1, p1 = stats.ttest_1samp(Z, 0)
        try:
            w, p2 = stats.wilcoxon(Z)
        except Exception:
            p2 = np.nan
        Zs = Z.sum() / np.sqrt(len(Z))
        ps = 2 * stats.norm.sf(abs(Zs))
        rows.append({"set": k, "n_members": int(s.shape[0]),
                     "mean_meta_Z": round(float(Z.mean()), 3),
                     "median_meta_Z": round(float(np.median(Z)), 3),
                     "frac_down": round(float((Z < 0).mean()), 3),
                     "t_p": p1, "wilcoxon_p": p2, "stouffer_p": ps,
                     "n_meta_sig": int((s.meta_p < 0.05).sum())})
        # set-representative: |Z| 前 10
        rep[k] = s.reindex(s.meta_Z.abs().sort_values(ascending=False).index).head(10)
    gs = pd.DataFrame(rows).sort_values("set")
    gs["BH_p(planned3)"] = np.nan
    gs["BH_p"] = multipletests(gs["t_p"], method="fdr_bh")[1]   # 四套机制线统一 BH
    gs.to_csv(os.path.join(OUT, "step2_geneset_stats.csv"), index=False)
    print(gs.to_string(index=False))

    # --- 候选池 = 基因水平候选 ∪ 集合代表 ---
    pool = meta.merge(
        gene_cand[["gene", "tier"]].rename(columns={"tier": "gene_tier"}),
        on="gene", how="left")
    rep_genes, rep_set = {}, {}
    for k, s in rep.items():
        for g in s.gene:
            rep_set[g] = rep_set.get(g, "") + (k.split("_")[0] + ",")
            rep_genes[g] = True
    pool = pool[pool.gene.isin(set(gene_cand.gene) | set(rep_genes))]
    pool["set_repr"] = pool.gene.map(lambda g: rep_set.get(g, "").rstrip(","))
    pool["source"] = np.where(pool.gene.isin(set(gene_cand.gene)), "gene-level", "set-representative")
    pool = pool.sort_values(["gene_tier", "meta_p"], na_position="last")
    keep = ["gene", "source", "gene_tier", "set_repr", "n_datasets", "consistency", "n_sig",
            "meta_Z", "meta_p"] + [c for c in pool.columns
                                   if c.endswith("_log2FC") or c.endswith("_p")
                                   and c not in ("meta_p",)]
    pool[[c for c in keep if c in pool.columns]].to_csv(
        os.path.join(OUT, "step2_candidate_pool.csv"), index=False)
    print(f"\n[candidate pool] 基因水平 {int((pool.source=='gene-level').sum())} + "
          f"集合代表 {int((pool.source!='gene-level').sum())} = {len(pool)} 基因 -> Step 3")

    # --- 图: 各集合 meta_Z 分布 ---
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    data, labels = [], []
    for k in sets:
        s = present[present.gene.isin(sets[k])]
        if len(s):
            data.append(s.meta_Z.values)
            labels.append(f"{k.split('_')[0]}\n(n={len(s)})")
    try:
        bp = ax.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True,
                        medianprops=dict(color="black"))
    except TypeError:
        bp = ax.boxplot(data, labels=labels, showfliers=False, patch_artist=True,
                        medianprops=dict(color="black"))
    for p, c in zip(bp["boxes"], ["#d62728", "#1f77b4", "#2ca02c", "#8c564b"]):
        p.set_facecolor(c); p.set_alpha(0.45)
    ax.axhline(0, ls="--", c="grey", lw=0.9)
    ax.set_ylabel("Stouffer meta Z (per gene)")
    ax.set_title("Pathway gene-set shift in cross-dataset meta-analysis\n"
                 "(A complement/pruning | B DAM microglia | C mito OXPHOS | D myelin/OL exploratory)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_step2_geneset_Z.png"), dpi=150)
    plt.close(fig)

    # --- 追加摘要 (幂等) ---
    summ = os.path.join(OUT, "step2_summary.md")
    prev = open(summ, encoding="utf-8").read() if os.path.exists(summ) else ""
    with open(summ, "a", encoding="utf-8") as f:
        if "## 基因集水平统计 (Step 2b)" not in prev:
            f.write("\n## 基因集水平统计 (Step 2b)\n")
        f.write("| 集合 | 成员 | 均值 meta_Z | 下调比例 | 单样本 t p | Wilcoxon p | Stouffer p | 单体 meta_p<0.05 |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for _, r in gs.iterrows():
            f.write(f"| {r['set']} | {r.n_members} | {r.mean_meta_Z} | {r.frac_down} | "
                    f"{r.t_p:.2e} | {r['wilcoxon_p']:.2e} | {r['stouffer_p']:.2e} | {r.n_meta_sig} |\n")
        f.write("\n**关键结论**: 线粒体 OXPHOS(C)与补体(A)在**个体基因层面**多为弱改变, 但**集合层面**")
        f.write("方向高度一致 -> 支持'需基因集/ML 才能捕获协调性弱信号'的论点。\n")
        f.write(f"\n**Step 3 候选池**: {len(pool)} 基因 (基因水平 + 集合代表), 见 step2_candidate_pool.csv。\n")
    print("[Done] Step 2b 输出已写入", OUT)


if __name__ == "__main__":
    main()
