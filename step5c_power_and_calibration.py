"""
Step 5c — 单细胞层"无显著"的校准与功效分析
============================================
要回答的问题:
  step5b/5b2 显示四条机制线在 9 种细胞类型上 BH 后无一显著。这究竟是
  (a) 集合真的没变化, 还是 (b) 3v3 设计功效不足、什么都测不出来?
判定方法(两条独立证据):
  1) 全基因组校准: 对每个细胞类型的 pseudobulk DEG 做 BH。若"全基因组范围内"
     也几乎无显著基因(接近零), 说明是该设计的功效极限, 而非本集合特异地没变;
     若全基因组有数百显著基因而机制线全无, 则说明集合确实没变 —— 这是强结论。
     同时看 p 值分布是否偏离均匀(偏离=有真信号但被样本量卡住)。
  2) 最小可检测效应(MDE): 用观测到的逐样本得分 SD, 反算 80% 功效下 n=3v3 能
     检出的最小 delta_Z, 并与 bulk 层的实际效应量对比 —— 直接量化"够不够用"。

输入: sc_out/sc_pseudobulkDEG_*.csv, sc_out/step5b_pseudobulk_*.csv, sc_out/step5b_geneset_stats.csv
输出: sc_out/step5c_*.csv, fig_sc5c_*.png, step5c_summary.md
"""
import os, glob, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "sc_out")


def log(*a):
    print(*a, flush=True)


def mde_n3(alpha=0.05, power=0.80):
    """两组 n=3, Welch t(近似 df=3.4) 下的标准化最小可检测效应 (delta/SD),
    再乘 SD 得原始尺度 MDE。用非中心 t 的近似解。"""
    from scipy.stats import nct, t as tdist
    df = 3.4                                   # Welch df for 3 vs 3 with equal var
    crit = tdist.ppf(1 - alpha / 2, df)
    lo, hi = 0, 60.0
    for _ in range(200):                       # 二分求非中心参数
        mid = (lo + hi) / 2
        p = 1 - nct.cdf(crit, df, mid) + nct.cdf(-crit, df, mid)
        if p < power:
            lo = mid
        else:
            hi = mid
    ncp = (lo + hi) / 2
    return ncp / np.sqrt(1.5)                  # d = ncp / sqrt(n/2), n=3


def main():
    # ---------- 1) 全基因组校准 ----------
    rows, pvals = [], {}
    for fp in sorted(glob.glob(os.path.join(WORK, "sc_pseudobulkDEG_*.csv"))):
        ct = os.path.basename(fp).replace("sc_pseudobulkDEG_", "").replace(".csv", "")
        d = pd.read_csv(fp)
        d["p"] = pd.to_numeric(d["p"], errors="coerce")
        d = d[np.isfinite(d["p"])]
        if len(d) < 100:
            continue
        bh = multipletests(d["p"].values, method="fdr_bh")[1]
        d["padj"] = bh
        pvals[ct] = d["p"].values
        # p 值分布偏离均匀的检验 (仅取 p<0.5 以外的整体 KS 检验对 U(0,1))
        ks = stats.kstest(d["p"].values, "uniform")
        rows.append({
            "celltype": ct, "n_genes": len(d),
            "padj_lt_0.05": int((d.padj < 0.05).sum()),
            "padj_lt_0.10": int((d.padj < 0.10).sum()),
            "p_lt_0.05": int((d.p < 0.05).sum()),
            "p_lt_0.01": int((d.p < 0.01).sum()),
            "exp_p_lt_0.05": round(0.05 * len(d), 1),
            "KS_D": round(ks.statistic, 4), "KS_p": ks.pvalue,
            "min_p": float(d.p.min()),
            "top_gene": d.loc[d.p.idxmin(), "gene"] if len(d) else "",
            "top_log2FC": float(d.loc[d.p.idxmin(), "log2FC"]) if len(d) else np.nan,
        })
    cal = pd.DataFrame(rows).sort_values("padj_lt_0.05", ascending=False)
    cal.round(5).to_csv(os.path.join(WORK, "step5c_genomewide_calibration.csv"), index=False)
    log("[全基因组校准] 每个细胞类型 pseudobulk DEG 的 BH 结果:")
    log(cal.to_string(index=False))

    # ---------- 2) 最小可检测效应 vs bulk 实际效应 ----------
    d_star = mde_n3()
    log(f"\n[MDE] n=3v3, alpha=0.05, power=0.80 -> 标准化 MDE d = {d_star:.2f} "
        f"(即 delta/SD 需达到 {d_star:.2f} 才能被检出)")
    log("[对照] bulk 层 Stouffer 均值 Z 的量级: D -1.48 / B -0.91 / C -0.50 / A -0.29")

    gs = pd.read_csv(os.path.join(WORK, "step5b_geneset_stats.csv"))
    plan = ["A_complement_synapse_pruning", "B_DAM_microglia",
            "C_mitochondrial_OXPHOS", "D_myelin_oligodendrocyte"]
    out = []
    for ct in gs.celltype.unique():
        fp = os.path.join(WORK, f"step5b_pseudobulk_{ct}.csv")
        if not os.path.exists(fp):
            # 该文件存的是基因 x 样本 的 pseudobulk, 需按基因集重算 SD
            continue
        pb = pd.read_csv(fp, index_col=0)
        for sname in plan:
            sub = gs[(gs.celltype == ct) & (gs.geneset == sname)]
            if sub.empty:
                continue
            # 用同一 z-均值 定义重算逐样本得分 SD
            pw = json.load(open(os.path.join(HERE, "step2_pathway_sets.json"), encoding="utf-8"))
            gl = [g for g in pw[sname] if g in pb.columns]
            if len(gl) < 3:
                continue
            s = pb[gl]
            z = (s - s.mean(0)) / s.std(0, ddof=1).replace(0, np.nan)
            sc_ = z.fillna(0.0).mean(1)             # 逐样本得分
            sd = float(sc_.std(ddof=1))
            obs = float(sub.delta_Z.iloc[0])
            out.append({"celltype": ct, "geneset": sname, "sd_score": round(sd, 3),
                        "observed_delta_Z": round(obs, 3),
                        "MDE_delta_Z(80%)": round(d_star * sd, 3),
                        "observed_over_MDE": round(abs(obs) / (d_star * sd), 3),
                        "p": float(sub.p.iloc[0])})
    mde = pd.DataFrame(out)
    mde.round(4).to_csv(os.path.join(WORK, "step5c_mde.csv"), index=False)
    log("\n[MDE vs 观测] (observed/MDE < 1 表示设计无法检出该量级):")
    log(mde.sort_values("observed_over_MDE", ascending=False).head(20).to_string(index=False))
    log(f"\n[汇总] 32 个 (集合 x 细胞类型) 组合中, 观测效应超过 MDE 的有 "
        f"{int((mde.observed_over_MDE >= 1).sum())} 个")

    # ---------- 图 ----------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    x = np.arange(len(pvals)); 
    names = list(pvals)
    for i, (ct, pv) in enumerate(pvals.items()):
        ax.hist(pv, bins=20, range=(0, 1), histtype="step", lw=1.2,
                label=f"{ct} (padj<0.05: {int((multipletests(pv, method='fdr_bh')[1]<0.05).sum())})")
    ax.axhline(len(pvals) and 0.05 * np.mean([len(v) for v in pvals.values()]) / 20,
               ls="--", c="black", lw=1, label="期望(均匀)")
    ax.set_xlabel("pseudobulk DEG p value"); ax.set_ylabel("基因数 / bin")
    ax.set_title("p-value distribution per cell type (uniform => no detectable signal)")
    ax.legend(fontsize=6)

    ax = axes[1]
    if len(mde):
        for sname, color in zip(plan, ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]):
            sub = mde[mde.geneset == sname]
            ax.scatter(sub["MDE_delta_Z(80%)"], sub["observed_delta_Z"], s=42,
                       c=color, label=sname.split("_")[0], edgecolors="k", linewidths=0.4)
        lim = max(mde["MDE_delta_Z(80%)"].max(), mde["observed_delta_Z"].abs().max()) * 1.1
        ax.plot([-lim, lim], [-lim, lim], ls="--", c="grey", lw=1)
        ax.axhline(0, ls=":", c="grey", lw=0.8); ax.axvline(0, ls=":", c="grey", lw=0.8)
        ax.set_xlim(0, lim)
        ax.set_xlabel("MDE delta_Z at 80% power (n=3v3)")
        ax.set_ylabel("observed delta_Z")
        ax.set_title("Observed effect vs detectable effect (all below the diagonal = underpowered)")
        ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc5c_calibration.png"), dpi=150)
    plt.close(fig)

    with open(os.path.join(WORK, "step5c_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Step 5c — 单细胞层\"无显著\"的校准与功效分析\n\n")
        f.write(f"- n=3v3, alpha=0.05, power=0.80 的标准化最小可检测效应 d = **{d_star:.2f}**\n")
        f.write("- 判定逻辑: 若全基因组 BH 显著基因≈0 且 p 值分布均匀 -> 是**功效极限**;\n")
        f.write("  若全基因组有大量显著基因而机制线全无 -> 是**集合确实没变**(强结论)。\n\n")
        f.write("## 全基因组校准\n" + cal.to_string(index=False) + "\n\n")
        f.write("## MDE vs 观测\n" + mde.round(3).to_string(index=False) + "\n")
    log("[Done] Step 5c 输出写入", WORK)


if __name__ == "__main__":
    main()
