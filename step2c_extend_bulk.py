"""
Step 2c — 扩展 bulk 元分析: 并入 GSE199318 (分选海马星形胶质细胞 RNA-seq)
========================================================================
背景: Step 5c 结论要求"先把 bulk 从 3 套扩到更多"以检验 4 机制线下调是否稳健。
      候选 6 个 GSE 经实测筛选 (见 step2c_候选筛选.md):
        - GSE95426 (海马 POCD 组织 6v6, 定制 Agilent 阵列): 平台无基因符号列,
          GPL22782.annot.gz 在沙箱 FTP 全路径 404 -> 无法并入基因级 meta -> 排除
        - GSE165798 (PND 海马, circRNA 芯片): 平台表仅 circRNA ID, 无 mRNA 符号 -> 排除
        - GSE303920 (海马, 标称 bulk): RAW 实为 10x scRNA (barcodes/features/matrix.mtx) -> 排除
        - GSE316433 (PBMC), GSE234493 (snRNA): 非海马 bulk -> 排除
        - **GSE199318 (分选海马星形胶质细胞, RNA-seq, 直接带 Gene Symbol+计数): 唯一干净可用**
      故本步把 GSE199318 作为"星形胶质细胞分辨"的第 4 个数据集并入, 同时保留 3 套
      组织层作基线对照, 重跑 Stouffer 元分析 + 四机制线基因集统计, 量化"加数据集是否改变结论"。

GSE199318 分组: SEV(sevoflurane 麻醉对照)=Control; LA(laparotomy 剖腹手术=PND)=Post。
"""
import os, re, json, tarfile, gzip, warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
OUT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(OUT, "raw_ext")

# 复用 step2 的 DEG / 元分析 / 通路集定义
from step2_bulk_integrate_pathway import (
    one_deg, meta_analysis, PATHWAY_DEF, NDUF_PREFIX,
    load_gse276942, load_gse215410, load_gse174412,
)
np.random.seed(42)


# ============================================================
# 载入 GSE199318 (分选海马星形胶质细胞, RNA-seq 计数)
# ============================================================
def load_gse199318():
    tf = tarfile.open(os.path.join(RAW, "GSE199318_RAW.tar"))
    cols = {}
    for m in tf.getmembers():
        if not m.name.endswith(".txt.gz"):
            continue
        base = os.path.basename(m.name).replace(".txt.gz", "")   # GSM..._SEV1 / _LA3
        mm = re.findall(r"(SEV\d|LA\d)", base)
        if not mm:
            continue
        tag = mm[0]                                              # SEV1 / LA3 (保留独立样本)
        f = tf.extractfile(m)
        df = pd.read_csv(f, sep="\t", compression="gzip",
                         names=["GeneID", "Symbol", "Count"])
        counts = df["Count"].apply(pd.to_numeric, errors="coerce")
        sym = df["Symbol"].astype(str).str.strip().str.strip("'\"")
        s = pd.Series(counts.values, index=sym.values)
        s = s[~s.index.isin(["nan", "-", "", "NA"])]
        s = s[s > 0]
        cols[tag] = s
    # 合并到共同符号, 缺失补 0
    mat = pd.DataFrame(cols)
    mat = mat.fillna(0.0)
    cpm = mat / mat.sum(axis=0) * 1e6
    expr = np.log2(cpm + 1.0)
    grp = pd.Series({t: ("Control" if t.startswith("SEV") else "Post") for t in cols})
    return expr, grp


def geneset_stats(meta, pathway_hit):
    present = meta[meta.n_datasets >= 2]
    rows = []
    for k, genes in pathway_hit.items():
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
                     "frac_down": round(float((Z < 0).mean()), 3),
                     "t_p": p1, "wilcoxon_p": p2, "stouffer_p": ps,
                     "n_meta_sig": int((s.meta_p < 0.05).sum())})
    gs = pd.DataFrame(rows).sort_values("set").reset_index(drop=True)
    gs["BH_p"] = multipletests(gs["t_p"], method="fdr_bh")[1]
    return gs


def main():
    # --- 载入 3 套组织层 ---
    e1, g1 = load_gse276942()
    e2, g2, ens2sym = load_gse215410()
    e3, g3 = load_gse174412(ens2sym)
    # --- 载入 GSE199318 星形胶质 ---
    e4, g4 = load_gse199318()
    print(f"[load] GSE276942{e1.shape} GSE215410{e2.shape} GSE174412{e3.shape} "
          f"GSE199318_astro{e4.shape}")
    print(f"[GSE199318] grp={g4.to_dict()}  n_post={int((g4=='Post').sum())} "
          f"n_ctrl={int((g4=='Control').sum())}")

    universe = set(e1.index) | set(e2.index) | set(e3.index) | set(e4.index)
    nduf = sorted(g for g in universe if g.startswith(NDUF_PREFIX))
    pathway_hit = {}
    for k, lst in PATHWAY_DEF.items():
        s = sorted(set(lst) & universe)
        if k.startswith("C_"):
            s = sorted(set(s) | set(nduf))
        pathway_hit[k] = s
    print("[pathway] " + ", ".join(f"{k}={len(v)}" for k, v in pathway_hit.items()))

    # --- 基线: 3 套组织层 ---
    deg3 = pd.concat([one_deg(e1, g1, "GSE276942"),
                      one_deg(e2, g2, "GSE215410"),
                      one_deg(e3, g3, "GSE174412")], ignore_index=True)
    meta3 = meta_analysis(deg3)
    gs3 = geneset_stats(meta3, pathway_hit)

    # --- 扩展: 3 套 + GSE199318 星形胶质 ---
    deg4 = pd.concat([deg3, one_deg(e4, g4, "GSE199318_astro")], ignore_index=True)
    meta4 = meta_analysis(deg4)
    gs4 = geneset_stats(meta4, pathway_hit)

    # --- 输出 ---
    meta4.to_csv(os.path.join(OUT, "step2c_meta_DEG.csv"), index=False)
    deg4.to_csv(os.path.join(OUT, "step2c_per_dataset_DEG.csv"), index=False)
    gs3.to_csv(os.path.join(OUT, "step2c_geneset_stats_3ds.csv"), index=False)
    gs4.to_csv(os.path.join(OUT, "step2c_geneset_stats_4ds.csv"), index=False)

    # --- 比较表 ---
    cmp = gs3[["set", "n_members", "mean_meta_Z", "frac_down", "t_p", "stouffer_p", "BH_p"]].merge(
        gs4[["set", "mean_meta_Z", "frac_down", "t_p", "stouffer_p", "BH_p"]],
        on="set", suffixes=("_3ds", "_4ds"))
    cmp["delta_mean_Z"] = (cmp["mean_meta_Z_4ds"] - cmp["mean_meta_Z_3ds"]).round(3)
    cmp.to_csv(os.path.join(OUT, "step2c_geneset_compare.csv"), index=False)
    print("\n==== 四机制线 基因集统计: 3 套 vs 4 套 ====")
    print(cmp.to_string(index=False))

    # --- 图: 4 套元分析 meta_Z 分布 (扩展前后) ---
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, meta, title in [(axes[0], meta3, "3 套组织层"),
                            (axes[1], meta4, "3 套 + GSE199318 星形胶质")]:
        present = meta[meta.n_datasets >= 2]
        data, labels = [], []
        for k in PATHWAY_DEF:
            s = present[present.gene.isin(pathway_hit[k])]
            if len(s):
                data.append(s.meta_Z.values); labels.append(f"{k.split('_')[0]}\n(n={len(s)})")
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
        ax.set_title(f"四机制线基因集位移\n({title})")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_step2c_geneset_Z.png"), dpi=150)
    plt.close(fig)

    # --- 摘要 ---
    with open(os.path.join(OUT, "step2c_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Step 2c 摘要 — 扩展 bulk 元分析 (并入 GSE199318 星形胶质)\n\n")
        f.write("## 候选筛选结果 (6 个 GSE 实测)\n")
        f.write("- GSE95426 海马 POCD 组织 6v6 (定制 Agilent 阵列): 平台无基因符号列, "
                "GPL22782.annot.gz 沙箱 FTP 全路径 404 -> **排除**\n")
        f.write("- GSE165798 PND 海马 (circRNA 芯片): 平台表仅 circRNA ID 无 mRNA 符号 -> **排除**\n")
        f.write("- GSE303920 海马 (标称 bulk): RAW 实为 10x scRNA -> **排除**\n")
        f.write("- GSE316433 PBMC / GSE234493 snRNA: 非海马 bulk -> **排除**\n")
        f.write("- **GSE199318 分选海马星形胶质 RNA-seq (Gene Symbol+计数): 唯一干净可用** -> 并入\n\n")
        f.write("## 纳入数据集\n")
        f.write(f"- 原 3 套组织层: GSE276942, GSE215410, GSE174412\n")
        f.write(f"- 新增: GSE199318 星形胶质 (SEV=麻醉对照 n=3, LA=剖腹手术 PND n=3)\n\n")
        f.write("## 四机制线基因集统计 (扩展前后)\n")
        f.write(cmp.to_string(index=False) + "\n\n")
        f.write("## 关键判读\n")
        f.write("- 见 step2c_geneset_compare.csv 与正文。\n")
    print("[Done] Step 2c 输出已写入", OUT)


if __name__ == "__main__":
    main()
