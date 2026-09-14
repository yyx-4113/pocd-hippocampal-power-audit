"""
Step 1a — GSE276942 (POCD 时序 bulk RNA-seq) 时序 DEG + Mfuzz 式软聚类轨迹分析
输入: GSE276942_gene_expression_filter.csv.gz  (processed/normalized expression, 16885 genes x 17 samples)
输出:
  - GSE276942_DEG_by_timepoint.csv        分时间点 vs Control 的差异基因 (log2FC, p, BH-padj)
  - GSE276942_trajectory_clusters.csv     每个候选基因 -> 软聚类簇 + 成员度
  - fig_trajectory_centroids.png          各簇标准化轨迹中心线
  - fig_trajectory_heatmap.png            候选基因轨迹热图 (按簇排序)
  - step1_summary.md                       文本摘要
说明:
  - 表达值为 GEO processed/normalized 矩阵(线性尺度, 疑似 TPM/FPKM 类), 非 raw counts,
    故 DEG 为探索性; 晚时间点(1M-4M)各 n=1, 仅入轨迹, 不做统计推断.
  - 16885 基因经 BH 多重校正后极难有 padj<0.05 者(最显著 p~1e-4 量级), 这与
    "bulk 海马 POCD 为细微改变"一致; 下游用探索性候选集(p<0.05 且 |log2FC|>0.5, 未校正)继续,
    并在结果中明确标注为探索性.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

OUT = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(OUT, "GSE276942_gene_expression_filter.csv.gz")
GENE_COL = "gene_name"

# ---------- 1. 载入 ----------
df = pd.read_csv(DATA)
expr = df.drop(columns=[GENE_COL]).copy()
expr = expr.groupby(df[GENE_COL]).mean()          # 同名基因取均值(Entrez ID, 理论唯一)
genes = expr.index.values
X = expr.values.astype(float)                     # genes x 17
sample_titles = list(expr.columns)
assert X.shape[1] == 17, f"期望 17 个样本列, 实际 {X.shape[1]}"

# ---------- 2. 样本 -> 分组/时间点 映射 ----------
meta = {
    "CON_1": ("Control", "CON"), "CON_2": ("Control", "CON"),
    "CON_4": ("Control", "CON"), "CON_5": ("Control", "CON"),
    "24H_2": ("24H", "24H"), "24H_3": ("24H", "24H"), "24H_4": ("24H", "24H"),
    "3D_1": ("3D", "3D"), "3D_3": ("3D", "3D"), "3D_4": ("3D", "3D"),
    "7D_1": ("7D", "7D"), "7D_2": ("7D", "7D"), "7D_3": ("7D", "7D"),
    "1M_1": ("1M", "1M"), "1M_2": ("2M", "2M"),
    "1M_3": ("3M", "3M"), "1M_4": ("4M", "4M"),
}
missing = [c for c in sample_titles if c not in meta]
assert not missing, f"未映射列: {missing}"
groups = np.array([meta[c][0] for c in sample_titles])
TIMEORDER = ["CON", "24H", "3D", "7D", "1M", "2M", "3M", "4M"]
ctrl_mask = groups == "Control"
post_times = ["24H", "3D", "7D", "1M", "2M", "3M", "4M"]

Xlog = np.log2(X + 1)   # 对数变换, 稳健

# ---------- 3. 分时间点 DEG (vs Control) ----------
deg_rows = []
for tp in post_times:
    mask = groups == tp
    a = Xlog[:, ctrl_mask]      # control
    b = Xlog[:, mask]           # post
    log2fc = b.mean(1) - a.mean(1)
    if mask.sum() >= 2:
        t, p = stats.ttest_ind(b, a, axis=1, equal_var=False)   # Welch; 仅 n>=2 做检验
        p = np.where(np.isfinite(p), np.clip(p, 0.0, 1.0), 1.0)
    else:
        p = np.full(X.shape[0], np.nan)   # 晚时间点 n=1, 不做推断, 仅供轨迹
    deg_rows.append(pd.DataFrame({"gene": genes, "timepoint": tp,
                                   "log2FC": log2fc, "p": p}))
deg = pd.concat(deg_rows, ignore_index=True)


def _fdr(v):
    if v.notna().sum() == 0:
        return v
    ok = v.dropna()
    out = v.copy()
    out.loc[ok.index] = false_discovery_control(ok.values, method="bh")
    return out


deg["padj"] = deg.groupby("timepoint")["p"].transform(_fdr)
deg.to_csv(os.path.join(OUT, "GSE276942_DEG_by_timepoint.csv"), index=False)

bh_pass = deg[(deg.padj < 0.05) & (deg.log2FC.abs() > 0.5)]
bh_n = len(set(bh_pass["gene"]))
print(f"[DEG] 总表行数={len(deg)}")
print(f"[DEG] BH-padj<0.05 显著基因(并集)={bh_n}  "
      f"-> bulk POCD 信号经多重校正后极弱, 与'细微改变'一致")
print("[DEG] 下游采用探索性候选集: p<0.05 且 |log2FC|>0.5 (未校正)")

cand_mask = (deg.p < 0.05) & (deg.log2FC.abs() > 0.5) & deg.p.notna()
cand_genes = sorted(set(deg[cand_mask]["gene"]))
print(f"[DEG] 探索性候选基因(并集)={len(cand_genes)}")
for tp in post_times:
    sub = deg[cand_mask & (deg.timepoint == tp)]
    print(f"    {tp}: up={int((sub.log2FC > 0).sum())}, down={int((sub.log2FC < 0).sum())}")

de_genes = cand_genes
if len(de_genes) < 50:   # 兜底: 跨时间点高变异基因
    prof0 = np.column_stack([Xlog[:, groups == tp].mean(1) for tp in TIMEORDER])
    tv = np.nanstd(prof0, axis=1)
    de_genes = sorted(genes[np.argsort(tv)[-500:]])
    print(f"[DEG] 候选过少, 退化用跨时间点 top500 高变异基因, n={len(de_genes)}")

# ---------- 4. 轨迹软聚类 (Mfuzz 式 模糊 c-means) ----------
gi = [list(genes).index(g) for g in de_genes]
prof = np.column_stack([Xlog[gi, :][:, groups == tp].mean(1) for tp in TIMEORDER])
prof = np.nan_to_num(prof, nan=0.0)
row_sd = prof.std(1)
keep = row_sd > 0
if keep.sum() == 0:
    keep = np.ones(len(row_sd), dtype=bool)
prof_std = np.zeros_like(prof)
prof_std[keep] = (prof[keep] - prof[keep].mean(1, keepdims=True)) / (prof[keep].std(1, keepdims=True) + 1e-9)
prof_std[~keep] = 0.0
de_genes = [de_genes[i] for i in range(len(de_genes)) if keep[i]]
prof_std = prof_std[keep]


def fcm(X, k, m=1.3, max_iter=300, seed=42):
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    C = X[rng.choice(n, k, replace=False)].copy()
    for _ in range(max_iter):
        d = np.linalg.norm(X[:, None, :] - C[None, :, :], axis=2) ** 2
        d = np.fmax(d, 1e-12)
        U = (1.0 / d) ** (1.0 / (m - 1))
        U /= U.sum(1, keepdims=True)
        Um = U ** m
        Cnew = (Um.T @ X) / np.fmax(Um.sum(0, keepdims=True).T, 1e-12)
        if np.linalg.norm(Cnew - C) < 1e-4:
            C = Cnew
            break
        C = Cnew
    return U, C


K = 5
U, C = fcm(prof_std, K, m=1.3)
clusters = U.argmax(1)
mem = U.max(1)

pd.DataFrame({"gene": de_genes, "cluster": clusters,
              "max_membership": np.round(mem, 3)}).to_csv(
    os.path.join(OUT, "GSE276942_trajectory_clusters.csv"), index=False)
print(f"[Cluster] K={K}, 各簇: " +
      ", ".join(f"C{c+1}={(clusters == c).sum()}" for c in range(K)))

# ---------- 5. 图形 ----------
sns.set_style("whitegrid")
COLORS = sns.color_palette("tab10", K)

fig, ax = plt.subplots(figsize=(8, 5))
for c in range(K):
    ax.plot(range(len(TIMEORDER)), C[c], marker="o", color=COLORS[c],
            label=f"Cluster {c+1} (n={(clusters == c).sum()})", lw=2)
ax.axhline(0, color="grey", ls="--", lw=0.8)
ax.set_xticks(range(len(TIMEORDER)))
ax.set_xticklabels(TIMEORDER)
ax.set_xlabel("Post-surgery timepoint")
ax.set_ylabel("Standardized expression (z)")
ax.set_title("GSE276942 time-trajectory clusters (Mfuzz-style FCM, K=5)")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_trajectory_centroids.png"), dpi=150)
plt.close(fig)

order = []
for c in range(K):
    idx = np.where(clusters == c)[0]
    corr = np.array([np.corrcoef(prof_std[i], C[c])[0, 1] for i in idx])
    order.extend(idx[np.argsort(-corr)].tolist())
order = np.array(order)
fig, ax = plt.subplots(figsize=(7, 11))
sns.heatmap(prof_std[order], cmap="RdBu_r", center=0, vmin=-2, vmax=2,
            xticklabels=TIMEORDER, yticklabels=False, ax=ax,
            cbar_kws={"shrink": 0.5})
y = 0
for c in range(K):
    y += (clusters == c).sum()
    if c < K - 1:
        ax.axhline(len(order) - y, color="black", lw=0.8)
ax.set_title(f"Candidate gene trajectories (n={len(de_genes)}, K={K})")
ax.set_xlabel("Post-surgery timepoint")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_trajectory_heatmap.png"), dpi=150)
plt.close(fig)

# ---------- 6. 摘要 ----------
peak_tp = [TIMEORDER[int(np.argmax(np.abs(C[c])))] for c in range(K)]
with open(os.path.join(OUT, "step1_summary.md"), "w", encoding="utf-8") as f:
    f.write("# Step 1a 摘要 — GSE276942 时序 DEG + 轨迹聚类\n\n")
    f.write(f"- 基因数: {X.shape[0]}; 样本数: {X.shape[1]} (8 时间点)\n")
    f.write(f"- BH-padj<0.05 显著基因(并集): **{bh_n}** (经 16885 基因多重校正后极弱)\n")
    f.write(f"- 探索性候选基因(p<0.05 & |log2FC|>0.5, 未校正, 并集): **{len(de_genes)}**\n")
    for tp in post_times:
        sub = deg[cand_mask & (deg.timepoint == tp)]
        f.write(f"  - {tp}: up={int((sub.log2FC > 0).sum())}, "
                f"down={int((sub.log2FC < 0).sum())}\n")
    f.write(f"- 轨迹聚类 K={K}, 各簇: " +
            ", ".join(f"C{c+1}={(clusters == c).sum()}" for c in range(K)) + "\n")
    f.write("\n## 各簇代表基因 (top 10 by |z| peak)\n")
    for c in range(K):
        idx = np.where(clusters == c)[0]
        peak = np.abs(prof_std[idx]).max(1)
        top = idx[np.argsort(-peak)[:10]]
        f.write(f"\n### Cluster {c+1} (n={(clusters == c).sum()}, 峰值时间={peak_tp[c]})\n")
        f.write(", ".join(str(de_genes[i]) for i in top) + "\n")

print("[Done] 输出已写入", OUT)
