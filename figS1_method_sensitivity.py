"""
Supplementary Figure S1 — pseudobulk 检验选择敏感性 (审稿意见 P2) 与 Ttr 指纹
左: 每细胞类型 BH<0.05 基因数, Welch/lognorm vs limma-voom/TMM
右: Ttr 在各类细胞的 logFC vs 该类型 Ttr 表达水平 -> 非表达者变化更大 = ambient RNA 指纹
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))
SCO = os.path.join(OUT, "sc_out")
voom = pd.read_csv(os.path.join(SCO, "step5e_limma_voom_stats.csv"))
cal = pd.read_csv(os.path.join(SCO, "step5c_genomewide_calibration.csv"))
ORDER = ["Microglia", "Ependymal", "Astrocyte", "Oligodendrocyte", "OPC",
         "Neuron", "Endothelial", "VLMC", "Pericyte"]
v = voom.set_index("celltype").loc[ORDER]
c = cal.set_index("celltype").loc[ORDER]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 4.9),
                               gridspec_kw={"width_ratios": [1.35, 1]})

x = np.arange(len(ORDER)); w = 0.38
axL.bar(x - w / 2, c["padj_lt_0.05"], w, color="#b9c3cf", label="Welch t, lognorm pseudobulk")
axL.bar(x + w / 2, v["n_adjp_lt_0.05"], w, color="#4c72b0", label="limma-voom, TMM pseudobulk")
for i, ct in enumerate(ORDER):
    axL.text(i - w / 2, c.loc[ct, "padj_lt_0.05"] + 0.06, f"{int(c.loc[ct,'padj_lt_0.05'])}",
             ha="center", fontsize=7.5, color="#555")
    axL.text(i + w / 2, v.loc[ct, "n_adjp_lt_0.05"] + 0.06, f"{int(v.loc[ct,'n_adjp_lt_0.05'])}",
             ha="center", fontsize=7.5, color="#2b4c7e")
axL.set_xticks(x); axL.set_xticklabels(ORDER, rotation=40, ha="right", fontsize=8)
axL.set_ylabel("genes with BH-adjusted p<0.05")
axL.set_ylim(0, 5.6)
axL.set_title("Two pseudobulk routes agree qualitatively\n"
              "0 / 173,421 (Welch)  vs  11 / 72,204 (voom, 0.015%)", fontsize=10.5)
axL.legend(fontsize=8, loc="upper right", frameon=False)
axL.spines[["top", "right"]].set_visible(False)

# right: Ttr fingerprint
pts = []
for f in sorted(glob.glob(os.path.join(SCO, "pbcounts", "voom_*.csv"))):
    ct = os.path.basename(f)[5:-4]
    d = pd.read_csv(f, index_col=0)
    if "Ttr" in d.index:
        r = d.loc["Ttr"]
        pts.append((ct, float(r["AveExpr"]), float(r["logFC"]), float(r["adj.P.Val"])))
p = pd.DataFrame(pts, columns=["celltype", "AveExpr", "logFC", "adjP"])
epi = p[p.celltype == "Ependymal"]
oth = p[p.celltype != "Ependymal"]
axR.scatter(oth.AveExpr, oth.logFC, s=58, c="#c44e52", zorder=3,
            edgecolors="white", linewidths=0.7, label="non-expressing types (log2CPM 6.5–8.6)")
axR.scatter(epi.AveExpr, epi.logFC, s=95, c="#2ca02c", zorder=4, marker="D",
            edgecolors="white", linewidths=0.7, label="Ependymal (its true source, 16.3)")
axR.axhline(0, ls=":", c="grey", lw=0.9)
for _, r in p.iterrows():
    axR.annotate(r.celltype, (r.AveExpr, r.logFC), textcoords="offset points",
                 xytext=(6, 4), fontsize=7.4, color="#333")
axR.set_xlabel("Ttr expression in that cell type (mean log2 CPM)")
axR.set_ylabel("Ttr log2 fold change (Surgery vs Control)")
axR.set_title("Ttr: the largest 'signal' is an ambient-RNA artifact\n"
              "non-expressers shift further than the source cell type", fontsize=10.5)
axR.legend(fontsize=7.6, loc="lower right", frameon=False)
axR.spines[["top", "right"]].set_visible(False)

fig.suptitle("Supplementary Figure S1 — pseudobulk method sensitivity and the Ttr artifact",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(OUT, "fig_ApC_S1_method_sensitivity.png"), dpi=300)
plt.close(fig)
print("[Done] fig_ApC_S1_method_sensitivity.png")
print(p.to_string(index=False))
