"""
Fig 1 — 路线 A+C 方法学稿中心图: 组成 vs 状态 功效不对称
左: 细胞组成 (9 类, Control vs Surgery %)，BH 显著者标 *
右: 状态层转录组 (9 类各 19269 基因) FDR 显著基因数 = 0, 对照随机期望线
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))
comp = pd.read_csv(os.path.join(OUT, "sc_out", "sc_celltype_composition.csv"))
cal = pd.read_csv(os.path.join(OUT, "sc_out", "step5c_genomewide_calibration.csv"))

ORDER = ["Microglia", "Ependymal", "Astrocyte", "Endothelial", "Oligodendrocyte",
         "OPC", "Neuron", "VLMC", "Pericyte"]
comp = comp.set_index("celltype").loc[ORDER].reset_index()
sig = dict(zip(comp.celltype, comp.BH_p < 0.05))

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5.2))

# ---- 左: 组成 ----
y = np.arange(len(ORDER))
axL.barh(y - 0.2, comp.pct_ctrl, height=0.4, color="#9ecae1", label="Control")
axL.barh(y + 0.2, comp.pct_surg, height=0.4, color="#fb6a4a", label="Surgery (POCD)")
for i, r in comp.iterrows():
    if sig[r.celltype]:
        axL.text(max(r.pct_ctrl, r.pct_surg) + 0.6, i, "*", va="center",
                 fontsize=14, color="black", fontweight="bold")
axL.set_yticks(y); axL.set_yticklabels(ORDER, fontsize=9)
axL.set_xlabel("% of cells"); axL.set_title(
    "Composition is powered\n(cell-type proportions; * BH<0.05)", fontsize=11)
axL.legend(fontsize=8, loc="lower right")
axL.spines[["top", "right"]].set_visible(False)

# ---- 右: 状态层 FDR 显著基因数 ----
cal = cal.set_index("celltype").loc[ORDER].reset_index()
x = np.arange(len(ORDER))
axR.bar(x, cal["padj_lt_0.05"], color="#bdbdbd", width=0.6)
axR.axhline(963.5, ls="--", color="#d62728", lw=1.3,
            label="random expectation (p<0.05)")
axR.set_xticks(x); axR.set_xticklabels(ORDER, rotation=40, ha="right", fontsize=8)
axR.set_ylabel("FDR-significant genes (padj<0.05)")
axR.set_ylim(0, 1100)
axR.set_title("State is NOT powered\n(0 / 173,421 tests; test conservative)",
              fontsize=11)
axR.text(0.5, 0.92, "MDE d\u22483.1 required (n=3v3, 80% power)",
        transform=axR.transAxes, ha="center", fontsize=9, color="#d62728")
axR.legend(fontsize=8, loc="upper right")
axR.spines[["top", "right"]].set_visible(False)

fig.suptitle("Power asymmetry in a 3-vs-3 hippocampal scRNA-seq design (POCD)",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(OUT, "fig_ApC_power_asymmetry.png"), dpi=160)
plt.close(fig)
print("[Done] fig_ApC_power_asymmetry.png")
