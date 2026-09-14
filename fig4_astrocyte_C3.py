"""
Fig 4 (revised) — 星形胶质细胞型分辨: Set A 补体轴在 GSE199318 星形胶质内的方向
修订点 (审稿意见 P2): 原图未报 p/FDR。此处逐基因标注 p 与 BH, 并明确标注
n=3v3 下无一基因达到 BH<0.05 —— 方向与文献一致, 但仅属提示性。
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))
deg = pd.read_csv(os.path.join(OUT, "step2c_per_dataset_DEG.csv"))
a = deg[deg.dataset == "GSE199318_astro"].set_index("gene")
want = ["C3", "Cfb", "C4a", "C2", "C1ra"]
genes = [g for g in want if g in a.index]
a = a.loc[genes]
fc = a.log2FC.values
pv = a.p.values
bh = a.padj.values

fig, ax = plt.subplots(figsize=(7.0, 4.6))
colors = ["#d62728" if g == "C3" else ("#2ca02c" if v < 0 else "#888")
          for g, v in zip(genes, fc)]
bars = ax.bar(genes, fc, color=colors, width=0.6)
for i, (b, v) in enumerate(zip(bars, fc)):
    ax.text(b.get_x() + b.get_width() / 2, v + (0.12 if v >= 0 else -0.12),
            f"{v:+.2f}", ha="center", va="bottom" if v >= 0 else "top",
            fontsize=9, fontweight="bold" if abs(v) > 2 else "normal")
    ax.text(i, -5.55, f"p={pv[i]:.2f}\nBH={bh[i]:.2f}", ha="center", va="top",
            fontsize=7.4, color="#444")
ax.axhline(0, color="black", lw=0.8)
ax.set_ylim(-6.4, 4.2)
ax.set_ylabel("log2 fold change (laparotomy vs anesthesia control)")
ax.set_title("Astrocyte-resolved complement axis (GSE199318, sorted astrocytes)\n"
             "no gene reaches BH<0.05 at n=3 vs 3; C3 rises, classical cascade falls",
             fontsize=10.5)
ax.text(0.5, 0.975, "direction is consistent with published reactive-astrocyte C3 release, "
                    "but at this n it is suggestive only",
        transform=ax.transAxes, ha="center", va="top", fontsize=8, color="#555")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_ApC_astrocyte_C3.png"), dpi=300)
plt.close(fig)
print("[Done] fig_ApC_astrocyte_C3.png")
print("genes:", genes)
print("p:", [round(float(x), 4) for x in pv], "BH:", [round(float(x), 4) for x in bh])
