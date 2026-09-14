"""
Fig 2 — 路线 A+C 方法学稿: 最小可检出效应 (MDE) 随样本量曲线
x = 每组 n; y = MDE d (80% power, alpha=0.05, Welch)
叠加: bulk 实测效应量级带 (0.3-1.5), n=3 (现状, 红虚) vs n=6 (可检, 绿虚)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(OUT, "sc_out", "step5c_mde_by_samplesize.csv"))
n = df["每组n"].values
d = df["可检出标准化效应 d (80% power)"].values

fig, ax = plt.subplots(figsize=(7.2, 4.8))
ax.plot(n, d, "o-", color="#1f77b4", lw=2, label="MDE d (required to detect)")
# bulk 实测效应量级带
ax.axhspan(0.3, 1.5, color="#2ca02c", alpha=0.12, label="observed bulk effect size (0.3-1.5)")
ax.axhline(0.8, ls=":", color="#2ca02c", lw=1, label="Cohen 'large' d=0.8")
# 现状 vs 可检
ax.axvline(3, ls="--", color="#d62728", lw=1.5, label="n=3 (current 3v3): not powered")
ax.axvline(6, ls="--", color="#9467bd", lw=1.5, label="n=6: enters detectable zone")
# 标注关键 n 的 MDE
for nn, dd in [(3, 3.07), (6, 1.80), (8, 1.51)]:
    ax.annotate(f"n={nn}\nd={dd:.2f}", (nn, dd), textcoords="offset points",
                xytext=(6, 8), fontsize=8, color="#333")
ax.set_xlabel("samples per group (n)"); ax.set_ylabel("minimum detectable effect (Cohen d)")
ax.set_title("Statistical power of scRNA/pseudobulk POCD designs\n"
             "small n cannot resolve the subtle (d<1) effects POCD produces", fontsize=11)
ax.set_ylim(0, 6); ax.set_xlim(1.5, 21)
ax.legend(fontsize=7.5, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_ApC_mde_curve.png"), dpi=160)
plt.close(fig)
print("[Done] fig_ApC_mde_curve.png")
