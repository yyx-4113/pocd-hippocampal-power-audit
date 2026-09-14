"""
Step 5f — 修订版 Fig 1 / Fig 2 + required-n 表 (响应审稿意见 P1-1 / P1-2 / P1-3)
=================================================================================
P1-1: 组成检验的推断单元必须是动物, 不是细胞。Fig 1 左改为"逐动物比例 strip plot +
      动物级 Welch BH", 显示九类 BH 后无一显著。
P1-2: 原稿"n>=6 进入可检带"与自家 MDE 表冲突 (n=6 的 MDE d=1.80 > 效应带上界 1.5)。
      Fig 2 改用"效应量 -> 所需每组 n"真实曲线, 标 d=1.5 -> n≈9, d=1.0 -> n≈17。
P1-3: 3.1 vs 0.3-1.5 是 2-10 倍, 非"数量级"。
输出: fig_ApC_power_asymmetry.png, fig_ApC_mde_curve.png, step5f_required_n.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import nct
from scipy.optimize import brentq

OUT = os.path.dirname(os.path.abspath(__file__))
SCO = os.path.join(OUT, "sc_out")

# ---------- 0. 精确功率函数 (与稿件 MDE 方法一致: 非中心 t, df=2n-2) ----------
ALPHA, POWER = 0.05, 0.80


def power_at(n, d):
    """向量化非中心 t 功效; n 为标量或数组, d 为标量或数组 (广播)。"""
    n = np.asarray(n, dtype=float)
    d = np.asarray(d, dtype=float)
    df = 2 * n - 2
    tcrit = nct.ppf(1 - ALPHA / 2, df, 0)
    nc = d * np.sqrt(n / 2.0)
    return nct.sf(tcrit, df, nc) + nct.cdf(-tcrit, df, nc)


def mde(n):
    ds = np.linspace(0.01, 8.0, 1600)
    ps = power_at(n, ds)
    ok = np.where(ps >= POWER)[0]
    return float(ds[ok[0]]) if len(ok) else np.nan


def req_n(d):
    ns = np.arange(2, 601)
    ps = power_at(ns, d)
    ok = np.where(ps >= POWER)[0]
    return int(ns[ok[0]]) if len(ok) else np.nan


# sanity: 复现稿件 MDE 表
print("[sanity] MDE n=3,6,8,20:", [round(mde(n), 3) for n in (3, 6, 8, 20)])

rows = []
for d in [1.5, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2]:
    rows.append({"effect_size_d": d, "required_n_per_group": req_n(d),
                 "total_n": 2 * req_n(d)})
req = pd.DataFrame(rows)
req.to_csv(os.path.join(SCO, "step5f_required_n.csv"), index=False)
print(req.to_string(index=False))

# ---------- 1. Fig 1 (revised) ----------
comp = pd.read_csv(os.path.join(SCO, "step5d_composition_animallevel.csv"))
cal = pd.read_csv(os.path.join(SCO, "step5c_genomewide_calibration.csv"))
ORDER = ["Microglia", "Ependymal", "Astrocyte", "Oligodendrocyte", "OPC",
         "Neuron", "Endothelial", "VLMC", "Pericyte"]
comp = comp.set_index("celltype").loc[ORDER].reset_index()
cal = cal.set_index("celltype").loc[ORDER].reset_index()

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.2, 5.6),
                               gridspec_kw={"width_ratios": [1.25, 1]})

# left: per-animal strip plot (animal = unit of inference)
y = np.arange(len(ORDER))
rng = np.random.default_rng(0)
for i, r in comp.iterrows():
    ctrl = np.array([float(x) for x in str(r.per_animal_ctrl).split(";")])
    surg = np.array([float(x) for x in str(r.per_animal_surg).split(";")])
    axL.scatter(ctrl, np.full(len(ctrl), i - 0.18) + rng.uniform(-0.05, 0.05, len(ctrl)),
                s=46, c="#3b6fb6", zorder=3, edgecolors="white", linewidths=0.6)
    axL.scatter(surg, np.full(len(surg), i + 0.18) + rng.uniform(-0.05, 0.05, len(surg)),
                s=46, c="#d1495b", zorder=3, edgecolors="white", linewidths=0.6)
    axL.plot([ctrl.mean()] * 2, [i - 0.28, i - 0.08], c="#1f3f66", lw=2.2, zorder=4)
    axL.plot([surg.mean()] * 2, [i + 0.08, i + 0.28], c="#8c2331", lw=2.2, zorder=4)
    bh = float(r.BH_Welch_animal)
    axL.text(1.015, i, "ns" if bh >= 0.05 else f"BH={bh:.1e}", va="center",
             fontsize=8, color="#333", transform=axL.get_yaxis_transform())
axL.set_yticks(y); axL.set_yticklabels(ORDER, fontsize=9)
axL.invert_yaxis()
axL.set_xlabel("% of hippocampal cells (per animal)")
axL.set_xlim(0, 56)
axL.set_title("Composition, animal-level test\n"
              "0 of 9 cell types reach BH p<0.05 (n=3 vs 3)", fontsize=11)
axL.scatter([], [], s=46, c="#3b6fb6", label="Control (n=3)")
axL.scatter([], [], s=46, c="#d1495b", label="Surgery (n=3)")
axL.legend(fontsize=8, loc="lower right", frameon=False)
axL.spines[["top", "right"]].set_visible(False)

# right: state level - raw p<0.05 counts vs null expectation; adj = 0
x = np.arange(len(ORDER))
axR.bar(x, cal["p_lt_0.05"], color="#b9c3cf", width=0.62, zorder=2,
        label="genes with raw p<0.05")
axR.axhline(0.05 * 19269, ls="--", color="#d62728", lw=1.4, zorder=3,
            label="null expectation (963.5)")
axR.set_xticks(x); axR.set_xticklabels(ORDER, rotation=40, ha="right", fontsize=8)
axR.set_ylabel("genes with raw p<0.05 (per cell type)")
axR.set_ylim(0, 1500)
axR.set_title("State, animal-level pseudobulk\n"
              "0 / 173,421 tests reach BH p<0.05", fontsize=11)
axR.text(0.5, -0.34, "most types fall below the null line: the test is conservative",
         transform=axR.transAxes, ha="center", fontsize=8.2, color="#444")
axR.legend(fontsize=7.6, loc="upper left", frameon=False)
axR.spines[["top", "right"]].set_visible(False)

fig.suptitle("Both composition and cell-state inference are underpowered at n=3 vs 3 "
             "(GSE267933)", fontsize=12.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "fig_ApC_power_asymmetry.png"), dpi=300)
plt.close(fig)
print("[Done] fig_ApC_power_asymmetry.png")

# ---------- 2. Fig 2 (revised) ----------
ns = np.arange(2, 41)
mdes = np.array([mde(n) for n in ns])

fig, ax = plt.subplots(figsize=(7.6, 5.0))
ax.plot(ns, mdes, "-", color="#1f77b4", lw=2.2, label="minimum detectable effect (80% power)")
ax.axhspan(0.3, 1.5, color="#2ca02c", alpha=0.13,
           label="effect sizes POCD produces (d 0.3-1.5)")
for d0, lab, col in [(1.5, "d=1.5", "#2ca02c"), (1.0, "d=1.0", "#7f7f7f"), (0.5, "d=0.5", "#9467bd")]:
    n0 = req_n(d0)
    ax.axvline(n0, ls="--", lw=1.2, color=col)
    ax.annotate(f"{lab}\nneeds n={n0}", (n0, 5.2), fontsize=8, color=col,
                ha="center", va="top")
ax.axvline(3, ls="-", lw=1.6, color="#d62728")
ax.annotate("n=3 (current 3v3)\ndetectable only at d≥3.1", (3, 3.35), fontsize=8.2,
            color="#d62728", ha="left", va="bottom", xytext=(6, 0),
            textcoords="offset points")
for nn in (3, 9, 17):
    ax.plot([nn], [mde(nn)], "o", color="#1f77b4", ms=5, zorder=4)
ax.set_xlabel("samples per group (n)")
ax.set_ylabel("minimum detectable effect (Cohen's d)")
ax.set_title("Detectable-effect boundary of pseudobulk testing\n"
             "even the largest POCD effects (d≈1.5) need n≈9 per group", fontsize=11)
ax.set_xlim(2, 40); ax.set_ylim(0, 6)
ax.legend(fontsize=7.6, loc="upper right", frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_ApC_mde_curve.png"), dpi=300)
plt.close(fig)
print("[Done] fig_ApC_mde_curve.png")
print("[Done] step5f_required_n.csv")
