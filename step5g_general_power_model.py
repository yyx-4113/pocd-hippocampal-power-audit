"""
Step 5g — 通用效能模型 + Monte Carlo 独立验证 (回应 PLOS ONE 拒因: "single dataset, no independent validation, not a research finding")
=========================================================================================================================
设计目标 (呼应编辑三点拒因):
  (1) 让 "所需样本量" 从单一 POCD 数据集脱钩 -> 变成由设计参数 (动物数 n, 每动物细胞数 k, 动物间方差 tau2, 细胞内噪声 sigma2, 生物学效应 Delta) 决定的通用量。
  (2) 用 Monte Carlo 模拟 (不依赖 GSE267933) 独立验证解析功率 power_at(), 即 "independent validation"。
  (3) 给出一个可行动的建设性结论: 细胞深度 k 只能压低细胞内噪声 sigma2/k, 永远压不到动物间方差 tau2;
      因此 k 对动物级推断的功效存在由 tau2 设定的天花板, 过了这个天花板只能靠增加动物 n。-> 正面产出 "research finding" + 可复用计算器。

统计约定 (与稿件既有 MDE/required-n 表完全一致, 复用 step5f 的 power_at):
  - 两样本 t 检验, 双尾, ALPHA=0.05, POWER=0.80
  - df = 2n-2 (等方差/合并), 非中心 t
  - 动物级标准化效应 d_a = Delta / sigma_a, 其中 sigma_a = sqrt(tau2 + sigma2/k) 为逐动物 pseudobulk 的标准差
  - 因此既有表 (d=1.5->9, 1.0->17, 0.5->64, 0.3->176) 被本脚本逐项复现 (= 同源)

输出:
  sc_out_general/sim_power_validation.csv   稀疏网格的 解析功效 vs 模拟功效 (验证)
  sc_out_general/required_n_general.csv     d_a -> required n (与 step5f 表同格式)
  sc_out_general/design_map.csv             (tau2, sigma2, Delta, k) -> d_a -> required n (设计映射, 通用)
  sc_out_general/fig_general_power_model.png  图: 左=功效~n (多 d_a, 含模拟点); 右=d_a~k 天花板效应
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import nct

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sc_out_general")
os.makedirs(OUT, exist_ok=True)

ALPHA, POWER = 0.05, 0.80
SEED = 20261027
rng = np.random.default_rng(SEED)

# ---------- 0. 精确功率函数 (与稿件既有方法逐字一致) ----------
def power_at(n, d):
    """向量化非中心 t 功效; n 标量或数组, d 标量或数组 (广播)。"""
    n = np.asarray(n, dtype=float)
    d = np.asarray(d, dtype=float)
    df = 2 * n - 2
    tcrit = nct.ppf(1 - ALPHA / 2, df, 0)
    nc = d * np.sqrt(n / 2.0)
    return nct.sf(tcrit, df, nc) + nct.cdf(-tcrit, df, nc)

def req_n(d):
    for n in range(2, 601):
        if power_at(n, d) >= POWER:
            return int(n)
    return int(np.nan)

# sanity: 必须复现稿件既有 required-n 表 (同源校验)
print("[sanity] required n for d in {1.5,1.0,0.5,0.3}:",
      {d: req_n(d) for d in (1.5, 1.0, 0.5, 0.3)})

# ---------- 1. Monte Carlo 模拟 (向量化 over reps) ----------
def sim_power(n, k, tau2, sigma2, delta, reps=500):
    """模拟动物级两样本 t 检验 (合并方差, 与 power_at 的 df=2n-2 一致)。
    返回经验功效 = P(|t| > tcrit)。不依赖任何真实数据集。"""
    # 逐动物均值: (reps, n)  control 中心 0, treatment 中心 delta, 均带动物间方差 tau2
    bC = rng.normal(0.0, np.sqrt(tau2), size=(reps, n))
    bT = rng.normal(delta, np.sqrt(tau2), size=(reps, n))
    # 细胞级噪声: (reps, 2, n, k)
    noise = rng.normal(0.0, np.sqrt(sigma2), size=(reps, 2, n, k))
    meanC = bC + noise[:, 0].mean(axis=-1)   # (reps, n) 逐元素相加, 不可加 None 维 (否则广播成 reps×reps×n)
    meanT = bT + noise[:, 1].mean(axis=-1)   # (reps, n)
    xbar = meanC.mean(axis=1); ybar = meanT.mean(axis=1)
    s2x = ((meanC - xbar[:, None]) ** 2).sum(axis=1) / (n - 1)
    s2y = ((meanT - ybar[:, None]) ** 2).sum(axis=1) / (n - 1)
    sp = np.sqrt(((n - 1) * s2x + (n - 1) * s2y) / (2 * n - 2))
    t = (ybar - xbar) / (sp * np.sqrt(2.0 / n))
    tcrit = nct.ppf(1 - ALPHA / 2, 2 * n - 2, 0)
    return float((np.abs(t) > tcrit).mean())

# ---------- 2. 验证网格: 固定 d_a, 扫 n, 比较 解析 vs 模拟 ----------
print("\n[validation] analytic vs Monte Carlo power (fixed d_a) ...")
val_rows = []
for d_a in (0.3, 0.5, 1.0, 1.5):
    # 选一组 tau2/sigma2/k 使得 d_a = delta/sqrt(tau2+sigma2/k) 恰好等于目标 d_a
    tau2, sigma2, k = 1.0, 1.0, 200
    sigma_a = np.sqrt(tau2 + sigma2 / k)
    delta = d_a * sigma_a
    for n in (3, 6, 9, 17, 30):
        ap = float(power_at(n, d_a))
        sp_ = sim_power(n, k, tau2, sigma2, delta, reps=500)
        val_rows.append({"d_a": d_a, "n": n, "k": k, "tau2": tau2, "sigma2": sigma2,
                         "delta": round(delta, 4), "analytic_power": round(ap, 4),
                         "sim_power": round(sp_, 4)})
        print(f"  d_a={d_a} n={n:2d}: analytic={ap:.3f} sim={sp_:.3f}")
val = pd.DataFrame(val_rows)
val.to_csv(os.path.join(OUT, "sim_power_validation.csv"), index=False)

# ---------- 3. required-n 表 (与既有表同格式, 同源) ----------
print("\n[required-n] general table (d_a -> n) ...")
rn_rows = []
for d in [1.5, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2]:
    rn_rows.append({"effect_size_d": d, "required_n_per_group": req_n(d),
                    "total_n": 2 * req_n(d)})
rn = pd.DataFrame(rn_rows)
rn.to_csv(os.path.join(OUT, "required_n_general.csv"), index=False)
print(rn.to_string(index=False))

# ---------- 4. 设计映射: (tau2, sigma2, Delta, k) -> d_a -> required n ----------
print("\n[design map] (tau2, sigma2, Delta, k) -> d_a -> required n ...")
dm_rows = []
tau2_grid = [0.25, 1.0]
sigma2_grid = [1.0, 4.0]
Delta_grid = [0.5, 1.0, 1.5]          # 生物学效应 (动物均值真实偏移)
k_grid = [50, 200, 1000, 10000]
for tau2 in tau2_grid:
    for sigma2 in sigma2_grid:
        for delta in Delta_grid:
            d_cell = delta / np.sqrt(tau2 + sigma2)   # 细胞级标准化效应
            for k in k_grid:
                sigma_a = np.sqrt(tau2 + sigma2 / k)
                d_a = delta / sigma_a
                dm_rows.append({"tau2": tau2, "sigma2": sigma2, "Delta": delta,
                                "d_cell": round(d_cell, 3), "k_cells_per_animal": k,
                                "sigma_a": round(sigma_a, 3), "d_a_animal": round(d_a, 3),
                                "required_n_per_group": req_n(d_a)})
dm = pd.DataFrame(dm_rows)
dm.to_csv(os.path.join(OUT, "design_map.csv"), index=False)

# ---------- 5. 图 ----------
print("\n[figure] plotting ...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 左: 功效 ~ n, 多 d_a, 叠加模拟验证点
ax = axes[0]
ns = np.arange(2, 41)
for d_a, c in zip((0.3, 0.5, 1.0, 1.5), ("#888", "#555", "#222", "black")):
    ax.plot(ns, power_at(ns, d_a), color=c, label=f"d_a={d_a}")
# 模拟点
for _, r in val.iterrows():
    ax.scatter(r["n"], r["sim_power"], color="red", s=18, zorder=5)
ax.axhline(POWER, ls="--", color="gray", lw=1)
ax.text(2.2, 0.82, "80% power", color="gray", fontsize=9)
ax.set_xlabel("Animals per group (n)"); ax.set_ylabel("Power (animal-level Welch/pooled t)")
ax.set_title("Power vs n — matches the manuscript required-n table")
ax.set_ylim(0, 1.02); ax.legend(fontsize=9, title="animal-level d_a")

# 右: d_a ~ k 天花板效应 (固定生物学效应 Delta, 不同 tau2/sigma2)
ax = axes[1]
kk = np.array([10, 30, 50, 100, 200, 500, 1000, 5000, 20000], dtype=float)
for tau2, sigma2, delta, ls in [
    (0.25, 4.0, 1.0, "-"), (1.0, 1.0, 1.0, "--"), (0.25, 1.0, 1.0, ":")]:
    d_a_of_k = delta / np.sqrt(tau2 + sigma2 / kk)
    ax.plot(kk, d_a_of_k, ls=ls,
            label=f"tau2={tau2}, sigma2={sigma2}, Delta={delta}")
    ax.axhline(delta / np.sqrt(tau2), ls=":", color="gray", lw=0.8)  # k->inf 天花板
ax.set_xscale("log")
ax.set_xlabel("Cells per animal (k, log)")
ax.set_ylabel("Resulting animal-level d_a")
ax.set_title("Cell depth k has a ceiling set by animal variance tau2")
ax.legend(fontsize=8)

fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_general_power_model.png"), dpi=300)
print("[done] outputs in", OUT)
