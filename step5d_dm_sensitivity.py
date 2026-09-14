# -*- coding: utf-8 -*-
"""m5 — Dirichlet-multinomial (DM) sensitivity for the animal-level composition claim.

Composition claim (step5d): across the 9 hippocampal cell types, no per-type proportion
differs significantly between Control and Surgery at the animal level (microglia Welch
p=0.47 BH; ependymal p=0.34 BH; 0/9 significant after BH).

Reviewer m3/m5 noted the 9 proportions are compositional (sum to 1), so the 9 Welch
tests are correlated and BH is an approximation. This script fits an overdispersed
Dirichlet-Multinomial at the PROPORTION level and asks: under a proper compositional
null with the observed between-animal overdispersion, does the observed animal-level
min-BH p (0.34) stay central? If yes, the 'none significant' conclusion is robust to
compositional overdispersion.

Data: sc_out/step5d_composition_animallevel.csv (per-animal proportions, %). We rebuild
the 6x9 proportion matrix from per_animal_ctrl / per_animal_surg.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "sc_out", "step5d_composition_animallevel.csv")
OUT = os.path.join(HERE, "sc_out", "step5d_dm_sensitivity.csv")
SEED = 20260914
RNG = np.random.default_rng(SEED)
B = 2000
CTRL, SURG = [0, 1, 2], [3, 4, 5]


def load_props():
    df = pd.read_csv(CSV)
    # per_animal_ctrl/surg are 3 animal values PER TYPE (row = cell type).
    # Rebuild the 6 x 9 matrix: rows = animals (0-2 Control, 3-5 Surgery), cols = cell types.
    ctrl = np.array([list(map(float, s.split(";"))) for s in df["per_animal_ctrl"]]) / 100.0
    surg = np.array([list(map(float, s.split(";"))) for s in df["per_animal_surg"]]) / 100.0
    P = np.vstack([ctrl.T, surg.T])          # (6, 9)
    assert P.shape == (6, len(df)), f"bad shape {P.shape}"
    return P, df["celltype"].tolist()


def animal_stats(P, ctrl, surg):
    ps = [ttest_ind(P[ctrl, j], P[surg, j], equal_var=False)[1] for j in range(P.shape[1])]
    _, bh, _, _ = multipletests(ps, method="fdr_bh")
    return float(bh.min()), ps


def fit_dm(P):
    mu = P.mean(0)
    s2 = P.var(0, ddof=1)
    a0 = max(np.mean(mu * (1 - mu)) / np.mean(s2) - 1.0, 1e-3)
    return a0, a0 * mu, mu, s2


def simulate(alpha, b=B):
    k = alpha.shape[0]
    minbh = np.empty(b)
    any_sig = np.empty(b)
    animals = np.arange(6)
    for r in range(b):
        Psim = RNG.dirichlet(alpha) if False else np.vstack(
            [RNG.dirichlet(alpha) for _ in range(6)])
        pick = RNG.choice(animals, 3, replace=False)
        c = list(pick); s = [a for a in animals if a not in pick]
        mb, ps = animal_stats(Psim, c, s)
        minbh[r] = mb
        any_sig[r] = float(multipletests(ps, method="fdr_bh")[1].min() < 0.05)
    return minbh, any_sig


def main():
    P, types = load_props()
    a0, alpha, mu, s2 = fit_dm(P)
    rho = 1.0 / (a0 + 1.0)
    print(f"[fit] alpha_0 = {a0:.3f}  rho = {rho:.4f}")
    print(f"[fit] mu = {np.round(mu,4)}")
    print(f"[fit] s2 = {np.round(s2,6)}")

    obs_minbh, obs_ps = animal_stats(P, CTRL, SURG)
    print(f"[observed] animal-level min BH p = {obs_minbh:.4f}")
    for t, p in zip(types, obs_ps):
        print(f"    {t:14s} Welch p = {p:.4f}")

    rows = []
    for g in [0.25, 0.5, 1.0, 2.0, 5.0]:
        sim, any_sig = simulate(g * alpha)
        med = float(np.median(sim))
        emp = float((sim <= obs_minbh).mean())
        frac_sig = float(any_sig.mean())
        print(f"[grid rho x{g:<4}] alpha_0={g*a0:8.2f}  median minBH={med:.3f}  "
              f"empP(<=0.34)={emp:.3f}  frac sims w/ any BH<0.05={frac_sig:.3f}")
        rows.append({"rho_mult": g, "alpha_0": g * a0, "median_minBH_p": med,
                     "emp_p_le_obs": emp, "frac_any_BH_lt_0.05": frac_sig})

    sim_fit, any_fit = simulate(alpha)
    print(f"\n[PRIMARY fitted DM rho={rho:.4f}] median minBH={np.median(sim_fit):.3f}  "
          f"empP(<=0.34)={(sim_fit<=obs_minbh).mean():.3f}  "
          f"10-90 pct={np.percentile(sim_fit,10):.3f}-{np.percentile(sim_fit,90):.3f}  "
          f"frac any BH<0.05={any_fit.mean():.3f}")

    pd.DataFrame(rows).to_csv(OUT, index=False)
    print("SAVED:", OUT)


if __name__ == "__main__":
    main()
