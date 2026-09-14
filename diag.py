import pandas as pd, numpy as np
from scipy import stats

df = pd.read_csv("GSE276942_gene_expression_filter.csv.gz")
expr = df.drop(columns=["gene_name"]).groupby(df["gene_name"]).mean()
genes = expr.index.values
X = expr.values.astype(float)
print("shape", X.shape, "min", float(X.min()), "max", float(X.max()),
      "median", float(np.median(X)), "mean", float(X.mean()))
print("per-sample mean:", {c: round(float(expr[c].mean()), 1) for c in expr.columns})

meta = {
    "CON_1": ("Control", "CON"), "CON_2": ("Control", "CON"),
    "CON_4": ("Control", "CON"), "CON_5": ("Control", "CON"),
    "24H_2": ("24H", "24H"), "24H_3": ("24H", "24H"), "24H_4": ("24H", "24H"),
    "3D_1": ("3D", "3D"), "3D_3": ("3D", "3D"), "3D_4": ("3D", "3D"),
    "7D_1": ("7D", "7D"), "7D_2": ("7D", "7D"), "7D_3": ("7D", "7D"),
    "1M_1": ("1M", "1M"), "1M_2": ("2M", "2M"), "1M_3": ("3M", "3M"), "1M_4": ("4M", "4M"),
}
groups = np.array([meta[c][0] for c in expr.columns])
print("group counts:", {g: int((groups == g).sum()) for g in ["Control", "24H", "3D", "7D", "1M", "2M", "3M", "4M"]})

Xlog = np.log2(X + 1)
ctrl = groups == "Control"
for tp in ["24H", "3D", "7D"]:
    m = groups == tp
    a = Xlog[:, ctrl]
    b = Xlog[:, m]
    lfc = b.mean(1) - a.mean(1)
    t, p = stats.ttest_ind(b, a, axis=1, equal_var=False)
    p = np.where(np.isfinite(p), np.clip(p, 0, 1), 1.0)
    print(f"{tp}: n_post={int(m.sum())} |lfc|>0.5={int((np.abs(lfc)>0.5).sum())} "
          f"p<0.05={int((p<0.05).sum())} pmin={float(np.nanmin(p)):.2e} lfc_median={float(np.median(lfc)):.3f}")

# overall gene variance in log2 space
rowvar = Xlog.var(1)
print("log2 row-var: median", float(np.median(rowvar)), "max", float(rowvar.max()))
