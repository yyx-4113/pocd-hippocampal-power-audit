"""
Step 5e (a) — 导出全基因 counts pseudobulk (供 R limma-voom 敏感性分析)
=====================================================================
动机 (审稿意见 P2):
  "pseudobulk 检验选择未论证" —— 主分析用 Welch t (log1p 均值), 审稿人要求补 limma-voom。
  但 layers['counts'] 仅含 2000 HVG, 无法做全基因 voom。故从原始 10x counts 重建:
    GSE267933.h5ad (20684 细胞 x 27998 基因, 原始整数 counts)
    + sc_out/GSE267933_processed.h5ad 的 QC 后 barcode 与 celltype 注释
  -> 逐 (细胞类型 x 样本) 求和 counts -> genes x 6 样本矩阵 CSV

输出: sc_out/pbcounts/pbcounts_<celltype>.csv  (行=基因, 列=6 样本; 另存 sample_meta.csv)
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "sc_out")
OUT = os.path.join(WORK, "pbcounts")
os.makedirs(OUT, exist_ok=True)


def log(*a):
    print(*a, flush=True)


def main():
    log("[load] raw counts GSE267933.h5ad")
    raw = ad.read_h5ad(os.path.join(HERE, "GSE267933.h5ad"))
    log(f"  raw: {raw.n_obs} cells x {raw.n_vars} genes; X dtype={raw.X.dtype}")
    is_int = sp.issparse(raw.X)
    log(f"  X sparse={is_int}")
    # 确认整数 counts
    sample_vals = raw.X[:5].data[:8] if sp.issparse(raw.X) else raw.X[:5, :8]
    log(f"  value sample: {sample_vals}")

    log("[load] processed obs (QC + celltype)")
    proc = ad.read_h5ad(os.path.join(WORK, "GSE267933_processed.h5ad"), backed="r")
    obs = proc.obs[["sample", "group", "celltype"]].copy()
    obs["sample"] = obs["sample"].astype(str)
    log(f"  processed: {proc.n_obs} cells; celltypes={sorted(obs.celltype.unique())}")

    # 对齐细胞 (processed obs_names 是 QC 后保留的 barcode)
    keep = obs.index.intersection(raw.obs_names)
    log(f"  overlap cells: {len(keep)} / processed {proc.n_obs} / raw {raw.n_obs}")
    assert len(keep) > 15000, "barcode 对齐异常"

    sub = raw[keep].copy()
    sub.obs["celltype"] = obs.loc[keep, "celltype"].values
    sub.obs["sample"] = obs.loc[keep, "sample"].values
    sub.obs["group"] = obs.loc[keep, "group"].values

    samples = sorted(sub.obs["sample"].unique())
    grp = {s: sub.obs.loc[sub.obs["sample"] == s, "group"].iloc[0] for s in samples}
    pd.DataFrame({"sample": samples, "group": [grp[s] for s in samples]}).to_csv(
        os.path.join(OUT, "sample_meta.csv"), index=False)
    log(f"  samples: {[(s, grp[s]) for s in samples]}")

    celltypes = [c for c in sub.obs.celltype.value_counts().index
                 if (sub.obs.celltype.values == c).sum() >= 100]
    log(f"  celltypes (>=100 cells): {celltypes}")

    genes = pd.Index(sub.var_names)
    for ct in celltypes:
        mask = (sub.obs.celltype.values == ct)
        Xc = sub.X[mask]
        cellsamp = sub.obs["sample"].values[mask]
        cols = []
        for s in samples:
            m = (cellsamp == s)
            if m.sum() == 0:
                cols.append(np.zeros(len(genes), dtype=np.int64))
            else:
                v = np.asarray(Xc[m].sum(axis=0)).ravel()
                cols.append(np.rint(v).astype(np.int64))
        df = pd.DataFrame(np.vstack(cols).T, index=genes, columns=samples)
        # 过滤全零基因
        df = df[df.sum(axis=1) > 0]
        p = os.path.join(OUT, f"pbcounts_{ct}.csv")
        df.to_csv(p)
        log(f"  [{ct}] {df.shape[0]} genes x {df.shape[1]} samples -> {os.path.basename(p)} "
            f"(library sizes: {df.sum(axis=0).tolist()})")

    log("[Done] pseudobulk counts exported ->", OUT)


if __name__ == "__main__":
    main()
