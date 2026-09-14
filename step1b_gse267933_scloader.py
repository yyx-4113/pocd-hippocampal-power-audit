"""
Step 1b — GSE267933 (海马 scRNA-seq, POCD) 单细胞骨架加载
输入: GSE267933_RAW.tar (suppl, 6 样本 10x filtered_feature_bc_matrix: C1/C2/C3=Control, S1/S2/S3=Surgery)
流程: 解压 -> 逐样本读 10x (scipy mmread + 手动解析) -> 合并 -> 基础 QC -> 存 .h5ad
输出:
  - GSE267933.h5ad                 合并后的 AnnData (obs: sample/group)
  - GSE267933_qc_summary.md        n_cells/n_genes/mito% 概览
  - fig_sc_qc_violin.png          QC 指标小提琴图(按样本)
说明: 骨架(Step 5 再做 Harmony 去批次/聚类/注释). 18 月龄雄 C57BL/6 小鼠海马.
"""
import os
import gzip
import tarfile
import numpy as np
import pandas as pd
import scipy.io
import scanpy as sc

OUT = os.path.dirname(os.path.abspath(__file__))
TAR = os.path.join(OUT, "GSE267933_RAW.tar")
RAW = os.path.join(OUT, "GSE267933_RAW")

# ---------- 1. 解压 ----------
if not os.path.isdir(RAW) or not os.listdir(RAW):
    os.makedirs(RAW, exist_ok=True)
    with tarfile.open(TAR, "r") as tf:
        tf.extractall(RAW)
    print("[extract] done ->", RAW)

# ---------- 2. 逐样本读 10x (扁平文件, 手动解析) ----------
prefixes = sorted({f.rsplit("_matrix", 1)[0]
                   for f in os.listdir(RAW) if f.endswith("_matrix.mtx.gz")})
print("[load] 发现样本前缀:", prefixes)


def read_sample(base):
    mtx = scipy.io.mmread(os.path.join(RAW, base + "_matrix.mtx.gz")).T.tocsr()
    with gzip.open(os.path.join(RAW, base + "_barcodes.tsv.gz"), "rt") as f:
        barcodes = [l.rstrip("\n").split("\t")[0] for l in f if l.strip()]
    genes = []
    with gzip.open(os.path.join(RAW, base + "_features.tsv.gz"), "rt") as f:
        for l in f:
            p = l.rstrip("\n").split("\t")
            genes.append(p[1] if len(p) > 1 else p[0])
    ad = sc.AnnData(mtx)
    ad.obs_names = [f"{base}-{b}" for b in barcodes]
    ad.var_names = pd.Index(genes)
    ad.var_names_make_unique()
    tag = base.split("_")[-1]                       # C1 / S1 ...
    ad.obs["sample"] = base
    ad.obs["gsm"] = base.split("_")[0]
    ad.obs["group"] = "Control" if tag.startswith("C") else "Surgery"
    print(f"  {base}: cells={ad.n_obs}, genes={ad.n_vars}, group={ad.obs['group'].values[0]}")
    return ad


adatas = [read_sample(p) for p in prefixes]

# ---------- 3. 合并 ----------
adata = sc.concat(adatas, label="sample", index_unique="-")
print(f"[concat] 总 cells={adata.n_obs}, genes={adata.n_vars}")

# ---------- 4. 基础 QC ----------
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                           log1p=False, inplace=True)
print("[qc] mito% median:", round(float(np.median(adata.obs["pct_counts_mt"])), 2))

# ---------- 5. 保存 ----------
adata.write(os.path.join(OUT, "GSE267933.h5ad"))
print("[save] GSE267933.h5ad")

# ---------- 6. 摘要 + 图 ----------
summ = adata.obs.groupby("sample").agg(
    group=("group", "first"),
    cells=("n_genes_by_counts", "size"),
    genes_median=("n_genes_by_counts", "median"),
    counts_median=("total_counts", "median"),
    mito_pct_median=("pct_counts_mt", "median"),
).round(1)
print(summ)

with open(os.path.join(OUT, "GSE267933_qc_summary.md"), "w", encoding="utf-8") as f:
    f.write("# Step 1b 摘要 — GSE267933 单细胞骨架 QC\n\n")
    f.write(f"- 总细胞数: {adata.n_obs}; 基因数: {adata.n_vars}\n")
    f.write(f"- mito% 中位: {float(np.median(adata.obs['pct_counts_mt'])):.2f}\n\n")
    f.write(summ.to_markdown() + "\n")

sc.settings.set_figure_params(dpi=120)
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
              groupby="sample", jitter=0.4, multi_panel=True,
              save="_sc_qc_violin.png", show=False)
print("[Done] Step 1b 完成")
