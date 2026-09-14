"""
Step 5-6 — GSE267933 术后海马单细胞图谱: 定位 + 深入分析
================================================================
Step 5 定位: QC -> 归一化 -> HVG -> PCA -> Harmony(去批次 by sample) -> 邻接 -> UMAP -> Leiden
            -> marker 打分自动注释细胞类型 -> 特征基因定位(DotPlot/FeaturePlot/Violin)
Step 6 深入: 细胞比例变化 / 亚群再聚类(少突谱系、小胶质稳态↔DAM) / 拟时序(DPT) /
            配体-受体轴(补体-趋化) / pseudobulk 逐细胞类型 DEG

输入: GSE267933.h5ad (20684 细胞 x 27998 基因; obs: sample/gsm/group)
输出: sc_out/ 下 h5ad / csv / 图 / step5_summary.md
注意: 3v3 设计不做个体级统计推断; 统计放在 pseudobulk(逐样本聚合并) 与 bulk 层;
      单细胞层定位为"定位 + 方向验证"。
"""
import os, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
OUT = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(OUT, "sc_out")
os.makedirs(WORK, exist_ok=True)
SEED = 42
np.random.seed(SEED)

FEATURE_GENES = ["Fos", "Gltp", "Junb", "Runx2", "Sgk1", "Sla2"]

MARKERS = {
    "Neuron":        ["Rbfox3", "Syt1", "Snap25", "Syp", "Elavl2", "Meg3", "Syt2"],
    "CA1":           ["Wfs1", "Fibcd1", "Amigo2", "Grin2b", "Ppp3ca"],
    "CA3":           ["Cpne4", "Grp", "Nptx2", "Trh", "Wnt7b"],
    "DG":            ["Prox1", "Calb2", "Dsp", "C1ql2", "Nrp2"],
    "Microglia":     ["Cx3cr1", "P2ry12", "Tmem119", "Hexb", "Csf1r", "Aif1", "Siglech"],
    "Astrocyte":     ["Gfap", "Aqp4", "Slc1a2", "Aldh1l1", "Gja1", "S100b"],
    "Oligodendrocyte": ["Mbp", "Plp1", "Mog", "Mag", "Cnp", "Mal", "Mobp", "Opalin"],
    "OPC":           ["Pdgfra", "Cspg4", "Lhfpl3", "Sox10", "Olig1"],
    "Endothelial":   ["Cldn5", "Pecam1", "Flt1", "Slco1a4"],
    "Ependymal":     ["Ttr", "Hdc", "Foxj1", "Ccdc153"],
    "VLMC":          ["Dcn", "Col1a1", "Col3a1", "Slc6a13"],
    "Pericyte":      ["Pdgfrb", "Rgs5", "Kcnj8", "Abcc9"],
}


def log(*a):
    print(*a, flush=True)


def main():
    log("[load] GSE267933.h5ad")
    ad = sc.read_h5ad(os.path.join(OUT, "GSE267933.h5ad"))
    ad.layers["counts"] = ad.X.copy()

    # ---------- QC ----------
    n0 = ad.n_obs
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    ad = ad[ad.obs.pct_counts_mt < 20].copy()
    ad.obs["group"] = ad.obs["group"].astype(str)
    ad.obs["sample"] = ad.obs["sample"].astype(str)
    log(f"[QC] {n0} -> {ad.n_obs} 细胞, {ad.n_vars} 基因; group={dict(ad.obs.group.value_counts())}")

    # ---------- 归一化 / HVG / PCA ----------
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.layers["lognorm"] = ad.X.copy()
    sc.pp.highly_variable_genes(ad, n_top_genes=2000, batch_key="sample")
    ad.raw = ad
    ad = ad[:, ad.var.highly_variable].copy()
    sc.pp.scale(ad, max_value=10)
    sc.tl.pca(ad, n_comps=50, svd_solver="arpack", random_state=SEED)

    # ---------- Harmony 去批次 (直接调 harmonypy, 规避 scanpy 包装器 torch 兼容问题) ----------
    rep, Z = "X_pca", None
    try:
        import harmonypy, torch
        ho = harmonypy.run_harmony(ad.obsm["X_pca"], ad.obs, ["sample"],
                                   max_iter_harmony=20, random_state=SEED)
        Z = ho.Z_corr
        Z = Z.cpu().detach().numpy() if isinstance(Z, torch.Tensor) else np.asarray(Z, dtype=np.float32)
        Z = np.asarray(Z, dtype=np.float32)
    except Exception as ex:
        log("[Harmony] 直接调用失败:", type(ex).__name__, ex)
    if Z is None:
        try:
            sc.external.pp.harmony_integrate(ad, "sample", random_state=SEED)
            rep = "X_pca_harmony" if "X_pca_harmony" in ad.obsm else "X_pca"
        except Exception as ex2:
            log("[Harmony] 包装器亦失败:", type(ex2).__name__)
    else:
        if Z.shape[0] != ad.n_obs:
            Z = Z.T
        ad.obsm["X_pca_harmony"] = Z
        rep = "X_pca_harmony"
    log(f"[Harmony] rep={rep}" + (f" shape={ad.obsm[rep].shape}" if rep != "X_pca" else " (未校正)"))
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30, use_rep=rep, random_state=SEED)
    sc.tl.umap(ad, random_state=SEED)
    sc.tl.leiden(ad, resolution=0.6, key_added="leiden", random_state=SEED, flavor="igraph",
                 n_iterations=2, directed=False)
    log(f"[Leiden] {ad.obs.leiden.nunique()} 簇")

    # ---------- 细胞类型注释 (marker 打分) ----------
    sig = {k: [g for g in v if g in ad.raw.var_names] for k, v in MARKERS.items()}
    for k, v in sig.items():
        sc.tl.score_genes(ad, v, score_name=f"sc_{k}", use_raw=True)
    score_cols = [f"sc_{k}" for k in sig]
    cl_mean = ad.obs.groupby("leiden")[score_cols].mean()
    cl_mean.columns = [c.replace("sc_", "") for c in cl_mean.columns]
    cl2type = cl_mean.idxmax(1)
    # 神经元亚型簇若总分低于泛神经元, 归并到 Neuron
    ad.obs["celltype"] = ad.obs.leiden.map(cl2type).astype(str)
    ad.obs["celltype"] = ad.obs.celltype.replace(
        {"CA1": "Neuron_CA1", "CA3": "Neuron_CA3", "DG": "Neuron_DG"})
    log("[annot] 簇->类型:\n" + cl_mean.round(2).to_string())
    log("[annot] 细胞类型分布:\n" + ad.obs.celltype.value_counts().to_string())

    # marcador 表 (仅用 HVG 加速; 注释已由 marker score 完成)
    sc.tl.rank_genes_groups(ad, "leiden", method="wilcoxon", use_raw=False)
    mk = sc.get.rank_genes_groups_df(ad, None).groupby("group").head(10)
    mk.to_csv(os.path.join(WORK, "sc_cluster_markers.csv"), index=False)
    cl_mean.to_csv(os.path.join(WORK, "sc_cluster_typescore.csv"))

    # ---------- 图: UMAP ----------
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.4))
    sc.pl.umap(ad, color="leiden", ax=axes[0], show=False, legend_loc="on data",
               legend_fontsize=7, title="Leiden")
    sc.pl.umap(ad, color="celltype", ax=axes[1], show=False, legend_loc="on data",
               legend_fontsize=6, title="Cell type")
    sc.pl.umap(ad, color="group", ax=axes[2], show=False, title="Group")
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_umap.png"), dpi=150)
    plt.close(fig)

    # ---------- 图: marker DotPlot ----------
    keys = [k for k in sig if k in ad.obs.celltype.unique() or True]
    dp_genes = [g for k in ["Neuron","CA1","CA3","DG","Microglia","Astrocyte",
                            "Oligodendrocyte","OPC","Endothelial","Ependymal"]
                for g in sig.get(k, [])]
    dp_genes = sorted(set(dp_genes))
    try:
        dp = sc.pl.dotplot(ad, dp_genes, groupby="celltype", use_raw=True, show=False,
                           standard_scale="var", return_fig=True)
        dp.savefig(os.path.join(WORK, "fig_sc_marker_dotplot.png"), dpi=150, bbox_inches="tight")
    except Exception as ex:
        log("[dotplot markers] 跳过:", ex)

    # ---------- 细胞比例 (Control vs Surgery) ----------
    ct = pd.crosstab(ad.obs.celltype, ad.obs.group)
    prop = ct.div(ct.sum(0), axis=1)
    comp_rows = []
    for c in ct.index:
        cc = ct.loc[c, "Control"]; sc_ = ct.loc[c, "Surgery"]
        tc = ct["Control"].sum(); ts = ct["Surgery"].sum()
        try:
            _, p = stats.fisher_exact([[cc, tc - cc], [sc_, ts - sc_]])
        except Exception:
            p = np.nan
        comp_rows.append({"celltype": c, "n_ctrl": int(cc), "n_surg": int(sc_),
                          "pct_ctrl": round(100 * cc / tc, 2), "pct_surg": round(100 * sc_ / ts, 2),
                          "fisher_p": p})
    comp = pd.DataFrame(comp_rows).sort_values("pct_surg", ascending=False)
    comp["BH_p"] = multipletests(comp.fisher_p.fillna(1), method="fdr_bh")[1]
    comp.to_csv(os.path.join(WORK, "sc_celltype_composition.csv"), index=False)
    log("[composition]\n" + comp.to_string(index=False))

    fig, ax = plt.subplots(figsize=(8, 4.2))
    x = np.arange(len(prop)); w = 0.38
    ax.bar(x - w / 2, 100 * prop["Control"].reindex(prop.index), w, label="Control", color="#4c72b0")
    ax.bar(x + w / 2, 100 * prop["Surgery"].reindex(prop.index), w, label="Surgery", color="#c44e52")
    ax.set_xticks(x); ax.set_xticklabels(prop.index, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("% of cells"); ax.set_title("Cell-type composition"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_composition.png"), dpi=150)
    plt.close(fig)

    # ---------- 特征基因定位 ----------
    # 注意: 不把基因名写进 ad.obs (会与 raw.var_names 冲突导致 dotplot 报错);
    #       统一用 ad.raw.X (全基因 lognorm) 自行聚合。
    present_feat = [g for g in FEATURE_GENES if g in ad.raw.var_names]
    RV = list(ad.raw.var_names)
    gidx = {g: i for i, g in enumerate(RV)}
    fidx = [gidx[g] for g in present_feat]
    Xf = ad.raw.X                                            # cells x 19269, 稀疏
    F = np.asarray(Xf[:, fidx].todense(), dtype=float)        # cells x k
    Fdf = pd.DataFrame(F, columns=present_feat)
    Fdf["celltype"] = ad.obs.celltype.values
    ct_meanf = Fdf.groupby("celltype")[present_feat].mean()
    pct_df = pd.DataFrame((F > 0).astype(float), columns=present_feat)
    pct_df["celltype"] = ad.obs.celltype.values
    ct_pctf = pct_df.groupby("celltype")[present_feat].mean()
    ct_meanf.to_csv(os.path.join(WORK, "sc_featuregene_meancount_by_celltype.csv"))
    ct_pctf.to_csv(os.path.join(WORK, "sc_featuregene_pct_by_celltype.csv"))
    # 气泡图: 颜色 = 平均表达(z per gene), 大小 = 表达细胞比例
    zc = (ct_meanf - ct_meanf.mean(0)) / ct_meanf.std(0).replace(0, np.nan)
    zc = zc.fillna(0)
    fig, ax = plt.subplots(figsize=(1.0 * len(present_feat) + 3.4, 0.42 * len(zc) + 2.2))
    for i, ct_ in enumerate(zc.index):
        for j, g in enumerate(present_feat):
            ax.scatter(j, i, s=20 + 380 * ct_pctf.loc[ct_, g], c=[zc.loc[ct_, g]],
                       cmap="Reds", vmin=-2, vmax=2, edgecolors="grey", linewidths=0.4)
    ax.set_xticks(range(len(present_feat))); ax.set_xticklabels(present_feat, rotation=0)
    ax.set_yticks(range(len(zc))); ax.set_yticklabels(zc.index, fontsize=8)
    ax.set_xlim(-0.6, len(present_feat) - 0.4); ax.set_ylim(-0.6, len(zc) - 0.4)
    ax.invert_yaxis(); ax.set_title("Feature genes by cell type (color=scaled mean, size=% expressing)")
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_feature_dotplot.png"), dpi=150)
    plt.close(fig)
    log("[feature genes] celltype 平均表达:\n" + ct_meanf.round(3).to_string())
    # FeaturePlot
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    for i, g in enumerate(present_feat[:6]):
        sc.pl.umap(ad, color=g, ax=axes[i // 3][i % 3], show=False, use_raw=True,
                   title=g, cmap="Reds")
    fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_feature_umap.png"), dpi=150)
    plt.close(fig)

    # ---------- pseudobulk: 逐细胞类型 DEG (Surgery vs Control) ----------
    # 用 ad.raw.X (全基因 lognorm) 聚合成逐样本伪 bulk —— 与 genes=raw.var_names 严格对齐
    genes = np.asarray(ad.raw.var_names)

    def pseudobulk(adx, ct_name):
        mask = (adx.obs.celltype.values == ct_name)
        if mask.sum() < 20:
            return None, None
        Xsub = adx.raw.X[mask]
        samples = adx.obs["sample"].values[mask]
        groups = adx.obs["group"].values[mask]
        rows, meta = [], []
        for s in pd.unique(samples):
            m = (samples == s)
            rows.append(np.asarray(Xsub[m].mean(0)).ravel())
            meta.append({"sample": s, "group": groups[m][0], "n": int(m.sum())})
        if len(rows) < 4:
            return None, None
        return pd.DataFrame(rows, index=[x["sample"] for x in meta]), pd.DataFrame(meta).set_index("sample")

    pb_rows = []
    for ct_name in ad.obs.celltype.unique():
        X, m = pseudobulk(ad, ct_name)
        if X is None or m.group.nunique() < 2:
            continue
        a = X[(m.group == "Surgery").values].values
        b = X[(m.group == "Control").values].values
        if len(a) < 2 or len(b) < 2:
            continue
        t, p = stats.ttest_ind(a, b, axis=0, equal_var=False)
        p = np.where(np.isfinite(p), p, 1.0)
        lfc = a.mean(0) - b.mean(0)
        for g in FEATURE_GENES:
            if g in gidx:
                gi = gidx[g]
                pb_rows.append({"celltype": ct_name, "gene": g, "log2FC_pb": lfc[gi],
                                "p": p[gi], "n_ctrl": len(b), "n_surg": len(a)})
        # 全基因 pseudobulk DEG 存档
        pd.DataFrame({"gene": genes, "log2FC": lfc, "p": p}).to_csv(
            os.path.join(WORK, f"sc_pseudobulkDEG_{ct_name}.csv"), index=False)
    if pb_rows:
        pbd = pd.DataFrame(pb_rows)
        pbd.to_csv(os.path.join(WORK, "sc_pseudobulk_featuregenes.csv"), index=False)
        piv = pbd.pivot_table(index="gene", columns="celltype", values="log2FC_pb")
        fig, ax = plt.subplots(figsize=(max(5, 1.1 * piv.shape[1] + 2), 3.6))
        sns.heatmap(piv, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax,
                    annot_kws={"size": 7}, cbar_kws={"label": "pseudobulk log2FC"})
        ax.set_title("Feature genes: pseudobulk log2FC (Surgery vs Control, per celltype)")
        fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_pseudobulk_feature.png"), dpi=150)
        plt.close(fig)
        log("[pseudobulk] 特征基因 log2FC:\n" + piv.round(2).to_string())

    # ---------- 少突胶质谱系再聚类 ----------
    ol_labels = [c for c in ad.obs.celltype.unique()
                 if c.startswith("Oligodendrocyte") or c.startswith("OPC")]
    ol = ad[ad.obs.celltype.isin(ol_labels)].copy()
    if ol.n_obs > 200:
        sc.pp.neighbors(ol, n_neighbors=15, n_pcs=30, use_rep=rep, random_state=SEED)
        sc.tl.umap(ol, random_state=SEED)
        sc.tl.leiden(ol, resolution=0.5, key_added="ol_sub", random_state=SEED,
                     flavor="igraph", n_iterations=2, directed=False)
        # 成熟度打分 (OPC -> 成熟髓鞘)
        m_opc = [g for g in ["Pdgfra","Cspg4","Lhfpl3","Sox10"] if g in ol.raw.var_names]
        m_mye = [g for g in ["Mbp","Plp1","Mog","Mal","Opalin","Cnp"] if g in ol.raw.var_names]
        sc.tl.score_genes(ol, m_opc, score_name="OPC_score", use_raw=True)
        sc.tl.score_genes(ol, m_mye, score_name="Myelin_score", use_raw=True)
        fig, axes = plt.subplots(1, 4, figsize=(23, 5))
        sc.pl.umap(ol, color="ol_sub", ax=axes[0], show=False, legend_loc="on data", title="OL subclusters")
        sc.pl.umap(ol, color="Myelin_score", ax=axes[1], show=False, title="Myelin score")
        sc.pl.umap(ol, color="OPC_score", ax=axes[2], show=False, title="OPC score")
        sc.pl.umap(ol, color="group", ax=axes[3], show=False, title="Group")
        fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_ol_subclusters.png"), dpi=150)
        plt.close(fig)
        # 特征基因沿 OL 亚群
        fig, axes = plt.subplots(1, len(present_feat), figsize=(4 * len(present_feat), 4))
        for i, g in enumerate(present_feat):
            sc.pl.umap(ol, color=g, ax=axes[i], show=False, use_raw=True, title=g, cmap="Reds")
        fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_ol_featuregenes.png"), dpi=150)
        plt.close(fig)
        olsub = pd.crosstab(ol.obs.ol_sub, ol.obs.group)
        olsub.to_csv(os.path.join(WORK, "sc_ol_subcluster_composition.csv"))
        log("[OL subclusters]\n" + olsub.to_string())

        # DPT 拟时序 (自 OPC 样簇起)
        try:
            sc.tl.diffmap(ol, random_state=SEED)
            root = ol.obs.groupby("ol_sub")["OPC_score"].mean().idxmax()
            ol.uns["iroot"] = int(np.where(ol.obs.ol_sub.values == root)[0][0])
            sc.tl.dpt(ol)
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
            sc.pl.umap(ol, color="dpt_pseudotime", ax=axes[0], show=False, cmap="viridis",
                       title="DPT pseudotime (OL lineage)")
            sc.pl.umap(ol, color="Myelin_score", ax=axes[1], show=False, title="Myelin score")
            fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_ol_dpt.png"), dpi=150)
            plt.close(fig)
            log(f"[DPT] 根簇={root}")
        except Exception as ex:
            log("[DPT] 跳过:", type(ex).__name__, ex)

    # ---------- 小胶质 稳态 ↔ DAM ----------
    mg_labels = [c for c in ad.obs.celltype.unique() if c.startswith("Microglia")]
    mg = ad[ad.obs.celltype.isin(mg_labels)].copy()
    if mg.n_obs > 100:
        dam = [g for g in ["Trem2","Apoe","Cst7","Lpl","Clec7a","Itgax","Spp1","Gpnmb","Cd9"]
               if g in mg.raw.var_names]
        hom = [g for g in ["P2ry12","Cx3cr1","Tmem119","Siglech","Hexb"] if g in mg.raw.var_names]
        sc.tl.score_genes(mg, dam, score_name="DAM_score", use_raw=True)
        sc.tl.score_genes(mg, hom, score_name="Homeo_score", use_raw=True)
        fig, axes = plt.subplots(1, 4, figsize=(23, 5))
        sc.pl.umap(mg, color="DAM_score", ax=axes[0], show=False, title="DAM score")
        sc.pl.umap(mg, color="Homeo_score", ax=axes[1], show=False, title="Homeostatic score")
        sc.pl.umap(mg, color="group", ax=axes[2], show=False, title="Group")
        sc.pl.umap(mg, color="Sla2" if "Sla2" in mg.raw.var_names else "Cx3cr1",
                   ax=axes[3], show=False, use_raw=True, title="Sla2")
        fig.tight_layout(); fig.savefig(os.path.join(WORK, "fig_sc_microglia_DAM.png"), dpi=150)
        plt.close(fig)
        damsc = mg.obs.groupby("group")[["DAM_score", "Homeo_score"]].mean()
        damsc.to_csv(os.path.join(WORK, "sc_microglia_scores.csv"))
        log("[microglia scores]\n" + damsc.round(3).to_string())

    # ---------- 配体-受体轴 (补体-趋化) ----------
    lr_pairs = [("C3", "C3ar1"), ("Cx3cl1", "Cx3cr1"), ("Gas6", "Mertk"),
                ("Spp1", "Cd44"), ("C1qa", "Lrp1"), ("Apoe", "Trem2"), ("Tyrobp", "Trem2")]
    ct_mean = {}
    for k in MARKERS:
        cells = (ad.obs.celltype.str.startswith(k.split("_")[0])).values
        if cells.sum() < 20:
            continue
        ct_mean[k] = np.asarray(ad.raw.X[cells].mean(0)).ravel()   # 全基因 lognorm
    lr_rows = []
    for lig, rec in lr_pairs:
        if lig not in gidx or rec not in gidx:
            continue
        li, ri = gidx[lig], gidx[rec]
        for sk, sv in ct_mean.items():
            for rk, rv in ct_mean.items():
                if sk == rk:
                    continue
                lr_rows.append({"ligand": lig, "receptor": rec, "sender": sk, "receiver": rk,
                                "score": float(sv[li] * rv[ri]),
                                "lig_expr": float(sv[li]), "rec_expr": float(rv[ri])})
    if lr_rows:
        lr = pd.DataFrame(lr_rows).sort_values("score", ascending=False)
        lr.to_csv(os.path.join(WORK, "sc_lr_axis.csv"), index=False)
        top = lr.head(15)
        log("[LR axis] top:\n" + top[["ligand","receptor","sender","receiver","score"]].to_string(index=False))

    # ---------- 保存 ----------
    ad.write(os.path.join(WORK, "GSE267933_processed.h5ad"))

    # ---------- 摘要 ----------
    with open(os.path.join(WORK, "step5_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Step 5-6 摘要 — GSE267933 单细胞定位与深入分析\n\n")
        f.write(f"- QC 后: {ad.n_obs} 细胞 × {ad.n_vars} HVG (原始 {n0} 细胞)\n")
        f.write(f"- 整合: {rep} (Harmony by sample)\n")
        f.write(f"- Leiden {ad.obs.leiden.nunique()} 簇 -> {ad.obs.celltype.nunique()} 种细胞类型\n\n")
        f.write("## 细胞类型分布\n" + ad.obs.celltype.value_counts().to_string() + "\n\n")
        f.write("## 细胞比例 (Control vs Surgery)\n" + comp.to_string(index=False) + "\n\n")
        f.write("## 特征基因定位 (pseudobulk log2FC, Surgery vs Control)\n")
        if pb_rows:
            f.write(piv.round(3).to_string() + "\n")
        f.write("\n注: 3v3 设计, 单细胞层仅作定位与方向验证; 统计推断依赖 bulk 层与 pseudobulk(逐样本聚合)。\n")
    log("[Done] Step 5-6 输出写入", WORK)


if __name__ == "__main__":
    main()
