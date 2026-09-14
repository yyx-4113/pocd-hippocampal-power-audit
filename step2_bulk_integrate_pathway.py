"""
Step 2 — bulk 多数据集整合 + 三套通路基因集接入 -> 候选基因
================================================================
输入:
  - GSE276942_gene_expression_filter.csv.gz   (小鼠时序 bulk, 17 样本, 线性尺度, 首列为基因符号)
  - GSE215410_genes.rawcount.anno.xlsx        (小鼠 bulk, 4 样本: C1,C2 vs P1,P2; 基因级 raw counts; 含 Symbol/Entrez 注释)
  - GSE174412_RAW/*.txt.gz                    (小鼠 bulk, 6 样本: CON1-3 vs PND1-3; 首列 Ensembl, log2 归一化)
  (GSE178995 仅 SRA 原始数据, GEO 无处理矩阵 -> 本步不可用, 见 step2_summary)

产出:
  - step2_per_dataset_DEG.csv        每个数据集内 手术组 vs 对照组 DEG
  - step2_meta_DEG.csv               跨数据集 Stouffer 元分析 (Z, p, BH, 方向一致性)
  - step2_combat_merged.csv          ComBat 去批次后的合并表达矩阵 (z 尺度)
  - step2_combat_DEG.csv             ComBat 合并矩阵上的 DEG (交叉验证)
  - step2_pathway_sets.json          三套通路基因集 (实际命中矩阵基因的部分)
  - step2_candidate_genes.csv        【核心产出】通路 ∩ 元分析 候选基因 (喂给 Step 3 双模型)
  - step2_summary.md                 文本摘要
  - fig_step2_*.png                  图

方法说明:
  - bulk POCD/PND 信号细微(见 Step1: BH 显著=0), 故本步以"探索性 + 跨数据集一致性"为主,
    元分析采用 Stouffer 加权 Z(权重=sqrt(n)), 强调"多数据集方向一致"而非单集强显著.
  - ComBat(scanpy.pp.combat) 需保留生物学协变量 group, 否则会把组间差异当批次移除;
    整合前对每个数据集做"基因内 z-score", 消除平台量纲差异.
  - 所有阈值在摘要中显式标注为探索性.
"""
import os, re, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
OUT = os.path.dirname(os.path.abspath(__file__))
np.random.seed(42)


# ============================================================
# 0. 三套通路基因集 (mouse symbols)
# ============================================================
SET_A = [
    # 补体级联
    "C1qa","C1qb","C1qc","C1r","C1ra","C1rb","C1s","C1sa","C1sb","C2","C3","C3ar1",
    "C4a","C4b","C4bpa","C4bpb","C5","C5ar1","C5ar2","C6","C7","C8a","C8b","C8g","C9",
    "Cfb","Cfd","Cfh","Cfi","Cfp","Cr1l","Cr2","Masp1","Masp2","Mbl1","Mbl2","Fcnb",
    "Serping1","Vtn","C1qtnf1","C1qtnf3","C1qtnf5",
    # 吞噬受体/识别
    "Itgam","Itgb2","Itgav","Itgb3","Itgb5","Cd68","Axl","Mertk","Gas6","Pros1",
    "Trem2","Tyrobp","Lrp1","Cd47","Sirpa","Megf10","Jedi1","Elmo1","Dock1","Rac1",
    "Crk","Crkl","Gulp1","Clu","Apoe","Lpl","Sorl1",
    # 突触机器/可塑性
    "Syt1","Syt2","Syt4","Syt11","Syn1","Syn2","Syn3","Syp","Snap25","Snap23",
    "Stx1a","Stx1b","Vamp2","Vamp7","Dlg1","Dlg2","Dlg4","Gria1","Gria2","Gria3",
    "Grin1","Grin2a","Grin2b","Grin2c","Grin3a","Grm1","Grm2","Grm3","Grm5",
    "Homer1","Homer2","Camk2a","Camk2b","Camk4","Bdnf","Ngf","Ntf3","Ntf4",
    "Arc","Egr1","Egr2","Fos","Fosb","Fosl2","Jun","Junb","Jund","Npas4","Creb1",
    "Mtor","Ppp1r1b","Dlgap1","Dlgap2","Dlgap3","Shank1","Shank2","Shank3",
    "Nlgn1","Nlgn2","Nlgn3","Nrxn1","Nrxn2","Nrxn3","Cntn1","Cntn2","Cntn4",
    "Cntn5","Cntn6","Ncam1","Ncam2","Sparcl1","Thbs1","Thbs2","Thbs4","Gpc4",
    "Slc17a6","Slc17a7","Gad1","Gad2",
    # 修剪/免疫-神经元互作
    "Cx3cr1","Cx3cl1","Spp1","Cd200","Cd200r1","Adgrg1","Mfge8","Gpnmb","Serpine2",
    "Plat","Plaur","Plg","Pla2g7","Mmp2","Mmp9","Timp1","Ctsb","Ctsd","Ctss","Ctsz",
    "Ctsc","Ctsl","Lgmn","Cst3","Cst7",
]

SET_B = [
    "Trem2","Tyrobp","Apoe","Cst7","Lpl","Clec7a","Itgax","Spp1","Gpnmb","Cd9","Cd63",
    "Lag3","Ank","Axl","Lilrb4","Msr1","Fabp5","Ctsb","Ctsd","Ctss","Ctsz","Ctsa","Ctsc",
    "Ctsl","B2m","H2-D1","H2-K1","H2-Aa","H2-Ab1","H2-Eb1","H2-DMb1","Cd74","Pld3",
    "Serinc3","Fth1","Ftl1","Lyz2","Gns","Gpx3","Tnf","Il1b","Il6","Il18","Ccl3","Ccl4",
    "Ccl6","Ccl9","Ccl12","Cxcl10","Cxcl16","Nos2","Ptgs2","Socs3","Cd86","Cd83","Cd80",
    "Fcgr1","Fcgr2b","Fcgr3","Tlr2","Tlr4","Tlr7","Tlr9","Nfkbia","Nfkb1","Nfkb2","Rela",
    "Relb","Stat1","Stat3","Stat6","Irf1","Irf3","Irf5","Irf7","Irf8","Spi1","Runx1",
    "Runx2","Cebpb","Mafb","P2ry12","Cx3cr1","Tmem119","Siglech","Sall1","Hexb","Gpr34",
    "Tgfbr1","Tgfb1","Csf1r","Csf1","Il34","Csf2rb","Plcg2","Syk","Inpp5d","Ptpn11",
    "Aif1","C1qa","C1qb","C1qc","C3","Serping1","Gfap","Vim","Lcn2","Steap4","Hmox1",
    "Mt1","Mt2","Gpr56","Gpr84","Ms4a7","Ms4a6c","Cd52","Ctla2a","Ctla2b","Ccl5",
]

SET_C = [
    # Complex II
    "Sdha","Sdhb","Sdhc","Sdhd","Sdhaf1","Sdhaf2","Sdhaf3","Sdhaf4",
    # Complex III
    "Uqcc1","Uqcc2","Uqcc3","Uqcr10","Uqcr11","Uqcrb","Uqcrc1","Uqcrc2","Uqcrfs1",
    "Uqcrh","Uqcrq","Cyc1",
    # Complex IV
    "Cox4i1","Cox4i2","Cox5a","Cox5b","Cox6a1","Cox6a2","Cox6b1","Cox6b2","Cox6c",
    "Cox7a1","Cox7a2","Cox7a2l","Cox7b","Cox7b2","Cox7c","Cox8a","Cox8b","Cox8c",
    "Cox10","Cox11","Cox15","Cox17","Cox20","Sco1","Sco2","Surf1",
    # Complex V
    "Atp5a1","Atp5b","Atp5c1","Atp5d","Atp5e","Atp5f1","Atp5f1a","Atp5f1b","Atp5g1",
    "Atp5g2","Atp5g3","Atp5h","Atp5j","Atp5j2","Atp5k","Atp5l","Atp5o","Atp5s","Atp5md",
    "Atpif1",
    # 生物发生/动力学/线粒体自噬
    "Ppargc1a","Ppargc1b","Tfam","Nrf1","Nfe2l2","Esrra","Esrrb","Tomm20","Tomm22",
    "Tomm40","Tomm70a","Timm8a1","Timm8b","Timm9","Timm10","Timm13","Timm17a","Timm17b",
    "Timm22","Timm23","Timm44","Timm50","Vdac1","Vdac2","Vdac3","Mfn1","Mfn2","Opa1",
    "Dnm1l","Fis1","Pink1","Prkn","Bnip3","Bnip3l","Fundc1","Map1lc3a","Map1lc3b",
    "Sqstm1","Optn","Slc25a4","Slc25a5","Slc25a3","Slc25a1","Slc25a10","Slc25a11",
    "Slc25a12","Slc25a13","Aifm1","Cycs","Apaf1","Casp9","Casp3","Htra2","Oma1","Yme1l1",
    "Spg7","Afg3l2","Clpp","Lonp1","Phb","Phb2","Cs","Aco2","Idh3a","Idh2","Ogdh","Dlst",
    "Sucla2","Suclg1","Fh1","Mdh2","Mpc1","Mpc2","Srebf1","Srebf2","Hk1","Hk2","Pfkm",
    "Pkm","Ldha","Pdk1","Pdk2","Pdk4",
]
# ---- 第四条机制线: 髓鞘/少突胶质 (Step 2 数据驱动发现全表最强信号, 经用户确认纳入) ----
SET_D = [
    # 致密髓鞘结构蛋白
    "Mbp","Plp1","Mog","Mag","Mal","Cnp","Cldn11","Mobp","Opalin","Pmp22","Pmp2",
    "Ncmap","Tspan2","Ermn","Prx","Gjc2","Gjb1","Sirt2","Galc","Ugt8a","Fa2h",
    # 少突胶质谱系/转录因子
    "Sox10","Olig1","Olig2","Myrf","Nkx2-2","St18","Zfp24","Zfp191","Creb3l2",
    "Cspg4","Pdgfra","Cntn2","Srebf1","Srebf2","Elovl1","Elovl7","Cers2","Aspa",
    # 少突胶质富集的脂质/胆固醇合成
    "Hmgcr","Hmgcs1","Idi1","Fdft1","Sqle","Lss","Mvd","Mvk","Dhcr7","Dhcr24",
    "Nsdhl","Msmo1","Ebp","Cyp51",
    # 髓鞘相关/OL 代谢与信号 (含 Step2 全表 top 命中)
    "Bcas1","Myo1d","Sla2","Gltp","Tmem63a","Prr5l","Prr18","Lpar1","Gpr37",
    "Gpr17","Gpr37l1","Adgrg1","Enpp6","Bcas3","Sgk1","Kctd17","Sh3tc2","Fig4",
    "Mtmr2","Ndrg1","Tppp","Wasf3","Tmcc3",
]

# Complex I 用前缀扩展 (Ndufa/b/c/s/v/ab 家族, 成员众多)
NDUF_PREFIX = ("Ndufa", "Ndufab", "Ndufb", "Ndufc", "Ndufs", "Ndufv")

PATHWAY_DEF = {
    "A_complement_synapse_pruning": SET_A,
    "B_DAM_microglia": SET_B,
    "C_mitochondrial_OXPHOS": SET_C,
    "D_myelin_oligodendrocyte": SET_D,
}


# ============================================================
# 1. 载入与标准化
# ============================================================
def log2(x):
    return np.log2(x + 1.0)


def collapse_dups(df, idcol):
    """同名列取均值"""
    return df.groupby(idcol).mean()


def load_gse276942():
    df = pd.read_csv(os.path.join(OUT, "GSE276942_gene_expression_filter.csv.gz"))
    idc = df.columns[0]                      # gene_name (符号为主)
    mat = collapse_dups(df.drop(columns=[idc]).apply(pd.to_numeric, errors="coerce")
                        .assign(_g=df[idc].values), "_g")
    expr = log2(mat)                          # 线性 -> log2
    meta = {
        "CON_1": ("Control", "CON"), "CON_2": ("Control", "CON"),
        "CON_4": ("Control", "CON"), "CON_5": ("Control", "CON"),
        "24H_2": ("Post", "24H"), "24H_3": ("Post", "24H"), "24H_4": ("Post", "24H"),
        "3D_1": ("Post", "3D"), "3D_3": ("Post", "3D"), "3D_4": ("Post", "3D"),
        "7D_1": ("Post", "7D"), "7D_2": ("Post", "7D"), "7D_3": ("Post", "7D"),
        "1M_1": ("Post", "1M"), "1M_2": ("Post", "2M"),
        "1M_3": ("Post", "3M"), "1M_4": ("Post", "4M"),
    }
    grp = pd.Series({c: meta[c][0] for c in expr.columns})
    return expr, grp


def load_gse215410():
    xl = pd.ExcelFile(os.path.join(OUT, "GSE215410_genes.rawcount.anno.xlsx"))
    df = xl.parse(xl.sheet_names[0])
    sym = df["Symbol"].astype(str)
    counts = df[["C1", "C2", "P1", "P2"]].apply(pd.to_numeric, errors="coerce")
    counts.index = sym.values
    counts = counts[~counts.index.isin(["nan", "-", ""])]
    counts = collapse_dups(counts, counts.index)
    # CPM -> log2
    cpm = counts / counts.sum(axis=0) * 1e6
    expr = log2(cpm)
    grp = pd.Series({"C1": "Control", "C2": "Control", "P1": "Post", "P2": "Post"})
    # Ensembl -> Symbol 映射 (供 GSE174412 使用)
    ens2sym = dict(zip(df["#id"].astype(str), sym))
    return expr, grp, ens2sym


def load_gse174412(ens2sym):
    d = os.path.join(OUT, "GSE174412_RAW")
    files = sorted(os.listdir(d))
    cols = {}
    for f in files:
        s = f.replace(".txt.gz", "")
        name = s.split("_", 1)[1]                       # PND1 / CON1 ...
        t = pd.read_csv(os.path.join(d, f), sep="\t", index_col=0)
        cols[name] = t.iloc[:, 0]
    expr = pd.DataFrame(cols)                            # index = Ensembl
    expr.index = [ens2sym.get(i, np.nan) for i in expr.index]
    expr = expr[expr.index.notna()]
    expr = expr[~expr.index.isin(["nan", "-", ""])]
    expr = collapse_dups(expr, expr.index)
    grp = pd.Series({c: ("Control" if c.startswith("CON") else "Post")
                     for c in expr.columns})
    return expr, grp


# ============================================================
# 2. 单数据集 DEG (Welch t + BH)
# ============================================================
def one_deg(expr, grp, name):
    g = grp.reindex(expr.columns)
    a = expr.loc[:, (g == "Post").values]
    b = expr.loc[:, (g == "Control").values]
    log2fc = a.mean(1) - b.mean(1)
    t, p = stats.ttest_ind(a.values, b.values, axis=1, equal_var=False)
    p = np.where(np.isfinite(p), np.clip(p, 0, 1), 1.0)
    padj = multipletests(p, method="fdr_bh")[1]
    return pd.DataFrame({"gene": expr.index, "dataset": name,
                         "log2FC": log2fc.values, "p": p, "padj": padj,
                         "n_post": a.shape[1], "n_ctrl": b.shape[1]})


# ============================================================
# 3. Stouffer 元分析
# ============================================================
def meta_analysis(deg_all):
    rows = []
    for gene, sub in deg_all.groupby("gene"):
        sub = sub.dropna(subset=["p", "log2FC"])
        if len(sub) < 2:
            continue
        z = np.sign(sub.log2FC.values) * stats.norm.isf(np.clip(sub.p.values, 1e-300, 1) / 2)
        w = np.sqrt(sub.n_post.values)
        Z = np.sum(w * z) / np.sqrt(np.sum(w ** 2))
        pm = 2 * stats.norm.sf(abs(Z))
        d = np.sign(sub.log2FC.values)
        rows.append({
            "gene": gene, "n_datasets": len(sub), "meta_Z": Z, "meta_p": pm,
            "consistency": max((d > 0).sum(), (d < 0).sum()) / len(d),  # 方向一致比例
            "n_sig": int((sub.p < 0.05).sum()),
            "mean_log2FC": sub.log2FC.mean(),
            **{f"{r.dataset}_log2FC": r.log2FC for _, r in sub.iterrows()},
            **{f"{r.dataset}_p": r.p for _, r in sub.iterrows()},
        })
    meta = pd.DataFrame(rows)
    meta["meta_padj"] = multipletests(meta.meta_p, method="fdr_bh")[1]
    return meta.sort_values("meta_p").reset_index(drop=True)


# ============================================================
# 4. ComBat 合并
# ============================================================
def param_combat(Y, batch, X, eb=True):
    """
    参数化 ComBat (Johnson 2007), 保留生物学协变量 X (含截距 + group).
    Y: G x N (genes x samples); batch: (N,) 标签; X: N x p 设计矩阵.
    位置参数 gamma 做经验贝叶斯收缩; 尺度参数 delta 采用非 EB(稳健).
    返回校正后的 G x N.
    """
    Y = np.asarray(Y, float)
    G, N = Y.shape
    uniq = list(dict.fromkeys(batch))
    idxs = [np.where(np.asarray(batch) == b)[0] for b in uniq]
    n_b = np.array([len(ix) for ix in idxs], float)

    B = np.linalg.lstsq(X, Y.T, rcond=None)[0]        # p x G
    Yhat = (X @ B).T                                   # G x N
    df = max(N - X.shape[1], 1)
    sigma = np.sqrt(((Y - Yhat) ** 2).sum(1) / df)     # G
    sigma[sigma == 0] = 1e-8
    Z = (Y - Yhat) / sigma[:, None]

    gamma = np.zeros((G, len(uniq)))
    delta2 = np.zeros((G, len(uniq)))
    for i, ix in enumerate(idxs):
        gamma[:, i] = Z[:, ix].mean(1)
        if n_b[i] > 1:
            delta2[:, i] = ((Z[:, ix] - gamma[:, i:i + 1]) ** 2).sum(1) / (n_b[i] - 1)
        else:
            delta2[:, i] = 1.0

    gamma_star = gamma.copy()
    if eb:
        for i in range(len(uniq)):
            g_i = gamma[:, i]
            g_bar = g_i.mean()
            tau2 = max(g_i.var(ddof=1), 0.0)
            d_bar = delta2[:, i].mean()
            den = n_b[i] * tau2 + d_bar
            gamma_star[:, i] = (n_b[i] * tau2 * g_i + d_bar * g_bar) / (den if den else 1.0)

    Yc = Y.copy()
    for i, ix in enumerate(idxs):
        Yc[:, ix] = (Z[:, ix] - gamma_star[:, i:i + 1]) / np.sqrt(delta2[:, i:i + 1]) \
                    * sigma[:, None] + Yhat[:, ix]
    return Yc


def combat_merge(exps, grps):
    common = None
    for e in exps.values():
        common = set(e.index) if common is None else (common & set(e.index))
    common = sorted(common)
    print(f"[ComBat] 三数据集共同基因 = {len(common)}")
    blocks, obs = [], []
    for name, e in exps.items():
        sub = e.loc[common]
        z = sub.sub(sub.mean(1), axis=0).div(sub.std(1).replace(0, np.nan), axis=0)
        blocks.append(z.fillna(0.0).T)          # samples x genes
        obs += [(name, grps[name].reindex(sub.columns)[c]) for c in sub.columns]
    X = pd.concat(blocks)                       # samples x genes
    obs = pd.DataFrame(obs, columns=["dataset", "group"], index=X.index)
    # 设计矩阵: 截距 + group(Post)
    D = np.column_stack([np.ones(len(obs)), (obs.group == "Post").astype(float).values])
    Yt = X.values.T                             # genes x samples
    corrected, used = X.copy(), "manual"
    try:
        from inmoose.pycombat import pycombat
        cov = pd.DataFrame({"group": (obs.group == "Post").astype(int).values}, index=obs.index)
        Yc = pycombat(pd.DataFrame(Yt, index=X.columns, columns=X.index),
                      obs.dataset.tolist(), covar_mod=cov)
        corrected = pd.DataFrame(Yc.values.T, index=X.index, columns=X.columns)
        used = "inmoose.pycombat(covar_mod=group)"
    except Exception as ex:
        print("[ComBat] inmoose 不可用, 使用内置参数化 ComBat:", type(ex).__name__, ex)
        Yc = param_combat(Yt, obs.dataset.values, D, eb=True)
        corrected = pd.DataFrame(Yc.T, index=X.index, columns=X.columns)
        used = "内置参数化 ComBat (位置 EB, 保留 group)"
    # 合并矩阵上的 DEG (样本 x 基因 -> 按列/基因检验 -> axis=0)
    a = corrected[(obs.group == "Post").values]
    b = corrected[(obs.group == "Control").values]
    t, p = stats.ttest_ind(a.values, b.values, axis=0, equal_var=False)
    p = np.where(np.isfinite(p), np.clip(p, 0, 1), 1.0)
    cdeg = pd.DataFrame({"gene": corrected.columns,
                         "log2FC_z": a.mean().values - b.mean().values,
                         "p": p, "padj": multipletests(p, method="fdr_bh")[1]})
    return corrected, obs, cdeg, used


# ============================================================
# MAIN
# ============================================================
def main():
    summaries = []

    # --- 载入 ---
    e1, g1 = load_gse276942()
    e2, g2, ens2sym = load_gse215410()
    e3, g3 = load_gse174412(ens2sym)
    print(f"[load] GSE276942 {e1.shape} | GSE215410 {e2.shape} | GSE174412 {e3.shape}")

    # --- 通路集命中情况 ---
    universe = set(e1.index) | set(e2.index) | set(e3.index)
    nduf = sorted(g for g in universe if g.startswith(NDUF_PREFIX))
    pathway_hit = {}
    for k, lst in PATHWAY_DEF.items():
        s = sorted(set(lst) & universe)
        if k.startswith("C_"):
            s = sorted(set(s) | set(nduf))
        pathway_hit[k] = s
    pu = sorted(set().union(*pathway_hit.values()))
    print("[pathway] 命中: " + ", ".join(f"{k}={len(v)}" for k, v in pathway_hit.items()) +
          f" | 并集={len(pu)}")

    # --- 单数据集 DEG ---
    deg_all = pd.concat([one_deg(e1, g1, "GSE276942"),
                         one_deg(e2, g2, "GSE215410"),
                         one_deg(e3, g3, "GSE174412")], ignore_index=True)
    deg_all.to_csv(os.path.join(OUT, "step2_per_dataset_DEG.csv"), index=False)
    for name in ["GSE276942", "GSE215410", "GSE174412"]:
        s = deg_all[deg_all.dataset == name]
        summaries.append(f"- {name}: n_post vs n_ctrl = "
                         f"{int(s.n_post.iloc[0])}v{int(s.n_ctrl.iloc[0])}; "
                         f"BH<0.05 = {int((s.padj<0.05).sum())}, p<0.05 = {int((s.p<0.05).sum())}")

    # --- 元分析 ---
    meta = meta_analysis(deg_all)
    meta.to_csv(os.path.join(OUT, "step2_meta_DEG.csv"), index=False)
    print(f"[meta] 基因数={len(meta)}; meta_p<0.05 = {int((meta.meta_p<0.05).sum())}; "
          f"meta_padj<0.05 = {int((meta.meta_padj<0.05).sum())}")

    # --- ComBat ---
    corr, obs, cdeg, used = combat_merge(
        {"GSE276942": e1, "GSE215410": e2, "GSE174412": e3},
        {"GSE276942": g1, "GSE215410": g2, "GSE174412": g3})
    corr.to_csv(os.path.join(OUT, "step2_combat_merged.csv"))
    cdeg.to_csv(os.path.join(OUT, "step2_combat_DEG.csv"), index=False)
    print(f"[ComBat] {used}; BH<0.05 = {int((cdeg.padj<0.05).sum())}, "
          f"p<0.05 = {int((cdeg.p<0.05).sum())}")

    # --- 候选基因: 通路 ∩ (元分析显著 或 跨集一致) ---
    pw = meta[meta.gene.isin(pu)].copy()
    pw["in_pathway"] = pw.gene.map(lambda g: "+".join(
        k.split("_")[0] for k, v in pathway_hit.items() if g in v))
    # 分层
    tier1 = pw[(pw.meta_p < 0.05) & (pw.n_datasets >= 2) & (pw.consistency >= 0.99)]
    tier2 = pw[(pw.meta_p < 0.05) & (pw.n_datasets >= 2) & (pw.consistency >= 0.66)]
    tier3 = pw[(pw.n_sig >= 2) | ((pw.meta_p < 0.05) & (pw.n_datasets >= 2))]
    cand = pw[(pw.meta_p < 0.05) & (pw.n_datasets >= 2) &
              ((pw.consistency >= 0.66) | (pw.n_sig >= 2))].copy()
    cand = cand.sort_values("meta_p")
    cand["tier"] = np.where(cand.index.isin(tier1.index), "T1_meta&全一致",
                            np.where(cand.index.isin(tier2.index), "T2_meta&多数一致",
                                     "T3_多集名义显著"))
    cand.to_csv(os.path.join(OUT, "step2_candidate_genes.csv"), index=False)
    with open(os.path.join(OUT, "step2_pathway_sets.json"), "w", encoding="utf-8") as f:
        json.dump(pathway_hit, f, ensure_ascii=False, indent=1)
    print(f"[candidate] 通路∩元分析候选基因 = {len(cand)} "
          f"(T1={int((cand.tier.str.startswith('T1')).sum())}, "
          f"T2={int((cand.tier.str.startswith('T2')).sum())}, "
          f"T3={int((cand.tier.str.startswith('T3')).sum())})")
    print("  top20:", ", ".join(cand.head(20).gene.tolist()))

    # --- 图: 候选基因跨数据集 log2FC 热图 ---
    if len(cand) >= 3:
        top = cand.head(min(30, len(cand)))
        cols = [f"{d}_log2FC" for d in ["GSE276942", "GSE215410", "GSE174412"]
                if f"{d}_log2FC" in top.columns]
        M = top.set_index("gene")[cols]
        M.columns = [c.replace("_log2FC", "") for c in cols]
        fig, ax = plt.subplots(figsize=(5.2, max(3, 0.28 * len(M) + 1.5)))
        sns.heatmap(M, cmap="RdBu_r", center=0, ax=ax, annot=True, fmt=".2f",
                    annot_kws={"size": 6}, cbar_kws={"label": "log2FC"})
        ax.set_title(f"Top {len(M)} candidate genes: per-dataset log2FC")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "fig_step2_candidate_log2FC.png"), dpi=150)
        plt.close(fig)

    # --- 图: 通路成员 meta_Z 火山 ---
    fig, ax = plt.subplots(figsize=(7, 5))
    base = meta[meta.n_datasets >= 2]
    ax.scatter(base.meta_Z, -np.log10(base.meta_p), s=6, c="lightgrey", label="all genes")
    hit = base[base.gene.isin(pu)]
    cA = hit[hit.gene.isin(pathway_hit["A_complement_synapse_pruning"])]
    cB = hit[hit.gene.isin(pathway_hit["B_DAM_microglia"])]
    cC = hit[hit.gene.isin(pathway_hit["C_mitochondrial_OXPHOS"])]
    cD = hit[hit.gene.isin(pathway_hit["D_myelin_oligodendrocyte"])]
    ax.scatter(cA.meta_Z, -np.log10(cA.meta_p), s=14, c="#d62728", label="A complement/pruning")
    ax.scatter(cB.meta_Z, -np.log10(cB.meta_p), s=14, c="#1f77b4", label="B DAM microglia")
    ax.scatter(cC.meta_Z, -np.log10(cC.meta_p), s=14, c="#2ca02c", label="C mito OXPHOS")
    ax.scatter(cD.meta_Z, -np.log10(cD.meta_p), s=14, c="#8c564b", label="D myelin/OL")
    ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.8)
    ax.set_xlabel("Stouffer meta Z"); ax.set_ylabel("-log10(meta p)")
    ax.set_title("Pathway genes (A/B/C/D) in cross-dataset meta-analysis")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_step2_meta_volcano.png"), dpi=150)
    plt.close(fig)

    # --- 摘要 ---
    with open(os.path.join(OUT, "step2_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Step 2 摘要 — bulk 整合 + 通路基因集接入\n\n")
        f.write("## 数据可用性修正\n")
        f.write("- **GSE178995 不可用**: 其 4 样本 `Sample_type=SRA`, GEO 仅有注释文件 "
                "`GSE178995_S_C.anno.txt.gz`, 无处理表达矩阵(原始 fastq 在 SRA SRP325799)。"
                "方案原计划三集合并, 实际整合 **GSE276942 + GSE215410 + GSE174412**。\n")
        f.write("- 另修正: GSE215410 样本为 C1,C2(对照) vs P1,P2(POCD) = **2v2**; "
                "GSE174412 = CON1-3 vs PND1-3 = 3v3。\n\n")
        f.write("## 单数据集 DEG\n" + "\n".join(summaries) + "\n\n")
        f.write(f"## 跨数据集元分析 (Stouffer 加权 Z)\n")
        f.write(f"- 可比基因(≥2 数据集): {len(meta)}\n")
        f.write(f"- meta_p<0.05: {int((meta.meta_p<0.05).sum())}; "
                f"meta_padj(BH)<0.05: {int((meta.meta_padj<0.05).sum())}\n\n")
        f.write(f"## ComBat 合并\n- 方法: {used}\n")
        f.write(f"- 共同基因: {corr.shape[1]}; 样本: {corr.shape[0]}\n")
        f.write(f"- 合并矩阵 DEG: BH<0.05={int((cdeg.padj<0.05).sum())}, "
                f"p<0.05={int((cdeg.p<0.05).sum())}\n\n")
        f.write("## 四套通路基因集 (命中本步表达矩阵)\n")
        for k, v in pathway_hit.items():
            f.write(f"- {k}: {len(v)}\n")
        f.write(f"- 并集: {len(pu)}\n\n")
        f.write(f"## 候选基因 (通路 ∩ 元分析)\n- 共 {len(cand)} 个\n")
        for t in ["T1", "T2", "T3"]:
            g = cand[cand.tier.str.startswith(t)].gene.tolist()
            f.write(f"- {t} (n={len(g)}): {', '.join(g[:40])}\n")
    print("[Done] Step 2 输出已写入", OUT)


if __name__ == "__main__":
    main()
