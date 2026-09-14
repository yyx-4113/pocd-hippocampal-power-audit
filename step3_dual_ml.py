"""
Step 3 — 机器学习双模型锁靶 + 表达/ROC 双重确认
================================================================
模型一: LASSO logistic (L1), 10 折 CV, lambda.1se 规则
模型二: Boruta (RandomForest, ntree=1000) + MeanDecreaseGini 排序
交集 = 特征基因; 若为空 -> 取"两法均入 Top10"。
加分: XGBoost 重要性 (+ SHAP 若可用)。

训练集:
  - 发现集 GSE276942 (n=17, 4 对照 / 13 术后)
  - 整合矩阵 ComBat (n=27, 9 对照 / 18 术后) —— 提高样本量的敏感性分析
ROC:
  - 多基因 logistic: 发现集 10 折 CV AUC(mean±SD); 验证集 GSE215410 / GSE174412
  - 单基因 ROC AUC (方向取 max(auc,1-auc))

输入: step2_candidate_pool.csv (+ step2_combat_merged.csv)
输出: step3_feature_genes.csv / step3_lasso_cv.csv / step3_boruta_importance.csv /
      step3_roc_summary.csv / fig_step3_*.png / step3_summary.md
注意: 样本数小(n=17~27), 结果均为**探索性**, 摘要中显式标注。
"""
import os, warnings, importlib.util
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
OUT = os.path.dirname(os.path.abspath(__file__))
SEED = 42


def _load_step2():
    spec = importlib.util.spec_from_file_location(
        "s2", os.path.join(OUT, "step2_bulk_integrate_pathway.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------------- LASSO (lambda.1se) ----------------
def lasso_1se(X, y, Cs=np.logspace(-3, 1, 40), cv=10, seed=SEED):
    n_min = int(np.bincount(y).min())
    cv = min(cv, n_min)
    skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=seed)
    means, ses = [], []
    for C in Cs:
        pipe = make_pipeline(StandardScaler(), LogisticRegression(
            penalty="l1", solver="liblinear", C=C, max_iter=5000))
        s = cross_val_score(pipe, X, y, cv=skf, scoring="roc_auc")
        means.append(s.mean()); ses.append(s.std(ddof=1) / np.sqrt(cv))
    means, ses = np.array(means), np.array(ses)
    best = int(means.argmax())
    thr = means[best] - ses[best]
    ok = np.where(means >= thr)[0]
    c_1se = Cs[ok[0]] if len(ok) else Cs[best]
    pipe = make_pipeline(StandardScaler(), LogisticRegression(
        penalty="l1", solver="liblinear", C=c_1se, max_iter=5000)).fit(X, y)
    coef = pipe.named_steps["logisticregression"].coef_[0]
    sel = np.where(np.abs(coef) > 1e-10)[0]
    if len(sel) == 0:                       # 过于稀疏 -> 退化到最优 C
        c_1se = Cs[best]
        pipe = make_pipeline(StandardScaler(), LogisticRegression(
            penalty="l1", solver="liblinear", C=c_1se, max_iter=5000)).fit(X, y)
        coef = pipe.named_steps["logisticregression"].coef_[0]
        sel = np.where(np.abs(coef) > 1e-10)[0]
    cvdf = pd.DataFrame({"C": Cs, "lambda": 1 / Cs, "cv_auc": means, "se": ses})
    return c_1se, sel, coef, cvdf


# ---------------- LASSO bootstrap 稳定性选择 ----------------
def lasso_stability(X, y, C, B=200, seed=SEED):
    """分层 bootstrap 重抽样 B 次, 统计每个特征被 LASSO 选中的频率 (0~1)。"""
    rng = np.random.default_rng(seed)
    n, p = X.shape
    freq = np.zeros(p)
    ok = 0
    for _ in range(B):
        ix = []
        for cls in np.unique(y):
            ci = np.where(y == cls)[0]
            ix += list(rng.choice(ci, size=len(ci), replace=True))
        ix = np.array(ix)
        if len(np.unique(y[ix])) < 2:
            continue
        try:
            pipe = make_pipeline(StandardScaler(), LogisticRegression(
                penalty="l1", solver="liblinear", C=C, max_iter=5000)).fit(X[ix], y[ix])
            coef = pipe.named_steps["logisticregression"].coef_[0]
            freq += (np.abs(coef) > 1e-10)
            ok += 1
        except Exception:
            continue
    return freq / max(ok, 1)


# ---------------- Boruta ----------------
def run_boruta(X, y, seed=SEED, max_iter=500):
    from boruta import BorutaPy
    rf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, random_state=seed,
                                class_weight="balanced")
    bor = BorutaPy(rf, n_estimators="auto", max_iter=max_iter, random_state=seed, verbose=0)
    bor.fit(X, y)
    return bor.support_, bor.ranking_, bor


def main():
    s2 = _load_step2()
    pool = pd.read_csv(os.path.join(OUT, "step2_candidate_pool.csv"))
    genes = [g for g in pool.gene.tolist()]
    print(f"[pool] {len(genes)} 候选基因")

    # ---- 载入各数据集表达 ----
    e1, g1 = s2.load_gse276942()
    e2, g2, ens = s2.load_gse215410()
    e3, g3 = s2.load_gse174412(ens)
    merg = pd.read_csv(os.path.join(OUT, "step2_combat_merged.csv"), index_col=0)

    def build(expr, grp, gset):
        gs = [g for g in gset if g in expr.index]
        X = expr.loc[gs].T
        y = np.where(grp.reindex(X.index).values == "Post", 1, 0)
        return X, y, gs

    def grp_from_index(idx):
        post = [("CON_" not in i) and (not i.startswith("CON"))
                and i not in ("C1", "C2") for i in idx]
        return pd.Series(np.where(post, "Post", "Control"), index=idx)

    X1, y1, gs1 = build(e1, g1, genes)
    X2, y2, gs2 = build(e2, g2, genes)
    X3, y3, gs3 = build(e3, g3, genes)
    Xm, ym, gsm = build(merg.T, grp_from_index(merg.index), genes)
    print(f"[matrices] GSE276942 {X1.shape} | merged {Xm.shape} | "
          f"GSE215410 {X2.shape} | GSE174412 {X3.shape}")

    results, lasso_cvs, feat_tables = {}, {}, {}

    for tag, X, y in [("discovery_GSE276942", X1, y1), ("integrated_ComBat", Xm, ym)]:
        print(f"\n===== {tag}: n={len(y)}, 特征={X.shape[1]} =====")
        feats = list(X.columns)
        # LASSO
        c, sel, coef, cvdf = lasso_1se(X.values, y)
        lasso_set = set(np.array(feats)[sel])
        lasso_cvs[tag] = cvdf.assign(c_1se=c)
        print(f"[LASSO] C*={c:.3f} (λ={1/c:.2f}); 选中 {len(lasso_set)} 基因; "
              f"CV-AUC={cvdf.cv_auc.max():.3f}")
        # RF importance
        rf = RandomForestClassifier(n_estimators=1000, n_jobs=-1, random_state=SEED,
                                    class_weight="balanced").fit(X.values, y)
        imp = pd.Series(rf.feature_importances_, index=feats).sort_values(ascending=False)
        top10 = set(imp.head(10).index)
        # Boruta
        try:
            supp, rank, bor = run_boruta(X.values, y)
            boruta_set = set(np.array(feats)[supp])
            bor_status = "BorutaPy"
        except Exception as ex:
            print("[Boruta] BorutaPy 失败, 用 RF-Top10 近似:", type(ex).__name__)
            boruta_set, bor_status = top10, "RF-Top10(近似)"
        print(f"[Boruta] 确认 {len(boruta_set)} 基因 ({bor_status})")
        inter = lasso_set & boruta_set
        if len(inter) == 0:
            inter = lasso_set & top10
            note = "交集空->LASSO∩RF-Top10"
        else:
            note = "LASSO∩Boruta"
        print(f"[特征基因] {note} = {len(inter)}: {sorted(inter)}")
        # Bootstrap 稳定性
        stab = lasso_stability(X.values, y, c, B=200)
        stable = set(np.array(feats)[stab >= 0.6])
        print(f"[stability] LASSO 选中频率>=0.6: {len(stable)} 基因 -> {sorted(stable)}")
        coef_map = dict(zip(np.array(feats), coef))
        feat_tables[tag] = pd.DataFrame({
            "gene": feats, "lasso_selected": [f in lasso_set for f in feats],
            "lasso_coef": [coef_map.get(f, 0.0) for f in feats],
            "lasso_stability": stab,
            "rf_gini": imp.reindex(feats).values,
            "rf_rank": imp.rank(ascending=False).reindex(feats).values,
            "boruta_selected": [f in boruta_set for f in feats],
            "is_feature": [f in inter for f in feats],
        })
        results[tag] = {"features": sorted(inter), "lasso": lasso_set,
                        "boruta": boruta_set, "stable": stable, "note": note, "c": c,
                        "bor_status": bor_status}

    # ---- 最终特征基因 + 稳定性/重要性汇总表 ----
    f_disc = set(results["discovery_GSE276942"]["features"])
    f_int = set(results["integrated_ComBat"]["features"])
    consensus = sorted(f_disc & f_int)
    final = consensus if len(consensus) >= 3 else sorted(f_disc | f_int)
    src = "两训练集共识" if len(consensus) >= 3 else "两训练集并集"
    print(f"\n[FINAL] 特征基因 {len(final)} ({src}): {final}")
    print(f"[注意] 两训练集特征交集={len(consensus)} -> "
          f"{'稳定' if len(consensus)>=3 else '不稳定(小样本), 特征选择需谨慎解读'}")
    robust_core = sorted(results["integrated_ComBat"]["stable"] |
                         results["discovery_GSE276942"]["stable"])
    print(f"[稳健核心] LASSO bootstrap 频率>=0.6 的基因: {robust_core}")

    sd = feat_tables["discovery_GSE276942"].set_index("gene")
    si = feat_tables["integrated_ComBat"].set_index("gene")
    allg = sorted(set(f_disc) | set(f_int) |
                  results["discovery_GSE276942"]["stable"] |
                  results["integrated_ComBat"]["stable"])
    rows = []
    for g in allg:
        rows.append({
            "gene": g,
            "is_final": g in final,
            "robust_core": g in robust_core,
            "discovery_feature": g in f_disc, "integrated_feature": g in f_int,
            "discovery_lasso_stab": round(float(sd.lasso_stability.get(g, np.nan)), 3),
            "integrated_lasso_stab": round(float(si.lasso_stability.get(g, np.nan)), 3),
            "discovery_rf_rank": int(sd.rf_rank.get(g, -1)),
            "integrated_rf_rank": int(si.rf_rank.get(g, -1)),
            "discovery_boruta": bool(sd.boruta_selected.get(g, False)),
            "integrated_boruta": bool(si.boruta_selected.get(g, False)),
        })
    ftab = pd.DataFrame(rows).sort_values(
        ["robust_core", "is_final", "integrated_lasso_stab"], ascending=[False, False, False])
    ftab.to_csv(os.path.join(OUT, "step3_feature_genes.csv"), index=False)
    for tag, t in feat_tables.items():
        t.to_csv(os.path.join(OUT, f"step3_features_{tag}.csv"), index=False)
    for tag, c in lasso_cvs.items():
        c.to_csv(os.path.join(OUT, f"step3_lasso_cv_{tag}.csv"), index=False)

    # ---- XGBoost 重要性 (加分) ----
    xgb_imp = None
    try:
        from xgboost import XGBClassifier
        xg = XGBClassifier(n_estimators=300, max_depth=2, learning_rate=0.05,
                           subsample=0.9, colsample_bytree=0.9, reg_lambda=1.0,
                           eval_metric="logloss", random_state=SEED).fit(X1.values, y1)
        xgb_imp = pd.Series(xg.feature_importances_, index=gs1).sort_values(ascending=False)
        xgb_imp.to_csv(os.path.join(OUT, "step3_xgboost_importance.csv"))
        print("[XGBoost] gain 重要性 top10:", ", ".join(xgb_imp.head(10).index))
    except Exception as ex:
        print("[XGBoost] 跳过:", type(ex).__name__, ex)

    # ---- ROC ----
    roc_rows = []
    skf = StratifiedKFold(n_splits=min(5, int(np.bincount(y1).min())), shuffle=True,
                          random_state=SEED)
    if len(final) >= 1:
        for tag, X, y, gs in [("discovery_GSE276942", X1, y1, gs1),
                              ("integrated_ComBat", Xm, ym, gsm),
                              ("validation_GSE215410", X2, y2, gs2),
                              ("validation_GSE174412", X3, y3, gs3)]:
            present = [g for g in final if g in X.columns]
            if not present:
                continue
            # 多基因 logistic
            Xs = X[present]
            if y.min() != y.max():
                if tag.startswith("discovery"):
                    s = cross_val_score(make_pipeline(StandardScaler(), LogisticRegression(
                        max_iter=5000)), Xs, y, cv=skf, scoring="roc_auc")
                    roc_rows.append({"dataset": tag, "type": "multi_gene_logistic_CV",
                                     "n": len(y), "AUC": s.mean(), "SD": s.std(ddof=1)})
                else:
                    pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)).fit(Xs, y)
                    auc = roc_auc_score(y, pipe.predict_proba(Xs)[:, 1])
                    roc_rows.append({"dataset": tag, "type": "multi_gene_logistic_fit",
                                     "n": len(y), "AUC": auc, "SD": np.nan})
                # 单基因
                for g in present:
                    v = Xs[g].values
                    if np.std(v) == 0:
                        continue
                    a = roc_auc_score(y, v)
                    roc_rows.append({"dataset": tag, "type": "single_gene", "gene": g,
                                     "n": len(y), "AUC": max(a, 1 - a), "SD": np.nan})
    roc = pd.DataFrame(roc_rows)
    roc.to_csv(os.path.join(OUT, "step3_roc_summary.csv"), index=False)
    print("\n[ROC] 多基因:")
    print(roc[roc.type.str.startswith("multi")].to_string(index=False))
    if len(roc[roc.type == "single_gene"]):
        piv = roc[roc.type == "single_gene"].pivot_table(
            index="gene", columns="dataset", values="AUC")
        print("\n[ROC] 单基因 AUC (行=特征基因):")
        print(piv.round(3).to_string())

    # ---- 重复交叉验证 (更稳 AUC) + 跨数据集 train->test (更诚实的泛化评估) ----
    from sklearn.model_selection import RepeatedStratifiedKFold
    rep, cross_rows = [], []
    rskf = RepeatedStratifiedKFold(n_splits=min(5, int(np.bincount(y1).min())),
                                   n_repeats=20, random_state=SEED)
    for tag, X, y in [("discovery_GSE276942", X1, y1), ("integrated_ComBat", Xm, ym)]:
        pr = [g for g in final if g in X.columns]
        if pr and 0 < y.sum() < len(y):
            s = cross_val_score(make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)),
                                X[pr], y, cv=rskf, scoring="roc_auc")
            rep.append({"dataset": tag, "n": len(y), "n_feat": len(pr),
                        "repeatedCV_AUC_mean": round(s.mean(), 3),
                        "repeatedCV_AUC_sd": round(s.std(ddof=1), 3)})
    for tr_tag, Xtr, ytr in [("discovery_GSE276942", X1, y1), ("integrated_ComBat", Xm, ym)]:
        for te_tag, Xte, yte in [("GSE215410", X2, y2), ("GSE174412", X3, y3),
                                 ("discovery_GSE276942", X1, y1), ("integrated_ComBat", Xm, ym)]:
            pr = [g for g in final if g in Xtr.columns and g in Xte.columns]
            if not pr or yte.sum() in (0, len(yte)):
                continue
            p = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)).fit(Xtr[pr], ytr)
            cross_rows.append({"train": tr_tag, "test": te_tag, "n_test": len(yte),
                               "n_feat": len(pr),
                               "AUC": round(roc_auc_score(yte, p.predict_proba(Xte[pr])[:, 1]), 3)})
    ev = pd.DataFrame(rep + cross_rows)
    ev.to_csv(os.path.join(OUT, "step3_evaluation.csv"), index=False)
    print("\n[重复CV] 5 折 x 20 次:")
    if rep:
        print(pd.DataFrame(rep).to_string(index=False))
    print("\n[跨数据集 train->test AUC]:")
    if cross_rows:
        print(pd.DataFrame(cross_rows).to_string(index=False))

    # ---- 图 ----
    # 1) LASSO CV 曲线
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (tag, cvdf) in zip(axes, lasso_cvs.items()):
        ax.errorbar(np.log10(cvdf.C), cvdf.cv_auc, yerr=cvdf.se, marker="o", ms=3, lw=1)
        ax.axvline(np.log10(cvdf.c_1se.iloc[0]), color="red", ls="--",
                   label=f"λ.1se (C={cvdf.c_1se.iloc[0]:.3f})")
        ax.set_xlabel("log10(C)"); ax.set_ylabel("10-fold CV AUC")
        ax.set_title(tag); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_step3_lasso_cv.png"), dpi=150)
    plt.close(fig)

    # 2) 特征基因跨数据集表达箱线图
    if len(final) >= 1:
        ncol = min(5, len(final)); nrow = int(np.ceil(len(final) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.4 * ncol, 2.6 * nrow), squeeze=False)
        for i, g in enumerate(final):
            ax = axes[i // ncol][i % ncol]
            if g in X1.columns:
                d = pd.DataFrame({"expr": X1[g].values, "group": np.where(y1 == 1, "Post", "Control")})
                sns.boxplot(data=d, x="group", y="expr", ax=ax, palette=["#4c72b0", "#c44e52"],
                            width=0.6, fliersize=2)
                sns.stripplot(data=d, x="group", y="expr", ax=ax, color="k", size=3, alpha=0.6)
            ax.set_title(g, fontsize=9); ax.set_xlabel(""); ax.set_ylabel("log2 expr", fontsize=8)
        for j in range(len(final), nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        fig.suptitle("Feature genes expression (GSE276942) — exploratory", fontsize=11)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_step3_feature_boxplots.png"), dpi=150)
        plt.close(fig)

    # 3) ROC 曲线 (发现集 CV + 验证集)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    if len(final) >= 2:
        Xs = X1[final]
        proba = cross_val_predict(make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)),
                                  Xs, y1, cv=skf, method="predict_proba")[:, 1]
        fpr, tpr, _ = roc_curve(y1, proba)
        ax.plot(fpr, tpr, color="#c44e52", lw=2,
                label=f"multi-gene (GSE276942, CV AUC={roc_auc_score(y1, proba):.2f})")
        for tag, X, yv, col in [("GSE215410", X2, y2, "#55a868"),
                                ("GSE174412", X3, y3, "#8172b2")]:
            pr = [g for g in final if g in X.columns]
            if len(pr) >= 2 and 0 < yv.sum() < len(yv):
                p = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000)).fit(X[pr], yv)
                sc = p.predict_proba(X[pr])[:, 1]
                fpr, tpr, _ = roc_curve(yv, sc)
                ax.plot(fpr, tpr, lw=1.6, color=col,
                        label=f"{tag} (AUC={roc_auc_score(yv, sc):.2f}, n={len(yv)})")
    ax.plot([0, 1], [0, 1], "grey", ls="--", lw=0.8)
    ax.set_xlabel("1 - Specificity"); ax.set_ylabel("Sensitivity")
    ax.set_title("Multi-gene logistic ROC"); ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_step3_roc.png"), dpi=150)
    plt.close(fig)

    # ---- 摘要 ----
    with open(os.path.join(OUT, "step3_summary.md"), "w", encoding="utf-8") as f:
        f.write("# Step 3 摘要 — 双机器学习锁靶 + ROC\n\n")
        f.write(f"候选池: {len(genes)} 基因 (Step 2: 四机制线 A/B/C/D)\n\n")
        f.write("## 双模型结果\n")
        for tag, r in results.items():
            f.write(f"### {tag} (LASSO C*={r['c']:.3f}, λ*={1/r['c']:.2f}; Boruta={r['bor_status']})\n")
            f.write(f"- LASSO 选中 ({len(r['lasso'])}): {sorted(r['lasso'])}\n")
            f.write(f"- Boruta 确认 ({len(r['boruta'])}): {sorted(r['boruta'])}\n")
            f.write(f"- 交集 ({r['note']}) = {r['features']}\n")
            f.write(f"- LASSO bootstrap 稳定性>=0.6 ({len(r['stable'])}): {sorted(r['stable'])}\n\n")
        f.write(f"## 特征基因汇总 ({src}, n={len(final)})\n{final}\n")
        f.write(f"- **稳健核心** (LASSO bootstrap 频率>=0.6): {robust_core}\n")
        f.write(f"- 两训练集特征交集 = {len(consensus)} 个 -> "
                f"{'稳定' if len(consensus)>=3 else '**不稳定**(小样本, 特征选择须谨慎)'}\n\n")
        f.write("## 特征基因稳定性/重要性 (step3_feature_genes.csv)\n")
        f.write(ftab.to_string(index=False) + "\n\n")
        f.write("## 重复交叉验证 (5 折 x 20 次)\n")
        if rep:
            f.write(pd.DataFrame(rep).to_string(index=False) + "\n")
        f.write("\n## 跨数据集 train -> test AUC\n")
        if cross_rows:
            f.write(pd.DataFrame(cross_rows).to_string(index=False) + "\n")
        f.write("\n## 单基因 ROC AUC\n")
        if len(roc[roc.type == "single_gene"]):
            f.write(roc[roc.type == "single_gene"].pivot_table(
                index="gene", columns="dataset", values="AUC").round(3).to_string() + "\n")
        f.write("\n## ⚠️ 关键局限 (必须写入稿件 Limitation)\n")
        f.write("1. **样本量极小**: 发现集 n=17, 整合矩阵 n=27, 验证集 n=4/6 -> 所有 ML 结果"
                "均为**探索性**, AUC 的高值部分来自过拟合。\n")
        f.write("2. **选择偏倚/双重浸用(double dipping)**: 候选基因本身由这些数据集的 DEG/元分析"
                "选出, 再用同一批数据训练+评估, AUC 必然偏高; 真实泛化须依赖**独立外部队列**"
                "(本方案 Step 8 临床血样或未来数据集)。\n")
        f.write("3. **特征选择不稳定**: 两个训练集选出的特征集交叠很小, 说明在当前样本量下"
                "特征选择不可复现; 报告应以\"机制方向\"而非\"精确基因列表\"为主。\n")
    print("[Done] Step 3 输出已写入", OUT)


if __name__ == "__main__":
    main()
