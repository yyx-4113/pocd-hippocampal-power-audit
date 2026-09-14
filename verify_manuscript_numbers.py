"""
稿件数值审计 — 逐条核对 Appendix 审计追踪表与正文引用的数字
================================================================
用途: 把"每个数字都能追溯"从承诺变成可执行检查。本脚本**只读**，不改任何产物；
      每一条检查给出【稿件声称值 / 实测值 / 来源文件 / 通过与否】，末行汇总，
      若有任一条不通过，退出码为 1（可直接接进 CI 或投稿前自检）。

覆盖范围: 稿件 Appendix 审计追踪表全部行 + Abstract / Results 引用的关键数字。

用法: python verify_manuscript_numbers.py
"""
import os
import re
import gzip
import sys

import numpy as np
import pandas as pd
from scipy.stats import nct

HERE = os.path.dirname(os.path.abspath(__file__))
SCO = os.path.join(HERE, "sc_out")
ALPHA, POWER = 0.05, 0.80

RESULTS = []          # (claim, expected, actual, ok, source)


def add(claim, expected, actual, ok, source):
    RESULTS.append((claim, str(expected), str(actual), bool(ok), source))


def num_ok(expected, actual, rtol=1e-3, atol=0.0):
    try:
        return bool(np.isclose(float(expected), float(actual), rtol=rtol, atol=atol))
    except (TypeError, ValueError):
        return False


def rounded_ok(expected, actual, nd):
    """按稿件给出的位数比较（如稿件写 0.14, 实测 0.1412 -> True）。"""
    try:
        return round(float(expected), nd) == round(float(actual), nd)
    except (TypeError, ValueError):
        return False


def grep_md(path, pattern, cast=float, group=1):
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(pattern, txt)
    if not m:
        raise SystemExit("[FAIL] 无法从 %s 解析: %s" % (os.path.basename(path), pattern))
    return cast(m.group(group))


def mde_noncentral(n):
    """非中心 t 最小可检测效应（与稿件 Methods 同法）。"""
    ds = np.linspace(0.01, 12, 6000)
    df = 2 * n - 2
    tcrit = nct.ppf(1 - ALPHA / 2, df, 0)
    nc = ds * np.sqrt(n / 2.0)
    pw = nct.sf(tcrit, df, nc) + nct.cdf(-tcrit, df, nc)
    ok = np.where(pw >= POWER)[0]
    return float(ds[ok[0]]) if len(ok) else np.nan


def main():
    # ================= 0. 单细胞对象的规模与 QC =================
    tot = grep_md(os.path.join(HERE, "GSE267933_qc_summary.md"),
                  r"总细胞数:\s*(\d+)", int)
    add("Cells before QC", 20684, tot, tot == 20684, "GSE267933_qc_summary.md")

    qc = grep_md(os.path.join(SCO, "step5_summary.md"), r"QC 后:\s*(\d+)", int)
    add("Cells retained after QC", 18328, qc, qc == 18328, "sc_out/step5_summary.md")

    nt = grep_md(os.path.join(SCO, "step5_summary.md"), r"->\s*(\d+)\s*种细胞类型", int)
    add("Cell types annotated", 9, nt, nt == 9, "sc_out/step5_summary.md")

    ps = pd.read_csv(os.path.join(SCO, "step5d_per_sample_cells.csv"), index_col=0)
    lo, hi = int(ps.n_cells_qc.min()), int(ps.n_cells_qc.max())
    add("Per-sample QC cells range", "2,384-3,680", "%d-%d" % (lo, hi),
        (lo, hi) == (2384, 3680), "sc_out/step5d_per_sample_cells.csv")

    # ================= 1. 组成：动物级 vs 池化细胞级 =================
    d5 = pd.read_csv(os.path.join(SCO, "step5d_composition_animallevel.csv")).set_index("celltype")
    n_sig = int((d5.BH_Welch_animal < 0.05).sum())
    add("Animal-level composition: cell types with BH<0.05", "0 of 9", "%d of %d" % (n_sig, len(d5)),
        n_sig == 0, "sc_out/step5d_composition_animallevel.csv")

    min_bh = float(d5.BH_Welch_animal.min())
    add("Smallest animal-level BH p", "0.34", round(min_bh, 2), min_bh == min_bh and
        rounded_ok(0.34, min_bh, 2), "sc_out/step5d_composition_animallevel.csv")

    mg = d5.loc["Microglia"]
    add("Microglia, animal level (pct)", "39.5 -> 48.0",
        "%.1f -> %.1f" % (mg.pct_ctrl_animalmean, mg.pct_surg_animalmean),
        rounded_ok(39.5, mg.pct_ctrl_animalmean, 1) and rounded_ok(48.0, mg.pct_surg_animalmean, 1),
        "sc_out/step5d_composition_animallevel.csv")
    add("Microglia, animal level (p)", "0.14", round(mg.Welch_p_animal, 2),
        rounded_ok(0.14, mg.Welch_p_animal, 2), "sc_out/step5d_composition_animallevel.csv")
    add("Microglia, animal level (BH)", "0.47", round(mg.BH_Welch_animal, 2),
        rounded_ok(0.47, mg.BH_Welch_animal, 2), "sc_out/step5d_composition_animallevel.csv")
    add("Microglia, animal level (Cohen d)", "1.50", mg.Cohen_d_animal,
        rounded_ok(1.50, mg.Cohen_d_animal, 2), "sc_out/step5d_composition_animallevel.csv")
    add("Microglia, required n per group", 9, int(mg.req_n_per_group_exact),
        int(mg.req_n_per_group_exact) == 9, "sc_out/step5d_composition_animallevel.csv")

    ep = d5.loc["Ependymal"]
    add("Ependymal, animal level (pct)", "10.7 -> 3.6",
        "%.1f -> %.1f" % (ep.pct_ctrl_animalmean, ep.pct_surg_animalmean),
        rounded_ok(10.7, ep.pct_ctrl_animalmean, 1) and rounded_ok(3.6, ep.pct_surg_animalmean, 1),
        "sc_out/step5d_composition_animallevel.csv")
    add("Ependymal, animal level (p / BH)", "0.038 / 0.34",
        "%.3f / %.2f" % (ep.Welch_p_animal, ep.BH_Welch_animal),
        rounded_ok(0.038, ep.Welch_p_animal, 3) and rounded_ok(0.34, ep.BH_Welch_animal, 2),
        "sc_out/step5d_composition_animallevel.csv")

    c5 = pd.read_csv(os.path.join(SCO, "sc_celltype_composition.csv")).set_index("celltype")
    add("Pooled-cell Fisher BH, microglia", "2.7e-27", "%.1e" % c5.loc["Microglia", "BH_p"],
        num_ok(2.747e-27, c5.loc["Microglia", "BH_p"], rtol=1e-2),
        "sc_out/sc_celltype_composition.csv")
    add("Pooled-cell Fisher BH, ependymal", "1.0e-87", "%.1e" % c5.loc["Ependymal", "BH_p"],
        num_ok(1.021e-87, c5.loc["Ependymal", "BH_p"], rtol=1e-2),
        "sc_out/sc_celltype_composition.csv")
    ratio = float(np.log10(c5.loc["Microglia", "BH_p"]) - np.log10(mg.BH_Welch_animal))
    add("Shift from inferential unit alone (orders of magnitude)", ">25",
        "%.1f" % abs(ratio), abs(ratio) > 25, "both rows above")

    # ================= 2. 状态：两条 pseudobulk 路线 =================
    cal = pd.read_csv(os.path.join(SCO, "step5c_genomewide_calibration.csv"))
    n_tests = int(cal.n_genes.sum())
    n_hits = int(cal["padj_lt_0.05"].sum())
    add("Welch pseudobulk: tests reaching FDR<0.05", "0 of 173,421",
        "%d of %s" % (n_hits, format(n_tests, ",")),
        n_hits == 0 and n_tests == 173421, "sc_out/step5c_genomewide_calibration.csv")

    below = int((cal["p_lt_0.05"] < cal["exp_p_lt_0.05"]).sum())
    add("Cell types with fewer raw p<0.05 than null expectation", "7 of 9",
        "%d of %d" % (below, len(cal)), below == 7, "sc_out/step5c_genomewide_calibration.csv")

    vo = pd.read_csv(os.path.join(SCO, "step5e_limma_voom_stats.csv"))
    vn, vt = int(vo["n_adjp_lt_0.05"].sum()), int(vo["n_genes_tested"].sum())
    add("limma-voom: tests reaching adjP<0.05", "11 of 72,204",
        "%d of %s" % (vn, format(vt, ",")), vn == 11 and vt == 72204,
        "sc_out/step5e_limma_voom_stats.csv")
    add("limma-voom: share significant", "0.015%", "%.3f%%" % (100.0 * vn / vt),
        abs(100.0 * vn / vt - 0.015) < 0.002, "sc_out/step5e_limma_voom_stats.csv")

    sig = pd.read_csv(os.path.join(SCO, "step5e_voom_significant_genes.csv"))
    add("Ttr among the 11 voom hits (ambient RNA)", "5 of 11",
        "%d of %d" % (int((sig.gene == "Ttr").sum()), len(sig)),
        int((sig.gene == "Ttr").sum()) == 5, "sc_out/step5e_voom_significant_genes.csv")

    ve = pd.read_csv(os.path.join(SCO, "pbcounts", "voom_Ependymal.csv"), index_col=0)
    vn_ = pd.read_csv(os.path.join(SCO, "pbcounts", "voom_Neuron.csv"), index_col=0)
    add("Ttr logFC, ependymal source", -3.06, round(float(ve.loc["Ttr", "logFC"]), 2),
        rounded_ok(-3.06, ve.loc["Ttr", "logFC"], 2), "sc_out/pbcounts/voom_Ependymal.csv")
    add("Ttr logFC, neuron (barely expresses)", -5.70, round(float(vn_.loc["Ttr", "logFC"]), 2),
        rounded_ok(-5.70, vn_.loc["Ttr", "logFC"], 2), "sc_out/pbcounts/voom_Neuron.csv")
    add("Ttr expression, ependymal (log2CPM)", "16.3", round(float(ve.loc["Ttr", "AveExpr"]), 1),
        rounded_ok(16.3, ve.loc["Ttr", "AveExpr"], 1), "sc_out/pbcounts/voom_Ependymal.csv")

    # ================= 3. 功效边界 =================
    mde3 = mde_noncentral(3)
    add("MDE at n=3 per group (non-central t)", "3.07", round(mde3, 2),
        rounded_ok(3.07, mde3, 2), "computed here, same method as Methods")

    rq = pd.read_csv(os.path.join(SCO, "step5f_required_n.csv")).set_index("effect_size_d")
    for d, n in [(1.5, 9), (1.0, 17), (0.5, 64), (0.3, 176)]:
        add("Required n per group at d=%s" % d, n, int(rq.loc[d, "required_n_per_group"]),
            int(rq.loc[d, "required_n_per_group"]) == n, "sc_out/step5f_required_n.csv")

    # ================= 4. 同数据两条目 =================
    g5 = pd.read_csv(os.path.join(SCO, "step5g_same_data_two_entries.csv"))
    with gzip.open(os.path.join(HERE, "GSE289098_barcodes.tsv.gz"), "rt") as f:
        n98 = sum(1 for line in f if line.strip())
    add("GSE289098 barcodes", 20684, n98, n98 == 20684, "GSE289098_barcodes.tsv.gz")
    add("Per-sample barcode sets identical (6/6)", "6 of 6",
        "%d of %d" % (int(g5.barcode_set_identical.sum()), len(g5)),
        bool(g5.barcode_set_identical.all()) and len(g5) == 6,
        "sc_out/step5g_same_data_two_entries.csv")
    add("Our per-sample totals match the deposited tallies", "match",
        "match" if (g5.our_n_cells == g5.deposited_n_cells).all() else "mismatch",
        bool((g5.our_n_cells == g5.deposited_n_cells).all()),
        "sc_out/step5g_same_data_two_entries.csv")

    # ================= 5. GSE178995 不可用 =================
    txt = gzip.open(os.path.join(HERE, "GSE178995_series_matrix.txt.gz"),
                    "rt", errors="replace").read()
    n_sra = len(re.findall(r'Sample_type\t"[^"]*"', txt)) and txt.count('"SRA"')
    n_rows_zero = txt.count('!Sample_data_row_count\t"0"')
    add("GSE178995: samples typed as SRA with no data rows", "4 SRA / 0 rows",
        "%d SRA / %d zero-row" % (n_sra, n_rows_zero),
        n_sra == 4 and n_rows_zero >= 1, "GSE178995_series_matrix.txt.gz")

    # ================= 6. bulk 四机制线元分析 =================
    s3 = pd.read_csv(os.path.join(HERE, "step2c_geneset_stats_3ds.csv")).set_index("set")
    s4 = pd.read_csv(os.path.join(HERE, "step2c_geneset_stats_4ds.csv")).set_index("set")
    dz = (s3.mean_meta_Z - s4.mean_meta_Z).abs().max()
    add("Bulk meta-Z change, 3 -> 4 datasets", "<= 0.013", round(float(dz), 3),
        float(dz) <= 0.013 + 1e-9, "step2c_geneset_stats_3ds/4ds.csv")
    strong = ["B_DAM_microglia", "C_mitochondrial_OXPHOS", "D_myelin_oligodendrocyte"]
    add("Bulk lines B/C/D BH<1e-13 at four datasets", "3 of 3",
        "%d of %d" % (int((s4.loc[strong, "BH_p"] < 1e-13).sum()), len(strong)),
        bool((s4.loc[strong, "BH_p"] < 1e-13).all()), "step2c_geneset_stats_4ds.csv")
    add("Complement set A BH, 4 datasets", "1.9e-03",
        "%.2e" % s4.loc["A_complement_synapse_pruning", "BH_p"],
        abs(s4.loc["A_complement_synapse_pruning", "BH_p"] - 1.93e-3) < 1e-4,
        "step2c_geneset_stats_4ds.csv")
    add("Complement set A BH, 3 datasets", "9.5e-04",
        "%.2e" % s3.loc["A_complement_synapse_pruning", "BH_p"],
        abs(s3.loc["A_complement_synapse_pruning", "BH_p"] - 9.55e-4) < 1e-5,
        "step2c_geneset_stats_3ds.csv")
    add("Bulk four lines all remain BH<0.05 at four datasets", "4 of 4",
        "%d of %d" % (int((s4.BH_p < 0.05).sum()), len(s4)),
        bool((s4.BH_p < 0.05).all()), "step2c_geneset_stats_4ds.csv")

    pdg = pd.read_csv(os.path.join(HERE, "step2c_per_dataset_DEG.csv"))
    ast = pdg[pdg.dataset == "GSE199318_astro"].set_index("gene")
    add("Astrocytic C3 (GSE199318) log2FC / p / BH", "+2.45 / 0.25 / 0.66",
        "%+.2f / %.2f / %.2f" % (ast.loc["C3", "log2FC"], ast.loc["C3", "p"], ast.loc["C3", "padj"]),
        rounded_ok(2.45, ast.loc["C3", "log2FC"], 2) and rounded_ok(0.25, ast.loc["C3", "p"], 2)
        and rounded_ok(0.66, ast.loc["C3", "padj"], 2), "step2c_per_dataset_DEG.csv")
    add("Astrocytic C1ra (GSE199318) log2FC", "-4.78", round(float(ast.loc["C1ra", "log2FC"]), 2),
        rounded_ok(-4.78, ast.loc["C1ra", "log2FC"], 2), "step2c_per_dataset_DEG.csv")

    # ================= 汇总 =================
    w = max(len(r[0]) for r in RESULTS)
    print("=" * (w + 46))
    print("稿件数值审计 — %d 条" % len(RESULTS))
    print("=" * (w + 46))
    for claim, exp, act, ok, src in RESULTS:
        print("[%s] %-*s  稿件=%-16s 实测=%-16s  %s"
              % ("通过" if ok else "未通过", w, claim, exp, act, src))
    bad = [r for r in RESULTS if not r[3]]
    print("=" * (w + 46))
    print("通过 %d / %d" % (len(RESULTS) - len(bad), len(RESULTS)))
    if bad:
        print("\n未通过条目：")
        for claim, exp, act, _, src in bad:
            print("  - %s：稿件 %s vs 实测 %s（%s）" % (claim, exp, act, src))
    print("=" * (w + 46))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
