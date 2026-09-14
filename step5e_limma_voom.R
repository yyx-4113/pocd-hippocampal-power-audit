# Step 5e (b) — pseudobulk limma-voom 敏感性分析 (审稿意见 P2)
# =====================================================================
# 输入: sc_out/pbcounts/pbcounts_<celltype>.csv (gene x 6 samples, integer counts)
#       sc_out/pbcounts/sample_meta.csv
# 流程: filterByExpr 风格低表达过滤 -> (可选 TMM) -> voom -> lmFit(~group) -> eBayes
# 输出: sc_out/step5e_limma_voom_stats.csv
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- args[1]
files <- list.files(OUTDIR, pattern = "^pbcounts_.*\\.csv$", full.names = TRUE)
meta <- read.csv(file.path(OUTDIR, "sample_meta.csv"), stringsAsFactors = FALSE)
grp <- factor(meta$group, levels = c("Control", "Surgery"))
design <- model.matrix(~grp)

has_edgeR <- requireNamespace("edgeR", quietly = TRUE)
cat("edgeR available:", has_edgeR, "\n")

res <- list()
for (f in files) {
  ct <- sub("^pbcounts_(.*)\\.csv$", "\\1", basename(f))
  m <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  storage.mode(m) <- "numeric"
  # 低表达过滤 (edgeR::filterByExpr 的稳健近似: 至少 3 个样本 counts>=10)
  keep <- rowSums(m >= 10) >= 3
  m2 <- m[keep, , drop = FALSE]
  if (nrow(m2) < 500) { cat("skip", ct, "too few genes\n"); next }

  libsize <- colSums(m2)
  d <- list(counts = m2, lib.size = libsize, samples = colnames(m2))
  if (has_edgeR) {
    d <- edgeR::calcNormFactors(edgeR::DGEList(counts = m2))
    v <- voom(d, design)
  } else {
    v <- voom(m2, design, lib.size = libsize)
  }
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

  res[[ct]] <- data.frame(
    celltype = ct,
    n_genes_tested = nrow(tt),
    n_genes_input = nrow(m),
    n_p_lt_0.05 = sum(tt$P.Value < 0.05),
    exp_p_lt_0.05 = round(0.05 * nrow(tt), 1),
    n_adjp_lt_0.05 = sum(tt$adj.P.Val < 0.05),
    min_p = min(tt$P.Value),
    min_adjp = min(tt$adj.P.Val),
    top_gene = rownames(tt)[1],
    top_logFC = round(tt$logFC[1], 3),
    stringsAsFactors = FALSE
  )
  cat(sprintf("%-16s tested=%6d  p<0.05=%4d (exp %.1f)  adjp<0.05=%3d  min_adjp=%.3g  top=%s(%+.2f)\n",
              ct, nrow(tt), sum(tt$P.Value < 0.05), 0.05 * nrow(tt),
              sum(tt$adj.P.Val < 0.05), min(tt$adj.P.Val), rownames(tt)[1], tt$logFC[1]))
  write.csv(tt, file.path(OUTDIR, paste0("voom_", ct, ".csv")))
}

out <- do.call(rbind, res)
write.csv(out, file.path(dirname(OUTDIR), "step5e_limma_voom_stats.csv"), row.names = FALSE)
cat("\n=== TOTAL ===\n")
cat("total genes tested:", sum(out$n_genes_tested),
    " total adjp<0.05:", sum(out$n_adjp_lt_0.05), "\n")
cat("saved step5e_limma_voom_stats.csv\n")
