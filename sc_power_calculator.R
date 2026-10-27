#!/usr/bin/env Rscript
# =============================================================================
# sc_power_calculator.R
# Reusable sample-size calculator for ANIMAL-LEVEL inference in single-cell
# RNA-seq (and similar pseudoreplicated -omics designs).
#
# Reproduces the exact power model used in the manuscript
# (two-sample t-test, two-sided, alpha = 0.05, power = 0.80, df = 2n - 2,
# non-central t), and adds the variance-decomposition layer:
#
#   per-animal (pseudobulk) SD   sigma_a = sqrt(tau2 + sigma2 / k)
#   animal-level Cohen's d        d_a     = Delta / sigma_a
#   required n per group          = smallest n with power(d_a, n) >= target
#
# where tau2   = between-animal variance (biological, immune to cell depth)
#       sigma2 = within-animal / cell-level noise variance
#       k      = cells profiled per animal
#       Delta  = true treatment shift in animal means (biological effect)
#
# KEY MESSAGE the tool encodes:
#   Increasing k shrinks sigma2/k but can NEVER shrink tau2. Therefore cell
#   depth has a ceiling (set by tau2) on how much it helps animal-level power;
#   beyond that ceiling only more ANIMALS (n) help.
#
# No external packages required (base R only). Shiny wrapper optional.
# =============================================================================

power_at <- function(n, d, alpha = 0.05) {
  # Non-central t power for two-sample t-test, df = 2n - 2.
  n <- as.numeric(n); d <- as.numeric(d)
  df  <- 2 * n - 2
  tcrit <- qt(1 - alpha / 2, df)
  nc   <- d * sqrt(n / 2)
  1 - pt(tcrit, df, nc) + pt(-tcrit, df, nc)
}

required_n <- function(d, power = 0.80, alpha = 0.05, n_max = 2000) {
  # Smallest n per group achieving >= power.
  d <- as.numeric(d)
  for (n in 2:n_max) {
    if (power_at(n, d, alpha) >= power) return(n)
  }
  NA
}

# ---- variance decomposition -> animal-level effect size -> required n ----
sigma_a_of <- function(tau2, sigma2, k) sqrt(tau2 + sigma2 / k)
d_animal   <- function(tau2, sigma2, k, Delta) Delta / sigma_a_of(tau2, sigma2, k)

required_n_from_design <- function(tau2, sigma2, k, Delta, power = 0.80, alpha = 0.05) {
  da <- d_animal(tau2, sigma2, k, Delta)
  list(tau2 = tau2, sigma2 = sigma2, k = k, Delta = Delta,
       sigma_a = sigma_a_of(tau2, sigma2, k), d_animal = da,
       required_n_per_group = required_n(da, power, alpha),
       total_n = 2 * required_n(da, power, alpha))
}

# ============================ self-validation =============================
# Must reproduce the manuscript required-n table: d = 1.5->9, 1.0->17, 0.5->64, 0.3->176
cat("== cross-check vs Python power_at (manuscript table) ==\n")
for (d in c(1.5, 1.0, 0.5, 0.3)) {
  cat(sprintf("  d_a=%4.1f -> n=%d (expect %s)\n",
              d, required_n(d),
              ifelse(d==1.5,"9",ifelse(d==1.0,"17",ifelse(d==0.5,"64","176")))))
}

cat("\n== example: design-driven required n ==\n")
# Same biology (Delta=1.0); cell noise dominates (tau2=0.25, sigma2=4).
# Watch required n DROP as k grows, then plateau (tau2 ceiling).
for (k in c(10, 50, 200, 1000, 10000)) {
  r <- required_n_from_design(tau2 = 0.25, sigma2 = 4.0, k = k, Delta = 1.0)
  cat(sprintf("  k=%6d -> sigma_a=%.3f, d_a=%.3f, n/group=%d\n",
              k, r$sigma_a, r$d_animal, r$required_n_per_group))
}

if (Sys.getenv("SHINY_RUN") == "1" && requireNamespace("shiny", quietly = TRUE)) {
  library(shiny)
  ui <- fluidPage(
    titlePanel("scRNA-seq animal-level sample-size calculator"),
    sidebarLayout(
      sidebarPanel(
        numericInput("tau2", "Between-animal variance (tau2)", 1, min = 0.001),
        numericInput("sigma2", "Within-animal / cell noise (sigma2)", 1, min = 0.001),
        numericInput("k", "Cells per animal (k)", 200, min = 1),
        numericInput("delta", "True treatment shift (Delta)", 1),
        numericInput("power", "Target power", 0.80, min = 0.5, max = 0.99),
        numericInput("alpha", "Alpha", 0.05, min = 0.001, max = 0.2)
      ),
      mainPanel(verbatimTextOutput("out"))
    )
  )
  server <- function(input, output) {
    output$out <- renderPrint({
      r <- required_n_from_design(input$tau2, input$sigma2, input$k,
                                 input$delta, input$power, input$alpha)
      cat("Per-animal SD (sigma_a)      :", round(r$sigma_a, 4), "\n")
      cat("Animal-level Cohen's d (d_a) :", round(r$d_animal, 4), "\n")
      cat("Required n per group         :", r$required_n_per_group, "\n")
      cat("Total n (both groups)        :", r$total_n, "\n")
    })
  }
  shinyApp(ui, server)
}
