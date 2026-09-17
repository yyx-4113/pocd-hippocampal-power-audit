# pocd-hippocampal-power-audit — reproducibility bundle

Version-controlled companion repository for the manuscript:

> Yang Y. *Underpowered at every level: an animal-level power audit of 3-vs-3
> hippocampal single-cell RNA-seq in postoperative neurocognitive dysfunction.*
> Anesthesia & Analgesia (submitted).

The manuscript asks a single question about a single public dataset: what can a
3-vs-3 hippocampal scRNA-seq design in a POCD model actually support? This
bundle contains the complete pipeline behind the answer: the GEO inventory,
every statistical recomputation, the figure sources, the limma-voom
sensitivity analysis, and the power calibration. It is released under the MIT
licence. A persistent DOI (Zenodo) will be minted on acceptance.

## The claim in one paragraph

Cell-type composition was originally compared with a pooled-cell Fisher test.
That test treats each of ~18,328 cells as an independent observation, but the
experimental unit is the animal (n = 3 per group). Re-testing the same animals
at animal level removes the significance entirely: 0 of 9 cell types reach
BH p<0.05 (largest effects: ependymal 10.7% to 3.6%, BH p=0.34; microglia
39.5% to 48.0%, BH p=0.47), while the pooled-cell test returns BH
p=2.7e-27 and 1.0e-87 for those same two comparisons. Cell state is limited in
the same way: 0 of 173,421 gene-level Welch tests reach FDR<0.05, and 11 of
72,204 limma-voom tests do, five of them one transcript (Ttr) that shifts
further in cell types that barely express it than in its ependymal source, the
signature of ambient RNA. At n=3/group the minimum detectable effect is
d=3.07; detecting d=1.5 needs n=9 per group and d=0.5 needs n=64.

## What is in here

```
step1b_gse267933_scloader.py        GSE267933 RAW tar -> AnnData (fetch + 10x assembly)
step5_sc_analysis.py                QC -> lognorm -> HVG -> PCA -> Harmony -> UMAP -> Leiden
                                    -> marker-score annotation -> per-celltype pseudobulk
step5d_composition_animallevel.py   per-animal composition, Welch/MWU, BH, Cohen's d, required n
step5b_sc_pseudobulk_geneset.py     four mechanism gene sets (complement / DAM / OXPHOS / myelin)
step5b2_sensitivity.py              gene-set results under two pseudobulk constructions
step5c_power_and_calibration.py     genome-wide BH calibration + MDE by sample size
step5e_export_pseudobulk_counts.py  integer pseudobulk counts per cell type (for limma-voom)
step5e_limma_voom.R                 limma-voom + edgeR TMM sensitivity (R 4.4.3)
step5e_report_numbers.py            Welch-vs-voom comparison table + Ttr ambient-RNA panel
step5d_dm_sensitivity.py            Dirichlet-multinomial overdispersion sensitivity for the composition null (m5)
step5e_rank_correlation.py          limma-voom vs Welch gene-ranking correlation (Spearman rho / Pearson r)
step5h_effect_cis.py                Welch / Cohen's d 95% CIs + Welch-Satterthwaite worst-case MDE band
step5f_figs_revised.py              Fig 1 / Fig 2 + required-n table
fig3_same_data_two_entries.py       Fig 3 (GSE267933 vs GSE289098, barcode correspondence)
fig4_astrocyte_C3.py                Fig 4 (astrocytic complement, p / BH per gene)
figS1_method_sensitivity.py         Fig S1 (method sensitivity + Ttr fingerprint)
verify_manuscript_numbers.py        read-only audit: every number quoted in the manuscript,
                                    checked against its source file (exit 1 on any mismatch)

step1_gse276942_timeseries.py       bulk temporal trajectory (GSE276942, 8 time points)
step2_bulk_integrate_pathway.py     per-dataset DEG + Stouffer meta-analysis + pathway sets
step2b_geneset_stats.py             gene-set level statistics (meta_Z, Wilcoxon, Stouffer)
step2c_extend_bulk.py               four-dataset meta-analysis (adds GSE178995-eligible sets)
step3_dual_ml.py                    dual-model feature selection + cross-dataset validation

download_raw.py, download_matrix.py, list_ftp.py, inspect_matrix.py, dl_parse_matrix.py
                                    GEO retrieval and format inspection (HTTPS route)

sc_out/                             every numerical output quoted in the manuscript
sc_out/pbcounts/                    pseudobulk counts, sample metadata, limma-voom fits
figures/, fig_*.png                 regenerated figures (the five manuscript figures are
                                    fig_ApC_*.png)
matrix_tmp/, GSE*_*.xml, *.gz, GSE174412_RAW/
                                    cached GEO records and matrices used to build the tables
```

Numbers cited in the manuscript map to files as follows.

| Manuscript figure / table | Source file |
|---|---|
| Fig. 1, left (animal-level composition) | `sc_out/step5d_composition_animallevel.csv` |
| Fig. 1, right (raw p<0.05 vs null) | `sc_out/step5c_genomewide_calibration.csv` |
| Fig. 2 (MDE curve, required n) | `sc_out/step5f_required_n.csv`, `sc_out/step5c_mde.csv` |
| Fig. 3 (two deposits, one experiment) | `sc_out/step5g_same_data_two_entries.csv`, `GSE289098_barcodes.tsv.gz` vs `GSE267933_RAW/` |
| Per-animal cells after QC (2,384-3,680) | `sc_out/step5d_per_sample_cells.csv` |
| Fig. 4 (astrocytic complement) | `sc_out/pbcounts/voom_Astrocyte.csv` |
| Fig. S1 (Welch vs voom, Ttr) | `sc_out/step5e_method_comparison.csv`, `sc_out/step5e_limma_voom_stats.csv` |
| Fig. S1 (voom vs Welch ranking correlation) | `sc_out/step5e_rank_correlation.csv` |
| Composition null under overdispersion (DM sensitivity) | `sc_out/step5d_dm_sensitivity.csv` |
| Gene-set statistics (4 mechanism lines) | `sc_out/step5b_geneset_stats.csv`, `sc_out/step5b2_sensitivity.csv` |
| Bulk DEG meta-analysis | `step2c_meta_DEG.csv`, `step2c_per_dataset_DEG.csv` |
| Feature-gene panel | `step3_feature_genes.csv`, `step3_roc_summary.csv` |

## Data provenance

All data are public GEO deposits; no new animals were used. Retrieval date for
every cached record in this bundle: **14 September 2026**, over HTTPS from
`https://ftp.ncbi.nlm.nih.gov/geo/series/`.

| Accession | Role | Samples |
|---|---|---|
| GSE267933 | hippocampal scRNA-seq, 18-month male C57BL/6 (3 control / 3 surgery) | 6 |
| GSE289098 | second deposit of the same experiment (see note below) | 6 |
| GSE276942 | hippocampal bulk RNA-seq time course (control, 24 h, 3 d, 7 d, 1-4 months) | 17 |
| GSE215410 | hippocampal bulk RNA-seq, POCD 2-vs-2 | 4 |
| GSE174412 | hippocampal bulk RNA-seq, PND 3-vs-3 | 6 |
| GSE178995 | annotated but **no processed matrix** deposited (raw reads only) | - |
| GSE199318 | candidate validation dataset (6 samples) | 6 |

Two large objects are **not** archived in the repository because GEO supplies
them and they are regenerable:

* `GSE267933_RAW.tar` (134 MB) - `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267933/suppl/GSE267933_RAW.tar`
* `GSE267933.h5ad` (549 MB) and `sc_out/GSE267933_processed.h5ad` (715 MB) - rebuilt by `step1b_gse267933_scloader.py` then `step5_sc_analysis.py`

## Reproduce

Commands are run from the repository root, in this order. Steps 0-1 need only
CPU and a few GB of RAM; step 2 (Harmony + Leiden on 18,328 cells) takes a few
minutes.

```bash
# 0. fetch the raw 10x matrices (134 MB)
python download_raw.py

# 1. assemble the AnnData object
python step1b_gse267933_scloader.py

# 2. QC, integration, annotation, per-celltype pseudobulk
python step5_sc_analysis.py

# 3. the manuscript's central table: composition at the level of the animal
python step5d_composition_animallevel.py

# 3b. Dirichlet-multinomial overdispersion sensitivity for the composition null
python step5d_dm_sensitivity.py

# 4. gene-set statistics and construction-method sensitivity
python step5b_sc_pseudobulk_geneset.py
python step5b2_sensitivity.py

# 5. genome-wide calibration and the detectable-effect boundary
python step5c_power_and_calibration.py

# 6. limma-voom sensitivity (R 4.4.3)
python step5e_export_pseudobulk_counts.py
Rscript step5e_limma_voom.R sc_out/pbcounts
python step5e_report_numbers.py

# 6b. limma-voom vs Welch gene-ranking correlation (P2-A)
python step5e_rank_correlation.py

# 7. figures
python step5f_figs_revised.py
python fig3_same_data_two_entries.py
python fig4_astrocyte_C3.py
python figS1_method_sensitivity.py

# 8. self-check: every manuscript number against its source file
python verify_manuscript_numbers.py

# 9. Welch / Cohen's d 95% CIs and the Welch-worst-case MDE band (Fig. 2)
python step5h_effect_cis.py
```

`step5d_composition_animallevel.py` regenerates the archived
`sc_out/step5d_composition_animallevel.csv` **bit-for-bit** (verified: all
columns identical to 1e-15 relative, including the superseded
`req_n_per_group` normal approximation that was kept for traceability). The
canonical effect-size-to-required-n column is `req_n_per_group_exact`.

The bulk branch (Steps 1-3, `step1_*`, `step2*`, `step3_*`) is independent of
the single-cell branch and only feeds the manuscript's background and the
feature-gene panel; run it in file-name order.

## Notes and caveats

* **`req_n_per_group` is superseded.** It is the early normal approximation
  `2(z(1-a/2)+z(1-b))^2 / d^2`. Use `req_n_per_group_exact` (non-central t,
  df = 2n-2), which is what the manuscript reports.
* **Cohen's d in `step5d_...csv`** uses the equal-n average-variance form
  `(mean_surg - mean_ctrl) / sqrt((var_surg + var_ctrl)/2)` computed on
  unrounded per-animal percentages, so that the archived values reproduce
  exactly. Group variances use ddof = 1.
* **GSE289098 is not an independent replication.** It contains the same 20,684
  cells as GSE267933, and the correspondence is stronger than a cell count: for
  all six samples the 10x 16-mer barcode **set** is identical between the two
  deposits (3,171-4,083 distinct barcodes each, matched one-to-one), and the
  depositor's own protocol labels them in the same order. `fig3_same_data_two_entries.py`
  derives the correspondence from the deposited barcode files and writes it to
  `sc_out/step5g_same_data_two_entries.csv`, so it can be re-checked directly.
* **Ambient RNA.** The 11 limma-voom hits are dominated by `Ttr`, which shifts
  more in cell types that barely express it than in ependymal cells, its true
  source. This is reported as an ambient-RNA signature, not as biology, and no
  decontamination (e.g. SoupX/CellBender) was applied; that is stated as a
  limitation in the manuscript.
* **GSE178995 cannot be used** for expression-level meta-analysis: the series
  annotates four samples as `Sample_type = SRA` and deposits no processed
  matrix.

## Software

| Component | Version |
|---|---|
| Python | 3.13.14 |
| scanpy / anndata | 1.12.4 / 0.13.3 |
| numpy / pandas / scipy | 2.4.6 / 3.0.5 / 1.18.1 |
| statsmodels / scikit-learn | 0.15.0 / 1.9.0 |
| harmonypy | 0.2.0 |
| matplotlib / seaborn | 3.11.1 / 0.13.2 |
| R / limma / edgeR | 4.4.3 / 3.62.2 / 4.4.2 |

Harmony was called through `harmonypy.run_harmony` directly rather than via the
scanpy wrapper, because the wrapper is incompatible with harmonypy 0.2.0
(torch-based) and silently falls back to uncorrected PCA.

## Licence and identifiers

MIT (see `LICENSE`). Full analysis-pipeline citation metadata is in
`CITATION.cff`. A Zenodo DOI will be added here on acceptance; until then,
please cite the repository URL and the commit hash.

## Author verification

See `author_verification_statement.md` for what was checked, by what method,
and on what date; a Word copy for reviewers is included alongside it as
`author_verification_statement.docx`. Section 11 of that statement records
compliance with the Editorial Office's technical check on manuscript
AA-D-26-01935 (all-English submission files, continuous page numbers in the
manuscript, and no author contribution statement in the manuscript file).
