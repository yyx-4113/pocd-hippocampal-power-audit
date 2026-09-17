# Author verification statement

**Manuscript:** "Underpowered at the animal level for composition and cell state: a power audit of 3-vs-3 hippocampal single-cell RNA-seq in postoperative neurocognitive dysfunction" (Original Laboratory Research Report, *Anesthesia & Analgesia*, submitted).

**Prepared by:** Yongxin Yang (author)
**Method of verification:** every externally checkable claim in the manuscript was re-derived from the primary source with the scripts shipped in this repository, or re-checked against PubMed / the publisher record by identifier. No number was carried over from an earlier draft without re-derivation.
**Date of verification:** 14 September 2026; re-verified 15 September 2026 against the submitted revision, and re-checked 17 September 2026 against the editorial-check revision.

## What was checked

### 1. Cell-type composition at the level of the animal (the manuscript's central analysis)
`step5d_composition_animallevel.py` was written as a clean reimplementation and its
output compared against the table the manuscript was drafted from. Per-animal
proportions, Welch p, Mann-Whitney U p, both BH columns, Cohen's d and the exact
required-n column are **identical to 1e-15 relative** across all nine cell types.
Result: **0 of 9 cell types reach BH p<0.05 at animal level** (smallest BH = 0.34,
ependymal). The largest effect is ependymal (10.7% to 3.6%, Welch p = 0.038) and
microglia (39.5% to 48.0%, p = 0.14; per-animal values 42.8/33.2/42.4% control and
41.4/49.7/53.1% surgery; d = 1.50, requiring n = 9 per group at 80% power).

A reviewer can re-run this file and diff it against `sc_out/step5d_composition_animallevel.csv`.
The column `req_n_per_group` is the superseded normal approximation retained for
traceability; `req_n_per_group_exact` is the value the manuscript reports.

### 2. Composition at the pooled-cell level (for comparison only)
`sc_out/sc_celltype_composition.csv` holds the pooled-cell Fisher test that the
literature uses. For the same six animals it returns BH p = 2.7e-27 (microglia) and
1.0e-87 (ependymal) against BH p = 0.47 and 0.34 at animal level. Both numbers were
recomputed from `sc_out/GSE267933_processed.h5ad` and are unchanged.

### 3. Cell state by two independent pseudobulk routes
* Welch on log1p pseudobulk means: **0 of 173,421** tests reach FDR<0.05
  (`sc_out/step5c_genomewide_calibration.csv`, 19,269 genes x 9 cell types; the
  smallest raw p is 1e-5).
* limma-voom with edgeR TMM on integer pseudobulk counts (R 4.4.3, limma 3.62.2,
  edgeR 4.4.2): **11 of 72,204** tests reach adjP<0.05
  (`sc_out/step5e_limma_voom_stats.csv`). The same 11 genes are obtained with and
  without TMM, so the result does not depend on the normalization choice.

### 4. Ambient-RNA attribution
Five of the eleven limma-voom hits are the single transcript `Ttr`. Its fold change
is larger in cell types that barely express it (neuron, -5.70) than in its true
source, ependymal cells (+16.3 log2CPM; -3.06). This pattern, not directed biology,
is what the manuscript reports; no decontamination was applied and that is stated
as a limitation.

### 5. Genome-wide false-positive calibration
For each of the nine cell types the number of genes with raw p<0.05 was enumerated
and compared with the null expectation (963.5 of 19,269). Seven of nine cell types
produce fewer raw hits than expected under the null, i.e. the tests are conservative
rather than anti-conservative. This is why the manuscript does not describe the
single-cell negatives as biological negatives.

### 6. GSE267933 and GSE289098 are the same experiment
Both series report 20,684 cells, and per-sample barcode sets correspond 1:1
(the deposited barcode lists are shipped in this repository,
`sc_out/step5g_same_data_two_entries.csv`, so the correspondence can be checked
directly). The manuscript therefore treats the field as having one such experiment,
not two.

### 7. References
All 19 entries were verified against PubMed or the publisher DOI (see the manuscript's
DOI verification note). Corrections recorded there: [1] is *Anesth Analg*
2018;127(2):496-505, [2] is *Glia* 2013;61(1):71-90, [3] is *Ann Neurol* 2011;70(6):986-995,
[4] author list and DOI corrected, and the foundational C3 paper is Lian H *et al.*,
*Neuron* 2015;85(1):101-115 ([5]). Entries [17] to [19] were added during revision and
each was resolved through Crossref before insertion ([18] was additionally corrected to
its full citation, *J Neuroinflammation* 2024;21(1):3216). Result: **19/19 resolve to a
real record that supports the sentence citing it**. The author's companion commentary on
sleep-deprivation pseudoreplication is disclosed in the manuscript's reference note as
related work under review, not as a primary numbered source.

### 8. Data availability
GSE178995 is annotated as four `Sample_type = SRA` samples with no processed matrix
deposited; it therefore cannot enter expression-level meta-analysis and is excluded
explicitly. All other datasets are public GEO deposits retrieved over HTTPS.

### 9. Four a-priori mechanism gene sets
The four gene sets (A complement/synaptic pruning; B broader glial/microglial
activation including DAM; C mitochondrial OXPHOS; D myelin/oligodendrocyte) are
listed with their curated members in `sc_out/mechanism_genesets.csv` (A:173, B:119,
C:149, D:76). The exact analyzable members (curated members present in the integrated
bulk matrices) are provided in `sc_out/step2_pathway_sets.json` (A:163, B:118, C:191,
D:75). The bulk Stouffer meta-analysis uses these matrix-present members.

### 10. Self-reported counts
Text length (3,949 words including the figure list, 3,611 without it; both within the
A&A 4,000-word limit), structured abstract (389 words), figure/table count (4 main
figures + 1 appendix table = 5 main items, ≤6; 3 online supplementary figures S1-S3 and
Supplementary Table S1) and reference count (19) were recomputed from the manuscript
source and the count scripts, not carried over.

### 11. Editorial check compliance (AA-D-26-01935)
The Editorial Office returned the submission for four changes, all of which have been
made. The corresponding author's name is held in Editorial Manager in the English form
Yongxin Yang (given name Yongxin, family name Yang). Every submitted file is in English:
a Chinese revision log that had been rendered into the cover letter was removed, and the
filename of every uploaded item is English. The manuscript file carries continuous page
numbers in the footer. The author contribution statement was removed from the manuscript
file at the Editorial Office's request; the author is the sole contributor, who conceived
the study, performed all reanalyses, wrote the manuscript and approved the final version.
Each of these is now gated by `verify_submission_docx.py` (no non-English characters, a
PAGE field present in the footer, and no author contribution statement in the manuscript).

## Repository and persistent identifier
The analysis pipeline, all intermediate tables, the cached GEO records and matrices,
the limma-voom fits, the four a-priori mechanism gene sets
(`sc_out/mechanism_genesets.csv`, `sc_out/step2_pathway_sets.json`) and the
figure-generation scripts are assembled as a version-controlled bundle under the MIT
licence and are public at https://github.com/yyx-4113/pocd-hippocampal-power-audit
(release **v1.2.0**). Large re-downloadable objects (`GSE267933_RAW.tar`, `.h5ad`
objects) are excluded and regenerated by the scripts; the retrieval route is documented
in `README.md`. A persistent Zenodo DOI will be added to this file and to the
manuscript's data availability statement on acceptance.

**Signed:** Yongxin Yang
**Date:** 14 September 2026 (re-verified 15 September 2026; editorial-check revision 17 September 2026)
