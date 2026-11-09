# Cluster calls

Hand-reviewed decisions for each Harmony cluster, one table per subset and
resolution: `cluster_calls_<subset>_res<res>_reviewed.tsv`. They are kept
here, under version control, because they are human decisions and not
pipeline output. Every edit shows up as a diff.

For res 0.6 the first seven columns are the first draft (2026-11-09). The
review added `celltype_label`, `reviewer_call`, `review` (agree / verify /
disagree) and `reviewer_note`. The draft columns are never edited; the
reviewer columns win where the two disagree. `check_result` records checks
run afterwards on existing v3 outputs (scDblFinder scores, nCount, dataset
enrichment over baseline). The reasoning and order of the cleanup are in the
Claude Doc "PDAC atlas: cluster-call cleanup plan (B, stromal, T/NK, res 0.6)".

Evidence comes from `rebuild_v3/step4_<subset>_res<res>_decision.tsv`,
`..._top20.tsv` and `step4b_<subset>_res<res>_gate.tsv`. The gate verdict is
evidence, not the decision. It passes pericytes/vSMC and Schwann cells as
fibroblast (HPCA `Tissue_stem_cells` maps to fibroblast), and it marks NK
clusters for reassignment because only `Tcell` counts as native in TNK.

The planned step 05 reads these tables. It must stop if a table's clusters
or `n_cells` do not match the object, because that means 03 was rerun and
the cluster numbers have changed.

| Column | Values |
|---|---|
| `cluster` | cluster id at that resolution |
| `n_cells` | cells in the cluster when the call was made (consistency check) |
| `call` | `keep` or `drop` |
| `lineage` | what the cluster is. Kept: `fibroblast`, `mural`, `schwann`, `Tcell`, `NK`, `Bcell`, `plasma`. Dropped rows record what was found instead (`Tcell`, `endothelial`, `RBC_lowquality`, `pDC`, ...) |
| `qc_flags` | `;`-separated, empty if none: `ambient_acinar`, `ambient_myeloid`, `ambient_fibroblast`, `ambient_epithelial`, `ambient_RBC`, `low_complexity` (median nFeature < 1000), `single_dataset` (> 80% one dataset), `stress_heatshock`, `metallothionein`, `batch_markers` (most top markers STRUCTURAL) |
| `confidence` | `clear`, or `judgement` (a reasonable reviewer could decide otherwise) |
| `evidence` | markers and the gate/panel numbers behind the call |

Policy behind the first draft (res 0.6, 2026-11-09):

- Mural and Schwann cells stay in the stromal atlas as their own lineages.
  Only `fibroblast` is used for CAF programme scoring.
- RBC and low-quality clusters are dropped.

The review changed two of the draft's policies:

- Misrouted cells are routed to the right subset, not dropped. The draft
  dropped them.
- Ambient-RNA and doublet problems are fixed per cell (DecontX, scDblFinder
  scores), not by keeping or dropping the whole cluster. Low-quality clusters
  with no cell identity are dropped.
