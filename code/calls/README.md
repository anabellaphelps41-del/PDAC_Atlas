# Cluster calls

Two rounds so far:

- `cluster_calls_{TNK,stromal,Bcell}_res0.6_reviewed.tsv`: v3 clusters, reviewed;
  applied by `rebuild_02a_routing.R` to build v4.
- `cluster_calls_{TNK,stromal,Bcell,myeloid,endothelial}_v4run1_res0.6.tsv`:
  draft calls on the first v4 clustering (`rebuild_v4/run1_snapshot/`). They
  were never reviewed row by row, and are kept for the record. The columns are
  the reviewed columns plus `v3_origin` and `followup`.
- `cluster_calls_{TNK,stromal,Bcell,myeloid,endothelial}_v4_res0.6.tsv`:
  calls on v4 run 2 (wider stress mask), reviewed and applied by step 05 on
  2026-11-10. The columns are:
  - `call` (keep / drop / route) and `route_to`;
  - `confidence`: `decided` for rows that carry a user decision from run 1;
  - `run1_origin`: which run-1 clusters the cells came from;
  - `samples`: samples with at least 5% of the cluster, and the top three;
  - `cell_flags`: the share of cells carrying each run-1 flag;
  - the reviewer columns;
  - `pool_label` (last column, added 2026-11-10): the shared label for
    clusters that the review pools for composition, because they are the same
    cell type split by sample, depth or ambient RNA. Empty means the
    cluster's own `celltype_label`. Step 05 writes it to `celltype_pool_v4`,
    next to the fine `celltype_v4`. Pools: TNK c8+c15; Bcell c0+c4 and
    c1+c9+c11; myeloid c0+c4 and c6+c7; endothelial c3+c4+c5 and c1+c6;
    stromal c6+c8.

  Myeloid c14 is split in `subcluster_res0.6` (14a cDC1 / 14b mregDC) by
  `rebuild_04c_subcluster.R`.

  Reviewer columns are filled for the rows decided on 2026-11-09
  (myeloid c16, TNK c5, stromal c7, endothelial c11/c12) and, on
  2026-11-10, for the 27 `judgement` rows and the 8 `decided` rows. The
  2026-11-10 review was delegated to Claude; its notes start with
  "2026-11-10 (Claude review, delegated)". Where it relabelled a cluster,
  `celltype_label` holds the new label and the note ends with the old one.
  The user confirmed those 35 rows on 2026-11-10 (the notes end with
  "Confirmed by the user 2026-11-10."). The 49 `clear` rows were accepted as
  drafted on the same day without a row-by-row review: `reviewer_call` copies
  `call`, and the note says so.
- `v4_run2_decisions.tsv`: the machine-readable form of those run-2
  decisions, in the run-1 format plus `n_cells`:
  - drops: myeloid c16, endothelial c11 and c12. DecontX is not rerun for
    the endothelial ones;
  - a whole-cluster `stress_overlay` flag with `exclude_from_de_lr` on
    TNK c5. Cells that already carry a flag keep it.

  The 2026-11-10 review added 11 rows:
  - drops: TNK c17 (T/B doublets), stromal c11 (CNV+ osteogenic cells of
    P19 and tumour-fibroblast doublets of GSE154778:P01), stromal c12 and
    myeloid c18 (low-quality droplets);
  - flags with `exclude_from_de_lr`: TNK c15 (`ambient_RBC`), Bcell c8
    (`ambient_mixed`), myeloid c12 (`few_donors`, ADJ2/P06), endothelial c4
    (`mt_high_sample`), and extensions of the run-1 flags to the whole of
    TNK c12 and Bcell c5 (the TNK c5 rule: unflagged cells as stressed as
    flagged ones);
  - flag without exclusion: endothelial c5 (`patient_specific`).

  Stromal c7 gets no cluster flag: its cells keep their run-1 flags only.
  These apply to the run-2 clusters as they are, without re-clustering.
  `rebuild_05_apply_calls.R` reads the file. It stops if `n_cells` or the
  action disagrees with the call table, because a rerun of 03 renumbers the
  clusters.
- `subcluster_labels_v4_res0.6.tsv`: names for the `subcluster_res0.6`
  groups made by `rebuild_04c_subcluster.R` (myeloid 14a cDC1, 14b mregDC,
  14u unassigned), with their `n_cells`. Step 05 uses them in place of the
  cluster's own label.

Step 05 takes `reviewer_call` where it is filled, and the draft `call`
otherwise. It only runs once every row is reviewed, unless
`PDAC_DRAFT_CALLS=1` is set. The user decided on 2026-11-10 to apply the
calls to the run-2 clusters without re-clustering. A re-clustering test
after the drops gave kept-cell ARI 0.84-0.93, with boundary shifts between
neighbouring states, so proportions between such states are approximate.
- `v4_run1_decisions.tsv`: the user's decisions on specific v4 run-1 clusters
  (2026-11-09). `rebuild_02a_routing.R` turns them into per-cell actions
  (`cell_calls_v4run1.tsv`):
  - one drop: myeloid c0;
  - `cell_flag` / `exclude_from_de_lr` on the stress-overlay, patient-specific,
    few-donor and high-contamination clusters.

  The flags follow each cell into whatever cluster it joins after the
  re-clustering with the wider stress mask (v4 run 2). Cells with
  `exclude_from_de_lr` stay out of DE and ligand-receptor work.

Hand-reviewed decisions for each Harmony cluster, one table per subset and
resolution: `cluster_calls_<subset>_res<res>_reviewed.tsv`. They are kept
here, under version control, because they are human decisions and not
pipeline output. Every edit shows up as a diff.

For res 0.6 the first seven columns are the first draft (2026-11-09). The
review added `celltype_label`, `reviewer_call`, `review` (agree / verify /
disagree) and `reviewer_note`. The draft columns are never edited; the
reviewer columns win where the two disagree. `check_result` records checks
run afterwards on existing v3 outputs (scDblFinder scores, nCount, dataset
enrichment over baseline). `action`, `route_to`, `followup` and
`lineage_reviewed` are what `rebuild_02a_routing.R` acts on. They follow
`reviewer_call` unless a later decision overrides it, in which case
`check_result` says why. Example: TNK 16 was routed to myeloid in the review,
then dropped after DecontX showed it is mostly red-cell / low-quality droplets.
The reasoning and order of the cleanup are in the
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
| `lineage` | what the cluster is. Kept: `fibroblast`, `mural`, `schwann`, `Tcell`, `NK`, `ILC`, `Bcell`, `plasma`, `myeloid`, `mast`, `neutrophil`, `pDC`, `endothelial`. Dropped rows record what was found instead (`Tcell`, `endothelial`, `RBC_lowquality`, `pDC`, ...) |
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
