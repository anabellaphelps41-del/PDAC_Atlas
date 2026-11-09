# Cluster calls

Two rounds so far:

- `cluster_calls_{TNK,stromal,Bcell}_res0.6_reviewed.tsv`: v3 clusters, reviewed;
  applied by `rebuild_02a_routing.R` to build v4.
- `cluster_calls_{TNK,stromal,Bcell,myeloid,endothelial}_v4run1_res0.6.tsv`:
  draft calls on the first v4 clustering (`rebuild_v4/run1_snapshot/`). They
  were never reviewed row by row, and are kept for the record. The columns are
  the reviewed columns plus `v3_origin` and `followup`.
- `cluster_calls_{TNK,stromal,Bcell,myeloid,endothelial}_v4_res0.6.tsv`: DRAFT
  calls on v4 run 2 (wider stress mask), awaiting review. The columns are:
  - `call` (keep / drop / route) and `route_to`;
  - `confidence`: `decided` for rows that carry a user decision from run 1;
  - `run1_origin`: which run-1 clusters the cells came from;
  - `samples`: samples with at least 5% of the cluster, and the top three;
  - `cell_flags`: the share of cells carrying each run-1 flag;
  - empty reviewer columns.

  Myeloid c14 is split in `subcluster_res0.6` (14a cDC1 / 14b mregDC) by
  `rebuild_04c_subcluster.R`.

  Reviewer columns are filled only for the rows decided on 2026-11-09
  (myeloid c16, TNK c5, stromal c7, endothelial c11/c12); the rest are
  still drafts.
- `v4_run2_decisions.tsv`: the machine-readable form of those run-2
  decisions, in the run-1 format plus `n_cells`:
  - drops: myeloid c16, endothelial c11 and c12. DecontX is not rerun for
    the endothelial ones;
  - a whole-cluster `stress_overlay` flag with `exclude_from_de_lr` on
    TNK c5. Cells that already carry a flag keep it.

  Stromal c7 gets no cluster flag: its cells keep their run-1 flags only.
  These apply to the run-2 clusters as they are, without re-clustering.
  Nothing reads the file yet. Step 05 will, and it must stop if
  `n_cells` no longer matches, because a rerun of 03 renumbers the clusters.
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
