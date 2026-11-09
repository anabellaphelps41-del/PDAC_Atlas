# PDAC atlas rebuild — code

Everything that produces `PDAC_Atlas/rebuild_v3/` and `PDAC_Atlas/rebuild_v4/`
lives here. The data does not: `.gitignore` tracks code and docs only.

## Pipeline

Steps 00m–01c are atlas-level and unchanged since v3; their outputs stay in
`rebuild_v3/`. v4 (from 2026-11-09) starts at 02a and writes to `rebuild_v4/`.
The v3 pass of 02–04b (three subsets, raw counts) is frozen in `rebuild_v3/`,
because the reviewed cluster calls in [`calls/`](calls/) refer to its res-0.6
cluster numbers.

| Step | Script | Reads | Writes |
|---|---|---|---|
| 00m | `rebuild/rebuild_00m_manifest.R` | `inspect/atlas_metadata.csv.gz` | v3: `sample_manifest.tsv` — one row per sample, `include_in_atlas` + `exclude_reason` |
| 01 | `rebuild/rebuild_01_harmonize.R` | original 52 GB atlas, `ref/hgnc_complete_set.txt` | v3: `atlas_harmonized_counts_v3.qs2` (counts only), `feature_map_v3.tsv`, `feature_codetect_v3.tsv` |
| 01c | `rebuild/rebuild_01c_doublets.R` | v3 atlas, manifest | v3: `doublets_scDblFinder_v3.tsv`, `doublets_by_sample_v3.tsv` |
| 02a | `rebuild/rebuild_02a_routing.R` | `calls/*_reviewed.tsv`, v3 `step4b_*_percell.tsv` | v4: `cell_routing.tsv` — keep / drop / route per cell |
| 01d | `rebuild/rebuild_01d_decontx.R` | v3 atlas, manifest, doublets, routing | v4: `decontx_counts.qs2`, `decontx_percell.tsv` (contamination, QC floor) |
| 02 | `rebuild/rebuild_02_subset.R` | v3 atlas, manifest, doublets, v4 DecontX + routing | v4: `subset_{TNK,stromal,Bcell,myeloid,endothelial}.qs2` |
| 02b | `rebuild/rebuild_02b_batchgenes.R` | subsets | `exclude_batchgenes.rds`, `batchgene_detection_*.tsv` |
| 03 | `rebuild/rebuild_03_pipeline.R <subset>` | subset, batch genes | `subset_<subset>_harmony.qs2`, QC pdf, cluster composition |
| 04 | `rebuild/rebuild_04_markers.R <subset> [res]` | harmony object | markers, decision table (**re-saves its input**) |
| 04b | `rebuild/rebuild_04b_singler.R <subset> [res]` | harmony object, `ref/celldex_*.rds` | SingleR lineage gate (**re-saves its input**) |
| 04c | `rebuild/rebuild_04c_subcluster.R <subset> <res> <cluster> <markersA> <markersB>` | harmony object | `subcluster_res<res>` column, marker-defined (**re-saves its input**) |
| 05 | `rebuild/rebuild_05_apply_calls.R <subset> [res]` | harmony object (read only), `calls/cluster_calls_<subset>_v4_res<res>.tsv`, `calls/v4_run2_decisions.tsv`, `calls/subcluster_labels_v4_res<res>.tsv` | `subset_<subset>_final.qs2` (kept cells, `celltype_v4` / `lineage_v4`, flags), `step5_*_cells.tsv`, `step5_*_summary.tsv` |
| 00 | `rebuild/rebuild_00_provenance_audit.R` | subsets | metadata audit tables |

02a runs before 01d because DecontX takes its cell-type labels from the
reviewed lineages. Every step writes a `PROVENANCE_<step>.txt` next to its
outputs. Directories are set in `rebuild/io.R` (`PDAC_UPSTREAM_DIR` for the v3
inputs, `PDAC_REBUILD_DIR` for the v4 outputs).

```bash
bash code/slurm/run_v3.sh              # 00m -> 01 -> 01c (atlas level)
bash code/slurm/run_v4.sh              # 02a -> 01d -> 02 -> 02b -> (03 -> 04 -> 04b) x 5
bash code/slurm/run_v4.sh 03           # resume from step 03
RES=0.8 bash code/slurm/run_v4.sh 04   # markers + gate at another resolution
bash code/slurm/run_v4.sh 05           # apply the reviewed calls (after review only)
PDAC_DRAFT_CALLS=1 bash code/slurm/run_v4.sh 05   # draft: unreviewed rows use the draft call
```

Step 05 is never chained after 04b, because the call tables are reviewed in
between. It stops in four cases:
- the call table's clusters or `n_cells` don't match the object (03 was
  rerun);
- a cluster has no `reviewer_call` (unless it is a draft run);
- a call is `route`, which needs a rerun from 02a;
- `v4_run2_decisions.tsv` disagrees with the table.

A draft run writes to `rebuild_v4/step5_draft/`, never next to the reviewed
outputs.

04 and 04b re-save the harmony object, so run_v4.sh chains them per subset;
never run two of them on the same subset at once. Logs:
`PDAC_Atlas/logs/<jobname>_<jobid>.{out,err}`.

SingleR references are cached once on a node with internet:
`Rscript code/setup_singler_refs.R` → `PDAC_Atlas/ref/`.

## Decisions encoded here

- **Which cells**: `sample_manifest.tsv` is the only place inclusion is
  decided. Kept: primary tumour, adjacent normal, CRA001160 normal pancreas.
  Excluded: peripheral blood, metastases (GSE205013 liver, GSE154778 MET),
  non-PDAC tumours, GSE211644 (sorted CD8-biased T cells), GSE202051.
- **Keys**: sample and patient names repeat across datasets. Join on
  `sample_key` / `patient_key` (`dataset:id`), never on bare IDs.
  `patient_id_source` says whether a patient key came from metadata, was
  inferred from the sample name, or is a one-sample proxy.
- **Gene symbols**: HGNC `prev_symbol` (unique) and Ensembl ID only, never
  `alias_symbol`. Two rows are only summed if they are never both detected in
  one dataset; rejected merges are listed in `feature_codetect_v3.tsv`.
- **Doublets**: scDblFinder per `sample_key` on the whole atlas, before any
  compartment split.
- **Ambient RNA (v4)**: DecontX per `sample_key` on the whole atlas, before
  the split, so that every contamination source is in the run. There are no
  raw droplet matrices, so it runs without a background profile, on all five
  datasets (GSE205013/GSE212966 were also SoupX'd by their authors).
  Cell-type labels are lineages from the reviewed calls (B and plasma
  separate, T and NK separate, fibroblast and mural separate), then level3,
  then level2. With level2 `B/plasma` as one label, B cells came out ~98%
  "contamination". The `decontX` assay is the default from step 02 on; `RNA`
  keeps raw counts and the QC metrics.
- **QC floor (v4)**: cells more than 3 MADs below the median log10 nFeature
  or nCount of their dataset × compartment fail. No percent.mt filter.
- **Reviewed calls (v4)**: `calls/*_reviewed.tsv` decide, per v3 res-0.6
  cluster, keep / drop / route to another subset. Misrouted cells are moved,
  not deleted. Myeloid (with mast) and endothelial subsets exist from v4 to
  receive them. `Neutrophil-like` cells are no longer removed. The `action`
  column can override `reviewer_call` when later evidence changes a decision;
  `check_result` then says why. Example: TNK 16 was first routed to myeloid
  as neutrophils, then dropped. After DecontX it showed a median 21.9%
  hemoglobin and 0.67 contamination, and sharing the neutrophil label had
  made hemoglobin native to the GSE205013 neutrophil profile.
- **Batch key**: Harmony on `dataset_id` (`PDAC_HARMONY_KEY` to change). Step
  03 also runs Harmony on `sample_key` as a comparison and reports mixing,
  same-tissue neighbours and ARI.
- **Batch genes**: 02b masks STRUCTURAL (absent from a dataset's reference)
  and GRADED (≥10× detection spread) genes from HVG selection. A per-subset
  `PROTECT` list in `rebuild_02b_batchgenes.R` keeps lineage markers that are
  graded only because of the small GSE154778 cohort (stromal mural markers
  such as DES/MUSTN1; T/NK CD40LG, PTGDR, IL18RAP, ADGRG1, MYBL1), or in the
  v4 myeloid/endothelial subsets because cell-type composition differs between
  datasets (neutrophil, mast, pDC, DC, lymphatic/arterial/venous/tip markers). Cell-cycle
  genes stay masked on purpose. Step 04 tags batch genes in its marker tables
  (`batch_class`, and `(S)`/`(G)` in the decision table) instead of dropping
  them.
- **Storage**: `/ix` reads at ~27 MB/s, so file size is load time. Seurat
  objects are saved as `qs2` (zstd, multithreaded) through `rebuild/io.R`:
  ~6× smaller than uncompressed `.rds` and loaded in seconds rather than
  minutes. `read_obj()` still reads a legacy `.rds` of the same name.
  `save_obj()` writes to a temp file and renames, so a crash mid-write cannot
  corrupt the existing copy. Step 1 stores counts only.

## Layout

```
code/
  rebuild/   pipeline scripts (above)
  slurm/     one sbatch per step + run_v3.sh (atlas level) + run_v4.sh
  calls/     reviewed per-cluster decisions (human input to 02a)
  legacy/    superseded scripts kept for the record (v1 inspect/subset,
             01b collapse fix, harmonization diagnostics)
  fetch_atlas.sh, setup_singler_refs.R
```
