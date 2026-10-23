# PDAC atlas rebuild — code

Everything that produces `PDAC_Atlas/rebuild_v3/` lives here. The data does
not: `.gitignore` tracks code and docs only.

## Pipeline (v3)

| Step | Script | Reads | Writes |
|---|---|---|---|
| 00m | `rebuild/rebuild_00m_manifest.R` | `inspect/atlas_metadata.csv.gz` | `sample_manifest.tsv` — one row per sample, `include_in_atlas` + `exclude_reason` |
| 01 | `rebuild/rebuild_01_harmonize.R` | original 52 GB atlas, `ref/hgnc_complete_set.txt` | `atlas_harmonized_counts_v3.rds` (counts only), `feature_map_v3.tsv`, `feature_codetect_v3.tsv` |
| 01c | `rebuild/rebuild_01c_doublets.R` | v3 atlas, manifest | `doublets_scDblFinder_v3.tsv`, `doublets_by_sample_v3.tsv` |
| 02 | `rebuild/rebuild_02_subset.R` | v3 atlas, manifest, doublets | `subset_{TNK,stromal,Bcell}.rds` |
| 02b | `rebuild/rebuild_02b_batchgenes.R` | subsets | `exclude_batchgenes.rds`, `batchgene_detection_*.tsv` |
| 03 | `rebuild/rebuild_03_pipeline.R <subset>` | subset, batch genes | `subset_<subset>_harmony.rds`, QC pdf, cluster composition |
| 04 | `rebuild/rebuild_04_markers.R <subset> [res]` | harmony object | markers, decision table (**re-saves its input**) |
| 04b | `rebuild/rebuild_04b_singler.R <subset> [res]` | harmony object, `ref/celldex_*.rds` | SingleR lineage gate (**re-saves its input**) |
| 00 | `rebuild/rebuild_00_provenance_audit.R` | subsets | metadata audit tables |

Every step writes a `PROVENANCE_<step>.txt` next to its outputs.

Submit the chain (00m → 01 → 01c → 02 → 02b → 03 ×3) with SLURM dependencies:

```bash
bash code/slurm/run_v3.sh              # everything
bash code/slurm/run_v3.sh 02           # resume from step 02
AFTER=<jobid> bash code/slurm/run_v3.sh 01c
```

Logs: `PDAC_Atlas/logs/<jobname>_<jobid>.{out,err}`. Output directory defaults
to `PDAC_Atlas/rebuild_v3`; set `PDAC_REBUILD_DIR` to redirect a run.

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
- **Storage**: objects are saved uncompressed for load speed; the `/ix`
  filesystem compresses them ~2.5× on disk. Step 1 stores counts only.

## Layout

```
code/
  rebuild/   pipeline scripts (above)
  slurm/     one sbatch per step + run_v3.sh
  legacy/    superseded scripts kept for the record (v1 inspect/subset,
             01b collapse fix, harmonization diagnostics)
  fetch_atlas.sh, setup_singler_refs.R
```
