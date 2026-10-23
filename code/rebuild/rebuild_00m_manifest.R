#!/usr/bin/env Rscript
# =============================================================================
# rebuild_00m_manifest.R   --   one row per sample, the single source of truth
#                               for which cells may enter the atlas
#
# Input : inspect/atlas_metadata.csv.gz   (all 212 columns, from 01_inspect_atlas.R)
# Output: <DIR>/sample_manifest.tsv
#         <DIR>/PROVENANCE_00m_manifest.txt
#
# Why this exists: the per-cell metadata carries sample context in several
# inconsistent columns (tissue has 9 spellings for 5 concepts,
# master_primary_or_metastatic is "Unknown" for every GSE205013 cell including
# the liver biopsies, patient_id is empty for three datasets). Step 2 used to
# subset on cell type alone, so peripheral blood, liver metastases and non-PDAC
# tumours went straight into the compartments. Every downstream exclusion now
# reads include_in_atlas from this table, and the reason is recorded per sample.
#
# Sample and patient names are NOT unique across datasets (P01-P10 exist in both
# GSE154778 and GSE205013). Always join on sample_key / patient_key.
# =============================================================================

suppressPackageStartupMessages(library(data.table))
options(width = 220, datatable.print.class = FALSE)

ROOT <- "/ix/tbruno/anp795/PDAC_Atlas"
META <- file.path(ROOT, "inspect", "atlas_metadata.csv.gz")
DIR  <- Sys.getenv("PDAC_REBUILD_DIR", file.path(ROOT, "rebuild_v3"))
OUT  <- file.path(DIR, "sample_manifest.tsv")
PROV <- file.path(DIR, "PROVENANCE_00m_manifest.txt")
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

# ============================ inclusion policy ================================
# Edit here, rerun, and everything downstream follows.
KEEP_TISSUE <- c("primary_tumor", "adjacent_normal", "normal_pancreas")
EXCLUDE_DATASETS <- c(
  GSE202051_main      = "GSE202051 dropped at step 1 (snRNA-seq; original rationale not recorded)",
  GSE202051_010nuc    = "GSE202051 dropped at step 1 (nuclei side sample)",
  GSE202051_010orgCRT = "GSE202051 dropped at step 1 (organoid)",
  GSE211644           = "GSE211644: sorted, CD8-biased T cell dataset; dominated 9 TNK clusters at res 0.6"
)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# ============================ 1. load ========================================
cols <- c("dataset_id", "sample_id", "patient_id", "tissue", "condition",
          "master_diagnosis", "master_treatment_status", "master_sequencing_modality")
msg("reading ", META)
md <- fread(cmd = paste("zcat", shQuote(META)), select = cols, na.strings = c("NA", ""))
msg("cells: ", format(nrow(md), big.mark = ","))

first_val <- function(x) { x <- unique(x[!is.na(x)]); if (length(x)) paste(sort(x), collapse = "|") else NA_character_ }
s <- md[, .(n_cells                 = .N,
            patient_id_raw          = first_val(patient_id),
            tissue_raw              = first_val(tissue),
            condition_raw           = first_val(condition),
            master_diagnosis        = first_val(master_diagnosis),
            treatment_status        = first_val(master_treatment_status),
            sequencing_modality     = first_val(master_sequencing_modality)),
        by = .(dataset_id, sample_id)]
setorder(s, dataset_id, sample_id)
s[, sample_key := paste(dataset_id, sample_id, sep = ":")]

multi <- s[grepl("|", tissue_raw, fixed = TRUE) | grepl("|", condition_raw, fixed = TRUE) |
           grepl("|", patient_id_raw, fixed = TRUE)]
if (nrow(multi)) {
  print(multi)
  stop("some samples carry more than one tissue/condition/patient value; resolve before building the manifest")
}

# ============================ 2. tissue class ================================
tissue_map <- c(
  "Pancreatic tumor"         = "primary_tumor",
  "Primary pancreatic tumor" = "primary_tumor",
  "Tumor"                    = "primary_tumor",
  "Pancreas"                 = "primary_tumor",   # GSE205013 pancreas biopsies (Werba 2023); assumed tumour
  "Adjacent pancreas"        = "adjacent_normal",
  "Normal pancreas"          = "normal_pancreas",
  "Peripheral blood"         = "blood",
  "Liver"                    = "liver_metastasis",
  "Metastatic lesion"        = "metastasis_site_unspecified"
)
unmapped <- setdiff(na.omit(unique(s$tissue_raw)), names(tissue_map))
if (length(unmapped)) stop("tissue values with no mapping: ", paste(unmapped, collapse = ", "))
s[, tissue_class := unname(tissue_map[tissue_raw])]
s[is.na(tissue_class) & grepl("org", sample_id, ignore.case = TRUE), tissue_class := "organoid"]
s[is.na(tissue_class), tissue_class := "unknown"]

s[, primary_or_metastatic := fcase(
  tissue_class == "primary_tumor",                    "Primary",
  tissue_class %in% c("liver_metastasis",
                      "metastasis_site_unspecified"), "Metastatic",
  tissue_class %in% c("adjacent_normal",
                      "normal_pancreas"),             "Normal",
  tissue_class == "blood",                            "Blood",
  default = "Unknown")]

# ============================ 3. patient diagnosis ===========================
# condition is about the PATIENT for adjacent normal (a PDAC patient's tissue)
# but about the SAMPLE for CRA001160 normals (non-PDAC patients).
s[, patient_diagnosis := fcase(
  condition_raw %in% c("PDAC", "Adjacent normal"),  "PDAC",
  condition_raw == "Healthy",                       "healthy_donor",
  condition_raw == "Non-malignant pancreas",        paste0("non-PDAC: ", master_diagnosis),
  !is.na(condition_raw),                            condition_raw,
  default = "unknown")]

# ============================ 4. patient key =================================
# Prefer the metadata patient_id; otherwise infer from the sample name and say
# so. Inferred pairings are NOT verified against the source papers.
s[, `:=`(patient_key = NA_character_, patient_id_source = NA_character_)]
s[!is.na(patient_id_raw), `:=`(patient_key = paste(dataset_id, patient_id_raw, sep = ":"),
                               patient_id_source = "metadata")]

infer <- function(ds, rx, repl, note) {
  i <- s$dataset_id == ds & is.na(s$patient_key) & grepl(rx, s$sample_id)
  s[i, `:=`(patient_key = paste(ds, sub(rx, repl, sample_id), sep = ":"),
            patient_id_source = note)]
}
infer("GSE211644", "^(MDA[12])_[TU]([0-9]+)$", "\\1_\\2", "inferred_from_sample_name")
infer("GSE155698", "^PDAC_(TISSUE|PBMC)_([0-9]+)[AB]?$", "PDAC_\\2",
      "inferred_from_sample_name (tissue/PBMC pairing by number unverified)")
infer("GSE155698", "^Healthy_PBMC_([0-9]+)$", "Healthy_\\1", "inferred_from_sample_name")
infer("GSE154778", "^(P[0-9]+)$", "\\1", "inferred_from_sample_name")
s[is.na(patient_key), `:=`(patient_key = sample_key, patient_id_source = "sample_proxy")]

# ============================ 5. inclusion ===================================
s[, exclude_reason := ""]
add_reason <- function(i, why) {
  why <- rep_len(why, nrow(s))[i]
  s[i, exclude_reason := ifelse(nzchar(exclude_reason), paste(exclude_reason, why, sep = "; "), why)]
}
for (d in names(EXCLUDE_DATASETS)) add_reason(s$dataset_id == d, EXCLUDE_DATASETS[[d]])
add_reason(s$tissue_class == "blood", "peripheral blood")
add_reason(s$tissue_class %in% c("liver_metastasis", "metastasis_site_unspecified"), "metastasis")
add_reason(s$tissue_class %in% c("unknown", "organoid"), paste0("tissue: ", s$tissue_class))
add_reason(s$tissue_class %in% c("primary_tumor", "adjacent_normal") & s$patient_diagnosis != "PDAC",
           paste0("non-PDAC: ", s$patient_diagnosis))
add_reason(!s$tissue_class %in% c(KEEP_TISSUE, "blood", "liver_metastasis",
                                  "metastasis_site_unspecified", "unknown", "organoid"),
           paste0("tissue class not kept: ", s$tissue_class))
s[, include_in_atlas := !nzchar(exclude_reason)]

# ============================ 6. write + report ==============================
out_cols <- c("sample_key", "dataset_id", "sample_id", "n_cells", "patient_key", "patient_id_source",
              "tissue_class", "primary_or_metastatic", "patient_diagnosis", "treatment_status",
              "sequencing_modality", "include_in_atlas", "exclude_reason",
              "tissue_raw", "condition_raw", "patient_id_raw")
fwrite(s[, ..out_cols], OUT, sep = "\t", na = "NA")
msg("wrote ", OUT, " (", nrow(s), " samples)")

tab_cells <- function(by) s[, .(samples = .N, cells = sum(n_cells)), by = by][order(get(by[1]))]
inc <- s[, .(samples = .N, cells = sum(n_cells)), by = include_in_atlas]
by_ds <- dcast(s, dataset_id ~ include_in_atlas, value.var = "n_cells", fun.aggregate = sum)
setnames(by_ds, c("FALSE", "TRUE"), c("cells_excluded", "cells_included"), skip_absent = TRUE)
reasons <- s[!(include_in_atlas), .(samples = .N, cells = sum(n_cells)), by = exclude_reason][order(-cells)]
pat <- s[(include_in_atlas), .(patients = uniqueN(patient_key), samples = .N), by = .(dataset_id, patient_id_source)]

show <- function(x) capture.output(print(x, row.names = FALSE))
report <- c(
  "PROVENANCE - sample manifest",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_00m_manifest.R",
  paste("input    :", META),
  paste("output   :", OUT),
  "",
  "POLICY",
  paste("  tissue classes kept :", paste(KEEP_TISSUE, collapse = ", ")),
  "  datasets excluded   :",
  paste0("    ", names(EXCLUDE_DATASETS), " - ", EXCLUDE_DATASETS),
  "  also excluded       : non-PDAC tumours (tumour/adjacent tissue from non-PDAC patients).",
  "  normal_pancreas (CRA001160) comes from non-PDAC patients; kept as reference tissue.",
  "",
  "INCLUDED vs EXCLUDED", show(inc), "",
  "BY DATASET (cells)", show(by_ds), "",
  "EXCLUSION REASONS", show(reasons), "",
  "TISSUE CLASS", show(tab_cells(c("tissue_class", "include_in_atlas"))), "",
  "PATIENT KEYS IN INCLUDED SAMPLES", show(pat), "",
  "CAVEATS",
  "  - 'Pancreas' (GSE205013) is assumed to be primary tumour.",
  "  - Patient keys marked inferred_from_sample_name are pattern-based, not",
  "    checked against the source papers. sample_proxy means one sample = one",
  "    patient was assumed; patient-level statistics on those are sample-level.",
  "  - treatment_status is only populated for GSE211644 (and GSE202051, dropped)."
)
writeLines(report, PROV)
cat(report, sep = "\n")
msg("wrote ", PROV)
