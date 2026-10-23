#!/usr/bin/env Rscript
# =============================================================================
# rebuild_01_harmonize.R   (v3)   --   Step 1 of the PDAC atlas rebuild
#
# Drop GSE202051, harmonize gene symbols across datasets, sum counts.
#
# Input : combined_9_datasets_hierarchical_annotation_seurat.rds  (the 52 GB original)
#         ref/hgnc_complete_set.txt
# Output: <DIR>/atlas_harmonized_counts_v3.qs2   counts layer only, no data layer
#         <DIR>/feature_map_v3.tsv               every original feature -> final feature, rule, decision
#         <DIR>/feature_codetect_v3.tsv          merges rejected because both rows were detected in one dataset
#         <DIR>/PROVENANCE_01_harmonize.txt
#
# Replaces the v1 script (HGNC alias_symbol merges, first-match-wins) and
# rebuild_01b_fixcollapse.R. alias_symbol is never used: HGNC alias fields list
# unrelated genes (AMN/ABCD1, AR/AKR1B1, MIF/AMH).
#
# Rules, applied per original feature name, first match wins:
#   approved        already a current HGNC symbol                         -> unchanged
#   ensembl_id      "<SYM>--ENSGxxx" and ENSGxxx is in HGNC                -> HGNC symbol for that ENSG
#   ensembl_strip   "<SYM>--ENSGxxx", ENSG not in HGNC (alt haplotype etc) -> resolve <SYM>
#   prev_symbol     a prev_symbol of exactly one approved gene             -> that gene
#   dup_suffix      make.unique "<SYM>.N" whose <SYM> resolves             -> resolve <SYM>
#   unresolved      none of the above (clone IDs, lncRNAs, retired names)  -> unchanged
#
# Merge guard: two rows that end up with the same name are only summed if they
# behave like one gene split across references, i.e. never both detected
# (>= CODETECT_PCT of cells) in the same sample. A source row that is
# co-detected with its group's anchor keeps its original name and is logged.
#
# The data layer is NOT written. LogNormalize is per cell and cheap; step 2
# normalizes each subset. This halves the file.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")   # read_obj() / save_obj()

ROOT <- "/ix/tbruno/anp795/PDAC_Atlas"
IN   <- Sys.getenv("PDAC_ATLAS_IN", file.path(ROOT, "combined_9_datasets_hierarchical_annotation_seurat.rds"))
HGNC <- file.path(ROOT, "ref", "hgnc_complete_set.txt")
DIR  <- Sys.getenv("PDAC_REBUILD_DIR", file.path(ROOT, "rebuild_v3"))
OUT  <- file.path(DIR, "atlas_harmonized_counts_v3.qs2")
FMAP <- file.path(DIR, "feature_map_v3.tsv")
FCOD <- file.path(DIR, "feature_codetect_v3.tsv")
PROV <- file.path(DIR, "PROVENANCE_01_harmonize.txt")

DROP         <- c("GSE202051_main", "GSE202051_010nuc", "GSE202051_010orgCRT")
CODETECT_PCT <- 0.5
KEEP_META    <- c("orig.ident", "dataset_id", "dataset_accession", "sample_id", "patient_id",
                  "condition", "tissue", "original_barcode", "nCount_RNA", "nFeature_RNA",
                  "percent.mt", "percent.ribo",
                  "celltype_level1", "celltype_level2", "celltype_level3")

dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
options(width = 200)
msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = ""); flush.console() }

# ============================ 1. HGNC lookups ================================
msg("reading ", HGNC)
h <- read.delim(HGNC, quote = "", stringsAsFactors = FALSE, colClasses = "character")
h <- h[h$status == "Approved", ]
hgnc_date <- max(h$date_modified, na.rm = TRUE)
approved  <- h$symbol
ens_map   <- setNames(h$symbol, h$ensembl_gene_id)[nzchar(h$ensembl_gene_id)]

# Multi-valued HGNC fields are double-quoted ("ATPSB|ATP5B"). With quote = ""
# the quotes survive and stick to the first and last symbol -- the reason v1
# left ~250 prev_symbol pairs (ATP5B/ATP5F1B, SEPT7/SEPTIN7) unmerged.
pv  <- strsplit(gsub('"', "", h$prev_symbol, fixed = TRUE), "|", fixed = TRUE)
prev <- data.frame(old = trimws(unlist(pv)), new = rep(h$symbol, lengths(pv)), stringsAsFactors = FALSE)
prev <- unique(prev[nzchar(prev$old) & !(prev$old %in% approved), ])
amb  <- unique(prev$old[duplicated(prev$old)])
prev <- prev[!(prev$old %in% amb), ]
prev_map <- setNames(prev$new, prev$old)
msg("HGNC approved: ", length(approved), " | unique prev_symbol: ", length(prev_map),
    " | ambiguous prev_symbol (skipped): ", length(amb), " | file date_modified max: ", hgnc_date)

resolve <- function(x) {
  out  <- x
  rule <- rep("unresolved", length(x))
  a <- x %in% approved
  rule[a] <- "approved"
  p <- !a & x %in% names(prev_map)
  out[p] <- prev_map[x[p]]; rule[p] <- "prev_symbol"
  d <- !a & !p & grepl("^(.+)\\.[0-9]+$", x)
  if (any(d)) {
    base <- sub("^(.+)\\.[0-9]+$", "\\1", x[d])
    ok   <- base %in% approved | base %in% names(prev_map)
    tgt  <- ifelse(base %in% approved, base, prev_map[base])
    idx  <- which(d)[ok]
    out[idx] <- tgt[ok]; rule[idx] <- "dup_suffix"
  }
  list(target = out, rule = rule)
}

# ============================ 2. load + drop =================================
msg("loading ", IN, " (expect ~40 min)")
obj <- readRDS(IN)
msg("loaded: ", format(nrow(obj), big.mark = ","), " features x ", format(ncol(obj), big.mark = ","), " cells")

keep_cells <- colnames(obj)[!obj$dataset_id %in% DROP]
msg("dropping ", format(ncol(obj) - length(keep_cells), big.mark = ","), " cells from ",
    paste(DROP, collapse = ", "), " -> retained ", format(length(keep_cells), big.mark = ","))

md  <- obj@meta.data[keep_cells, intersect(KEEP_META, colnames(obj@meta.data))]
cts <- LayerData(obj, assay = "RNA", layer = "counts")[, keep_cells]
rm(obj); gc(verbose = FALSE)
msg("counts extracted: ", format(nrow(cts), big.mark = ","), " x ", format(ncol(cts), big.mark = ","),
    ", nnz ", format(length(cts@x), big.mark = ","))
if (!all(cts@x == round(cts@x))) msg("WARNING: counts layer contains non-integer values")

# ============================ 3. map features ================================
g  <- rownames(cts)
fm <- data.frame(feature = g, target = g, rule = "unresolved", stringsAsFactors = FALSE)

is_ens <- grepl("^(.+)--(ENSG[0-9]+)(\\.[0-9]+)?$", g)
ens_id <- sub("^(.+)--(ENSG[0-9]+)(\\.[0-9]+)?$", "\\2", g[is_ens])
ens_sy <- sub("^(.+)--(ENSG[0-9]+)(\\.[0-9]+)?$", "\\1", g[is_ens])

r <- resolve(g[!is_ens])
fm$target[!is_ens] <- r$target
fm$rule[!is_ens]   <- r$rule

in_hgnc <- ens_id %in% names(ens_map)
r2 <- resolve(ens_sy)
fm$target[is_ens] <- ifelse(in_hgnc, ens_map[ens_id], r2$target)
fm$rule[is_ens]   <- ifelse(in_hgnc, "ensembl_id", "ensembl_strip")

# Seurat replaces underscores on object creation; do it here so the map is final
fm$target <- gsub("_", "-", fm$target, fixed = TRUE)
msg("rule tally:"); print(table(fm$rule))

# ============================ 4. co-detection guard ==========================
# Checked per SAMPLE, not per dataset: one sample is one alignment run, but a
# dataset can mix reference builds. GSE154778:P01 was aligned to a newer
# reference (ATP5MC1, MICOS10, SEPTIN7) and its other 14 samples to an older
# one (ATP5G1, MINOS1, SEPT7). Per dataset those pairs look co-detected; per
# sample they are perfectly complementary, i.e. one gene.
ds  <- factor(md$dataset_id)
ind <- Matrix::sparseMatrix(i = seq_along(ds), j = as.integer(ds), x = 1,
                            dims = c(length(ds), nlevels(ds)), dimnames = list(NULL, levels(ds)))
n_ds <- Matrix::colSums(ind)
sk   <- factor(paste(md$dataset_id, md$sample_id, sep = ":"))
ind_s <- Matrix::sparseMatrix(i = seq_along(sk), j = as.integer(sk), x = 1,
                              dims = c(length(sk), nlevels(sk)), dimnames = list(NULL, levels(sk)))
n_s  <- Matrix::colSums(ind_s)

fm$decision <- "kept"
codetect_log <- list()
for (pass in 1:5) {
  grp   <- split(seq_len(nrow(fm)), fm$target)
  multi <- grp[lengths(grp) > 1]
  if (!length(multi)) break
  rows <- unique(unlist(multi))
  det  <- as.matrix((cts[rows, , drop = FALSE] > 0) %*% ind_s)
  det  <- sweep(det, 2, n_s, "/") * 100
  rownames(det) <- fm$feature[rows]

  rejected <- integer(0)
  for (tg in names(multi)) {
    m   <- multi[[tg]]
    anc <- m[fm$feature[m] == tg]
    if (!length(anc)) anc <- m[which.max(Matrix::rowSums(det[fm$feature[m], , drop = FALSE]))]
    anc <- anc[1]
    for (s in setdiff(m, anc)) {
      both <- det[fm$feature[s], ] >= CODETECT_PCT & det[fm$feature[anc], ] >= CODETECT_PCT
      if (any(both)) {
        rejected <- c(rejected, s)
        codetect_log[[length(codetect_log) + 1]] <- data.frame(
          pass = pass, source = fm$feature[s], anchor = fm$feature[anc], target = tg,
          rule = fm$rule[s], samples_both_detected = paste(names(both)[both], collapse = ";"),
          source_pct = paste(round(det[fm$feature[s], both], 1), collapse = ";"),
          anchor_pct = paste(round(det[fm$feature[anc], both], 1), collapse = ";"),
          stringsAsFactors = FALSE)
      }
    }
  }
  msg("co-detection pass ", pass, ": ", length(multi), " merge groups, ", length(rejected), " sources rejected")
  if (!length(rejected)) break
  fm$target[rejected]   <- gsub("_", "-", fm$feature[rejected], fixed = TRUE)
  fm$decision[rejected] <- "rejected_codetected"
}
grp <- split(seq_len(nrow(fm)), fm$target)
gsz <- lengths(grp)[fm$target]
fm$decision[fm$decision == "kept" & gsz > 1]                          <- "merged"
fm$decision[fm$decision == "kept" & gsz == 1 & fm$target != fm$feature] <- "renamed"
fm$group_size <- as.integer(gsz)
codetect <- if (length(codetect_log)) do.call(rbind, codetect_log) else
  data.frame(pass = integer(), source = character(), anchor = character(), target = character())

write.table(fm, FMAP, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(codetect, FCOD, sep = "\t", quote = FALSE, row.names = FALSE)
msg("wrote ", FMAP, " and ", FCOD)
msg("decision x rule:"); print(table(fm$rule, fm$decision))

# ============================ 5. aggregate ===================================
total_before <- sum(cts@x)
cells_before <- Matrix::colSums(cts)
f   <- factor(fm$target, levels = unique(fm$target))
agg <- Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                            dims = c(length(f), nlevels(f)), dimnames = list(NULL, levels(f)))
new_cts <- as(Matrix::t(agg) %*% cts, "CsparseMatrix")
colnames(new_cts) <- colnames(cts)
rm(cts, agg); gc(verbose = FALSE)

if (!isTRUE(all.equal(total_before, sum(new_cts@x)))) stop("TOTAL COUNTS CHANGED during aggregation")
if (!isTRUE(all.equal(unname(cells_before), unname(Matrix::colSums(new_cts))))) stop("PER-CELL TOTALS CHANGED")
if (anyDuplicated(rownames(new_cts))) stop("duplicate feature names after aggregation")
msg("verified: total and per-cell UMI identical; no duplicate names. features ",
    format(nrow(fm), big.mark = ","), " -> ", format(nrow(new_cts), big.mark = ","))

spot <- intersect(c("ATP5F1B", "SELENOM", "SEPTIN7", "ACKR1", "PSMA2", "MICOS10", "SNU13",
                    "HLA-A", "HLA-DRB1", "TRAC", "IGKC", "NCR1"), rownames(new_cts))
det_spot <- sweep(as.matrix((new_cts[spot, , drop = FALSE] > 0) %*% ind), 2, n_ds, "/") * 100
msg("spot check, % cells detecting per dataset (a 0 is a structural absence or an unmerged split):")
print(round(det_spot, 1))

# ============================ 6. build + save ================================
new <- CreateSeuratObject(counts = new_cts, meta.data = md, min.cells = 0, min.features = 0)
rm(new_cts); gc(verbose = FALSE)
d <- max(abs(new$nCount_RNA - md$nCount_RNA))
msg("max |nCount_RNA change| vs source metadata: ", d)

msg("writing ", OUT, " (qs2)")
save_obj(new, OUT)
msg("wrote ", sprintf("%.2f GB", file.size(OUT) / 2^30))

show <- function(x) capture.output(print(x))
writeLines(c(
  "PROVENANCE - step 1 harmonization (v3)",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_01_harmonize.R (v3)",
  paste("input    :", IN),
  paste("output   :", OUT),
  paste("HGNC     :", HGNC, "(latest date_modified", hgnc_date, ")"),
  "",
  "CELLS",
  paste0("  dropped datasets : ", paste(DROP, collapse = ", ")),
  paste0("  retained cells   : ", format(ncol(new), big.mark = ",")),
  "",
  "FEATURES",
  paste0("  ", format(nrow(fm), big.mark = ","), " -> ", format(nrow(new), big.mark = ",")),
  "  alias_symbol is never used. prev_symbol only when it names exactly one approved gene.",
  paste0("  merge guard: rows co-detected in >= ", CODETECT_PCT, "% of one sample's cells are not summed"),
  "", "  rule x decision:", show(table(fm$rule, fm$decision)), "",
  paste0("  co-detection rejections: ", nrow(codetect), " (see ", basename(FCOD), ")"),
  "",
  "SPOT CHECK (% cells detecting)", show(round(det_spot, 1)),
  "",
  "LAYERS",
  "  counts only. Normalize after subsetting (LogNormalize, scale.factor 1e4).",
  "",
  "VERIFIED",
  "  total UMI and per-cell UMI identical before/after aggregation; no duplicate names."
), PROV)
msg("wrote ", PROV)
msg("done")
sessionInfo()
