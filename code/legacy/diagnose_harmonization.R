#!/usr/bin/env Rscript
# diagnose_harmonization.R
#
# Step 1 collapsed gene symbols but the batch-gene screen shows known HGNC
# synonym pairs surviving as SEPARATE rows with complementary per-dataset
# detection (ATP5B/ATP5F1B, SELM/SELENOM, SEPT7/SEPTIN7, ...). This script
# enumerates every such pair still present, plus every make.unique() ".N"
# feature whose base symbol also exists, and writes a collapse table for the
# step-1 rerun to consume.
#
# Reads rownames + the HGNC file only. No count matrix is touched.
# Safe and fast on a login node.

DIR   <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"
HGNC  <- "/ix/tbruno/anp795/PDAC_Atlas/ref/hgnc_complete_set.txt"
# smallest subset - same 49,570-feature space as the full object, loads in seconds
SRC   <- file.path(DIR, "subset_Bcell.rds")

msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep=""); flush.console() }

suppressPackageStartupMessages(library(Seurat))

msg("reading feature names from ", SRC)
s <- readRDS(SRC); g <- rownames(s); rm(s); gc(verbose = FALSE)
msg("features: ", format(length(g), big.mark = ","))

# --------------------------------------------------------------------------
# BUG 1 - make.unique() ".N" suffixes that step 1 never stripped
#
# Only flag when removing the suffix yields a symbol that is ALSO present as
# its own row. This deliberately spares clone-based Ensembl IDs such as
# AL627171.2, whose base (AL627171) is not in the matrix.
# --------------------------------------------------------------------------
has_suffix <- grepl("\\.[0-9]+$", g)
base       <- sub("\\.[0-9]+$", "", g)
dot_hits   <- which(has_suffix & base %in% g & base != g)

dot_tbl <- data.frame(feature = g[dot_hits],
                      collapse_into = base[dot_hits],
                      reason = "make.unique .N suffix",
                      stringsAsFactors = FALSE)

msg("BUG 1  make.unique '.N' features whose base symbol also exists: ", nrow(dot_tbl))
if (nrow(dot_tbl)) print(utils::head(dot_tbl, 30))

# also report .N features whose base is NOT present - these are fine, but
# their base symbol was still never HGNC-checked
orphan <- g[has_suffix & !(base %in% g) & grepl("^[A-Z][A-Z0-9orf-]+\\.[0-9]+$", g) &
            !grepl("^[A-Z]{1,2}[0-9]{6}\\.[0-9]+$", g)]
msg("       .N features whose base is absent (review, likely benign): ", length(orphan))
if (length(orphan)) print(utils::head(orphan, 20))

# --------------------------------------------------------------------------
# BUG 2 - HGNC synonym pairs where BOTH members are still separate rows
# --------------------------------------------------------------------------
msg("loading HGNC")
h <- read.delim(HGNC, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
msg("HGNC rows: ", format(nrow(h), big.mark = ","), " | approved: ", sum(h$status == "Approved"))
h <- h[h$status == "Approved", ]

# prev_symbol and alias_symbol are PIPE-DELIMITED multi-value fields. A mapper
# that matches the whole field instead of splitting on "|" silently misses
# every gene that has more than one previous symbol - a strong candidate for
# what went wrong in step 1.
explode <- function(col) {
  v <- gsub('^"|"$', '', h[[col]])
  l <- strsplit(v, "\\|", perl = TRUE)
  n <- lengths(l)
  data.frame(old = trimws(unlist(l)),
             cur = rep(h$symbol, n),
             src = col,
             stringsAsFactors = FALSE)
}
map <- rbind(explode("prev_symbol"), explode("alias_symbol"))
map <- map[nzchar(map$old) & map$old != map$cur, ]
msg("synonym -> current mappings after pipe-splitting: ", format(nrow(map), big.mark = ","))

# An 'old' that is itself an approved symbol for a DIFFERENT gene must never be
# collapsed - that is a genuine ambiguity, not a rename.
approved <- unique(h$symbol)
map$old_is_approved <- map$old %in% approved

both <- map[map$cur %in% g & map$old %in% g, ]
both <- both[order(both$src, both$cur), ]

safe   <- both[both$src == "prev_symbol" & !both$old_is_approved, ]
review <- both[!(both$src == "prev_symbol" & !both$old_is_approved), ]

msg("BUG 2  synonym pairs with BOTH members present as separate rows: ", nrow(both))
msg("       high confidence (prev_symbol, old not an approved symbol): ", nrow(safe))
msg("       needs review (alias-based, or old is approved elsewhere):  ", nrow(review))
if (nrow(safe)) {
  msg("       high-confidence pairs (first 40):")
  print(utils::head(safe[, c("old", "cur", "src")], 40), row.names = FALSE)
}
if (nrow(review)) {
  msg("       review pairs (first 25) - DO NOT collapse these blindly:")
  print(utils::head(review[, c("old", "cur", "src", "old_is_approved")], 25), row.names = FALSE)
}

# --------------------------------------------------------------------------
# How much of the batch-gene screen does this explain?
# --------------------------------------------------------------------------
explained <- unique(c(dot_tbl$feature, dot_tbl$collapse_into, safe$old, safe$cur))
msg("")
msg("-- coverage of the structural-absence lists --")
for (nm in c("TNK", "stromal", "Bcell")) {
  tf <- file.path(DIR, paste0("batchgene_detection_", nm, ".tsv"))
  if (!file.exists(tf)) next
  tab <- read.delim(tf, check.names = FALSE, stringsAsFactors = FALSE)
  st  <- tab$gene[tab$class == "STRUCTURAL"]
  gr  <- tab$gene[tab$class == "GRADED"]
  msg(sprintf("%-8s structural %4d, of which %4d (%.0f%%) are unmerged synonyms/.N features",
              nm, length(st), sum(st %in% explained),
              100 * mean(st %in% explained)))
  msg(sprintf("%-8s graded     %4d, of which %4d (%.0f%%) are unmerged synonyms/.N features",
              "", length(gr), sum(gr %in% explained),
              100 * mean(gr %in% explained)))
  msg(sprintf("%-8s residual structural after a corrected rerun: ~%d",
              "", sum(!st %in% explained)))
}

# --------------------------------------------------------------------------
# Collapse table for the step-1 rerun
# --------------------------------------------------------------------------
collapse <- rbind(
  dot_tbl,
  data.frame(feature = safe$old, collapse_into = safe$cur,
             reason = "HGNC prev_symbol", stringsAsFactors = FALSE)
)
collapse <- collapse[!duplicated(collapse$feature), ]

write.table(collapse, file.path(DIR, "harmonization_fix_collapse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(review, file.path(DIR, "harmonization_fix_review.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

msg("")
msg("wrote ", file.path(DIR, "harmonization_fix_collapse.tsv"), "  (", nrow(collapse), " features to merge)")
msg("wrote ", file.path(DIR, "harmonization_fix_review.tsv"),   "  (", nrow(review),   " to eyeball first)")
msg("done")
