#!/usr/bin/env Rscript
# diagnose_harmonization_v2.R   -- supersedes diagnose_harmonization.R
#
# WHY v1 WAS WRONG
# v1 asked: "is there an HGNC synonym pair where BOTH the old symbol and the
# CURRENT symbol are present as separate rows?" That misses every gene renamed
# more than once, because the matrix can hold two DEAD symbols and not the live
# one. USMG5 -> ATP5MD -> ATP5MK is the case that exposed it: the matrix has
# USMG5 and ATP5MD, has no ATP5MK, and v1 saw nothing.
#
# CORRECT APPROACH
# Resolve every feature independently to its current HGNC symbol, then group by
# the resolved symbol. Any group with >1 member is a split gene. The target is
# the group key, whether or not it already exists as a row. HGNC's prev_symbol
# field is cumulative (ATP5MK lists USMG5 AND ATP5MD), so a single hop is enough
# once you stop requiring the target to be present.
#
# Guards kept from v1, both load-bearing:
#   - never resolve a symbol that is ITSELF currently approved (a live gene is
#     not a synonym for another live gene)
#   - prev_symbol only. alias_symbol is noise: it lists AMN under ABCD1, AR
#     under AKR1B1, MIF under AMH. Collapsing on aliases merges distinct genes.
#   - a symbol that is prev_symbol for TWO different current genes is ambiguous
#     and is never collapsed. STRA13 is the live example: historically used for
#     both BHLHE40 and CENPX.
#
# Reads feature names + HGNC only. Login node, well under a minute.

DIR  <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"
HGNC <- "/ix/tbruno/anp795/PDAC_Atlas/ref/hgnc_complete_set.txt"
SRC  <- file.path(DIR, "subset_Bcell.rds")

msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep=""); flush.console() }
suppressPackageStartupMessages(library(Seurat))

msg("reading feature names from ", SRC)
s <- readRDS(SRC); g <- rownames(s); rm(s); gc(verbose = FALSE)
msg("features: ", format(length(g), big.mark=","))

msg("loading HGNC")
h <- read.delim(HGNC, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
h <- h[h$status == "Approved", ]
approved <- unique(h$symbol)
msg("approved symbols: ", format(length(approved), big.mark=","))

# -- prev_symbol -> current, pipe-split, ambiguity removed ------------------
v <- gsub('^"|"$', '', h$prev_symbol)
l <- strsplit(v, "\\|", perl = TRUE)
pm <- data.frame(old = trimws(unlist(l)),
                 cur = rep(h$symbol, lengths(l)),
                 stringsAsFactors = FALSE)
pm <- pm[nzchar(pm$old) & pm$old != pm$cur, ]
pm <- unique(pm)

amb <- unique(pm$old[duplicated(pm$old)])
msg("prev_symbol mappings: ", format(nrow(pm), big.mark=","),
    " | ambiguous (map to >1 current gene, never collapsed): ", length(amb))
if (length(amb)) {
  show <- pm[pm$old %in% intersect(amb, g), ]
  if (nrow(show)) {
    msg("ambiguous symbols PRESENT in your matrix - these stay as-is:")
    print(show[order(show$old), ], row.names = FALSE)
  }
}
pm <- pm[!pm$old %in% amb, ]
# a live gene is never a synonym for another live gene
pm <- pm[!pm$old %in% approved, ]
resolver <- setNames(pm$cur, pm$old)
msg("usable prev_symbol resolutions: ", format(length(resolver), big.mark=","))

# -- resolve every feature ---------------------------------------------------
# 1. strip make.unique ".N", but ONLY when the base is a real symbol or is
#    itself present. This spares clone IDs like AL627171.2, whose base is not.
base <- sub("\\.[0-9]+$", "", g)
strip_ok <- grepl("\\.[0-9]+$", g) & (base %in% approved | base %in% g) & base != g
step1 <- ifelse(strip_ok, base, g)
msg("features with a strippable make.unique '.N' suffix: ", sum(strip_ok))

# 2. HGNC-resolve anything not already approved
need <- !(step1 %in% approved) & (step1 %in% names(resolver))
canon <- step1
canon[need] <- unname(resolver[step1[need]])
msg("features resolved to a newer HGNC symbol: ", sum(need))

# -- group and find splits ---------------------------------------------------
grp <- split(g, canon)
multi <- grp[lengths(grp) > 1]
msg("")
msg("SPLIT GENES (>1 feature resolving to the same current symbol): ", length(multi))
msg("features involved: ", sum(lengths(multi)),
    "  -> will merge away ", sum(lengths(multi)) - length(multi), " rows")

if (length(multi)) {
  tgt_present <- names(multi) %in% g
  msg("  groups whose target already exists as a row : ", sum(tgt_present))
  msg("  groups that are a pure rename (target absent): ", sum(!tgt_present))
  msg("")
  msg("first 50 groups:")
  o <- head(order(-lengths(multi), names(multi)), 50)
  for (i in o) cat(sprintf("   %-16s <- %s\n", names(multi)[i], paste(multi[[i]], collapse=", ")))
}

# sanity: did we catch the pairs that exposed the v1 bug?
msg("")
msg("-- confirmation on the pairs that exposed the v1 bug --")
for (pair in list(c("USMG5","ATP5MD"), c("SELM","SELENOM"), c("SF3B14","SF3B6"),
                  c("ATP5MPL"), c("HIST1H3D"), c("SEPTIN2"), c("RPSA2"), c("GPX1"))) {
  pr <- pair[pair %in% g]
  if (!length(pr)) next
  cat(sprintf("   %-22s -> resolves to %s\n",
              paste(pr, collapse=" + "),
              paste(unique(canon[match(pr, g)]), collapse=" / ")))
}

# -- write the collapse table ------------------------------------------------
collapse <- do.call(rbind, lapply(names(multi), function(k) {
  data.frame(feature = multi[[k]], collapse_into = k,
             reason = ifelse(multi[[k]] == k, "target (kept)", "resolved to current HGNC symbol"),
             stringsAsFactors = FALSE)
}))
collapse <- collapse[collapse$feature != collapse$collapse_into, ]

write.table(collapse, file.path(DIR, "harmonization_fix_collapse_v2.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
msg("")
msg("wrote ", file.path(DIR, "harmonization_fix_collapse_v2.tsv"),
    "  (", nrow(collapse), " features to merge)")

# -- how much of the structural lists does this now explain? -----------------
expl <- unique(c(collapse$feature, collapse$collapse_into))
msg("")
msg("-- coverage of the structural-absence lists --")
for (nm in c("TNK","stromal","Bcell")) {
  tf <- file.path(DIR, paste0("batchgene_detection_", nm, ".tsv"))
  if (!file.exists(tf)) next
  tab <- read.delim(tf, check.names=FALSE, stringsAsFactors=FALSE)
  st <- tab$gene[tab$class == "STRUCTURAL"]
  msg(sprintf("%-8s structural %4d, explained %4d (%2.0f%%), residual ~%d",
              nm, length(st), sum(st %in% expl), 100*mean(st %in% expl),
              sum(!st %in% expl)))
}

# what survives, and does it matter?
msg("")
msg("-- residual named (non clone-ID) genes still structurally absent, TNK --")
tf <- file.path(DIR, "batchgene_detection_TNK.tsv")
if (file.exists(tf)) {
  tab <- read.delim(tf, check.names=FALSE, stringsAsFactors=FALSE)
  res <- setdiff(tab$gene[tab$class == "STRUCTURAL"], expl)
  nc  <- grepl("\\.[0-9]+$|^LINC|-AS[0-9]+$|-DT$|-[IO]T[0-9]$|^MIR|^SNOR|^RN7|^RNU|^SNHG|orf[0-9]|P[0-9]+$",
               res)
  msg("residual: ", length(res), "  (noncoding-looking ", sum(nc),
      ", named ", sum(!nc), ")")
  print(head(res[!nc], 80))
}
msg("done")
