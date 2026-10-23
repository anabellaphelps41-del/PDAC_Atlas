#!/usr/bin/env Rscript
# diagnose_harmonization_v3.R  -- final pass, supersedes v1 and v2
#
# WHY A THIRD PASS
#   v1 required the current HGNC symbol to be present in the matrix, so it
#      missed genes renamed twice (USMG5 -> ATP5MD -> ATP5MK).
#   v2 fixed that but used prev_symbol only. HGNC files some renames under
#      alias_symbol instead:  SELENOM prev="" alias="SELM|SEPM"
#                             STK26   prev="" alias="MST4|MASK"
#   Neither can catch SF3B14/SF3B6: HGNC records the alias as "SF3B14a" while
#      the matrix holds "SF3B14". No exact match exists at any setting.
#
# APPROACH
#   prev_symbol links      -> accepted ungated (reliable across v1/v2)
#   alias_symbol links     -> candidates only, must pass a COMPLEMENTARITY GATE
#   normalized alias links -> same gate ("SF3B14a" -> "SF3B14": strip one
#                             trailing lowercase letter)
#
# THE GATE
#   For a candidate pair, using per-dataset detection rates from the
#   batchgene_detection_*.tsv tables:
#       overlap = sum(pmin(a,b)) / sum(pmax(a,b))
#   A gene split across two rows is present in disjoint dataset sets, so
#   overlap ~ 0.  SELM/SELENOM = 0.028.  Two unrelated genes both genuinely
#   expressed score near 1.  We also require each member to be exclusively
#   present (other member exactly 0) in at least one dataset - that is the
#   signature of a split and cannot be produced by a real expression difference
#   without one gene being wholly absent from a cohort.
#
#   This is what rejects the alias noise: AMN/ABCD1, AR/AKR1B1, MIF/AMH are all
#   expressed in the same datasets, so their overlap is high and they fail.
#
# Reads feature names, HGNC, and the detection TSVs. Login node, under a minute.

DIR  <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"
HGNC <- "/ix/tbruno/anp795/PDAC_Atlas/ref/hgnc_complete_set.txt"
SRC  <- file.path(DIR, "subset_Bcell.rds")

OVERLAP_MAX <- 0.15   # accept alias candidates below this
MIN_MAX_DET <- 5      # each member must reach 5% detection in some dataset

msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep=""); flush.console() }
suppressPackageStartupMessages(library(Seurat))

msg("reading feature names")
s <- readRDS(SRC); g <- rownames(s); rm(s); gc(verbose = FALSE)
msg("features: ", format(length(g), big.mark=","))

# --------------------------------------------------------------------------
# Detection tables -> one matrix per subset, for the gate
# --------------------------------------------------------------------------
det <- list()
for (nm in c("TNK","stromal","Bcell")) {
  tf <- file.path(DIR, paste0("batchgene_detection_", nm, ".tsv"))
  if (!file.exists(tf)) next
  tb <- read.delim(tf, check.names=FALSE, stringsAsFactors=FALSE)
  ds <- setdiff(colnames(tb), c("gene","max_pct","min_pct","class"))
  m  <- as.matrix(tb[, ds, drop=FALSE]); rownames(m) <- tb$gene
  det[[nm]] <- m
  msg("detection table ", nm, ": ", nrow(m), " genes x ", ncol(m), " datasets")
}
if (!length(det)) stop("no detection tables found - run rebuild_02b_batchgenes.R first")

# Returns list(pass=, overlap=, excl=, subset=) using the most informative
# subset in which both genes clear MIN_MAX_DET.
gate <- function(a, b) {
  best <- list(pass=FALSE, overlap=NA_real_, excl=FALSE, subset=NA_character_)
  for (nm in names(det)) {
    m <- det[[nm]]
    if (!(a %in% rownames(m) && b %in% rownames(m))) next
    x <- m[a, ]; y <- m[b, ]
    if (max(x) < MIN_MAX_DET || max(y) < MIN_MAX_DET) next
    ov   <- sum(pmin(x, y)) / max(sum(pmax(x, y)), 1e-9)
    excl <- any(x > 0 & y == 0) && any(y > 0 & x == 0)
    if (is.na(best$overlap) || ov < best$overlap)
      best <- list(pass = (ov <= OVERLAP_MAX && excl), overlap = ov,
                   excl = excl, subset = nm)
  }
  best
}

# --------------------------------------------------------------------------
# HGNC candidate maps
# --------------------------------------------------------------------------
msg("loading HGNC")
h <- read.delim(HGNC, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
h <- h[h$status == "Approved", ]
approved <- unique(h$symbol)

explode <- function(col, normalize = FALSE) {
  v <- gsub('^"|"$', '', h[[col]])
  l <- strsplit(v, "\\|", perl=TRUE)
  old <- trimws(unlist(l))
  if (normalize) old <- sub("[a-z]$", "", old)   # SF3B14a -> SF3B14
  d <- data.frame(old = old, cur = rep(h$symbol, lengths(l)),
                  stringsAsFactors = FALSE)
  d[nzchar(d$old) & d$old != d$cur, ]
}

clean <- function(d) {
  d <- unique(d)
  amb <- unique(d$old[duplicated(d$old)])       # maps to >1 current gene
  d <- d[!d$old %in% amb, ]
  d[!d$old %in% approved, ]                     # a live gene is not a synonym
}

pm <- clean(explode("prev_symbol"))
al <- clean(rbind(explode("alias_symbol"), explode("alias_symbol", normalize=TRUE)))
al <- al[!al$old %in% pm$old, ]                 # prev wins where both exist

msg("prev_symbol resolutions : ", nrow(pm))
msg("alias candidates        : ", nrow(al))

# --------------------------------------------------------------------------
# Apply the gate to alias candidates that would actually merge something
# --------------------------------------------------------------------------
cand <- al[al$old %in% g, ]
cand$target_present <- cand$cur %in% g
test <- cand[cand$target_present, ]
msg("alias candidates where BOTH members are rows (gate applies): ", nrow(test))

if (nrow(test)) {
  res <- lapply(seq_len(nrow(test)), function(i) gate(test$old[i], test$cur[i]))
  test$overlap  <- vapply(res, function(r) r$overlap, numeric(1))
  test$mutual   <- vapply(res, function(r) r$excl,    logical(1))
  test$subset   <- vapply(res, function(r) r$subset,  character(1))
  test$pass     <- vapply(res, function(r) r$pass,    logical(1))
  test$pass[is.na(test$pass)] <- FALSE

  msg("")
  msg("ACCEPTED alias pairs (overlap <= ", OVERLAP_MAX, ", mutually exclusive):")
  acc <- test[test$pass, c("old","cur","overlap","subset")]
  if (nrow(acc)) print(acc[order(acc$overlap), ], row.names=FALSE) else msg("   none")

  msg("")
  msg("REJECTED alias pairs, top 25 by overlap (these are the noise the gate exists for):")
  rej <- test[!test$pass & !is.na(test$overlap), c("old","cur","overlap","mutual")]
  if (nrow(rej)) print(head(rej[order(-rej$overlap), ], 25), row.names=FALSE)

  msg("")
  msg("untestable (one member never reaches ", MIN_MAX_DET, "% anywhere): ",
      sum(is.na(test$overlap)), " - left unmerged, cannot be verified")

  al_ok <- test[test$pass, c("old","cur")]
} else {
  al_ok <- al[0, c("old","cur")]
}

# pure renames (target not a row) merge nothing on their own; safe to include
al_ren <- cand[!cand$target_present, c("old","cur")]

resolver <- setNames(c(pm$cur, al_ok$cur, al_ren$cur),
                     c(pm$old, al_ok$old, al_ren$old))
msg("")
msg("total usable resolutions: ", length(resolver))

# --------------------------------------------------------------------------
# Resolve, group, write
# --------------------------------------------------------------------------
base <- sub("\\.[0-9]+$", "", g)
strip_ok <- grepl("\\.[0-9]+$", g) & (base %in% approved | base %in% g) & base != g
step1 <- ifelse(strip_ok, base, g)

need  <- !(step1 %in% approved) & (step1 %in% names(resolver))
canon <- step1; canon[need] <- unname(resolver[step1[need]])

grp   <- split(g, canon)
multi <- grp[lengths(grp) > 1]

msg("")
msg("SPLIT GENES: ", length(multi), " groups, ", sum(lengths(multi)), " features")
msg("rows merged away: ", sum(lengths(multi)) - length(multi))

for (chk in list(c("SELM","SELENOM"), c("SF3B14","SF3B6"), c("MST4","STK26"),
                 c("USMG5","ATP5MD"), c("ATP5B","ATP5F1B"))) {
  pr <- chk[chk %in% g]
  if (!length(pr)) next
  cat(sprintf("   %-20s -> %s\n", paste(pr, collapse=" + "),
              paste(unique(canon[match(pr, g)]), collapse=" / ")))
}

collapse <- do.call(rbind, lapply(names(multi), function(k)
  data.frame(feature = multi[[k]], collapse_into = k, stringsAsFactors = FALSE)))
collapse <- collapse[collapse$feature != collapse$collapse_into, ]
collapse$reason <- "resolved to current HGNC symbol"

write.table(collapse, file.path(DIR, "harmonization_fix_collapse_v3.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
msg("")
msg("wrote harmonization_fix_collapse_v3.tsv  (", nrow(collapse), " features to merge)")

expl <- unique(c(collapse$feature, collapse$collapse_into))
msg("")
msg("-- coverage --")
for (nm in names(det)) {
  tf <- file.path(DIR, paste0("batchgene_detection_", nm, ".tsv"))
  tb <- read.delim(tf, check.names=FALSE, stringsAsFactors=FALSE)
  st <- tb$gene[tb$class == "STRUCTURAL"]
  msg(sprintf("%-8s structural %4d, explained %4d (%2.0f%%), residual ~%d",
              nm, length(st), sum(st %in% expl), 100*mean(st %in% expl),
              sum(!st %in% expl)))
}
msg("done")
