#!/usr/bin/env Rscript
# marker_audit.R
# Reads only the TSVs written by rebuild_02b_batchgenes.R. No Seurat object is
# loaded, so this runs in seconds and is safe on a login node.
#
# Answers three questions:
#   1. WHICH datasets are responsible for the structural absences? If they
#      concentrate in one or two datasets, it is a known reference-version
#      problem. If they are scattered, marker-based annotation is shakier.
#   2. WHAT KIND of genes are they? Clone-ID / lncRNA / small-RNA absences are
#      reference-version noise and no one annotates with them. Named
#      protein-coding absences are the ones that matter.
#   3. Are any of the STEP 4 MARKER PANELS compromised? This is the one that
#      actually changes what you do next.

DIR <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"

panels <- list(
  TNK = list(
    "pan T/NK"      = c("PTPRC","CD3D","CD3E","CD3G","CD2","TRAC","TRBC1","TRBC2"),
    "CD4 / CD8"     = c("CD4","CD8A","CD8B"),
    "naive / memory"= c("CCR7","SELL","TCF7","LEF1","IL7R","S100A4"),
    "Treg"          = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF4","TNFRSF18"),
    "exhaustion"    = c("PDCD1","LAG3","HAVCR2","TIGIT","TOX","ENTPD1","CXCL13"),
    "cytotoxic"     = c("GZMA","GZMB","GZMK","GZMH","PRF1","NKG7","GNLY","KLRG1"),
    "NK"            = c("NCAM1","FCGR3A","KLRD1","KLRF1","NCR1","TYROBP","FCER1G"),
    "gamma-delta"   = c("TRDC","TRGC1","TRGC2","TRDV1","TRDV2","TRGV9"),
    "proliferating" = c("MKI67","TOP2A","STMN1")
  ),
  stromal = list(
    "pan fibroblast"= c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRB","VIM"),
    "myCAF"         = c("ACTA2","TAGLN","MYL9","POSTN","COL10A1","COL11A1","THBS2","CTHRC1"),
    "iCAF"          = c("IL6","CXCL12","CXCL14","PDGFRA","CFD","C3","C7","HAS1","LMNA","DPT"),
    "apCAF"         = c("CD74","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","SLPI","SAA3"),
    "stellate"      = c("RGS5","PDGFRB","DES","LRAT","PDGFRA","NGFR","CYGB"),
    "pericyte/vSMC" = c("RGS5","ACTA2","NOTCH3","MCAM","MYH11","CNN1")
  ),
  Bcell = list(
    "pan B"         = c("MS4A1","CD79A","CD79B","CD19","BANK1","CD37"),
    "naive"         = c("TCL1A","IGHD","FCER2","IL4R"),
    "memory"        = c("CD27","AIM2","TNFRSF13B","CD80"),
    "germinal ctr"  = c("BCL6","AICDA","RGS13","MEF2B","LMO2"),
    "plasma"        = c("MZB1","XBP1","SDC1","TNFRSF17","PRDM1","DERL3","JCHAIN"),
    "isotype"       = c("IGHM","IGHG1","IGHA1","IGKC","IGLC1")
  )
)

noncoding_rx <- paste(
  "^(AC|AL|AP|BX|CR|CU|FP|Z)[0-9]{5,6}\\.[0-9]+$",
  "^LINC[0-9]+$", "-AS[0-9]+$", "-DT$", "-IT[0-9]$", "-OT[0-9]$",
  "^MIR[0-9]", "^SNOR", "^RNU[0-9]", "^SCARNA", "^RNA5S", "^RN7S",
  "^LOC[0-9]+$", "orf[0-9]+$",
  sep = "|")

out <- character()
say <- function(...) { s <- paste0(...); cat(s, "\n"); out <<- c(out, s) }

for (nm in names(panels)) {
  tf <- file.path(DIR, paste0("batchgene_detection_", nm, ".tsv"))
  if (!file.exists(tf)) { say("MISSING: ", tf); next }
  tab <- read.delim(tf, check.names = FALSE, stringsAsFactors = FALSE)

  meta <- c("gene", "max_pct", "min_pct", "class")
  ds   <- setdiff(colnames(tab), meta)
  m    <- as.matrix(tab[, ds, drop = FALSE]); rownames(m) <- tab$gene

  say("")
  say("=========================================================")
  say("  ", nm, "   (", length(ds), " datasets screened: ",
      paste(ds, collapse = ", "), ")")
  say("=========================================================")

  # ---- 1. who is responsible for the zeros -------------------------------
  st <- tab$class == "STRUCTURAL"
  say("")
  say("-- structural absences: ", sum(st), " genes --")
  if (sum(st)) {
    who <- apply(m[st, , drop = FALSE], 1,
                 function(x) paste(ds[x == 0], collapse = " + "))
    tw <- sort(table(who), decreasing = TRUE)
    say("zero in which dataset(s):")
    for (i in seq_len(min(12, length(tw))))
      say(sprintf("   %6d  %s", tw[i], names(tw)[i]))
    say("READ: if one dataset dominates this list it is a reference-version")
    say("      problem in that dataset. If it is spread evenly, every marker")
    say("      needs checking individually.")

    # ---- 2. what kind of genes ------------------------------------------
    g  <- rownames(m)[st]
    nc <- grepl(noncoding_rx, g)
    say("")
    say(sprintf("   noncoding / clone-ID / small-RNA : %4d  (%.1f%%)",
                sum(nc), 100 * mean(nc)))
    say(sprintf("   named, likely protein-coding     : %4d  (%.1f%%)",
                sum(!nc), 100 * mean(!nc)))
    say("   -- top 30 named protein-coding absences by max detection --")
    coding <- g[!nc]
    if (length(coding)) {
      ord <- coding[order(-tab$max_pct[match(coding, tab$gene)])]
      print(round(m[head(ord, 30), , drop = FALSE], 1))
      out <<- c(out, capture.output(print(round(m[head(ord, 30), , drop = FALSE], 1))))
    }
  }

  # ---- 3. the panels that matter ----------------------------------------
  say("")
  say("-- STEP 4 MARKER PANEL AUDIT --")
  for (pn in names(panels[[nm]])) {
    mk <- panels[[nm]][[pn]]
    st_hit <- mk[mk %in% tab$gene[tab$class == "STRUCTURAL"]]
    gr_hit <- mk[mk %in% tab$gene[tab$class == "GRADED"]]
    absent <- mk[!mk %in% tab$gene]
    ok     <- setdiff(mk, c(st_hit, gr_hit, absent))

    flag <- if (length(st_hit) || length(absent)) "  <<< COMPROMISED" else ""
    say("")
    say("  [", pn, "]", flag)
    if (length(ok))     say("     usable      : ", paste(ok, collapse = ", "))
    if (length(gr_hit)) say("     GRADED      : ", paste(gr_hit, collapse = ", "),
                            "   (check tsv - may be ambient, may be real)")
    if (length(st_hit)) say("     STRUCTURAL  : ", paste(st_hit, collapse = ", "),
                            "   (UNUSABLE - absent from >=1 dataset)")
    if (length(absent)) say("     not in matrix: ", paste(absent, collapse = ", "))
    if (length(st_hit)) {
      print(round(m[st_hit, , drop = FALSE], 1))
      out <<- c(out, capture.output(print(round(m[st_hit, , drop = FALSE], 1))))
    }
  }

  say("")
  say("  A panel is only COMPROMISED if it has no usable members left.")
  say("  Losing one gene from a six-gene panel is survivable; losing the")
  say("  gene the whole call rests on is not.")
}

writeLines(out, file.path(DIR, "marker_audit.txt"))
cat("\nwrote ", file.path(DIR, "marker_audit.txt"), "\n", sep = "")
