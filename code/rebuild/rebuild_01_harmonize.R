# rebuild_01_harmonize.R
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
IN   <- "/ix/tbruno/anp795/PDAC_Atlas/combined_9_datasets_hierarchical_annotation_seurat.rds"
OUT  <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"
HGNC <- "/ix/tbruno/anp795/PDAC_Atlas/ref/hgnc_complete_set.txt"
DROP <- c("GSE202051_main", "GSE202051_010nuc", "GSE202051_010orgCRT")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(IN)
cat("loaded:", nrow(obj), "genes x", ncol(obj), "cells\n")
## --- 1. archive the dropped cells, then remove them ---
arch <- colnames(obj)[obj$dataset_id %in% DROP]
cat("archiving/dropping:", length(arch), "cells\n")
saveRDS(subset(obj, cells = arch),
        file.path(OUT, "ARCHIVE_GSE202051_nuclei_organoid.rds"), compress = FALSE)
obj <- subset(obj, cells = setdiff(colnames(obj), arch)); gc()
cat("retained:", ncol(obj), "cells\n")
## --- 2. build old-symbol -> current-symbol map ---
g <- rownames(obj)
cat("\nfeatures:", length(g), "| with --ENSG suffix:", sum(grepl("--ENSG", g)), "\n")
h <- read.delim(HGNC, sep = "\t", quote = "", stringsAsFactors = FALSE)
h <- h[h$status == "Approved", ]
old_list <- strsplit(paste(h$alias_symbol, h$prev_symbol, sep = "|"), "\\|")
alias <- data.frame(old = trimws(unlist(old_list)),
                    new = rep(h$symbol, lengths(old_list)),
                    stringsAsFactors = FALSE)
alias <- alias[nzchar(alias$old), ]
alias <- alias[!duplicated(alias$old), ]
alias <- alias[!(alias$old %in% h$symbol), ]   # never remap an already-current symbol
map <- setNames(alias$new, alias$old)
newname <- sub("--ENSG.*$", "", g)             # strip Ensembl disambiguation
hit <- newname %in% names(map)
cat("remapped by HGNC alias:", sum(hit), "\n")
newname[hit] <- map[newname[hit]]
cat("features collapsing into one another:", length(g) - length(unique(newname)), "\n")
dupes <- sort(table(newname[newname %in% newname[duplicated(newname)]]), decreasing = TRUE)
cat("\n--- top collapsed groups ---\n"); print(head(dupes, 40))
## --- 3. collapse counts by summing ---
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
grp <- factor(newname)
A <- sparseMatrix(i = as.integer(grp), j = seq_along(grp), x = 1,
                  dims = c(nlevels(grp), length(grp)),
                  dimnames = list(levels(grp), NULL))
collapsed <- as(A %*% counts, "dgCMatrix")
colnames(collapsed) <- colnames(counts)
rm(counts, A); gc()
cat("\ncollapsed:", nrow(collapsed), "genes x", ncol(collapsed), "cells\n")
## --- 4. rebuild clean object (reductions/graphs intentionally discarded) ---
keep <- c("orig.ident","dataset_id","dataset_accession","sample_id","patient_id",
          "condition","tissue","batch","nCount_RNA","nFeature_RNA",
          "percent.mt","percent.ribo","doublet_score","doublet_call",
          "master_stage","master_treatment_status","master_primary_or_metastatic",
          "master_sequencing_modality",
          "celltype_level1","celltype_level2","celltype_level3")
keep <- intersect(keep, colnames(obj@meta.data))
md <- obj@meta.data[, keep]
rm(obj); gc()
new <- CreateSeuratObject(counts = collapsed, meta.data = md)
rm(collapsed); gc()
new <- NormalizeData(new, normalization.method = "LogNormalize", scale.factor = 1e4)
saveRDS(new, file.path(OUT, "atlas_harmonized_counts.rds"), compress = FALSE)
cat("\ndone:", nrow(new), "genes x", ncol(new), "cells\n")
print(table(new$celltype_level2))
sessionInfo()