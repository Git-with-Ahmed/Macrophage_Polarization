# Map all gene IDs to mouse SYMBOL, then recompute common genes

library(tidyverse)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Detect ID type and map to SYMBOL
to_symbol <- function(ids) {
  ids <- as.character(ids)
  ens <- grepl("^ENSMUSG", ids, ignore.case = TRUE)
  ent <- grepl("^[0-9]+$", ids)
  sym <- !(ens | ent)
  
  out <- rep(NA_character_, length(ids))
  
  if (any(ens)) {
    keys <- sub("\\.\\d+$", "", ids[ens])  # strip Ensembl version if present
    out[ens] <- mapIds(org.Mm.eg.db, keys=keys, keytype="ENSEMBL",
                       column="SYMBOL", multiVals="first")
  }
  if (any(ent)) {
    out[ent] <- mapIds(org.Mm.eg.db, keys=ids[ent], keytype="ENTREZID",
                       column="SYMBOL", multiVals="first")
  }
  if (any(sym)) {
    out[sym] <- ids[sym]  # already symbols
  }
  
  # Clean up oddities
  out <- toupper(out)           # normalize case so joins match
  out <- na_if(out, "")         # drop empties
  out
}

# Harmonize one dataset: map to SYMBOL and de-duplicate rows
harmonize_one <- function(d) {
  counts <- d$counts
  counts$gene <- to_symbol(counts$gene)
  counts <- counts %>%
    filter(!is.na(gene)) %>%
    distinct(gene, .keep_all = TRUE) %>%   # keep first if duplicates
    arrange(gene)
  list(counts = counts, meta = d$meta)
}

# 1) Map all to SYMBOL
data_list_sym <- lapply(data_list, harmonize_one)

# 2) Recompute common genes
gene_sets <- lapply(data_list_sym, \(d) d$counts$gene)
common_genes <- Reduce(intersect, gene_sets)

length(common_genes)  # should now be > 0

# 3) Subset each dataset to common genes, ordered alphabetically identically
data_list_common <- lapply(
  data_list_sym,
  \(d) {
    counts_sub <- d$counts %>%
      filter(gene %in% common_genes) %>%
      arrange(gene)
    list(counts = counts_sub, meta = d$meta)
  }
)

# 4) Sanity: identical gene order across all elements
stopifnot(all(vapply(data_list_common,
                     \(d) identical(data_list_common[[1]]$counts$gene, d$counts$gene),
                     logical(1))))

rm(data_list,data_list_sym,gene_sets)

# source
src <- data_list_common[["GSE304646_counts"]]
ct  <- src$counts
mt  <- src$meta

# convenience
ctrl <- c("S1","S6","S11","S16","S21")
m1_3to5   <- c("S12","S17","S22")
m2a_3to5  <- c("S13","S18","S23")
m2c_3to5  <- c("S14","S19","S24")
tgfb <- c("S5","S10","S15","S20","S25")

# -------- A) update ORIGINAL: keep all controls + M1/M2a/M2c from blocks 3–5, drop M2c_TGFb --------
keep_meta_rows <- mt$new_name %in% c(ctrl, m1_3to5, m2a_3to5, m2c_3to5)
mt_keep <- mt[keep_meta_rows, , drop = FALSE]

# order to match your example
order_orig <- c("S1","S6","S11","S12","S13","S14",
                "S16","S17","S18","S19",
                "S21","S22","S23","S24")
mt_keep <- mt_keep[match(order_orig, mt_keep$new_name), ]
stopifnot(!any(is.na(mt_keep$new_name)))

ct_keep <- ct[, c("gene", mt_keep$new_name), drop = FALSE]

data_list_common[["GSE304646_counts"]] <- list(
  counts = ct_keep,
  meta   = mt_keep
)

# -------- B) create NEW: controls + M2c_TGFb (renamed to M2c) --------
order_new <- c("S1","S5","S6","S10","S11","S15","S16","S20","S21","S25")
mt_new <- mt[match(order_new, mt$new_name), c("dataset","new_name","treatment")]
mt_new$treatment[mt_new$new_name %in% tgfb] <- "M2c"     # relabel
mt_new$dataset <- "GSE304646_1_counts"

ct_new <- ct[, c("gene", order_new), drop = FALSE]

data_list_common[["GSE304646_1_counts"]] <- list(
  counts = ct_new,
  meta   = mt_new
)

# -------- sanity --------
stopifnot(
  identical(colnames(data_list_common[["GSE304646_counts"]]$counts)[-1],
            data_list_common[["GSE304646_counts"]]$meta$new_name),
  identical(colnames(data_list_common[["GSE304646_1_counts"]]$counts)[-1],
            data_list_common[["GSE304646_1_counts"]]$meta$new_name),
  !any(data_list_common[["GSE304646_counts"]]$meta$treatment == "M2c_TGFb"),
  all(data_list_common[["GSE304646_1_counts"]]$meta$treatment %in% c("Control","M2c"))
)
