# --- 1) Canonical label mapper (exact target case) ---
canon <- function(x) {
  y <- toupper(x)
  y <- ifelse(y == "M2C", "M2c", y)
  y <- ifelse(y == "M2A", "M2a", y)
  y <- ifelse(y == "M2B", "M2b", y)
  y <- ifelse(y == "M1",  "M1",  y)
  y
}

# --- 2) Uniquify per-dataset names but keep duplicates (M2c, M2c.1, etc.) ---
res_list_u <- lapply(res_list, function(per_ds) {
  names(per_ds) <- make.unique(names(per_ds), sep=".")
  per_ds
})

# Quick sanity: what classes exist after canon?
avail <- unlist(lapply(res_list_u, function(per_ds) canon(names(per_ds))))
table(avail)

suppressPackageStartupMessages(library(tidyverse))
FDR_MAX <- 0.05; MIN_LOGFC <- 0; TOP_N <- 200
classes <- c("M1","M2a","M2b","M2c")

norm_tt <- function(tt){
  df <- as.data.frame(tt)
  if ("genes" %in% names(df) && !"gene"%in% names(df)) df$gene <- df$genes
  if ("P.Value" %in% names(df))  names(df)[names(df)=="P.Value"]  <- "PValue"
  if ("adj.P.Val" %in% names(df)) names(df)[names(df)=="adj.P.Val"] <- "FDR"
  stopifnot(all(c("gene","logFC") %in% names(df)))
  if (!"FDR" %in% names(df)) df$FDR <- NA_real_
  df %>% mutate(gene = toupper(as.character(gene))) %>%
    filter(!is.na(gene), is.finite(logFC))
}

up_genes_one <- function(tt, fdr_max=FDR_MAX, min_logfc=MIN_LOGFC){
  norm_tt(tt) %>%
    filter(logFC > min_logfc, if(!all(is.na(FDR))) FDR <= fdr_max else TRUE) %>%
    distinct(gene) %>% pull(gene)
}

# collect all up-genes for a class across ALL entries whose canon(label)==class
collect_up_genes <- function(obj, cl){
  unlist(lapply(obj, function(per_ds){
    out <- character(0)
    for (i in seq_along(per_ds)) if (identical(canon(names(per_ds)[i]), cl))
      out <- c(out, up_genes_one(per_ds[[i]]))
    out
  }), use.names = FALSE)
}

up_by_class <- setNames(lapply(classes, \(cl) collect_up_genes(res_list_u, cl)), classes)

# build UNIQUE vs all-other classes
sig_union <- list(); sig_topN <- list(); sig_stats <- tibble()
for (cl in classes) {
  focal <- up_by_class[[cl]]
  if (!length(focal)) next
  other <- unique(unlist(up_by_class[setdiff(classes, cl)], use.names = FALSE))
  focal_unique <- focal[!(focal %in% other)]
  if (!length(focal_unique)) next
  
  votes <- sort(table(focal_unique), decreasing = TRUE)
  ranked <- names(votes)
  
  sig_union[[cl]] <- ranked
  sig_topN[[cl]]  <- head(ranked, TOP_N)
  
  sig_stats <- bind_rows(sig_stats, tibble(
    class = cl,
    unique_genes_all = length(ranked),
    topN = min(TOP_N, length(ranked)),
    votes_median = median(as.integer(votes)),
    votes_max    = max(as.integer(votes))
  ))
}

sig_stats