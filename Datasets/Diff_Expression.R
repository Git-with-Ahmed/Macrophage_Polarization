# ===== 0) Setup =====
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
})

# Input: data_list_common from your str() output
obj <- data_list_common

# Output dirs
OUT_DIR_DE   <- "DE_full"
OUT_DIR_PREV <- "DE_top"
dir.create(OUT_DIR_DE,   showWarnings = FALSE)
dir.create(OUT_DIR_PREV, showWarnings = FALSE)

# ===== 1) Small helpers =====
is_count_like <- function(df_counts) {
  M <- as.matrix(df_counts)
  is.numeric(M) &&
    all(M >= 0, na.rm = TRUE) &&
    all(abs(M - round(M)) < 1e-6, na.rm = TRUE)
}

maybe_log2 <- function(M) {
  x <- as.numeric(M); x <- x[is.finite(x)]
  if (!length(x)) return(M)
  if (max(x, na.rm = TRUE) > 100 || median(x, na.rm = TRUE) > 20) log2(M + 0.5) else M
}

normalize_tt <- function(tt) {
  df <- as.data.frame(tt)
  if (!"gene" %in% names(df)) {
    if (!is.null(rownames(df))) {
      df <- tibble::rownames_to_column(df, "gene")
    } else if ("genes" %in% names(df)) {
      df <- df %>% dplyr::rename(gene = genes)
    }
  }
  if ("P.Value" %in% names(df))  df <- df %>% dplyr::rename(PValue = P.Value)
  if ("adj.P.Val" %in% names(df)) df <- df %>% dplyr::rename(FDR = adj.P.Val)
  rownames(df) <- NULL
  df
}

# ===== 2) edgeR path (counts) =====
run_edger <- function(counts_df, meta_df, trt_name) {
  ctl <- meta_df$new_name[meta_df$treatment == "Control"]
  trt <- meta_df$new_name[meta_df$treatment == trt_name]
  stopifnot(length(ctl) > 0, length(trt) > 0)
  
  y <- DGEList(
    counts = as.matrix(counts_df[, c(ctl, trt), drop = FALSE]),
    genes  = counts_df$gene
  )
  keep <- filterByExpr(y, group = factor(c(rep("C", length(ctl)), rep("T", length(trt)))))
  if (!any(keep)) return(NULL)
  
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ group, data = data.frame(
    group = factor(c(rep("Control", length(ctl)), rep(trt_name, length(trt))))
  ))
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  tt <- edgeR::topTags(qlf, n = Inf)$table
  tt <- normalize_tt(tt) %>% arrange(FDR)
  tt
}

# ===== 3) limma path (already log-like) =====
run_limma <- function(expr_df, meta_df, trt_name) {
  ctl <- meta_df$new_name[meta_df$treatment == "Control"]
  trt <- meta_df$new_name[meta_df$treatment == trt_name]
  stopifnot(length(ctl) > 0, length(trt) > 0)
  
  E <- as.matrix(expr_df[, c(ctl, trt), drop = FALSE])
  mode(E) <- "numeric"
  E <- maybe_log2(E)
  rownames(E) <- expr_df$gene
  
  design <- model.matrix(~ group, data = data.frame(
    group = factor(c(rep("Control", length(ctl)), rep(trt_name, length(trt))))
  ))
  
  fit <- lmFit(E, design)
  fit <- eBayes(fit, trend = TRUE)
  
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P") %>%
    dplyr::rename(PValue = P.Value, FDR = adj.P.Val) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(FDR)
  tt
}

# ===== 4) Main loop: compute and save =====
res_list  <- list()
manifest  <- list()

for (ds_name in names(obj)) {
  d <- obj[[ds_name]]
  meta <- d$meta %>% arrange(new_name)
  
  # safety
  stopifnot(all(meta$new_name %in% names(d$counts)))
  mat <- d$counts %>% dplyr::select(all_of(meta$new_name))
  rownames(mat) <- d$counts$gene
  
  use_counts <- is_count_like(mat)
  trts <- setdiff(unique(meta$treatment), "Control")
  if (!length(trts)) next
  
  # per-dataset folder
  ds_dir_de   <- file.path(OUT_DIR_DE,   ds_name); dir.create(ds_dir_de,   showWarnings = FALSE)
  ds_dir_prev <- file.path(OUT_DIR_PREV, ds_name); dir.create(ds_dir_prev, showWarnings = FALSE)
  
  for (trt in trts) {
    tt <- if (use_counts) run_edger(d$counts, meta, trt) else run_limma(d$counts, meta, trt)
    if (is.null(tt)) {
      message("Skipped: ", ds_name, " | ", trt)
      next
    }
    
    # save full table
    out_full <- file.path(ds_dir_de, paste0("DE_", ds_name, "__", trt, "_vs_Control.csv"))
    readr::write_csv(tt, out_full)
    
    # save a small preview for quick browsing
    out_prev <- file.path(ds_dir_prev, paste0("Top200_", ds_name, "__", trt, ".csv"))
    readr::write_csv(dplyr::slice_head(tt, n = 200), out_prev)
    
    # keep in memory for downstream work
    if (is.null(res_list[[ds_name]])) res_list[[ds_name]] <- list()
    res_list[[ds_name]][[trt]] <- tt
    
    # manifest row
    manifest[[length(manifest) + 1]] <- tibble(
      dataset      = ds_name,
      treatment    = trt,
      method       = if (use_counts) "edgeR" else "limma",
      n_genes_full = nrow(tt),
      n_sig_FDR5   = sum(!is.na(tt$FDR) & tt$FDR <= 0.05),
      path_full    = out_full,
      path_top200  = out_prev
    )
  }
}

manifest_tbl <- dplyr::bind_rows(manifest)
readr::write_csv(manifest_tbl, file.path(OUT_DIR_DE, "DE_manifest.csv"))

# Quick check
print(manifest_tbl)
str(lapply(res_list, names))

res_list <- lapply(res_list, function(dataset) {
  lapply(dataset, function(df) {
    if ("genes" %in% names(df)) {
      df$gene <- df$genes
      df$genes <- NULL
    }
    df
  })
})


