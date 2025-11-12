suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
})

# ---------- Inputs ----------
signatures_raw <- sig_topN
validation     <- res_list
OUT_DIR <- "Validation_GSEA"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------- Helpers ----------
prep_sig <- function(v) toupper(unique(v[nzchar(v)]))

prep_ranks <- function(de_tbl) {
  de_tbl %>%
    transmute(gene = toupper(as.character(gene)),
              logFC = as.numeric(logFC)) %>%
    filter(!is.na(gene), is.finite(logFC)) %>%
    group_by(gene) %>% summarise(logFC = mean(logFC), .groups="drop") %>%
    arrange(desc(logFC)) %>%
    { setNames(.$logFC, .$gene) }
}

run_fgsea_one <- function(sig_genes, ranks, minSize = 10L) {
  ov <- intersect(names(ranks), sig_genes)
  if (length(ov) < minSize) return(NULL)
  res <- fgseaMultilevel(pathways = list(sig = sig_genes),
                         stats = ranks, scoreType = "std", nproc = 1)
  as_tibble(res)[, c("NES","pval","padj","size")]
}

# ---------- Normalize signatures ----------
signatures <- lapply(signatures_raw, prep_sig)

# ---------- Overlap diagnostics (FIXED INIT) ----------
overlap_report <- tibble(
  dataset  = character(),
  contrast = character(),
  pathway  = character(),
  overlap  = integer()
)

for (ds in names(validation)) {
  for (trt in names(validation[[ds]])) {
    ranks <- prep_ranks(validation[[ds]][[trt]])
    for (pw in names(signatures)) {
      ov <- sum(names(ranks) %in% signatures[[pw]])
      overlap_report <- add_row(
        overlap_report,
        dataset = ds, contrast = trt, pathway = pw, overlap = ov
      )
    }
  }
}

# ---------- Run all GSEA ----------
rows <- list()
for (ds in names(validation)) {
  for (trt in names(validation[[ds]])) {
    ranks <- prep_ranks(validation[[ds]][[trt]])
    for (pw in names(signatures)) {
      res <- run_fgsea_one(signatures[[pw]], ranks, minSize = 10L)
      if (is.null(res)) next
      rows[[length(rows)+1]] <- tibble(
        dataset = ds, contrast = trt, pathway = pw,
        NES = res$NES, pval = res$pval, padj = res$padj, size = res$size
      )
    }
  }
}

gsea_res <- if (length(rows)) bind_rows(rows) else tibble(
  dataset=character(), contrast=character(), pathway=character(),
  NES=double(), pval=double(), padj=double(), size=integer()
)

# ---------- Summaries + save ----------
if (nrow(gsea_res)) {
  gsea_summary <- gsea_res %>%
    group_by(pathway, contrast) %>%
    summarise(mean_NES = mean(NES), min_padj = min(padj), .groups="drop") %>%
    arrange(desc(mean_NES))
  readr::write_csv(gsea_res, file.path(OUT_DIR, "GSEA_all_results.csv"))
  readr::write_csv(gsea_summary, file.path(OUT_DIR, "GSEA_summary.csv"))
  readr::write_csv(overlap_report, file.path(OUT_DIR, "GSEA_overlap_report.csv"))
  print(gsea_summary)
} else {
  message("No valid GSEA rows. See overlap report for counts per pair.")
  readr::write_csv(overlap_report, file.path(OUT_DIR, "GSEA_overlap_report.csv"))
}
