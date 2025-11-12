# manually fix M2c_TGFb naming inside res_list
for (ds in names(res_list)) {
  names(res_list[[ds]]) <- gsub("^M2c_TGFb$", "M2c", names(res_list[[ds]]))
}


# ===== 1) Fix naming =====
for (ds in names(res_list)) {
  names(res_list[[ds]]) <- gsub("^M2c_TGFb$", "M2c", names(res_list[[ds]]))
}
str(lapply(res_list, names))


# ===== 2) Define validation selection rules =====
val_tbl <- tibble::tribble(
  ~class, ~dataset,          ~treatment,
  "M1",   "GSE123596_counts","M1",
  "M2a",  "GSE123596_counts","M2a",
  "M2b",  "GSE112081_counts","M2b",
  "M2c",  "GSE178999_counts","M2c"
)

print(val_tbl)


# ===== 3) Build validation_list =====
validation_list <- list()
for (i in seq_len(nrow(val_tbl))) {
  ds  <- val_tbl$dataset[i]
  trt <- val_tbl$treatment[i]
  if (!is.null(res_list[[ds]][[trt]])) {
    if (is.null(validation_list[[ds]])) validation_list[[ds]] <- list()
    validation_list[[ds]][[trt]] <- res_list[[ds]][[trt]]
  }
}
str(lapply(validation_list, names))


# ===== 4) Build training_list =====
training_list <- res_list
for (i in seq_len(nrow(val_tbl))) {
  ds  <- val_tbl$dataset[i]
  trt <- val_tbl$treatment[i]
  if (!is.null(training_list[[ds]][[trt]])) {
    training_list[[ds]][[trt]] <- NULL
    if (length(training_list[[ds]]) == 0L) training_list[[ds]] <- NULL
  }
}
str(lapply(training_list, names))
