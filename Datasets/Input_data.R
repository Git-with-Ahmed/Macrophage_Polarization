# 1) Locate files
base <- "C:/Users/Ahmed/OneDrive - Rutgers University/Birge Lab/Macrophage_Polarization/Datasets"
setwd(base)
files <- list.files(pattern = "_counts\\.csv$", full.names = TRUE)

# 2) Process each file → renamed table + metadata
process_one <- function(path) {
  ds <- tools::file_path_sans_ext(basename(path))
  x  <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(x)[1] <- "gene"
  
  orig <- colnames(x)[-1]                 # original sample/treatment names
  new  <- paste0("S", seq_along(orig))    # S1..Sx
  
  meta <- data.frame(
    dataset   = ds,
    new_name  = new,
    treatment = orig,
    stringsAsFactors = FALSE
  )
  
  colnames(x) <- c("gene", new)          # rename count columns
  
  list(counts = x, meta = meta)
}

data_list <- setNames(lapply(files, process_one), tools::file_path_sans_ext(basename(files)))