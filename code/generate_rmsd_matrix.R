# =====================================================================
# code/generate_rmsd_matrix.R
# Function: Native R execution pipeline for Foldseek or TM-align to 
#           generate an all-to-all structural RMSD matrix.
# =====================================================================

#' Execute structural alignment pipeline natively
#' 
#' @param engine String. Either "Foldseek" or "TM-align"
#' @param pdb_dir String. Directory containing .pdb or .cif files
#' @param output_csv String. Output file path for the CSV matrix
#' @param custom_exe_path String. Optional path to the executable
#' @return TRUE if successful, stops with error otherwise
run_rmsd_pipeline <- function(engine, pdb_dir, output_csv, custom_exe_path = "") {
  
  if (!dir.exists(pdb_dir)) {
    stop(paste("Target directory does not exist:", pdb_dir))
  }
  
  files <- list.files(pdb_dir, pattern = "\\.(pdb|cif)$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop(paste("No .pdb or .cif files found in directory:", pdb_dir))
  }
  
  n <- length(files)
  basenames <- tools::file_path_sans_ext(basename(files))
  
  rmsd_mat <- matrix(NA, nrow = n, ncol = n, dimnames = list(basenames, basenames))
  diag(rmsd_mat) <- 0
  
  # =========================================================
  # Core Mechanism: Smart Binary Search
  # =========================================================
  find_binary <- function(bin_names, custom_path = "") {

    if (custom_path != "" && file.exists(custom_path)) {
      return(custom_path)
    }
    
    for (b in bin_names) {
      res <- Sys.which(b)
      if (res != "") return(unname(res))
    }
    
    home_dir <- Sys.getenv("HOME")
    common_paths <- c(
      "/usr/local/bin",
      "/opt/homebrew/bin",
      "/usr/bin",
      file.path(home_dir, "miniconda3/bin"),
      file.path(home_dir, "anaconda3/bin"),
      file.path(home_dir, ".local/bin")
    )
    
    for (p in common_paths) {
      for (b in bin_names) {
        full_path <- file.path(p, b)
        if (file.exists(full_path)) return(full_path)
    
        full_path_exe <- file.path(p, paste0(b, ".exe"))
        if (file.exists(full_path_exe)) return(full_path_exe)
      }
    }
    
    return("") 
  }
  
  # =========================================================
  # 2. Foldseek Execution Pathway
  # =========================================================
  if (grepl("Foldseek", engine, ignore.case = TRUE)) {
    
    fs_bin <- find_binary(c("foldseek"), custom_exe_path)
    
    if (fs_bin == "") {
      stop("Foldseek executable not found automatically. Please locate the file using the 'Browse' button under Executable Path.")
    }
    
    tmp_dir <- tempfile("fs_tmp")
    out_tsv <- tempfile("fs_out", fileext = ".tsv")
    dir.create(tmp_dir, showWarnings = FALSE)
    message(sprintf("[INFO] Auto-detected Foldseek at: %s", fs_bin))
    message(sprintf("[INFO] Starting Foldseek fast all-vs-all on %d structures...", n))
    
    cmd_args <- c("easy-search", shQuote(pdb_dir), shQuote(pdb_dir), shQuote(out_tsv), shQuote(tmp_dir),
                  "--format-output", "query,target,rmsd", "-e", "inf")
    
    res <- system2(fs_bin, args = cmd_args, stdout = TRUE, stderr = TRUE)
    
    if (!file.exists(out_tsv) || file.info(out_tsv)$size == 0) {
      stop("Foldseek failed to produce valid alignments.\nLogs:\n", paste(res, collapse = "\n"))
    }
    
    fs_data <- read.table(out_tsv, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(fs_data) <- c("query", "target", "rmsd")
    
    for (i in seq_len(nrow(fs_data))) {
      q <- tools::file_path_sans_ext(basename(fs_data$query[i]))
      t <- tools::file_path_sans_ext(basename(fs_data$target[i]))
      if (q %in% basenames && t %in% basenames) {
        rmsd_mat[q, t] <- fs_data$rmsd[i]
      }
    }
    
    for (i in 1:n) {
      for (j in i:n) {
        if (is.na(rmsd_mat[i, j]) && !is.na(rmsd_mat[j, i])) rmsd_mat[i, j] <- rmsd_mat[j, i]
        if (!is.na(rmsd_mat[i, j]) && is.na(rmsd_mat[j, i])) rmsd_mat[j, i] <- rmsd_mat[i, j]
      }
    }
    
    unlink(tmp_dir, recursive = TRUE)
    unlink(out_tsv)
    
  # =========================================================
  # 3. TM-align Execution Pathway
  # =========================================================
  } else if (grepl("TM-align", engine, ignore.case = TRUE)) {
    
    tm_bin <- find_binary(c("TMalign", "TM-align"), custom_exe_path)
    
    if (tm_bin == "") {
      stop("TM-align executable not found automatically. Please locate the file using the 'Browse' button under Executable Path.")
    }
    
    message(sprintf("[INFO] Auto-detected TM-align at: %s", tm_bin))
    message(sprintf("[INFO] Starting TM-align all-vs-all on %d structures (Computationally intensive)...", n))
    
    total_pairs <- (n * (n - 1)) / 2
    counter <- 0
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        counter <- counter + 1
        
        if (counter %% max(1, floor(total_pairs/10)) == 0) {
          message(sprintf("  -> Progress: %d / %d pairs aligned (%.1f%%)...", 
                          counter, total_pairs, (counter/total_pairs)*100))
        }
        
        tm_out <- system2(tm_bin, args = c(shQuote(files[i]), shQuote(files[j])), stdout = TRUE)
        rmsd_line <- grep("RMSD=", tm_out, value = TRUE)
        
        if (length(rmsd_line) > 0) {
          rmsd_val <- as.numeric(sub(".*RMSD=\\s*([0-9.]+).*", "\\1", rmsd_line[1]))
          rmsd_mat[i, j] <- rmsd_val
          rmsd_mat[j, i] <- rmsd_val 
        }
      }
    }
  } else {
    stop("Unknown alignment engine selected. Please choose Foldseek or TM-align.")
  }
  
  # =========================================================
  # 4. Final Output Generation
  # =========================================================
  write.csv(rmsd_mat, file = output_csv, row.names = TRUE)
  message("\n[SUCCESS] Matrix computation completed normally.")
  message(sprintf("[INFO] Generated an RMSD matrix for %d structures.", n))
  
  return(TRUE)
}