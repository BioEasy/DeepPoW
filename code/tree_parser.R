# =====================================================================
# code/tree_parsing.R
# Function: Parse annotated time-scaled trees from BEAST (1.x/2.x) and TreeTime,
#           extracting evolutionary rates and branch lengths.
# =====================================================================

library(ape)
library(treeio)
library(dplyr)
library(tidytree)

parse_beast_tree <- function(file_path, origin = "beast") {
  tryCatch({
    if (origin == "treetime") {
      tree_obj <- tryCatch(
        treeio::read.beast(file_path), 
        error = function(e) treeio::read.newick(file_path, comment.char = "[")
      )
    } else {
      tree_obj <- treeio::read.beast(file_path)
    }
    tree_data <- tidytree::as_tibble(tree_obj)
    return(list(tree = tree_obj, data = tree_data))
  }, error = function(e) {
    stop(paste("Tree parsing failed! Please ensure the file is a valid annotated tree (BEAST or TreeTime). Details:", e$message))
  })
}

extract_node_features <- function(tree_data, rate_col = NULL) {
  available_cols <- colnames(tree_data)
  has_rate <- any(grepl("rate$|rate_median$|rate_mean$|mut_rate$|branch_rate$|clockrate$|mutationrate$|mu$", available_cols, ignore.case = TRUE))
  
  if (!has_rate && "mutations" %in% available_cols) {
    message("[INFO] Detected TreeTime 'mutations' string. Inferring sequence length and calculating per-site substitution rates...")
    
    all_muts <- unlist(tree_data$mutations)
    all_muts <- all_muts[!is.na(all_muts) & all_muts != "" & all_muts != "NA"]
    
    if (length(all_muts) > 0) {
      pos_nums <- suppressWarnings(as.numeric(gsub("[A-Za-z-]", "", all_muts)))
      seq_length <- max(pos_nums, na.rm = TRUE)
      if (is.infinite(seq_length) || is.na(seq_length) || seq_length <= 0) seq_length <- 1000 
    } else {
      seq_length <- 1000
    }
    
    message(sprintf("[INFO] Auto-inferred sequence length from mutations: %d bp", seq_length))
    
    bl <- as.numeric(tree_data$branch.length)
    
    mut_counts <- sapply(tree_data$mutations, function(m) {
      m_clean <- m[!is.na(m) & m != "" & m != "NA"]
      if (length(m_clean) == 0) return(0)
      if (length(m_clean) == 1 && grepl(",", m_clean[1])) {
        return(length(strsplit(as.character(m_clean[1]), ",")[[1]]))
      }
      return(length(m_clean))
    })
    
    tree_data$mut_rate <- (mut_counts / seq_length) / (bl + 1e-6)
    available_cols <- colnames(tree_data)
  }
  

  if (is.null(rate_col) || rate_col == "Auto-detect") {
    
    candidates <- grep("rate$|rate_median$|rate_mean$|mut_rate$|branch_rate$|clockrate$|mutationrate$|mu$", available_cols, value = TRUE, ignore.case = TRUE)
    
    if (length(candidates) > 0) {
      seq_rates <- candidates[!grepl("Region|Location|Trait|Host|State", candidates, ignore.case = TRUE)]
      
      if (length(seq_rates) > 0) {
        rate_col <- seq_rates[1]
      } else {
        rate_col <- candidates[1]
      }
    }
    
    if (is.null(rate_col)) {
      stop(paste0(
        "Auto-detect failed. \n",
        "Your tree contains the following annotations: [ ", paste(available_cols, collapse = ", "), " ]. \n",
        "None of these match standard rate suffix patterns. Please select the column manually from the dropdown."
      ))
    }
  } else {
    if (!(rate_col %in% available_cols)) {
      stop(paste0("Selected rate column '", rate_col, "' not found. Available columns are: [ ", paste(available_cols, collapse = ", "), " ]"))
    }
  }
  
  df <- tree_data %>%
    dplyr::select(node, branch.length, all_of(rate_col)) %>%
    dplyr::rename(rate = !!sym(rate_col)) %>%
    dplyr::filter(!is.na(branch.length), branch.length > 0)
  
  if (is.list(df$rate)) {
    df$rate <- sapply(df$rate, function(x) as.numeric(x[1]))
  }
  
  df$rate <- as.numeric(df$rate)
  df$branch.length <- as.numeric(df$branch.length)
  df <- df %>% dplyr::filter(!is.na(rate), rate > 1e-12)
  
  return(df)
}

filter_extremes <- function(df, min_time = 0.0, max_rate_q = 0.99) {
  rate_threshold <- quantile(df$rate, probs = max_rate_q, na.rm = TRUE)
  filtered_df <- df %>%
    dplyr::filter(branch.length >= min_time) %>%
    dplyr::filter(rate <= rate_threshold)
  return(filtered_df)
}