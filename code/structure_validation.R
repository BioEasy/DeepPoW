# =====================================================================
# code/structure_validation.R
# Function: Extract phylogenetic distances, align with external RMSD
# matrix, run Mantel Test, and generate publication-ready plots.
# =====================================================================

library(ape)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(treeio)

# =====================================================================
# 1. Resolve PoW-transformed tree
# =====================================================================

resolve_pow_tree <- function(input_tree = NULL,
                             cache_path = "cache/pow_transformed_tree.tree") {
  tree_obj <- NULL
  if (!is.null(input_tree)) {
    message("[DeepPoW] Using active session PoW-transformed tree object.")
    tree_obj <- input_tree
  } else if (file.exists(cache_path)) {
    message(paste("[DeepPoW] Found cached tree at:", cache_path,
                  ". Loading automatically..."))
    tree_obj <- load_pow_tree_robust(cache_path)
  } else {
    message("[DeepPoW] No active tree or cache detected. Please select your PoW-transformed tree file...")
    if (interactive()) {
      chosen_file <- file.choose()
      tree_obj <- load_pow_tree_robust(chosen_file)
    } else {
      stop("[DeepPoW Error] No valid PoW-transformed tree provided or found in cache.")
    }
  }
  return(tree_obj)
}

# =====================================================================
# 2. Robust tree loader
# =====================================================================

load_pow_tree_robust <- function(file_path) {
  tree_obj <- tryCatch({
    parsed <- treeio::read.beast(file_path)
    if (inherits(parsed, "treedata")) parsed@phylo else parsed
  }, error = function(e1) {
    tryCatch({
      res <- ape::read.nexus(file_path)
      if (inherits(res, "multiPhylo")) res[[1]] else res
    }, error = function(e2) {
      ape::read.tree(file_path)
    })
  })
  return(tree_obj)
}

# =====================================================================
# 3. Calculate pairwise phylogenetic distances
# =====================================================================

calculate_phylo_distance <- function(tree_obj) {
  if (inherits(tree_obj, "treedata")) tree_obj <- tree_obj@phylo
  if (!inherits(tree_obj, "phylo")) stop("Invalid tree format loaded.")

  tree_obj$tip.label <- trimws(gsub("['\"]", "", tree_obj$tip.label))
  dist_matrix <- ape::cophenetic.phylo(tree_obj)

  phylo_long <- as.data.frame(as.table(dist_matrix),
                              stringsAsFactors = FALSE) %>%
    dplyr::rename(Taxon1 = Var1, Taxon2 = Var2, Phylo_Dist = Freq) %>%
    dplyr::filter(as.character(Taxon1) != as.character(Taxon2)) %>%
    dplyr::mutate(
      Taxon1 = as.character(Taxon1),
      Taxon2 = as.character(Taxon2),
      Pair1 = pmin(Taxon1, Taxon2),
      Pair2 = pmax(Taxon1, Taxon2),
      Pair_Key = paste(Pair1, Pair2, sep = "::")
    ) %>%
    dplyr::distinct(Pair_Key, .keep_all = TRUE) %>%
    dplyr::select(Taxon1, Taxon2, Pair_Key, Phylo_Dist)

  return(phylo_long)
}

# =====================================================================
# 4. Taxon matching diagnostic table
# =====================================================================

check_taxon_matching <- function(tree_obj, rmsd_matrix_path) {
  tree_taxa <- character(0)

  if (!is.null(tree_obj)) {
    if (inherits(tree_obj, "treedata")) tree_obj <- tree_obj@phylo
    if (inherits(tree_obj, "phylo")) {
      tree_obj$tip.label <- trimws(gsub("['\"]", "", tree_obj$tip.label))
      tree_taxa <- sort(unique(as.character(tree_obj$tip.label)))
    }
  }

  matrix_taxa <- character(0)

  if (!is.null(rmsd_matrix_path) && file.exists(rmsd_matrix_path)) {
    rmsd_raw <- tryCatch(
      read.csv(rmsd_matrix_path, row.names = 1,
               check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    if (!is.null(rmsd_raw)) {
      clean_names <- function(nms) {
        nms <- as.character(nms)
        nms <- trimws(gsub("['\"]", "", nms))
        nms <- gsub("fold_|_model_0|_model_1", "", nms)
        nms
      }
      matrix_taxa <- sort(unique(clean_names(rownames(rmsd_raw))))
    }
  }

  if (length(tree_taxa) == 0 && length(matrix_taxa) == 0) {
    return(data.frame(
      Taxon_ID = "N/A",
      Status = "Please upload an RMSD matrix and ensure a PoW tree is loaded.",
      stringsAsFactors = FALSE
    ))
  }

  all_taxa <- sort(unique(c(tree_taxa, matrix_taxa)))

  status_vec <- rep("[Missing in Matrix]", length(all_taxa))
  status_vec[all_taxa %in% tree_taxa & all_taxa %in% matrix_taxa] <- "[Matched]"
  status_vec[!all_taxa %in% tree_taxa & all_taxa %in% matrix_taxa] <- "[Missing in Tree]"

  data.frame(
    Taxon_ID = as.character(all_taxa),
    Status = as.character(status_vec),
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# 5. Match phylogenetic distances with RMSD/structure matrix
# =====================================================================

match_structure_matrix <- function(phylo_df, rmsd_matrix_path) {
  if (is.null(rmsd_matrix_path) || !file.exists(rmsd_matrix_path)) {
    stop("[DeepPoW Error] RMSD matrix file path is invalid or missing.")
  }

  rmsd_raw <- tryCatch(
    read.csv(rmsd_matrix_path, row.names = 1,
             check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("[DeepPoW Error] Failed to read CSV matrix:", e$message))
    }
  )

  clean_names <- function(nms) {
    nms <- as.character(nms)
    nms <- trimws(gsub("['\"]", "", nms))
    nms <- gsub("fold_|_model_0|_model_1", "", nms)
    nms
  }

  rownames(rmsd_raw) <- clean_names(rownames(rmsd_raw))
  colnames(rmsd_raw) <- clean_names(colnames(rmsd_raw))

  rmsd_matrix <- as.matrix(rmsd_raw)
  suppressWarnings(storage.mode(rmsd_matrix) <- "numeric")

  rmsd_long <- as.data.frame(as.table(rmsd_matrix),
                             stringsAsFactors = FALSE) %>%
    dplyr::rename(Taxon1 = Var1, Taxon2 = Var2, Struct_Dist = Freq) %>%
    dplyr::filter(as.character(Taxon1) != as.character(Taxon2)) %>%
    dplyr::mutate(
      Taxon1 = as.character(Taxon1),
      Taxon2 = as.character(Taxon2),
      Struct_Dist = suppressWarnings(as.numeric(as.character(Struct_Dist)))
    ) %>%
    dplyr::filter(!is.na(Struct_Dist)) %>%
    dplyr::mutate(
      Pair1 = pmin(Taxon1, Taxon2),
      Pair2 = pmax(Taxon1, Taxon2),
      Pair_Key = paste(Pair1, Pair2, sep = "::")
    ) %>%
    dplyr::distinct(Pair_Key, .keep_all = TRUE) %>%
    dplyr::select(Pair_Key, Struct_Dist)

  matched_df <- phylo_df %>%
    dplyr::inner_join(rmsd_long, by = "Pair_Key")

  if (nrow(matched_df) < 3) {
    stop("[DeepPoW Error] Less than 3 overlapping pairs found! Mantel test requires at least 3 pairs. Check the Diagnostic Table.")
  }

  matched_df <- matched_df %>% dplyr::select(-Pair_Key)
  matched_df$Taxon1 <- as.character(matched_df$Taxon1)
  matched_df$Taxon2 <- as.character(matched_df$Taxon2)
  matched_df$Phylo_Dist <- as.numeric(matched_df$Phylo_Dist)
  matched_df$Struct_Dist <- as.numeric(matched_df$Struct_Dist)

  return(matched_df)
}

# =====================================================================
# 6. Run Mantel Test
# =====================================================================

run_mantel_test <- function(matched_df, permutations = 999) {
  taxa <- unique(c(as.character(matched_df$Taxon1),
                   as.character(matched_df$Taxon2)))
  taxa <- as.character(taxa)
  num_taxa <- length(taxa)

  phylo_mat <- matrix(0, nrow = num_taxa, ncol = num_taxa,
                      dimnames = list(taxa, taxa))
  struct_mat <- matrix(0, nrow = num_taxa, ncol = num_taxa,
                       dimnames = list(taxa, taxa))

  for (i in seq_len(nrow(matched_df))) {
    t1 <- as.character(matched_df$Taxon1[i])
    t2 <- as.character(matched_df$Taxon2[i])
    phylo_value <- as.numeric(matched_df$Phylo_Dist[i])
    struct_value <- as.numeric(matched_df$Struct_Dist[i])

    phylo_mat[t1, t2] <- phylo_mat[t2, t1] <- phylo_value
    struct_mat[t1, t2] <- struct_mat[t2, t1] <- struct_value
  }

  mantel_res <- vegan::mantel(
    as.dist(phylo_mat),
    as.dist(struct_mat),
    method = "pearson",
    permutations = as.integer(permutations)
  )

  r_val <- round(as.numeric(unname(mantel_res$statistic)), 4)
  p_val <- as.numeric(unname(mantel_res$signif))
  p_str <- if (p_val < 0.001) {
    "< 0.001"
  } else {
    paste0("= ", round(p_val, 4))
  }

  list(
    r_statistic = r_val,
    p_value = p_val,
    p_string = as.character(p_str),
    permutations = as.integer(permutations),
    n_pairs = as.integer(nrow(matched_df))
  )
}

# =====================================================================
# 7. Publication-ready Mantel plot
# =====================================================================

plot_mantel_publication <- function(matched_df,
                                    mantel_res,
                                    virus_name = "Sobemovirus") {
  if (is.null(matched_df) || nrow(matched_df) == 0 || is.null(mantel_res)) {
    return(
      ggplot() +
        annotate("text", x = 1, y = 1,
                 label = "No valid data to plot.",
                 size = 4, color = "red") +
        theme_void()
    )
  }

  plot_df <- matched_df
  plot_df$Phylo_Dist <- as.numeric(unname(plot_df$Phylo_Dist))
  plot_df$Struct_Dist <- as.numeric(unname(plot_df$Struct_Dist))

  plot_df <- plot_df[
    is.finite(plot_df$Phylo_Dist) &
      is.finite(plot_df$Struct_Dist),
    ,
    drop = FALSE
  ]

  r_val <- as.numeric(unname(mantel_res$r_statistic))
  p_val <- as.numeric(unname(mantel_res$p_value))

  p_text <- if (p_val < 0.001) {
    "p < 0.001"
  } else {
    paste0("p = ", format(round(p_val, 3), nsmall = 3))
  }

  annot_text <- paste0(
    "Mantel r = ", format(r_val, nsmall = 4),
    ", ", p_text,
    " (n = ", nrow(plot_df), " pairs)"
  )

  ggplot(plot_df, aes(x = Phylo_Dist, y = Struct_Dist)) +
    geom_point(
      color = "#2C3E50", 
      alpha = 0.6,          
      size = 2.2
    ) +
geom_smooth(aes(fill = "95% CI"), method = "lm", color = "#D9534F", alpha = 0.2) +
scale_fill_manual(name = "", values = c("95% CI" = "#D9534F")) +
theme(legend.position = "bottom") +
    labs(
      title = paste0(
        "Structural & Phylogenetic Congruence (",
        virus_name,
        ")"
      ),
      subtitle = annot_text,
      x = "PoW-Calibrated Phylogenetic Distance",
      y = "Protein Structural Distance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold", size = 16, hjust = 0.5,
        color = "black", margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        size = 12, hjust = 0.5, color = "#333333",
        margin = margin(b = 10)
      ),
      axis.title = element_text(
        face = "bold", size = 13, color = "black"
      ),
      axis.text = element_text(
        size = 11, color = "#4A4A4A"
      ),
      panel.grid.major = element_line(
        color = "#EBEBEB", linewidth = 0.7
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "#333333", fill = NA, linewidth = 0.7
      ),
      plot.margin = margin(12, 15, 12, 15)
    )
}