# =====================================================================
# code/pow_multitree.R
# Function: PoW transformation of sampled MCMC posterior trees based on Ghafari et al. (2021).
#           Features spline-based interpolation for faster computation, multicore parallelization,
#           preservation of the multiPhylo class, and support for both HKY85 and JC69 models.
# =====================================================================

library(ape)
library(nleqslv)
library(stringr)
library(parallel)

read_beast_log <- function(file_path) {
  df <- tryCatch({
    read.table(file_path, header = TRUE, sep = "\t", comment.char = "#",
               check.names = FALSE, fill = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(df) || ncol(df) <= 1) {
    df <- read.table(file_path, header = TRUE, comment.char = "#",
                     check.names = FALSE, fill = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
  }
  return(df)
}

find_freq_cols <- function(col_names) {
  f1 <- grep("frequencies1$|freqParameter\\.1$|freq\\.1$|freq_1$|freq_A$", col_names, ignore.case = TRUE)[1]
  f2 <- grep("frequencies2$|freqParameter\\.2$|freq\\.2$|freq_2$|freq_C$", col_names, ignore.case = TRUE)[1]
  f3 <- grep("frequencies3$|freqParameter\\.3$|freq\\.3$|freq_3$|freq_G$", col_names, ignore.case = TRUE)[1]
  f4 <- grep("frequencies4$|freqParameter\\.4$|freq\\.4$|freq_4$|freq_T$", col_names, ignore.case = TRUE)[1]
  
  if (!is.na(f1)) {
    if (is.na(f2)) f2 <- f1 + 1
    if (is.na(f3)) f3 <- f1 + 2
    if (is.na(f4)) f4 <- f1 + 3
  }
  return(c(f1, f2, f3, f4))
}

transform_single_tree <- function(el, importedtree_el, meanRate, pA, pC, pG, pT, kap, muMax_val) {
  pY <- pT + pC
  pR <- pA + pG
  delta <- 0.1
  betaMax <- muMax_val / (2 * ((pT * pC + pA * pG) * kap + pY * pR))
  
  if (is.na(betaMax) || is.infinite(betaMax) || betaMax <= 0) {
    stop(sprintf("betaMax calculation failed at tree %d.", el))
  }
  
  max_steps <- log10(betaMax) + 9
  steps <- seq(0, max_steps, delta)
  steps2 <- seq(1, length(steps), 1)
  
  lambda <- function(l) {
    sum(exp(l * (1/delta * steps + 1)) * 
          (2 * ((pT * pC + pA * pG) * kap + pY * pR)) * 10^(-9 + steps)) /
      sum(exp(l * steps2)) - meanRate
  }
  
  LAMBDA <- tryCatch(nleqslv(0, lambda)$x,
                     error = function(e) stop("Solving LAMBDA failed for tree ", el))
  proportions <- exp(LAMBDA * (1/delta * steps + 1)) / sum(exp(LAMBDA * steps2))
  ratespread <- 10^(-9 + steps)
  
  e2 <- function(t) exp(-ratespread * t)
  e3 <- function(t) exp(-(pR * kap + pY) * ratespread * t)
  e4 <- function(t) exp(-(pY * kap + pR) * ratespread * t)
  pTC <- function(t) pC + (pC * pR) / pY * e2(t) - pC / pY * e4(t)
  pAG <- function(t) pG + (pG * pY) / pR * e2(t) - pG / pR * e3(t)
  pTA <- function(t) pA * (1 - e2(t))
  pTG <- function(t) pG * (1 - e2(t))
  pCA <- function(t) pA * (1 - e2(t))
  pCG <- function(t) pG * (1 - e2(t))
  s1 <- function(t) sum(proportions * 2 * pT * pTC(t))
  s2 <- function(t) sum(proportions * 2 * pA * pAG(t))
  v <- function(t) sum(proportions * (2 * pT * pTA(t) + 2 * pT * pTG(t) + 2 * pC * pCA(t) + 2 * pC * pCG(t)))
  a1 <- function(t) -log(1 - (pY * s1(t)) / (2 * pT * pC) - v(t) / (2 * pY))
  a2 <- function(t) -log(1 - (pR * s2(t)) / (2 * pA * pG) - v(t) / (2 * pR))
  b <- function(t) -log(1 - v(t) / (2 * pY * pR))
  
  dist_func_scalar <- function(t) {
    if (t <= 0) return(0)
    val <- 2 * pT * pC / pY * (a1(t) - pR * b(t)) +
           (2 * pA * pG / pR) * (a2(t) - pY * b(t)) +
           2 * pY * pR * b(t)
    if (is.na(val) || is.nan(val)) return(0)
    return(val)
  }
  
  current_tree <- importedtree_el
  Ntip <- length(current_tree$tip.label)
  
  node_depths <- ape::node.depth.edgelength(current_tree)
  max_depth <- max(node_depths[1:Ntip])
  node_dist_from_tips <- pmax(max_depth - node_depths, 0)
  max_target_d <- max(node_dist_from_tips)
  
  # =========================================================
  # Spline Inversion
  # =========================================================
  sol_root <- tryCatch(nleqslv(100, function(t) dist_func_scalar(t) - max_target_d)$x, error = function(e) 1000)
  t_max_upper <- max(sol_root * 1.5, 100)
  
  t_mesh <- c(0, exp(seq(log(1e-5), log(t_max_upper), length.out = 300)))
  d_mesh <- sapply(t_mesh, dist_func_scalar)
  
  valid_idx <- !is.na(d_mesh) & !is.nan(d_mesh)
  t_mesh <- t_mesh[valid_idx]
  d_mesh <- d_mesh[valid_idx]
  
  dups <- duplicated(d_mesh)
  if (any(dups)) {
    t_mesh <- t_mesh[!dups]
    d_mesh <- d_mesh[!dups]
  }
  
  inv_spline <- approxfun(x = d_mesh, y = t_mesh, method = "linear", rule = 2)
  
  node_height_years <- inv_spline(node_dist_from_tips)
  node_height_years[node_dist_from_tips <= 1e-12] <- 0
  
  parent <- current_tree$edge[, 1]
  child  <- current_tree$edge[, 2]
  len_yr <- node_height_years[parent] - node_height_years[child]
  current_tree$edge.length <- pmax(len_yr, 1e-6)
  
  class(current_tree) <- "phylo"
  return(current_tree)
}

pow_transform <- function(rate_log_path, trees_path, hky_log_path = NULL,
                          muMax_val, sample_size = 100, sub_model = "HKY") {
  
  if (is.na(muMax_val) || muMax_val <= 0) stop("muMax value is invalid.")
  
  rate_df <- read_beast_log(rate_log_path)
  rate_col_idx <- grep("meanRate|rate\\.mean|clock\\.rate|clockRate|rate_median", colnames(rate_df), ignore.case = TRUE)[1]
  if (is.na(rate_col_idx)) stop("Could not identify a mean rate column in Rate Log.")
  
  burnin_rate <- max(1, floor(nrow(rate_df) * 0.1))
  rate_df <- rate_df[(burnin_rate + 1):nrow(rate_df), , drop = FALSE]
  if (nrow(rate_df) < sample_size) stop("Rate Log has too few post-burn-in lines.")
  
  sampled_idx_rate <- sample(1:nrow(rate_df), sample_size)
  rates <- suppressWarnings(as.numeric(as.character(rate_df[sampled_idx_rate, rate_col_idx])))
  
  if (sub_model == "JC69") {
    K     <- rep(1.0, sample_size)
    freq1 <- rep(0.25, sample_size)
    freq2 <- rep(0.25, sample_size)
    freq3 <- rep(0.25, sample_size)
    freq4 <- rep(0.25, sample_size)
  } else {
    if (is.null(hky_log_path) || !file.exists(hky_log_path)) {
      stop("HKY Log file is required when Substitution Model is set to HKY85.")
    }
    hky_df <- read_beast_log(hky_log_path)
    kappa_idx <- grep("^kappa$|kappa", colnames(hky_df), ignore.case = TRUE)[1]
    freq_indices <- find_freq_cols(colnames(hky_df))
    
    if (is.na(kappa_idx)) stop("Column 'kappa' not found in HKY Log.")
    if (any(is.na(freq_indices))) stop("Could not locate nucleotide frequency columns in HKY Log.")
    
    burnin_hky <- max(1, floor(nrow(hky_df) * 0.1))
    hky_df <- hky_df[(burnin_hky + 1):nrow(hky_df), , drop = FALSE]
    if (nrow(hky_df) < sample_size) stop("HKY Log has too few post-burn-in lines.")
    
    sampled_idx_hky <- sample(1:nrow(hky_df), sample_size)
    K     <- suppressWarnings(as.numeric(as.character(hky_df[sampled_idx_hky, kappa_idx])))
    freq1 <- suppressWarnings(as.numeric(as.character(hky_df[sampled_idx_hky, freq_indices[1]])))
    freq2 <- suppressWarnings(as.numeric(as.character(hky_df[sampled_idx_hky, freq_indices[2]])))
    freq3 <- suppressWarnings(as.numeric(as.character(hky_df[sampled_idx_hky, freq_indices[3]])))
    freq4 <- suppressWarnings(as.numeric(as.character(hky_df[sampled_idx_hky, freq_indices[4]])))
  }
  
  con2 <- readLines(trees_path, warn = FALSE)
  tree_start <- grep("^tree", con2, ignore.case = TRUE)[1]
  tree_lines <- con2[tree_start:(length(con2) - 1)]
  
  burnin_trees <- max(1, floor(length(tree_lines) * 0.1))
  tree_lines <- tree_lines[(burnin_trees + 1):length(tree_lines)]
  if (length(tree_lines) < sample_size) stop("Trees file after burn-in has too few trees.")
  
  sampled_trees <- sample(tree_lines, sample_size)
  tmp_nexus <- tempfile(fileext = ".trees")
  writeLines(c(con2[1:(tree_start - 1)], sampled_trees, "End;"), tmp_nexus)
  importedtree <- read.nexus(tmp_nexus, force.multi = TRUE)
  unlink(tmp_nexus)
  
  num_cores <- max(1, parallel::detectCores() - 1)
  
  if (.Platform$OS.type == "unix") {
    transformed_list <- parallel::mclapply(1:sample_size, function(el) {
      transform_single_tree(el, importedtree[[el]], rates[el], freq1[el], freq2[el], freq3[el], freq4[el], K[el], muMax_val)
    }, mc.cores = num_cores)
  } else {
    cl <- parallel::makeCluster(num_cores)
    parallel::clusterExport(cl, varlist = c("transform_single_tree", "importedtree", "rates", "freq1", "freq2", "freq3", "freq4", "K", "muMax_val"), envir = environment())
    parallel::clusterEvalQ(cl, { library(ape); library(nleqslv) })
    transformed_list <- parallel::parLapply(cl, 1:sample_size, function(el) {
      transform_single_tree(el, importedtree[[el]], rates[el], freq1[el], freq2[el], freq3[el], freq4[el], K[el], muMax_val)
    })
    parallel::stopCluster(cl)
  }
  
  class(transformed_list) <- "multiPhylo"
  return(transformed_list)
}