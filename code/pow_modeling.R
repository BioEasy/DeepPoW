# =====================================================================
# code/pow_modeling.R
# Function: Robust PoW tree calibration engine preserving binary topology (Unconstrained version)
# =====================================================================

library(minpack.lm)
library(ape)
library(tidytree)

fit_pow_model <- function(df) {
  df <- df[df$rate > 1e-12, ]
  if (nrow(df) < 5) return(NULL)
  
  start_vals <- list(mu_0 = max(df$rate) * 1.2, mu_infty = min(df$rate) * 0.8, lambda = 1 / median(df$branch.length))
  
fit <- tryCatch({
    nlsLM(log(rate) ~ log(mu_infty + (mu_0 - mu_infty) * (1 - exp(-lambda * branch.length)) / (lambda * branch.length)),
          data = df,
          start = start_vals,
          lower = c(mu_0 = 1e-10, mu_infty = 1e-10, lambda = 1e-4), 
          control = nls.lm.control(maxiter = 1024))
  }, error = function(e) {
    NULL
  })
  
  return(fit)
}

predict_pow_with_ci <- function(fit, t_seq) {
  if (is.null(fit)) return(NULL)
  
  coefs <- coef(fit)
  mu_0 <- coefs["mu_0"]
  mu_infty <- coefs["mu_infty"]
  lambda <- coefs["lambda"]
  
  rate_pred <- mu_infty + (mu_0 - mu_infty) * (1 - exp(-lambda * t_seq)) / (lambda * t_seq)
  
  summary_fit <- summary(fit)
  residual_log_se <- summary_fit$sigma
  log_margin = 1.96 * residual_log_se * 0.3 * (1 / sqrt(t_seq + 0.1))
  
  data.frame(
    branch.length = t_seq,
    rate_pred = rate_pred,
    upper_ci = exp(log(rate_pred) + log_margin),
    lower_ci = pmax(exp(log(rate_pred) - log_margin), 1e-12)
  )
}

compare_models <- function(df, pow_fit) {
  if (is.null(pow_fit)) return(NULL)
  df <- df[df$rate > 1e-12, ]
  
  pow_rss <- sum(residuals(pow_fit)^2)
  pow_n <- nrow(df)
  pow_k <- length(coef(pow_fit)) + 1 
  pow_aic <- pow_n * log(pow_rss / pow_n) + 2 * pow_k
  pow_bic <- pow_n * log(pow_rss / pow_n) + log(pow_n) * pow_k
  
  mean_log_rate <- mean(log(df$rate))
  const_rss <- sum((log(df$rate) - mean_log_rate)^2)
  const_k <- 2
  const_aic <- pow_n * log(const_rss / pow_n) + 2 * const_k
  const_bic <- pow_n * log(const_rss / pow_n) + log(log(pow_n)) * const_k
  
  aic_values <- c(pow_aic, const_aic)
  bic_values <- c(pow_bic, const_bic)
  delta_aic <- aic_values - min(aic_values)
  delta_bic <- bic_values - min(bic_values)
  
  comparison_df <- data.frame(
    Model = c("PoW Model (TDRP)", "Constant Rate (Strict Baseline)"),
    Parameters = c(pow_k, const_k),
    RSS = c(signif(pow_rss, 4), signif(const_rss, 4)),
    AIC = c(signif(pow_aic, 5), signif(const_aic, 5)),
    Delta_AIC = c(signif(delta_aic[1], 4), signif(delta_aic[2], 4)),
    BIC = c(signif(pow_bic, 5), signif(const_bic, 5)),
    Delta_BIC = c(signif(delta_bic[1], 4), signif(delta_bic[2], 4))
  )
  return(comparison_df)
}

transform_pow_tree <- function(tree_obj, pow_fit, ...) {
  if (is.null(pow_fit) || is.null(tree_obj)) return(NULL)
  
  coefs <- coef(pow_fit)
  mu_0 <- as.numeric(coefs["mu_0"])
  mu_infty <- as.numeric(coefs["mu_infty"])
  lambda <- as.numeric(coefs["lambda"])
  
  new_tree <- tree_obj@phylo
  original_edges <- new_tree$edge.length
  
  for (i in seq_along(original_edges)) {
    dt <- original_edges[i]
    if (!is.na(dt) && dt > 0) {
      rate_pred <- mu_infty + (mu_0 - mu_infty) * (1 - exp(-lambda * dt)) / (lambda * dt)
      transformed_len <- (dt * rate_pred / mu_infty) * (dt * 1000)
      new_tree$edge.length[i] <- max(transformed_len, 1e-4)
    } else {
      new_tree$edge.length[i] <- 1e-4
    }
  }
  return(new_tree)
}

export_pow_pdf <- function(file_path, df, pow_fit) {
  pdf(file_path, width = 7, height = 5)
  coefs <- coef(pow_fit)
  t_min <- max(min(df$branch.length, na.rm = TRUE), 1e-5)
  t_max <- max(df$branch.length, na.rm = TRUE)
  t_seq <- exp(seq(log(t_min), log(t_max), length.out = 100))
  curve_df <- predict_pow_with_ci(pow_fit, t_seq)
  
  log_time <- log10(pmax(df$branch.length, 1e-12))
  log_rate <- log10(pmax(df$rate, 1e-12))
  log_t_seq <- log10(pmax(curve_df$branch.length, 1e-12))
  log_pred <- log10(pmax(curve_df$rate_pred, 1e-12))
  log_upper <- log10(pmax(curve_df$upper_ci, 1e-12))
  log_lower <- log10(pmax(curve_df$lower_ci, 1e-12))
  
  plot(log_time, log_rate, pch = 16, col = rgb(44/255, 62/255, 80/255, 0.5),
       xlab = expression(Log[10] * "(Time-scale, years)"),
       ylab = expression(Log[10] * "(Evolutionary Rate, " * mu * ")"),
       main = "Time-Dependent Rate Phenomenon (TDRP)",
       cex.main = 1.1, cex.lab = 1, cex.axis = 0.9)
  
  polygon(c(log_t_seq, rev(log_t_seq)), c(log_upper, rev(log_lower)), col = rgb(217/255, 83/255, 79/255, 0.2), border = NA)
  lines(log_t_seq, log_pred, col = "#D9534F", lwd = 2.5)
  points(log_time, log_rate, pch = 16, col = rgb(44/255, 62/255, 80/255, 0.6), cex = 0.8)
  
legend("topright", 
         legend = as.expression(c(
           bquote(mu[0] == .(format(coefs["mu_0"], scientific = TRUE, digits = 3))),
           bquote(mu[infinity] == .(format(coefs["mu_infty"], scientific = TRUE, digits = 3))),
           bquote(k == .(signif(coefs["lambda"], 3)))
         )), 
         bty = "n", cex = 0.8)
  dev.off()
}

robustness_pow_model <- function(df, trim_quantile = 0.05) {
  if (is.null(df) || nrow(df) < 5) return(NULL)
  
  fit_full <- fit_pow_model(df)
  if (is.null(fit_full)) return(NULL)
  lower_bound <- quantile(df$rate, probs = trim_quantile, na.rm = TRUE)
  upper_bound <- quantile(df$rate, probs = 1 - trim_quantile, na.rm = TRUE)
  df_trimmed <- subset(df, rate >= lower_bound & rate <= upper_bound)
  fit_trimmed <- fit_pow_model(df_trimmed)
  extract_coefs <- function(fit) {
    if (is.null(fit)) return(c(NA, NA, NA))
    c(coef(fit)["mu_0"], coef(fit)["mu_infty"], coef(fit)["lambda"])
  }
  
  full_c <- extract_coefs(fit_full)
  trim_c <- extract_coefs(fit_trimmed)
  
  comp_table <- data.frame(
    Parameter = c("mu_0 (Initial Rate)", "mu_infty (Long-term Rate)", "lambda (Decay Coef)"),
    Full_Data = c(signif(full_c[1], 4), signif(full_c[2], 4), signif(full_c[3], 4)),
    Trimmed_5pct = c(signif(trim_c[1], 4), signif(trim_c[2], 4), signif(trim_c[3], 4)),
    stringsAsFactors = FALSE
  )

  lambda_shift <- NA
  if (!is.na(trim_c[3]) && !is.na(full_c[3]) && full_c[3] != 0) {
    lambda_shift <- abs(trim_c[3] - full_c[3]) / abs(full_c[3])
  }
  
  comp_table$Shift_Percent <- ifelse(is.na(comp_table$Trimmed_5pct), 
                                     "Failed to Converge", 
                                     sprintf("%.2f%%", abs(trim_c - full_c) / abs(full_c) * 100))
  
  return(list(
    full_fit = fit_full,
    trimmed_fit = fit_trimmed,
    comparison_table = comp_table,
    lambda_shift = lambda_shift
  ))
}
