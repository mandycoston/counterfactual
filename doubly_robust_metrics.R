
# This file contains functions for computing counterfactual metrics via doubly-robust estimation. 

# For axis formatting
scale_format <- function(x) {
  ifelse((x > 0) & (x < 1), sprintf("%.2f", x), sprintf("%.0f", x))
}

compute_FNR <- function(t, dat, pred = quo(mu_obs), treat = quo(A), pi = quo(pi), Y = quo(Y), mu_true = quo(mu_count)) {
  pred_label <- pull(dat, !!pred) >= t
  y_true <- (1 - pull(dat, !!treat)) / (1 - pull(dat, !!pi)) * (pull(dat, !!Y) - pull(dat, !!mu_true)) + pull(dat, !!mu_true)
  return(mean((1 - pred_label) * y_true) / mean(y_true))
}

compute_FPR <- function(t, dat, pred = quo(mu_obs), treat = quo(A), pi = quo(pi), Y = quo(Y), mu_true = quo(mu_count)) {
  pred_label <- pull(dat, !!pred) >= t
  y_true <- 1 - (1 - pull(dat, !!treat)) / (1 - pull(dat, !!pi)) * (pull(dat, !!Y) - pull(dat, !!mu_true)) - pull(dat, !!mu_true)
  return(mean((pred_label) * y_true) / mean(y_true))
}

# Compute doubly robust estimate of counterfactual precision
compute_precision <- function(t, dat, pred = quo(mu_obs), treat = quo(A), pi = quo(pi), Y = quo(Y), mu_true = quo(mu_count)) {
  pred_pos <- filter(dat, (!!pred) >= t) 
  nuis <- (1 - pull(pred_pos, !!treat)) / (1 - pull(pred_pos, !!pi)) * (pull(pred_pos, !!Y) - pull(pred_pos, !!mu_true)) + pull(pred_pos, !!mu_true)
  precision <- list("precision" = mean(nuis), "precision.lower" = mean(nuis) - 1.96 * sqrt(var(nuis) / nrow(pred_pos)), "precision.upper" = mean(nuis) + 1.96 * sqrt(var(nuis) / nrow(pred_pos)))
  return(precision)
}

# Function to create dataframe for plotting PR curve
compute_PR_df_dr <- function(t_arr, dat, group, obs_preds = quo(mu_obs), count_preds = quo(mu_count), treat = quo(A), pi = quo(pi), Y = quo(Y), mu_true = quo(mu_count)) {
  precision_obs <- sapply(t_arr, compute_precision, dat = dat, pred = obs_preds, treat, pi, Y, mu_true)
  precision_obs_low <- apply(precision_obs, 2, function(x) {
    x$precision.lower
  })
  precision_obs_high <- apply(precision_obs, 2, function(x) {
    x$precision.upper
  })
  precision_obs <- apply(precision_obs, 2, function(x) {
    x$precision
  })
  fnr_obs <- sapply(t_arr, compute_FNR, dat = dat, pred = obs_preds, treat, pi, Y, mu_true)
  df_pr_obs <- data.frame("Threshold" = t_arr, "Precision" = precision_obs, "Precision.lower" = precision_obs_low, "Precision.upper" = precision_obs_high, "Recall" = 1 - fnr_obs, "Method" = "Observational", "Group" = group)
  
  precision_count <- sapply(t_arr, compute_precision, dat = dat, pred = count_preds, treat, pi, Y, mu_true)
  precision_count_low <- apply(precision_count, 2, function(x) {
    x$precision.lower
  })
  precision_count_high <- apply(precision_count, 2, function(x) {
    x$precision.upper
  })
  precision_count <- apply(precision_count, 2, function(x) {
    x$precision
  })
  fnr_count <- sapply(t_arr, compute_FNR, dat = dat, pred = count_preds, treat, pi, Y, mu_true)
  df_pr_count <- data.frame("Threshold" = t_arr, "Precision" = precision_count, "Precision.lower" = precision_count_low, "Precision.upper" = precision_count_high, "Recall" = 1 - fnr_count, "Method" = "Counterfactual", "Group" = group)
  
  return(rbind(df_pr_count, df_pr_obs))
}

# Function creates dataframe for plotting ROC curves. The argument names are more generic than compute_PR_df_dr since we use this to also evaluate the post-processing model.
compute_ROC_df_dr <- function(t_arr, dat, group, method1 = quo(mu_obs), method2 = quo(mu_count), method1_name = "Observational", method2_name = "Counterfactual", treat = quo(A), pi = quo(pi), Y = quo(Y), mu_true = quo(mu_count)) {
  fpr_method1 <- sapply(t_arr, compute_FPR, dat = dat, pred = method1, treat, pi, Y, mu_true)
  fnr_method1 <- sapply(t_arr, compute_FNR, dat = dat, pred = method1, treat, pi, Y, mu_true)
  df_method1 <- data.frame("Threshold" = t_arr, "FPR" = fpr_method1, "Recall" = 1 - fnr_method1, "Method" = method1_name, "Group" = group)
  
  fpr_method2 <- sapply(t_arr, compute_FPR, dat = dat, pred = method2, treat, pi, Y, mu_true)
  fnr_method2 <- sapply(t_arr, compute_FNR, dat = dat, pred = method2, treat, pi, Y, mu_true)
  # fnr.caus.arr <- apply(fnr.caus.arr, 2, function(x) {x$fnr})
  df_method2 <- data.frame("Threshold" = t_arr, "FPR" = fpr_method2, "Recall" = 1 - fnr_method2, "Method" = method2_name, "Group" = group)
  
  return(rbind(df_method2, df_method1))
}

# Function to compute quantile-based calibration
# Num bins specifies how many bins to segment the risk scores into
compute_calibration_dr <- function(num_bins = 10, mu_true, pred, pi, Y, A, method_name, group) {
  steps <- seq(1 / num_bins, 1, 1 / num_bins)
  threshs <- sapply(steps, function(x) {
    unname(quantile(pred, probs = x))
  })
  calib <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(calib) <- c("Average score", "Rate", "Low", "High")
  i <- 1
  while (i <= num_bins) {
    prev <- ifelse(i > 1, threshs[i - 1], 0)
    sub <- (pred > prev & pred <= threshs[i])
    est <- ((1 - A[sub]) / (1 - pi[sub])) * (Y[sub] - mu_true[sub]) + mu_true[sub]
    calib <- rbind(calib, data.frame("Ventile risk score" = i, "Average score" = mean(pred[sub]), "Rate" = mean(est), "Low" = mean(est) - 1.96 * sqrt(var(est) / length(est)), "High" = mean(est) + 1.96 * sqrt(var(est) / length(est))))
    i <- i + 1
  }
  calib$Method <- method_name
  calib$Group <- group
  return(calib)
}

# wrapper function for computing calibration
# standardize order of arguments across functions i.e.
# compute_PR_df_dr <- function(t_arr, dat, group, obs_preds = quo(mu_obs), count_preds = quo(mu_count), treat = quo(A), pi = quo(pi), Y = quo(Y), mu_true = quo(mu_count)) {

compute_calib_df_dr <- function(num_bins, dat,obs_preds = quo(mu_obs),  count_preds = quo(mu_count), treat = quo(A),   pi = quo(pi), Y = quo(Y), mu_true = quo(mu_count),   attr = quo(attr), attr_name, attr_name_other) {
  all_count <- compute_calibration_dr(num_bins, pull(dat, (!!mu_true)), pull(dat, (!!count_preds)), pull(dat, (!!pi)), pull(dat, (!!Y)), pull(dat, (!!treat)), "Counterfactual", "All")
  attr_count <- compute_calibration_dr(num_bins, pull(filter(dat,(!! attr) == 1), (!!mu_true)), pull(filter(dat,(!! attr) == 1), (!!count_preds)), pull(filter(dat, !! attr == 1), (!!pi)), pull(filter(dat, !! attr == 1), (!!Y)), pull(filter(dat, !! attr == 1), (!!treat)), "Counterfactual", attr_name)
  attr_other_count <- compute_calibration_dr(num_bins, pull(filter(dat, !! attr == 0), (!!mu_true)), pull(filter(dat, !! attr == 0), (!!count_preds)), pull(filter(dat, !! attr == 0), (!!pi)), pull(filter(dat, !! attr == 0), (!!Y)), pull(filter(dat, !! attr == 0), (!!treat)), "Counterfactual", attr_name_other)
  
  all_obs <- compute_calibration_dr(num_bins, pull(dat, (!!mu_true)), pull(dat, (!!obs_preds)), pull(dat, (!!pi)), pull(dat, (!!Y)), pull(dat, (!!treat)), "Observational", "All")
  attr_obs <- compute_calibration_dr(num_bins, pull(filter(dat, !! attr == 1), (!!mu_true)), pull(filter(dat, !! attr == 1), (!!obs_preds)), pull(filter(dat, !! attr == 1), (!!pi)), pull(filter(dat, !! attr == 1), (!!Y)), pull(filter(dat, !! attr == 1), (!!treat)), "Observational", attr_name)
  attr_other_obs <- compute_calibration_dr(num_bins, pull(filter(dat, !! attr == 0), (!!mu_true)), pull(filter(dat, !! attr == 0), (!!obs_preds)), pull(filter(dat, !! attr == 0), (!!pi)), pull(filter(dat, !! attr == 0), (!!Y)), pull(filter(dat, !! attr == 0), (!!treat)), "Observational", attr_name_other)
  return(rbind(all_count, attr_count, attr_other_count, all_obs, attr_obs, attr_other_obs))
}


# Function to compute threshold based on proportion classified
determine_thresh <- function(base, risks) {
  return(unname(quantile(risks, probs = (1 - base))))
}
