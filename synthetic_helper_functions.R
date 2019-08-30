
create_ROC_df_synthetic <- function(t_vals, dat, name, score_name, adj_score_name, y0_name, method, adjustment) {
  FPR_arr <- sapply(t_vals, compute_fpr_obs, pred = dat[[score_name]], label = dat[[y0_name]])
  FNR_arr <- sapply(t_vals, compute_fnr_obs, pred = dat[[score_name]], label = dat[[y0_name]])
  df <- data.frame("Threshold" = t_vals, "FPR" = FPR_arr, "TPR" = 1 - FNR_arr, "Method" = method, "Group" = name)
  FPR_arr_adj <- sapply(t_vals, compute_fpr_obs, pred = dat[[adj_score_name]], label = dat[[y0_name]])
  FNR_arr_adj <- sapply(t_vals, compute_fnr_obs, pred = dat[[adj_score_name]], label = dat[[y0_name]])
  df_adj <- data.frame("Threshold" = t_vals, "FPR" = FPR_arr_adj, "TPR" = 1 - FNR_arr_adj, "Method" = adjustment, "Group" = name)
  return(rbind(df, df_adj))
}

process_post_synthetic <- function(dat, A_name = "A", B_name = "B", score_name = "s", adj_score_name = "eo_fair_pred", y0_name = "y0", y_name= "y", method = "Observational", adjustment = "Post Proc. (obs)") {
  dat_A <- filter(dat, A == 1)
  dat_B <- filter(dat, A == 0)
  
  t_vals <- seq(0.1, 0.95, 0.01)
  post_auc_A <- create_ROC_df_synthetic(t_vals, dat_A, A_name, score_name, adj_score_name, y0_name, method, adjustment)
  post_auc_B <- create_ROC_df_synthetic(t_vals, dat_B, B_name, score_name, adj_score_name, y0_name, method, adjustment)
  post_auc_c <- rbind(post_auc_A, post_auc_B)
  post_auc_c$Evaluation <- "Counterfactual"
  post_auc_A_obs <- create_ROC_df_synthetic(t_vals, dat_A, A_name, score_name, adj_score_name, y_name, method, adjustment)
  post_auc_B_obs <- create_ROC_df_synthetic(t_vals, dat_B, B_name, score_name, adj_score_name, y_name, method, adjustment)
  post_auc_c_obs <- rbind(post_auc_A_obs, post_auc_B_obs)
  post_auc_c_obs$Evaluation <- "Observational"
  
  return(rbind(post_auc_c, post_auc_c_obs))
}


# helper function to get the counterfactual and observational metrics
# returns a dataframe
compute_rates_synthetic <- function(t, pred_values, y0_values, y_values, group, method) {
  fpr <- compute_fpr_obs(t, pred = pred_values, label= y0_values)
  fnr <- compute_fnr_obs(t, pred = pred_values, label= y0_values)
  fpr_cost <- compute_fpr_cost_obs(pred = pred_values, label= y0_values)
  fnr_cost <- compute_fnr_cost_obs(pred = pred_values, label= y0_values)
  fpr_obs <- compute_fpr_obs(t, pred = pred_values, label= y_values)
  fnr_obs <- compute_fnr_obs(t, pred = pred_values, label = y_values)
  fnr_obs_cost <- compute_fnr_cost_obs(pred_values, y_values)
  fpr_obs_cost <- compute_fpr_cost_obs(pred_values, y_values)
  return(data.frame("Group" = group, "Model" = method, "cGFPR" = fpr_cost, "cGFNR" = fnr_cost, "oGFPR" = fpr_obs_cost, "oGFNR" = fnr_obs_cost))
}