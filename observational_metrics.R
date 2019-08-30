# Compute fpr
compute_fpr_obs <- function(t, pred, label) {
  pred_label <- ifelse(pred >= t, 1, 0)
  neg <- 1-label
  return(sum(pred_label*neg)/sum(neg))
}

# Compute fnr
compute_fnr_obs <- function(t, pred, label) {
  pred_label <- ifelse(pred >= t, 1, 0)
  return(sum((1-pred_label)*label)/sum(label))
}


# Compute fpr_cost
compute_fpr_cost_obs <- function(pred, label) {
  return(mean(pred[label==0]))
}

# Compute fnr_cost
compute_fnr_cost_obs <- function(pred, label) {
  return(1-mean(pred[label==1]))
}

# Compute Recall
compute_recall_obs <- function(t, pred, label) {
  pred_label <- ifelse(pred >= t, 1, 0)
  return(sum(pred_label*label)/sum(label))
}

compute_precision_thresh <- function(t, pred, label) {
  pred_label <- ifelse(pred >= t, 1, 0)
  return(sum(pred_label*label)/sum(pred_label))
}

compute_pr_df <- function(dat.eval, y, pred_caus, pred_batch, group) {
  #t_vals <- seq(0.1, 1, 0.01)
  t_vals <- seq(0, 1, 0.01)
  r.caus <- sapply(t_vals, compute_recall_obs, pred = pull(dat.eval,!! pred_caus), label = pull(dat.eval, !! y))
  r.batch <- sapply(t_vals, compute_recall_obs, pred = pull(dat.eval, !! pred_batch), label = pull(dat.eval, !! y))
  p.caus <- sapply(t_vals, compute_precision_thresh, pred = pull(dat.eval, !! pred_caus), label = pull(dat.eval, !! y))
  p.batch <- sapply(t_vals, compute_precision_thresh, pred = pull(dat.eval, !! pred_batch), label = pull(dat.eval,!! y))
  obs.all.caus <- data.frame("Recall" = r.caus, "Precision" = p.caus, "Threshold"= t_vals, "Method" = "Counterfactual")
  obs.all.batch <- data.frame("Recall" = r.batch, "Precision" = p.batch,"Threshold"= t_vals, "Method" = "Observational")
  obs.all <- rbind(obs.all.caus, obs.all.batch)
  obs.all$Group <- group
  return(obs.all)
}

compute_roc_df <- function(dat.eval, y, pred_caus, pred_batch, group) {
  #t_vals <- seq(0.1, 1, 0.01)
  t_vals <- seq(0, 1, 0.01)
  r.caus <- sapply(t_vals, compute_recall_obs, pred = pull(dat.eval,!! pred_caus), label = pull(dat.eval, !! y))
  r.batch <- sapply(t_vals, compute_recall_obs, pred = pull(dat.eval, !! pred_batch), label = pull(dat.eval, !! y))
  fpr.caus <- sapply(t_vals, compute_fpr_obs, pred = pull(dat.eval, !! pred_caus), label = pull(dat.eval, !! y))
  fpr.batch <- sapply(t_vals, compute_fpr_obs, pred = pull(dat.eval, !! pred_batch), label = pull(dat.eval,!! y))
  obs.all.caus <- data.frame("Recall" = r.caus, "FPR" = fpr.caus, "Threshold"= t_vals, "Method" = "Counterfactual")
  obs.all.batch <- data.frame("Recall" = r.batch, "FPR" = fpr.batch,"Threshold"= t_vals, "Method" = "Observational")
  obs.all <- rbind(obs.all.caus, obs.all.batch)
  obs.all$Group <- group
  return(obs.all)
}

# Calibration based on deciles (default) or ventiles if num_bins = 20
compute_calibration_obs <- function(num_bins = 10, mu_hat, Y) {
  steps <- seq(1/num_bins, 1, 1/num_bins)
  threshs <- sapply(steps, function(x) {unname(quantile(mu_hat, probs=x))}) 
  calib <- data.frame(matrix(nrow=0, ncol=4))
  colnames(calib) <- c("Average score", "Obs.rate", "Low", "High")
  i <- 1
  while(i <= num_bins) {
    prev <- ifelse(i >1, threshs[i-1], 0)
    sub <- (mu_hat > prev & mu_hat <= threshs[i])
    calib <- rbind(calib, data.frame("Ventile risk score" = i, "Average score" = mean(mu_hat[sub]), "Rate"= mean(Y[sub]), "Low" = mean(Y[sub]) - 1.96*sqrt(var(Y[sub])/length(Y[sub])), "High" =  mean(Y[sub]) + 1.96*sqrt(var(Y[sub])/length(Y[sub]))))
    i= i+1
  }
  return(calib)
}

# num_bins is a scalar. dat.eval is a dataframe. pred_caus/batch is the counterfactual/observational predictions as a quosure. y is a quosure of the outcome
# group is a string
compute_calib_df <- function(num_bins, dat.eval, pred_caus, pred_batch, y, group) {
  calib.batch <- compute_calibration_obs(num_bins, pull(dat.eval, !! pred_batch), pull(dat.eval, !! y))
  calib.batch$Method <- "Observational"
  calib.caus <- compute_calibration_obs(num_bins, pull(dat.eval, !! pred_caus), pull(dat.eval, !! y))
  calib.caus$Method <- "Counterfactual"
  calib <- rbind(calib.caus, calib.batch)
  calib$Group <- group
  return(calib)
}