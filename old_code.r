# old code for calculating acc
# acc_df <- est_df
# acc_df[, 1:10] <- 0
# names(acc_df) <- c(paste0("acc_coef_", paste0("beta", 0:4)),
#                    paste0("acc_sd_", paste0("beta", 0:4)),
#                    paste0("X_", c("col_level", "R2", "sample_size", "model")))
# 
# for (i in seq(1,
#               length(col_levels) *
#               length(R2s) *
#               length(sample_sizes) *
#               length(models) * reps - 1,
#               by = reps)) {
#   sub_coefs <- as.matrix(est_df[i:(i + reps - 1), 1:5])
#   sub_sds <- as.matrix(est_df[i:(i + reps - 1), 6:10])
#   sub_coef_means <- colMeans(sub_coefs)
#   est_true_sub_sds <- sqrt(colSums(
#     (t(t(sub_coefs) - sub_coef_means))^2
#   ) / (reps - 1))
#   acc_df[i:(i + reps - 1), 6:10] <- abs(t(t(sub_sds) - est_true_sub_sds))
# }
# 
# model_index <- 1
# for (model in models) {
#   sub_coefs <- as.matrix(est_df[est_df$X_model == model_index, 1:5])
#   acc_df[est_df$X_model == model_index, 1:5] <- abs(t(t(sub_coefs) - model))
#   model_index <- model_index + 1
# }