rm(list = ls())

library(MASS)
library(data.table)
library(ggplot2)
library(envalysis)

f_0 <- function(X1, X2) {
  return(rnorm(length(X1)))
}

f_1 <- function(X1, X2) {
  result <- X1 / (X2^2)
  result[which(is.nan(result))] <- 0
  result[which(is.infinite(result))] <- 0
  return(result)
}

f_2 <- function(X1, X2) {
  return(X1 * X2)
}

# f_3 <- function(X1, X2) {
#   return(X1 + X2)
# }

f_4 <- function(X1, X2) {
  return(X1^2 + X2^2)
}

f_5 <- function(X1, X2) {
  return(X1^2 + X2^3)
}

f_6 <- function(X1, X2) {
  return(X1^3 + X2^3)
}

f_7 <- function(X1, X2) {
  return(exp(X1 + X2))
}

f_8 <- function(X1, X2) {
  return(exp(X1))
}


# X3 = f(X1, X2) + noise
# without intercept
generate_samples <- function(noise_sd, sample_size, f_X1_X2) {
  # X1 and X2
  X1_X2 <- mvrnorm(n = sample_size, 
                        mu = rep(0, 2), 
                        Sigma = diag(1, 2, 2))
  X3 <- f_X1_X2(X1_X2[, 1], X1_X2[, 2]) + rnorm(sample_size, 
                                                mean = 0, sd = noise_sd)
  # standardize
  X3 <- (X3 - mean(X3)) / sd(X3)
  X <- cbind(X1_X2, X3)
  return(X)
}

calculate_condition_number <- function(X) {
  e <- eigen(t(X) %*% X)
  return(sqrt(max(e$val) / min(e$val)))
}

generate_y <- function(X, expected_beta, R2) {
  cov_matrix_X <- cov(X)
  s2 <- t(expected_beta) %*% 
    cov_matrix_X %*% 
    expected_beta
  f2 <- 1 / R2 - 1

  Y <- X %*% expected_beta + rnorm(nrow(X), mean = 0, sd = sqrt(s2 * f2))
  return(Y)
}

noise_level_factor <- seq(0.0, 3.0, by = 0.5)
R2_factor <- 0.5
sample_size_factor <- 100
model_factor <- list(c(1, 1, 0))
f_factor <- list(f_0, f_1, f_2, f_4, f_5, f_6, f_7, f_8)

# without intercept
generate_dataset_for_test_3_2 <- function(noise_levels, R2s, sample_sizes, models, fs,
                             reps = 100) {
  noise_level_index <- 1
  R2_index <- 1
  sample_size_index <- 1
  model_index <- 1
  f_index <- 1
  total_index <- 1
  
  results <- list(estimated_coefs = list(),
                  estimated_coef_sds = list(),
                  condition_number = list(),
                  X = list())
  for (noise_level in noise_levels) {
    for (R2 in R2s) {
      for (sample_size in sample_sizes) {
        for (model in models) {
          for (f in fs) {
            for (i in 1:reps) {
              cur_samples <- generate_samples(noise_sd = noise_level, 
                                              sample_size = sample_size,
                                              f_X1_X2 = f)
              
              cur_Y <- generate_y(cur_samples, 
                                  expected_beta = model, 
                                  R2 = R2)
              
              cur_model_summary <- summary(lm(cur_Y ~ cur_samples + 0))
              
              results[["estimated_coefs"]][[total_index]] <- cur_model_summary$coefficients[, 1]
              results[["estimated_coef_sds"]][[total_index]] <- cur_model_summary$coefficients[, 2]
              results[["condition_number"]][[total_index]] <- calculate_condition_number(cur_samples)
              results[["X"]][[total_index]] <- c(noise_level_factor = noise_level_index, 
                                                 R2_factor = R2_index, 
                                                 sample_size_factor = sample_size_index, 
                                                 model_factor = model_index,
                                                 f_factor = f_index)
              
              total_index <- total_index + 1
            }
            f_index <- f_index + 1
          }
          f_index <- 1
          model_index <- model_index + 1
        }
        model_index <- 1
        sample_size_index <- sample_size_index + 1
      }
      sample_size_index <- 1
      R2_index <- R2_index + 1
    }
    R2_index <- 1
    noise_level_index <- noise_level_index + 1
  }
  
  X <- do.call(rbind, results[["X"]])
  estimated_coefs <- do.call(rbind, results[["estimated_coefs"]])
  estimated_coef_sds <- do.call(rbind, results[["estimated_coef_sds"]])
  condition_numbers <- do.call(rbind, results[["condition_number"]])
  
  est_df <- data.frame(estimated_coefs = estimated_coefs,
                       estimated_coefsds = estimated_coef_sds,
                       condition_numbers = condition_numbers,
                       X = X)
  names(est_df) <- c(paste0("coef_", paste0("beta", 1:3)),
                     paste0("coef_sd_", paste0("beta", 1:3)), 
                     "condition_number",
                     paste0("X_", c("noise_level", "R2", "sample_size", "model", "f")))
  
  return(list(est_df = est_df,
              reps = reps))
}


ds <- generate_dataset_for_test_3_2(noise_levels = noise_level_factor,
                                     R2s = R2_factor,
                                     sample_sizes = sample_size_factor,
                                     models = model_factor,
                                     fs = f_factor,
                                     reps = 1000)

calculate_sd_and_condition_number <- function(est_df) {
  est_dt <- as.data.table(est_df)
  
  est_dt <- est_dt[,
                   .(est_sd_beta1 = sqrt(sum((coef_beta1 - mean(coef_beta1))^2) / (.N - 1)),
                     est_sd_beta2 = sqrt(sum((coef_beta2 - mean(coef_beta2))^2) / (.N - 1)),
                     est_sd_beta3 = sqrt(sum((coef_beta3 - mean(coef_beta3))^2) / (.N - 1)),
                     mean_condition_number = mean(condition_number)),
                   by = .(X_noise_level, X_R2, X_sample_size, X_model, X_f)]
  
  return(est_dt)
}

est_df <- ds$est_df
est_sd_ci_dt <- calculate_sd_and_condition_number(est_df)
est_sd_ci_dt$X_noise_level <- as.factor(est_sd_ci_dt$X_noise_level)
levels(est_sd_ci_dt$X_noise_level) <- as.character(noise_level_factor)
est_sd_ci_dt$X_f <- as.factor(est_sd_ci_dt$X_f)
levels(est_sd_ci_dt$X_f) <- c("X3", "X1 / (X2^2)", "X1 * X2", "X1^2 + X2^2",
                              "X1^2 + X2^3", "X1^3 + X2^3", "exp(X1 + X2)",
                              "exp(X1)")

# est_sd_ci_dt$log_est_sd_beta1 <- log(est_sd_ci_dt$est_sd_beta1 + 1)
# est_sd_ci_dt$log_est_sd_beta2 <- log(est_sd_ci_dt$est_sd_beta2 + 1)
# est_sd_ci_dt$log_est_sd_beta3 <- log(est_sd_ci_dt$est_sd_beta3 + 1)
# est_sd_ci_dt$log_mean_condition_number <- log(est_sd_ci_dt$mean_condition_number + 1)

# plot est sd of beta1
# second_y_axis_ratio <- max(est_sd_ci_dt$log_est_sd_beta1) / max(est_sd_ci_dt$log_mean_condition_number)
# ggplot(est_sd_ci_dt, aes(x = X_noise_level)) +
#   geom_col(mapping = aes(y = log_mean_condition_number * second_y_axis_ratio, fill = X_f), 
#            position = "dodge") +
#   geom_point(mapping = aes(y = log_est_sd_beta1, color = X_f)) +
#   geom_line(mapping = aes(y = log_est_sd_beta1, color = X_f, group = X_f)) +
#   scale_y_continuous(name = "log(Estimated Beta1 SD + 1)",
#                      sec.axis = sec_axis(~ . / second_y_axis_ratio,
#                                          name = "log(Condition Number + 1)")) +
#   scale_x_discrete(name = "Noise Level") +
#   scale_color_discrete(name = "X3 =") +
#   scale_fill_discrete(name = "X3 =") +
#   theme_publish()

# plot est sd of beta1
second_y_axis_ratio <- max(est_sd_ci_dt$est_sd_beta1) / max(est_sd_ci_dt$mean_condition_number)
ggplot(est_sd_ci_dt, aes(x = X_noise_level)) +
  geom_col(mapping = aes(y = mean_condition_number * second_y_axis_ratio / 3, fill = X_f),
           position = "dodge") +
  geom_point(mapping = aes(y = est_sd_beta1, color = X_f)) +
  geom_line(mapping = aes(y = est_sd_beta1, color = X_f, group = X_f)) +
  scale_y_continuous(name = expression(paste("Estimated ", hat(beta)[1], " SD")),
                     sec.axis = sec_axis(~ . / second_y_axis_ratio * 3,
                                         name = expression(paste("Condition Number of ", bold(X))))) +
  scale_x_discrete(name = "Noise Level") +
  scale_color_discrete(name = "X3 =") +
  scale_fill_discrete(name = "X3 =") +
  theme_publish(base_size = 15, base_linewidth = 1)


# plot est sd of beta2
second_y_axis_ratio <- max(est_sd_ci_dt$est_sd_beta2) / max(est_sd_ci_dt$mean_condition_number)
ggplot(est_sd_ci_dt, aes(x = X_noise_level)) +
  geom_col(mapping = aes(y = mean_condition_number * second_y_axis_ratio / 3, fill = X_f),
           position = "dodge") +
  geom_point(mapping = aes(y = est_sd_beta2, color = X_f)) +
  geom_line(mapping = aes(y = est_sd_beta2, color = X_f, group = X_f)) +
  scale_y_continuous(name = expression(paste("Estimated ", hat(beta)[2], " SD")),
                     sec.axis = sec_axis(~ . / second_y_axis_ratio * 3,
                                         name = expression(paste("Condition Number of ", bold(X))))) +
  scale_x_discrete(name = "Noise Level") +
  scale_color_discrete(name = "X3 =") +
  scale_fill_discrete(name = "X3 =") +
  theme_publish(base_size = 15, base_linewidth = 1)

# plot est sd of beta3
second_y_axis_ratio <- max(est_sd_ci_dt$est_sd_beta2) / max(est_sd_ci_dt$mean_condition_number)
ggplot(est_sd_ci_dt, aes(x = X_noise_level)) +
  geom_col(mapping = aes(y = mean_condition_number * second_y_axis_ratio / 3, fill = X_f),
           position = "dodge") +
  geom_point(mapping = aes(y = est_sd_beta3, color = X_f)) +
  geom_line(mapping = aes(y = est_sd_beta3, color = X_f, group = X_f)) +
  scale_y_continuous(name = expression(paste("Estimated ", hat(beta)[3], " SD")),
                     sec.axis = sec_axis(~ . / second_y_axis_ratio * 3,
                                         name = expression(paste("Condition Number of ", bold(X))))) +
  scale_x_discrete(name = "Noise Level") +
  scale_color_discrete(name = "X3 =") +
  scale_fill_discrete(name = "X3 =") +
  theme_publish(base_size = 15, base_linewidth = 1)


# labels = c(expression(X[3]), expression(X[1] %/% X[2]^2),
#            expression(X[1] %*% X[2]), expression(X[1]^2 + X[2]^2),
#            expression(X[1]^2 + X[2]^3)), expression(X[1]^3 + X[2]^3),
#            expression(e^{X[1] + X[2]}), expression(e^{X[1]})


