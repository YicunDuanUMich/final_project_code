rm(list = ls())

library(data.table)
library(ggplot2)
library(envalysis)

source("./data_generation.r")

generate_decay_samples <- function(decay, sample_size, num_variables) {
  result_list <- list()
  pre_variable <- NULL
  for (i in 1:num_variables) {
    if (is.null(pre_variable)) {
      pre_variable <- runif(sample_size)
      result_list[[i]] <- pre_variable
      next
    } 
    
    cur_variable <- pre_variable + rnorm(sample_size, mean = 0, sd = decay)
    result_list[[i]] <- cur_variable
    pre_variable <- cur_variable
  }
  
  result_matrix <- do.call(cbind, result_list)
  
  # standardize
  result_matrix <- apply(result_matrix, 2, function(x) {
    return((x - mean(x)) / sd(x))
  })
  
  return(result_matrix)
}

generate_non_linear_decay_samples <- function(decay_for_low_values, decay_for_high_values,
                                              sample_size, num_variables) {
  result_list <- list()
  pre_variable <- NULL
  for (i in 1:num_variables) {
    if (is.null(pre_variable)) {
      pre_variable <- runif(sample_size)
      result_list[[i]] <- pre_variable
      next
    } 
    cur_variable <- pre_variable
    pre_variable_low_values <- pre_variable[pre_variable < 0.5]
    pre_variable_high_values <- pre_variable[pre_variable >= 0.5]
    cur_variable[pre_variable < 0.5] <- pre_variable_low_values + 
                                        rnorm(length(pre_variable_low_values), 
                                              mean = 0, sd = decay_for_low_values)
    cur_variable[pre_variable >= 0.5] <- pre_variable_high_values + 
                                          rnorm(length(pre_variable_high_values),
                                                mean = 0, sd = decay_for_high_values)
    result_list[[i]] <- cur_variable
    pre_variable <- cur_variable
  }
  
  result_matrix <- do.call(cbind, result_list)
  
  # standardize
  result_matrix <- apply(result_matrix, 2, function(x) {
    return((x - mean(x)) / sd(x))
  })
  
  return(result_matrix)
}

generate_temporal_col_samples <- function(decay, sample_size, num_variables,
                                          second_stage_process) {
  if (sample_size %% 2 != 0) {
    stop("sample_size should be an even number")
  }
  
  first_stage_samples <- generate_decay_samples(decay, sample_size / 2, num_variables)
  
  # same decay
  if (second_stage_process == 1) {
    second_stage_samples <- generate_decay_samples(decay, sample_size / 2, num_variables)
  
    # low decay
  } else if (second_stage_process == 2) {
    second_stage_samples <- generate_decay_samples(decay / 2, sample_size / 2, num_variables)
    
    # high decay
  } else if (second_stage_process == 3) {
    second_stage_samples <- generate_decay_samples(decay * 2, sample_size / 2, num_variables)
    
    # non linear decay
  } else if (second_stage_process == 4) {
    second_stage_samples <- generate_non_linear_decay_samples(decay_for_low_values = decay,
                                                              decay_for_high_values = decay * 2,
                                                              sample_size = sample_size / 2,
                                                              num_variables = num_variables)
    
    # all independent
  } else if (second_stage_process == 5) {
    second_stage_samples <- generate_mvnorm_samples(cov_matrix = diag(1, num_variables, num_variables), 
                                                    sample_size = sample_size / 2)
  } else {
    stop("Invalid second_stage_process.")
  }
  
  result_samples <- rbind(first_stage_samples, second_stage_samples)
  
  return(result_samples)
}

decay_level_factor <- seq(0.1, 0.9, by = 0.1)
R2_factor <- 0.5
sample_size_factor <- 100
second_stage_process_factor_index <- 1:5
second_stage_process_factor <- c("same decay", "low decay", "high decay",
                                 "non linear decay", "all independent")

# please note we use the second_stage_process_factor_index rather than second_stage_process_factor.
# this is because transferring character vector to factor is error-prone
ds_factors <- list(decay_level_factor = decay_level_factor,
                   R2_factor = R2_factor,
                   sample_size_factor = sample_size_factor,
                   second_stage_process_factor = second_stage_process_factor_index)

cell_generate_function_for_test_4_1 <- function(cell_factors, reps) {
  decay_level <- cell_factors$decay_level_factor
  R2 <- cell_factors$R2_factor
  sample_size <- cell_factors$sample_size_factor
  second_stage_process <- cell_factors$second_stage_process_factor
  
  results <- list(estimated_coefs = list(),
                  estimated_coef_sds = list(),
                  condition_numbers = list(),
                  X = list())
  for (i in 1:reps) {
    cur_samples <- generate_temporal_col_samples(decay = decay_level,
                                                 sample_size = sample_size,
                                                 num_variables = 3,
                                                 second_stage_process = second_stage_process)
    
    cur_Y <- generate_y(cur_samples, 
                        expected_beta = c(1, 1, 1), 
                        R2 = R2,
                        has_intercept = FALSE)
    
    cur_model_summary <- summary(lm(cur_Y ~ cur_samples + 0))
    
    results[["estimated_coefs"]][[i]] <- cur_model_summary$coefficients[, 1]
    results[["estimated_coef_sds"]][[i]] <- cur_model_summary$coefficients[, 2]
    results[["condition_numbers"]][[i]] <- calculate_condition_number(cur_samples)
    
    results[["X"]][[i]] <- c(decay_level_factor = decay_level,
                             R2_factor = R2, 
                             sample_size_factor = sample_size, 
                             second_stage_process_factor = second_stage_process)
  }
  
  X <- do.call(rbind, results[["X"]])
  estimated_coefs <- do.call(rbind, results[["estimated_coefs"]])
  estimated_coef_sds <- do.call(rbind, results[["estimated_coef_sds"]])
  condition_numbers <- do.call(rbind, results[["condition_numbers"]])
  
  est_dt <- data.table(estimated_coefs = estimated_coefs,
                       estimated_coef_sds = estimated_coef_sds,
                       condition_numbers = condition_numbers,
                       X = X)
  
  names(est_dt) <- c(paste0("coef_", paste0("beta", 1:ncol(estimated_coefs))),
                     paste0("coef_sd_", paste0("beta", 1:ncol(estimated_coef_sds))), 
                     "condition_number",
                     paste0("X_", c("decay_level", 
                                    "R2", "sample_size", 
                                    "second_stage_process")))
  
  est_dt <- est_dt[,
                   .(est_sd_beta1 = sqrt(sum((coef_beta1 - mean(coef_beta1))^2) / (.N - 1)),
                     est_sd_beta2 = sqrt(sum((coef_beta2 - mean(coef_beta2))^2) / (.N - 1)),
                     est_sd_beta3 = sqrt(sum((coef_beta3 - mean(coef_beta3))^2) / (.N - 1)),
                     mean_condition_number = mean(condition_number)),
                   by = .(X_decay_level, X_R2, X_sample_size, X_second_stage_process)]
  
  return(est_dt)
}

est_result <- generate_dataset_general_version(ds_factors, 
                                               cell_generate_function_for_test_4_1, 
                                               reps = 1000, use_parallel = TRUE,
                                               help_functions = c("generate_decay_samples",
                                                                  "generate_non_linear_decay_samples",
                                                                  "generate_temporal_col_samples",
                                                                  "generate_mvnorm_samples",
                                                                  "generate_y",
                                                                  "calculate_condition_number"))

est_sd_ci_dt <- est_result$result_dt
est_sd_ci_dt$X_decay_level <- as.factor(est_sd_ci_dt$X_decay_level)
est_sd_ci_dt$X_second_stage_process <- as.factor(est_sd_ci_dt$X_second_stage_process)
levels(est_sd_ci_dt$X_second_stage_process) <- second_stage_process_factor


# plot est sd of beta1
second_y_axis_ratio <- max(est_sd_ci_dt$est_sd_beta1) / max(est_sd_ci_dt$mean_condition_number)
ggplot(est_sd_ci_dt, aes(x = X_decay_level)) +
  geom_col(mapping = aes(y = mean_condition_number * second_y_axis_ratio, fill = X_second_stage_process),
           position = "dodge") +
  geom_point(mapping = aes(y = est_sd_beta1, color = X_second_stage_process)) +
  geom_line(mapping = aes(y = est_sd_beta1, color = X_second_stage_process, group = X_second_stage_process)) +
  scale_y_continuous(name = "Estimated Beta1 SD",
                     sec.axis = sec_axis(~ . / second_y_axis_ratio,
                                         name = "Condition Number")) +
  scale_x_discrete(name = "Decay Level") +
  scale_color_discrete(name = "Second Stage Process") +
  scale_fill_discrete(name = "Second Stage Process") +
  theme_publish()

# plot est sd of beta2
second_y_axis_ratio <- max(est_sd_ci_dt$est_sd_beta2) / max(est_sd_ci_dt$mean_condition_number)
ggplot(est_sd_ci_dt, aes(x = X_decay_level)) +
  geom_col(mapping = aes(y = mean_condition_number * second_y_axis_ratio, fill = X_second_stage_process),
           position = "dodge") +
  geom_point(mapping = aes(y = est_sd_beta2, color = X_second_stage_process)) +
  geom_line(mapping = aes(y = est_sd_beta2, color = X_second_stage_process, group = X_second_stage_process)) +
  scale_y_continuous(name = "Estimated Beta2 SD",
                     sec.axis = sec_axis(~ . / second_y_axis_ratio,
                                         name = "Condition Number")) +
  scale_x_discrete(name = "Decay Level") +
  scale_color_discrete(name = "Second Stage Process") +
  scale_fill_discrete(name = "Second Stage Process") +
  theme_publish()


# plot est sd of beta3
second_y_axis_ratio <- max(est_sd_ci_dt$est_sd_beta3) / max(est_sd_ci_dt$mean_condition_number)
ggplot(est_sd_ci_dt, aes(x = X_decay_level)) +
  geom_col(mapping = aes(y = mean_condition_number * second_y_axis_ratio, fill = X_second_stage_process),
           position = "dodge") +
  geom_point(mapping = aes(y = est_sd_beta3, color = X_second_stage_process)) +
  geom_line(mapping = aes(y = est_sd_beta3, color = X_second_stage_process, group = X_second_stage_process)) +
  scale_y_continuous(name = "Estimated Beta3 SD",
                     sec.axis = sec_axis(~ . / second_y_axis_ratio,
                                         name = "Condition Number")) +
  scale_x_discrete(name = "Decay Level") +
  scale_color_discrete(name = "Second Stage Process") +
  scale_fill_discrete(name = "Second Stage Process") +
  theme_publish()





