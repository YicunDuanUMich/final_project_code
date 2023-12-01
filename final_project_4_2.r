rm(list = ls())

library(data.table)
library(ggplot2)
library(envalysis)

source("./data_generation.r")

generate_decay_samples <- function(decay, sample_size, num_variables, group_size) {
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
    if (i %% group_size == 0) {
      pre_variable <- NULL
    } else {
      pre_variable <- cur_variable
    }
  }
  
  result_matrix <- do.call(cbind, result_list)
  
  # standardize
  result_matrix <- apply(result_matrix, 2, function(x) {
    return((x - mean(x)) / sd(x))
  })
  
  return(result_matrix)
}

generate_non_linear_decay_samples <- function(decay_for_low_values, decay_for_high_values,
                                              sample_size, num_variables, group_size) {
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
    
    if (i %% group_size == 0) {
      pre_variable <- NULL
    } else {
      pre_variable <- cur_variable
    }
  }
  
  result_matrix <- do.call(cbind, result_list)
  
  # standardize
  result_matrix <- apply(result_matrix, 2, function(x) {
    return((x - mean(x)) / sd(x))
  })
  
  return(result_matrix)
}

decay_level_factor <- c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
SD_factor <- 0.5
sample_size_factor <- 1000
test_set_type_factor_index <- 1:5
test_set_type_factor <- c("same decay", "low decay", "high decay",
                      "non linear decay", "all independent")

ds_factors <- list(decay_level_factor = decay_level_factor,
                   SD_factor = SD_factor,
                   sample_size_factor = sample_size_factor,
                   test_set_type_factor = test_set_type_factor_index)

cell_generate_function_for_test_4_2 <- function(cell_factors, reps) {
  decay_level <- cell_factors$decay_level_factor
  SD <- cell_factors$SD_factor
  sample_size <- cell_factors$sample_size_factor
  test_set_type <- cell_factors$test_set_type_factor
  
  results <- list(RMSE = list(),
                  X = list())
  expected_beta <-  c(25,
                      0, 0, 1, 0, 0,
                      rep(0, 5),
                      rep(0, 5),
                      rep(0, 5),
                      0)
  for (i in 1:reps) {
    train_samples <- generate_decay_samples(decay = decay_level,
                                          sample_size = sample_size,
                                          num_variables = 21,
                                          group_size = 5)
    colnames(train_samples) <- c(paste0("A", 1:5),
                                 paste0("B", 1:5),
                                 paste0("C", 1:5),
                                 paste0("D", 1:5),
                                 "E")
    train_set_Y <- generate_y_according_to_sd(train_samples, 
                                          expected_beta = expected_beta, 
                                          sd = SD,
                                          has_intercept = TRUE)
    train_set <- data.frame(Y = train_set_Y, train_samples)
    cur_model <- lm(Y ~ ., data = train_set)
    training_RMSE <- rmse(train_set_Y, cur_model$fitted.values)
    
    # same decay
    if (test_set_type == 1) {
      test_set_samples <- generate_decay_samples(decay_level, sample_size, 
                                                 num_variables = 21, group_size = 5)
      
      # low decay
    } else if (test_set_type == 2) {
      test_set_samples <- generate_decay_samples(decay_level / 2, sample_size, 
                                                 num_variables = 21, group_size = 5)
      
      # high decay
    } else if (test_set_type == 3) {
      test_set_samples <- generate_decay_samples(decay_level * 2, sample_size, 
                                                 num_variables = 21, group_size = 5)
      
      # non linear decay
    } else if (test_set_type == 4) {
      test_set_samples <- generate_non_linear_decay_samples(decay_for_low_values = 0,
                                                            decay_for_high_values = decay_level * 20,
                                                            sample_size = sample_size,
                                                            num_variables = 21,
                                                            group_size = 5)
      
      # all independent
    } else if (test_set_type == 5) {
      test_set_samples <- generate_mvnorm_samples(cov_matrix = diag(1, 21, 21),
                                                      sample_size = sample_size)
      
      # test_set_samples <- cbind(runif(sample_size), runif(sample_size), runif(sample_size))
      # test_set_samples <- apply(test_set_samples, 2, function(x) {
      #   return((x - mean(x)) / sd(x))
      # })
    } else {
      stop("Invalid test_set_type.")
    }
    
    colnames(test_set_samples) <- c(paste0("A", 1:5),
                                 paste0("B", 1:5),
                                 paste0("C", 1:5),
                                 paste0("D", 1:5),
                                 "E")
    test_set_Y <- generate_y_according_to_sd(test_set_samples, 
                                             expected_beta = expected_beta, 
                                             sd = SD,
                                             has_intercept = TRUE)
    predicted_test_set_Y <- predict(cur_model, data.frame(test_set_samples))
    testing_RMSE <- rmse(test_set_Y, predicted_test_set_Y)
    
    results[["RMSE"]][[i]] <- c(training_RMSE = training_RMSE,
                                testing_RMSE = testing_RMSE)
    results[["condition_numbers"]][[i]] <- calculate_condition_number(train_samples)
    results[["X"]][[i]] <- c(decay_level_factor = decay_level,
                             SD_factor = SD, 
                             sample_size_factor = sample_size, 
                             test_set_factor = test_set_type)
  }
  
  X <- do.call(rbind, results[["X"]])
  RMSE <- do.call(rbind, results[["RMSE"]])
  condition_numbers <- do.call(rbind, results[["condition_numbers"]])
  
  est_dt <- data.table(RMSE = RMSE,
                       condition_numbers = condition_numbers,
                       X = X)
  
  names(est_dt) <- c("training_RMSE", "testing_RMSE",
                     "condition_number",
                     paste0("X_", c("decay_level", 
                                    "SD", "sample_size", 
                                    "test_set_type")))
  
  # special: no reduce
  
  return(est_dt)
}

est_result <- generate_dataset_general_version(ds_factors, 
                                               cell_generate_function_for_test_4_2, 
                                               reps = 100, use_parallel = TRUE,
                                               help_functions = c("generate_decay_samples",
                                                                  "generate_non_linear_decay_samples",
                                                                  "generate_mvnorm_samples",
                                                                  "generate_y_according_to_sd",
                                                                  "calculate_condition_number"))

est_rmse_ci_dt <- est_result$result_dt
est_rmse_ci_dt$X_decay_level <- as.factor(est_rmse_ci_dt$X_decay_level)
est_rmse_ci_dt$X_test_set_type <- as.factor(est_rmse_ci_dt$X_test_set_type)
levels(est_rmse_ci_dt$X_test_set_type) <- paste0("Testing Data: ", test_set_type_factor)


ggplot(est_rmse_ci_dt, aes(x = X_decay_level, y = testing_RMSE)) +
  geom_boxplot() +
  scale_y_continuous(name = "Testing RMSE", trans = "log") +
  scale_x_discrete(name = "Training Data Decay Level") +
  facet_grid(col = vars(X_test_set_type)) +
  theme_publish()



