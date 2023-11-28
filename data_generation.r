library(MASS)
library(parallel)

generate_mvnorm_samples <- function(cov_matrix, sample_size) {
  return(mvrnorm(n = sample_size, 
                 mu = rep(0, ncol(cov_matrix)), 
                 Sigma = cov_matrix))
}

calculate_condition_number <- function(X) {
  e <- eigen(t(X) %*% X)
  return(sqrt(max(e$val) / min(e$val)))
}

generate_y <- function(X, expected_beta, R2, has_intercept = TRUE) {
  if (has_intercept) {
    cov_matrix_X <- cov(cbind(1, X))
    s2 <- t(expected_beta) %*% 
      cov_matrix_X %*% 
      expected_beta
    f2 <- 1 / R2 - 1
    
    X_w_intercept <- cbind(1, X)
    Y <- X_w_intercept %*% expected_beta + rnorm(nrow(X_w_intercept), 
                                                 mean = 0, sd = sqrt(s2 * f2))
  } else {
    cov_matrix_X <- cov(X)
    s2 <- t(expected_beta) %*% 
      cov_matrix_X %*% 
      expected_beta
    f2 <- 1 / R2 - 1
    
    Y <- X %*% expected_beta + rnorm(nrow(X), mean = 0, sd = sqrt(s2 * f2))
  }
  
  return(Y)
}

generate_y_according_to_sd <- function(X, expected_beta, sd, has_intercept = TRUE) {
  if (has_intercept) {
    X_w_intercept <- cbind(1, X)
    Y <- X_w_intercept %*% expected_beta + rnorm(nrow(X_w_intercept), 
                                                 mean = 0, sd = sd)
  } else {
    Y <- X %*% expected_beta + rnorm(nrow(X), mean = 0, sd = sd)
  }
  
  return(Y)
}


generate_dataset <- function(col_levels, R2s, sample_sizes, models,
                                   reps = 100, has_intercept = TRUE) {
  col_level_index <- 1
  R2_index <- 1
  sample_size_index <- 1
  model_index <- 1
  total_index <- 1
  
  results <- list(estimated_coefs = list(),
                  estimated_coef_sds = list(),
                  X = list())
  for (col_level in col_levels) {
    for (R2 in R2s) {
      for (sample_size in sample_sizes) {
        for (model in models) {
          for (i in 1:reps) {
            cur_samples <- generate_mvnorm_samples(col_level, 
                                                   sample_size = sample_size)
            
            cur_Y <- generate_y(cur_samples, 
                                expected_beta = model, 
                                R2 = R2, has_intercept = has_intercept)
            
            if (has_intercept) {
              cur_model_summary <- summary(lm(cur_Y ~ cur_samples))
            } else {
              cur_model_summary <- summary(lm(cur_Y ~ cur_samples + 0))
            }
            
            results[["estimated_coefs"]][[total_index]] <- cur_model_summary$coefficients[, 1]
            results[["estimated_coef_sds"]][[total_index]] <- cur_model_summary$coefficients[, 2]
            results[["X"]][[total_index]] <- c(col_level_factor = col_level_index, 
                                               R2_factor = R2_index, 
                                               sample_size_factor = sample_size_index, 
                                               model_factor = model_index)
            
            total_index <- total_index + 1
          }
          model_index <- model_index + 1
        }
        model_index <- 1
        sample_size_index <- sample_size_index + 1
      }
      sample_size_index <- 1
      R2_index <- R2_index + 1
    }
    R2_index <- 1
    col_level_index <- col_level_index + 1
  }
  
  X <- do.call(rbind, results[["X"]])
  estimated_coefs <- do.call(rbind, results[["estimated_coefs"]])
  estimated_coef_sds <- do.call(rbind, results[["estimated_coef_sds"]])
  
  est_df <- data.frame(estimated_coefs = estimated_coefs,
                       estimated_coef_sds = estimated_coef_sds,
                       X = X)
  if (has_intercept) {
    names(est_df) <- c(paste0("coef_", paste0("beta", 0:(ncol(estimated_coefs) - 1))),
                       paste0("coef_sd_", paste0("beta", 0:(ncol(estimated_coef_sds) - 1))), 
                       paste0("X_", c("col_level", "R2", "sample_size", "model")))
  } else {
    names(est_df) <- c(paste0("coef_", paste0("beta", 1:ncol(estimated_coefs))),
                       paste0("coef_sd_", paste0("beta", 1:ncol(estimated_coef_sds))), 
                       paste0("X_", c("col_level", "R2", "sample_size", "model")))
  }
  
  
  return(list(est_df = est_df,
              reps = reps))
}

generate_dataset_general_version <- function(factors, cell_generate_function,
                                              reps = 100, use_parallel = TRUE,
                                             help_functions = c("generate_mvnorm_samples", 
                                                                "generate_y")) {
  # factors should be a list whose elements are numeric vectors
  check_factors <- sapply(factors, is.numeric)
  if (!all(check_factors)) {
    stop("`factors` should be a list whose elements are numeric vectors.")
  }
  
  factors_df <- expand.grid(ds_factors)
  factors_list <- split(factors_df, 1:nrow(factors_df))
  
  if (use_parallel) {
    cl <- makeCluster(detectCores())
    
    clusterEvalQ(cl, library(MASS))
    clusterEvalQ(cl, library(data.table))
    clusterEvalQ(cl, library(envalysis))
    for (help_function in help_functions) {
      clusterExport(cl, help_function)
    }
    # the cell_generate_function should generate `reps` estimation samples 
    # and reduce the samples to an one row data.table
    result_list <- parLapply(cl, factors_list, cell_generate_function, reps = reps)
    
    stopCluster(cl)
  } else {
    result_list <- lapply(factors_list, cell_generate_function, reps = reps)
  }
  
  result_dt <- do.call(rbind, result_list)
 
  return(list(result_dt = result_dt,
              reps = reps))
}
