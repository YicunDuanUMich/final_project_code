rm(list = ls())

library(data.table)
library(ggplot2)
library(envalysis)
library(grDevices)
library(RColorBrewer)

source(file = "./data_generation.r")

col_level_x_factor <- seq(-0.45, 0.45, by = 0.05)
col_level_y_factor <- seq(-0.45, 0.45, by = 0.05)

R2_factor <- 0.25
sample_size_factor <- 30
model_factor <- list(c(1, 1, 0))

generate_dataset_for_test_3_1 <- function(col_levels_x, col_levels_y, 
                                       R2s, sample_sizes, models,
                                       reps = 100, has_intercept = TRUE) {
  col_level_x_index <- 1
  col_level_y_index <- 1
  R2_index <- 1
  sample_size_index <- 1
  model_index <- 1
  total_index <- 1
  
  results <- list(estimated_coefs = list(),
                  estimated_coef_sds = list(),
                  X = list())
  for (col_level_x in col_levels_x) {
    for (col_level_y in col_levels_y) {
      for (R2 in R2s) {
        for (sample_size in sample_sizes) {
          for (model in models) {
            for (i in 1:reps) {
              a <- col_level_x #x1x3
              b <- col_level_y #x2x3
              col_level <- matrix(c(1.0, 0.0, a,
                                    0.0, 1.0, b,
                                    a, b, 1.0), nrow = 3)
              
              cur_samples <- generate_mvnorm_samples(col_level, 
                                                     sample_size = sample_size)
              
              cur_Y <- generate_y(cur_samples, 
                                  expected_beta = model, 
                                  R2 = R2,
                                  has_intercept = has_intercept)
              
              if (has_intercept) {
                cur_model_summary <- summary(lm(cur_Y ~ cur_samples))
              } else {
                cur_model_summary <- summary(lm(cur_Y ~ cur_samples + 0))
              }
              
              results[["estimated_coefs"]][[total_index]] <- cur_model_summary$coefficients[, 1]
              results[["estimated_coef_sds"]][[total_index]] <- cur_model_summary$coefficients[, 2]
              results[["X"]][[total_index]] <- c(col_level_x_factor = col_level_x_index,
                                                 col_level_y_factor = col_level_y_index,
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
      col_level_y_index <- col_level_y_index + 1
    }
    col_level_y_index <- 1
    col_level_x_index <- col_level_x_index + 1
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
                       paste0("X_", c("col_level_x", "col_level_y", "R2", "sample_size", "model")))
  } else {
    names(est_df) <- c(paste0("coef_", paste0("beta", 1:ncol(estimated_coefs))),
                       paste0("coef_sd_", paste0("beta", 1:ncol(estimated_coef_sds))), 
                       paste0("X_", c("col_level_x", "col_level_y", "R2", "sample_size", "model")))
  }
  
  return(list(est_df = est_df,
              reps = reps))
}

ds <- generate_dataset_for_test_3_1(col_level_x_factor, col_level_y_factor,
                                 R2_factor, sample_size_factor, model_factor,
                                 reps = 5000, has_intercept = FALSE)

calculate_sd <- function(est_df) {
  est_dt <- as.data.table(est_df)
  
  est_dt <- est_dt[,
                   .(est_sd_beta1 = sqrt(sum((coef_beta1 - mean(coef_beta1))^2) / (.N - 1)),
                     est_sd_beta2 = sqrt(sum((coef_beta2 - mean(coef_beta2))^2) / (.N - 1)),
                     est_sd_beta3 = sqrt(sum((coef_beta3 - mean(coef_beta3))^2) / (.N - 1))),
                   by = .(X_col_level_x, X_col_level_y, X_R2, X_sample_size, X_model)]
  
  return(est_dt)
}

est_sd_dt <- calculate_sd(ds$est_df)
est_sd_dt$X_col_level_x <- as.factor(est_sd_dt$X_col_level_x)
levels(est_sd_dt$X_col_level_x) <- as.character(col_level_x_factor)
est_sd_dt$X_col_level_y <- as.factor(est_sd_dt$X_col_level_y)
levels(est_sd_dt$X_col_level_y) <- as.character(col_level_y_factor)

# ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y)) +
#   geom_point(mapping = aes(size = est_sd_beta1)) +
#   scale_size_continuous(name = "Estimated Beta1 SD") +
#   scale_x_discrete(name = "Corr(X1, X3)") + 
#   scale_y_discrete(name = "Corr(X2, X3)") +
#   theme_publish()
# 
# ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y)) +
#   geom_point(mapping = aes(size = est_sd_beta2)) +
#   scale_size_continuous(name = "Estimated Beta2 SD") +
#   scale_x_discrete(name = "Corr(X1, X3)") + 
#   scale_y_discrete(name = "Corr(X2, X3)") +
#   theme_publish()
# 
# ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y)) +
#   geom_point(mapping = aes(size = est_sd_beta3)) +
#   scale_size_continuous(name = "Estimated Beta3 SD") +
#   scale_x_discrete(name = "Corr(X1, X3)") + 
#   scale_y_discrete(name = "Corr(X2, X3)") +
#   theme_publish()


colormap <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(32)

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta1)) +
  geom_tile(aes(fill = est_sd_beta1)) +
  scale_x_discrete(name = "Corr(X1, X3)") + 
  scale_y_discrete(name = "Corr(X2, X3)") +
  scale_fill_gradientn(colours = colormap, name = "Estimated Beta1 SD") +
  theme_publish()

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta2)) +
  geom_tile(aes(fill = est_sd_beta2)) +
  scale_x_discrete(name = "Corr(X1, X3)") + 
  scale_y_discrete(name = "Corr(X2, X3)") +
  scale_fill_gradientn(colours = colormap, name = "Estimated Beta2 SD") +
  theme_publish()

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta3)) +
  geom_tile(aes(fill = est_sd_beta3)) +
  scale_x_discrete(name = "Corr(X1, X3)") + 
  scale_y_discrete(name = "Corr(X2, X3)") +
  scale_fill_gradientn(colours = colormap, name = "Estimated Beta3 SD") +
  theme_publish()

