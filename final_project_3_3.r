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

model_alpha <- seq(0, pi / 2, by = pi / 8)
model_beta <- seq(0, pi / 2, by = pi / 8)

ds_factors <- list(col_level_x_factor = col_level_x_factor,
                   col_level_y_factor = col_level_y_factor,
                   R2_factor = R2_factor,
                   sample_size_factor = sample_size_factor,
                   model_alpha = model_alpha,
                   model_beta = model_beta)

cell_generate_function_for_test_3_3 <- function(cell_factors, reps) {
  col_level_x <- cell_factors$col_level_x_factor
  col_level_y <- cell_factors$col_level_y_factor
  R2 <- cell_factors$R2_factor
  sample_size <- cell_factors$sample_size_factor
  alpha <- cell_factors$model_alpha
  beta <- cell_factors$model_beta
  
  results <- list(estimated_coefs = list(),
                  estimated_coef_sds = list(),
                  X = list())
  for (i in 1:reps) {
    a <- col_level_x #x1x3
    b <- col_level_y #x2x3
    col_level <- matrix(c(1.0, 0.0, a,
                          0.0, 1.0, b,
                          a, b, 1.0), nrow = 3)
    
    cur_samples <- generate_mvnorm_samples(col_level, 
                                           sample_size = sample_size)
    
    model_x <- cos(beta) * sin(alpha)
    model_y <- cos(beta) * cos(alpha)
    model_z <- sin(beta)
    
    cur_Y <- generate_y(cur_samples, 
                        expected_beta = c(model_x, model_y, model_z), 
                        R2 = R2,
                        has_intercept = FALSE)
    
    cur_model_summary <- summary(lm(cur_Y ~ cur_samples + 0))
    
    results[["estimated_coefs"]][[i]] <- cur_model_summary$coefficients[, 1]
    results[["estimated_coef_sds"]][[i]] <- cur_model_summary$coefficients[, 2]
    results[["X"]][[i]] <- c(col_level_x_factor = col_level_x,
                             col_level_y_factor = col_level_y,
                             R2_factor = R2, 
                             sample_size_factor = sample_size, 
                             model_alpha_factor = alpha,
                             model_beta_factor = beta)
  }
  
  X <- do.call(rbind, results[["X"]])
  estimated_coefs <- do.call(rbind, results[["estimated_coefs"]])
  estimated_coef_sds <- do.call(rbind, results[["estimated_coef_sds"]])
  
  est_dt <- data.table(estimated_coefs = estimated_coefs,
                       estimated_coef_sds = estimated_coef_sds,
                       X = X)
  
  names(est_dt) <- c(paste0("coef_", paste0("beta", 1:ncol(estimated_coefs))),
                     paste0("coef_sd_", paste0("beta", 1:ncol(estimated_coef_sds))), 
                     paste0("X_", c("col_level_x", "col_level_y", 
                                    "R2", "sample_size", 
                                    "model_alpha", "model_beta")))
    
  est_dt <- est_dt[,
                   .(est_sd_beta1 = sqrt(sum((coef_beta1 - mean(coef_beta1))^2) / (.N - 1)),
                     est_sd_beta2 = sqrt(sum((coef_beta2 - mean(coef_beta2))^2) / (.N - 1)),
                     est_sd_beta3 = sqrt(sum((coef_beta3 - mean(coef_beta3))^2) / (.N - 1))),
                   by = .(X_col_level_x, X_col_level_y, X_R2, X_sample_size, 
                          X_model_alpha, X_model_beta)]
  
  return(est_dt)
}

est_result <- generate_dataset_general_version(ds_factors, 
                                       cell_generate_function_for_test_3_3, 
                                       reps = 1000, use_parallel = TRUE)

est_sd_dt <- est_result$result_dt
est_sd_dt$X_col_level_x <- as.factor(est_sd_dt$X_col_level_x)
est_sd_dt$X_col_level_y <- as.factor(est_sd_dt$X_col_level_y)
est_sd_dt$X_model_alpha <- as.factor(est_sd_dt$X_model_alpha)
levels(est_sd_dt$X_model_alpha) <- paste0("alpha = ", c("0", "pi / 8", "pi / 4", "3 pi / 8", "pi / 2"))
est_sd_dt$X_model_beta<- as.factor(est_sd_dt$X_model_beta)
levels(est_sd_dt$X_model_beta) <- paste0("beta = ", c("0", "pi / 8", "pi / 4", "3 pi / 8", "pi / 2"))


colormap <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(32)

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta1)) +
  geom_tile(aes(fill = est_sd_beta1)) +
  scale_x_discrete(name = "Corr(X1, X3)") + 
  scale_y_discrete(name = "Corr(X2, X3)") +
  scale_fill_gradientn(colours = colormap, name = "Estimated Beta1 SD") +
  facet_grid(vars(X_model_alpha), vars(X_model_beta)) +
  theme_publish()

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta2)) +
  geom_tile(aes(fill = est_sd_beta2)) +
  scale_x_discrete(name = "Corr(X1, X3)") + 
  scale_y_discrete(name = "Corr(X2, X3)") +
  scale_fill_gradientn(colours = colormap, name = "Estimated Beta2 SD") +
  facet_grid(vars(X_model_alpha), vars(X_model_beta)) +
  theme_publish()

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta3)) +
  geom_tile(aes(fill = est_sd_beta3)) +
  scale_x_discrete(name = "Corr(X1, X3)") + 
  scale_y_discrete(name = "Corr(X2, X3)") +
  scale_fill_gradientn(colours = colormap, name = "Estimated Beta3 SD") +
  facet_grid(vars(X_model_alpha), vars(X_model_beta)) +
  theme_publish()


