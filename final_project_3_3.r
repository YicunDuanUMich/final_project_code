rm(list = ls())

library(data.table)
library(ggplot2)
library(envalysis)
library(grDevices)
library(RColorBrewer)

source(file = "./data_generation.r")

col_level_x_factor <- seq(-0.45, 0.45, by = 0.05)
col_level_y_factor <- seq(-0.45, 0.45, by = 0.05)

R2_factor <- 0.5
sample_size_factor <- 100

model_psi <- seq(0, pi / 2, by = pi / 8)
model_theta <- seq(0, pi / 2, by = pi / 8)

ds_factors <- list(col_level_x_factor = col_level_x_factor,
                   col_level_y_factor = col_level_y_factor,
                   R2_factor = R2_factor,
                   sample_size_factor = sample_size_factor,
                   model_psi = model_psi,
                   model_theta = model_theta)

cell_generate_function_for_test_3_3 <- function(cell_factors, reps) {
  col_level_x <- cell_factors$col_level_x_factor
  col_level_y <- cell_factors$col_level_y_factor
  R2 <- cell_factors$R2_factor
  sample_size <- cell_factors$sample_size_factor
  psi <- cell_factors$model_psi
  theta <- cell_factors$model_theta
  
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
    
    model_x <- sin(psi) * cos(theta)
    model_y <- sin(psi) * sin(theta)
    model_z <- cos(psi)
    
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
                             model_psi_factor = psi,
                             model_theta_factor = theta)
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
                                    "model_psi", "model_theta")))
    
  est_dt <- est_dt[,
                   .(est_sd_beta1 = sqrt(sum((coef_beta1 - mean(coef_beta1))^2) / (.N - 1)),
                     est_sd_beta2 = sqrt(sum((coef_beta2 - mean(coef_beta2))^2) / (.N - 1)),
                     est_sd_beta3 = sqrt(sum((coef_beta3 - mean(coef_beta3))^2) / (.N - 1))),
                   by = .(X_col_level_x, X_col_level_y, X_R2, X_sample_size, 
                          X_model_psi, X_model_theta)]
  
  return(est_dt)
}

est_result <- generate_dataset_general_version(ds_factors, 
                                       cell_generate_function_for_test_3_3, 
                                       reps = 5000, use_parallel = TRUE)

est_sd_dt <- est_result$result_dt
est_sd_dt$X_col_level_x <- as.factor(est_sd_dt$X_col_level_x)
est_sd_dt$X_col_level_y <- as.factor(est_sd_dt$X_col_level_y)
est_sd_dt$X_model_psi <- as.factor(est_sd_dt$X_model_psi)
levels(est_sd_dt$X_model_psi) <- paste0("psi = ", c("0", "pi / 8", "pi / 4", "3 pi / 8", "pi / 2"))
est_sd_dt$X_model_theta<- as.factor(est_sd_dt$X_model_theta)
levels(est_sd_dt$X_model_theta) <- paste0("theta = ", c("0", "pi / 8", "pi / 4", "3 pi / 8", "pi / 2"))


colormap <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(32)

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta1)) +
  geom_tile(aes(fill = est_sd_beta1)) +
  scale_x_discrete(name = expression(rho[13])) + 
  scale_y_discrete(name = expression(rho[23])) +
  scale_fill_gradientn(colours = colormap, name = expression(paste("Estimated ", hat(beta)[1], " SD"))) +
  facet_grid(vars(X_model_psi), vars(X_model_theta)) +
  theme_publish(base_size = 10, base_linewidth = 1) +
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        legend.text = element_text(size = 12),
        legend.title=element_text(size = 14),
        axis.title=element_text(size = 15, face = "bold"))

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta2)) +
  geom_tile(aes(fill = est_sd_beta2)) +
  scale_x_discrete(name = expression(rho[13])) + 
  scale_y_discrete(name = expression(rho[23])) +
  scale_fill_gradientn(colours = colormap, name = expression(paste("Estimated ", hat(beta)[2], " SD"))) +
  facet_grid(vars(X_model_psi), vars(X_model_theta)) +
  theme_publish(base_size = 10, base_linewidth = 1) +
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        legend.text = element_text(size = 12),
        legend.title=element_text(size = 14),
        axis.title=element_text(size = 15, face = "bold"))

ggplot(est_sd_dt, aes(x = X_col_level_x, y = X_col_level_y, z = est_sd_beta3)) +
  geom_tile(aes(fill = est_sd_beta3)) +
  scale_x_discrete(name = expression(rho[13])) + 
  scale_y_discrete(name = expression(rho[23])) +
  scale_fill_gradientn(colours = colormap, name = expression(paste("Estimated ", hat(beta)[3], " SD"))) +
  facet_grid(vars(X_model_psi), vars(X_model_theta)) +
  theme_publish(base_size = 10, base_linewidth = 1) +
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(3, "cm"),
        legend.text = element_text(size = 12),
        legend.title=element_text(size = 14),
        axis.title=element_text(size = 15, face = "bold"))


