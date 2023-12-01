rm(list = ls())

library(data.table)
library(ggplot2)
library(envalysis)

source(file = "./data_generation.r")

col_level_factor <- list()
col_strengths <- seq(-0.80, 0.80, by = 0.1)
for (col_strength in col_strengths) {
  a <- col_strength
  col_level_factor[[as.character(col_strength)]] <- matrix(c(1.0, a,
                                                             a, 1.0), nrow = 2)
}

# test condition number
condition_numbers <- sapply(col_level_factor, function(x) {
  return(calculate_condition_number(generate_mvnorm_samples(x, 10000)))
})
condition_numbers <- data.frame(col_level = col_strengths, condition_number = condition_numbers)
condition_numbers$col_level <- as.factor(condition_numbers$col_level)
levels(condition_numbers$col_level) <- as.character(col_strengths)

R2_factor <- 0.5
sample_size_factor <- 100
model_factor <- list()
model_angles <- seq(0, pi / 2, by = pi / 8)
for (model_angle in model_angles) {
  model_factor[[as.character(model_angle)]] <- c(cos(model_angle), sin(model_angle))
}

ds <- generate_dataset(col_levels = col_level_factor,
                       R2s = R2_factor,
                       sample_sizes = sample_size_factor,
                       models = model_factor,
                       reps = 500,
                       has_intercept = FALSE)

calculate_sd <- function(est_df) {
  est_dt <- as.data.table(est_df)
  
  est_dt <- est_dt[,
                   .(est_sd_beta1 = sqrt(sum((coef_beta1 - mean(coef_beta1))^2) / (.N - 1)),
                     est_sd_beta2 = sqrt(sum((coef_beta2 - mean(coef_beta2))^2) / (.N - 1))),
                   by = .(X_col_level, X_R2, X_sample_size, X_model)]
  
  return(est_dt)
}

est_df <- ds$est_df
est_sd_dt <- calculate_sd(est_df)
est_sd_dt$X_col_level <- as.factor(est_sd_dt$X_col_level)
levels(est_sd_dt$X_col_level) <- as.character(col_strengths)
est_sd_dt$X_model <- as.factor(est_sd_dt$X_model)
levels(est_sd_dt$X_model) <- c("0", "pi / 8", "pi / 4", "3 pi / 8", 
                               "pi / 2")


# R2_facet_labeller <- function(value) {
#   return(data.frame("R2 = 0.25", "R2 = 0.5", "R2 = 0.75"))
# }

second_y_axis_ratio <- max(est_sd_dt$est_sd_beta1) / max(condition_numbers$condition_number)
ggplot() +
  geom_col(data = condition_numbers, 
           mapping = aes(x = col_level, y = condition_number * second_y_axis_ratio / 3),
           fill = "#E5FFCC", color = "black",
           position = "dodge") +
  geom_point(data = est_sd_dt, mapping = aes(x = X_col_level, y = est_sd_beta1, 
                                             color = X_model)) +
  geom_line(data = est_sd_dt, mapping = aes(x = X_col_level, y = est_sd_beta1,
                                            group = X_model, color = X_model)) + 
  scale_y_continuous(name = expression(paste("Estimated ", hat(beta)[1], " SD")),
                     sec.axis = sec_axis(~ . / second_y_axis_ratio * 3,
                                         name = expression(paste("Condition Number of ", bold(X))))) +
  scale_x_discrete(name = expression(paste(rho[12], " of correlation matrix"))) +
  scale_color_discrete(name = expression(paste("Angle of (" , beta[1], ", ",beta[2],")"))) + 
  theme_publish(base_size = 15, base_linewidth = 1)

second_y_axis_ratio <- max(est_sd_dt$est_sd_beta2) / max(condition_numbers$condition_number)
ggplot() +
  geom_col(data = condition_numbers, 
           mapping = aes(x = col_level, y = condition_number * second_y_axis_ratio / 3),
           fill = "#E5FFCC", color = "black",
           position = "dodge") +
  geom_point(data = est_sd_dt, mapping = aes(x = X_col_level, y = est_sd_beta2, 
                                             color = X_model)) +
  geom_line(data = est_sd_dt, mapping = aes(x = X_col_level, y = est_sd_beta2,
                                            group = X_model, color = X_model)) + 
  scale_y_continuous(name = expression(paste("Estimated ", hat(beta)[2], " SD")),
                     sec.axis = sec_axis(~ . / second_y_axis_ratio * 3,
                                         name = expression(paste("Condition Number of ", bold(X))))) +
  scale_x_discrete(name = expression(paste(rho[12], " of correlation matrix"))) +
  scale_color_discrete(name = expression(paste("Angle of (" , beta[1], ", ",beta[2],")"))) + 
  theme_publish(base_size = 15, base_linewidth = 1)

# # plot est sd of beta1
# ggplot(est_sd_dt, aes(x = X_col_level, y = est_sd_beta1, 
#                       color = X_model)) +
#   geom_point() +
#   geom_line(mapping = aes(group = X_model)) + 
#   facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
#   scale_y_continuous(name = "Beta1 SD") +
#   scale_x_discrete(name = "Collinearity Level") +
#   scale_color_discrete(name = "Angle of (Beta1, Beta2)") + 
#   theme_publish()
# 
# # plot est sd of beta2
# ggplot(est_sd_dt, aes(x = X_col_level, y = est_sd_beta2, 
#                       color = X_model)) +
#   geom_point() +
#   geom_line(mapping = aes(group = X_model)) + 
#   facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
#   scale_y_continuous(name = "Beta2 SD") +
#   scale_x_discrete(name = "Collinearity Level") +
#   scale_color_discrete(name = "Angle of (Beta1, Beta2)") + 
#   theme_publish()


