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
sapply(col_level_factor, function(x) {
  return(calculate_condition_number(generate_mvnorm_samples(x, 10000)))
})

R2_factor <- c(0.25, 0.5, 0.75)
sample_size_factor <- 100
model_factor <- list()
model_angles <- seq(pi / 12, pi / 2, by = pi / 12)
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
levels(est_sd_dt$X_model) <- c("pi / 12", "pi / 6", "pi / 4", 
                               "pi / 3", "5 pi / 12", "pi / 2")


R2_facet_labeller <- function(value) {
  return(data.frame("R2 = 0.25", "R2 = 0.5", "R2 = 0.75"))
}


# plot est sd of beta1
ggplot(est_sd_dt, aes(x = X_col_level, y = est_sd_beta1, 
                      color = X_model)) +
  geom_point() +
  geom_line(mapping = aes(group = X_model)) + 
  facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
  scale_y_continuous(name = "Beta1 SD") +
  scale_x_discrete(name = "Collinearity Level") +
  scale_color_discrete(name = "Angle of (Beta1, Beta2)") + 
  theme_publish()

# plot est sd of beta2
ggplot(est_sd_dt, aes(x = X_col_level, y = est_sd_beta2, 
                      color = X_model)) +
  geom_point() +
  geom_line(mapping = aes(group = X_model)) + 
  facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
  scale_y_continuous(name = "Beta2 SD") +
  scale_x_discrete(name = "Collinearity Level") +
  scale_color_discrete(name = "Angle of (Beta1, Beta2)") + 
  theme_publish()


