rm(list = ls())

library(data.table)
library(ggplot2)
library(envalysis)

source(file = "./data_generation.r")

col_level_1_matrix <- matrix(c(1.0, 0.5, 0.2, 0.0,
                               0.5, 1.0, 0.2, 0.0,
                               0.2, 0.2, 1.0, 0.0,
                               0.0, 0.0, 0.0, 1.0), nrow = 4)
col_level_2_matrix <- matrix(c(1.0, 0.65, 0.4, 0.0,
                               0.65, 1.0, 0.4, 0.0,
                               0.4, 0.4, 1.0, 0.0,
                               0.0, 0.0, 0.0, 1.0), nrow = 4)
col_level_3_matrix <- matrix(c(1.0, 0.8, 0.6, 0.0,
                               0.8, 1.0, 0.6, 0.0,
                               0.6, 0.6, 1.0, 0.0,
                               0.0, 0.0, 0.0, 1.0), nrow = 4)
col_level_4_matrix <- matrix(c(1.0, 0.95, 0.8, 0.0,
                               0.95, 1.0, 0.8, 0.0,
                               0.8, 0.8, 1.0, 0.0,
                               0.0, 0.0, 0.0, 1.0), nrow = 4)
col_level_factor <- list(col_level_1_matrix,
                         col_level_2_matrix,
                         col_level_3_matrix,
                         col_level_4_matrix)
R2_factor <- c(0.25, 0.5, 0.75)
sample_size_factor <- c(30, 100, 150, 200, 250, 300)
true_coef_1 <- c(2, 0.5, 0.265, 0.0, 0.25)
true_coef_2 <- c(2, 0.4, 0.4, 0.0, 0.25)
model_factor <- list(true_coef_1, true_coef_2)

generate_acc <- function(est_df, models, mean_method = TRUE) {
  # calculate acc
  est_dt <- as.data.table(est_df)

  if (mean_method) {
    acc_dt <- est_dt[,
                     .(acc_coef_beta0 = abs(coef_beta0 - models[[X_model]][1]),
                       acc_coef_beta1 = abs(coef_beta1 - models[[X_model]][2]),
                       acc_coef_beta2 = abs(coef_beta2 - models[[X_model]][3]),
                       acc_coef_beta3 = abs(coef_beta3 - models[[X_model]][4]),
                       acc_coef_beta4 = abs(coef_beta4 - models[[X_model]][5]),
                       acc_sd_beta0 = abs(coef_sd_beta0 - mean(coef_sd_beta0) ),
                       acc_sd_beta1 = abs(coef_sd_beta1 - mean(coef_sd_beta1) ),
                       acc_sd_beta2 = abs(coef_sd_beta2 - mean(coef_sd_beta2) ),
                       acc_sd_beta3 = abs(coef_sd_beta3 - mean(coef_sd_beta3) ),
                       acc_sd_beta4 = abs(coef_sd_beta4 - mean(coef_sd_beta4) )),
                     by = .(X_col_level, X_R2, X_sample_size, X_model)]
  } else {
    acc_dt <- est_dt[,
                     .(acc_coef_beta0 = abs(coef_beta0 - models[[X_model]][1]),
                       acc_coef_beta1 = abs(coef_beta1 - models[[X_model]][2]),
                       acc_coef_beta2 = abs(coef_beta2 - models[[X_model]][3]),
                       acc_coef_beta3 = abs(coef_beta3 - models[[X_model]][4]),
                       acc_coef_beta4 = abs(coef_beta4 - models[[X_model]][5]),
                       acc_sd_beta0 = abs(coef_sd_beta0 - sqrt(sum((coef_beta0 - mean(coef_beta0))^2) / (.N - 1)) ),
                       acc_sd_beta1 = abs(coef_sd_beta1 - sqrt(sum((coef_beta1 - mean(coef_beta1))^2) / (.N - 1)) ),
                       acc_sd_beta2 = abs(coef_sd_beta2 - sqrt(sum((coef_beta2 - mean(coef_beta2))^2) / (.N - 1)) ),
                       acc_sd_beta3 = abs(coef_sd_beta3 - sqrt(sum((coef_beta3 - mean(coef_beta3))^2) / (.N - 1)) ),
                       acc_sd_beta4 = abs(coef_sd_beta4 - sqrt(sum((coef_beta4 - mean(coef_beta4))^2) / (.N - 1)) )),
                     by = .(X_col_level, X_R2, X_sample_size, X_model)]
  }
 
  return(acc_dt)
}

ds <- generate_dataset(col_levels = col_level_factor,
                       R2s = R2_factor,
                       sample_sizes = sample_size_factor,
                       models = model_factor,
                       reps = 100)

acc_data <- generate_acc(ds$est_df, model_factor, mean_method = TRUE)
acc_data$X_col_level <- as.factor(acc_data$X_col_level)
levels(acc_data$X_col_level) <- paste("Level", 1:4)
acc_data$X_R2 <- as.factor(acc_data$X_R2)
levels(acc_data$X_R2) <- R2_factor
acc_data$X_sample_size <- as.factor(acc_data$X_sample_size)
levels(acc_data$X_sample_size) <- sample_size_factor
acc_data$X_model <- as.factor(acc_data$X_model)
levels(acc_data$X_model) <- paste("Model", 1:2)

anova_model_1 <- anova(
  lm(acc_coef_beta1 ~ X_col_level * X_R2 * X_sample_size * X_model, 
     data = acc_data))
anova_model_1
anova_model_2 <- anova(
  lm(acc_sd_beta1 ~ X_col_level * X_R2 * X_sample_size * X_model, 
     data = acc_data))
anova_model_2

calculate_mean_acc <- function(acc_data) {
  result_dt <- acc_data[,
                   .(mean_acc_coef_beta0 = mean(acc_coef_beta0),
                     mean_acc_coef_beta1 = mean(acc_coef_beta1),
                     mean_acc_coef_beta2 = mean(acc_coef_beta2),
                     mean_acc_coef_beta3 = mean(acc_coef_beta3),
                     mean_acc_coef_beta4 = mean(acc_coef_beta4),
                     mean_acc_sd_beta0 = mean(acc_sd_beta0),
                     mean_acc_sd_beta1 = mean(acc_sd_beta1),
                     mean_acc_sd_beta2 = mean(acc_sd_beta2),
                     mean_acc_sd_beta3 = mean(acc_sd_beta3),
                     mean_acc_sd_beta4 = mean(acc_sd_beta4)),
                   by = .(X_col_level, X_R2, X_sample_size)]
  
  return(result_dt)
}

mean_acc_data <- calculate_mean_acc(acc_data)

R2_facet_labeller <- function(value) {
  return(data.frame("R2 = 0.25", "R2 = 0.5", "R2 = 0.75"))
}


# plot mean acc coef beta1
ggplot(mean_acc_data, aes(x = X_sample_size, y = mean_acc_coef_beta1, 
                          color = X_col_level)) +
  geom_point() +
  geom_line(mapping = aes(group = X_col_level)) + 
  facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
  scale_y_continuous(name = "Mean Beta1 Coef Estimation Error") +
  scale_x_discrete(name = "Sample Size") +
  scale_color_discrete(name = "Collinearity Level") + 
  theme_publish()

# plot mean acc coef beta3
ggplot(mean_acc_data, aes(x = X_sample_size, y = mean_acc_coef_beta3, 
                          color = X_col_level)) +
  geom_point() +
  geom_line(mapping = aes(group = X_col_level)) + 
  facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
  scale_y_continuous(name = "Mean Beta3 Coef Estimation Error") +
  scale_x_discrete(name = "Sample Size") +
  scale_color_discrete(name = "Collinearity Level") + 
  theme_publish()


# plot mean acc coef beta4
ggplot(mean_acc_data, aes(x = X_sample_size, y = mean_acc_coef_beta4, 
                          color = X_col_level)) +
  geom_point() +
  geom_line(mapping = aes(group = X_col_level)) + 
  facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
  scale_y_continuous(name = "Mean Beta4 Coef Estimation Error") +
  scale_x_discrete(name = "Sample Size") +
  scale_color_discrete(name = "Collinearity Level") + 
  theme_publish()


# plot mean acc sd beta1
ggplot(mean_acc_data, aes(x = X_sample_size, y = mean_acc_sd_beta1, 
                          color = X_col_level)) +
  geom_point() +
  geom_line(mapping = aes(group = X_col_level)) + 
  facet_grid(cols = vars(X_R2), labeller = labeller(X_R2 = R2_facet_labeller)) +
  scale_y_continuous(name = "Mean Beta1 SD Estimation Error") +
  scale_x_discrete(name = "Sample Size") +
  scale_color_discrete(name = "Collinearity Level") + 
  theme_publish()
  




