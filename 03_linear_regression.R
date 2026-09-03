library("dplyr")
library("ggplot2")

list_se <- readRDS("list_se.rds")
final_matrices <- readRDS("final_matrices.rds")
changes <- readRDS("changes_data_cleaning.rds")
wang_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "H4c3", "H4c9", "Ifit1", "Igfbp4", "LOC118567921",
                          "Lrrc70", "Ms4a1", "Phlda3", "Ptprn", "Rps20", "Serpine2", "Thy1")
#Wang data

#read wang data
wang_training_data <- read.csv("wang_training_dataset.csv")
wang_training_data <- as.data.frame(t(wang_training_data))
colnames(wang_training_data) <- wang_training_data[1,]
wang_training_data <- wang_training_data[-1,]
wang_training_data <- type.convert(wang_training_data, as.is = TRUE)

#select genes
wang_training_data <- wang_training_data %>% dplyr::select(c("Dose", all_of(wang_signature_genes)))

#Linear Regression Model
linear_model <- lm(Dose ~ ., data = wang_training_data)
residuals <- linear_model$residuals
hist(residuals)
# Plot the residuals
qqnorm(residuals)
# Plot the Q-Q line
qqline(residuals)

#Plot the model line
# Extract fitted values
wang_training_data$pred <- predict(linear_model)

ggplot(wang_training_data, aes(x = pred, y = Dose)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "darkgreen", linetype = "dashed") +
  labs(
    x = "Predicted Dose",
    y = "Actual Dose",
    title = "Actual vs. Predicted Values"
  ) +
  theme_minimal()

#Residual vs predicted
# Extract residuals
wang_training_data$res <- residuals(linear_model)

ggplot(wang_training_data, aes(x = pred, y = res)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "orange") +
  labs(
    x = "Predicted Dose",
    y = "Residual (Actual - Predicted)",
    title = "Residuals vs. Predicted Values"
  ) +
  theme_minimal()

library(dotwhisker)

# Plots coefficients with confidence intervals
dwplot(linear_model) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Variable Impacts (Coefficients)")


#Performance metrics
# 1. Get predictions and actual values
actual <- wang_training_data$Dose
predicted <- predict(linear_model, newdata = wang_training_data)
residuals <- actual - predicted

# 2. Calculate metrics
rmse <- sqrt(mean(residuals^2))
mae <- mean(abs(residuals))

# R-squared (from summary or manually)
r_squared <- summary(linear_model)$r.squared
adj_r_squared <- summary(linear_model)$adj.r.squared

# Display results
cat("RMSE:", rmse, "\nMAE:", mae, "\nR2:", r_squared, "\nAdj R2:", adj_r_squared)

