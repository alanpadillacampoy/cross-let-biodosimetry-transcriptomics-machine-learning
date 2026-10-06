## Results from the reconstruction ----
model_tested <- as.data.frame(t(model_performance))
colnames(model_tested) <- as.character(model_tested[1,])
model_tested <- model_tested[-1,]
rownames(model_tested) <- sub("^seed_", "RS", rownames(model_tested))
model_tested <- as.data.frame(lapply(model_tested, as.numeric))

## Descriptive stats ---- 
mean <- model_tested %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) mean(x)))
sd <- model_tested %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) sd(x)))
median <- model_tested %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) median(x)))
iqr <- model_tested %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) IQR(x)))

## Confidence Intervals ----
CI_RMSE <- confidence_interval(model_tested, "RMSE")
CI_MAE <- confidence_interval(model_tested, "MAE")
CI_R2 <- confidence_interval(model_tested, "R2")
CI_MRE <- confidence_interval(model_tested, "MRE")
CI_REmax <- confidence_interval(model_tested, "REmax")

## Paired t-tests ----
# t.test(not_normalized$ RMSE[2:11], actb_normalized$RMSE[2:11], paired = TRUE)
# t.test(not_normalized$MAE[2:11], actb_normalized$MAE[2:11], paired = TRUE)
# t.test(not_normalized$RSquared[2:11], actb_normalized$RSquared[2:11], paired = TRUE)
# t.test(not_normalized$MRE[2:11], actb_normalized$MRE[2:11], paired = TRUE)
# t.test(not_normalized$REmax[2:11], actb_normalized$REmax[2:11], paired = TRUE)

## Seed metrics ----
metrics <- data.frame(
  Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
  Original = c(0.97, 0.785, 0.949, 11.138, 28.687),
  Mean = c(mean$RMSE, mean$MAE, mean$R2, mean$MRE, mean$REmax),
  SD = c(sd$RMSE, sd$MAE, sd$R2, sd$MRE, sd$REmax),
  Median = c(median$RMSE, median$MAE, median$R2, median$MRE, median$REmax),
  CI_lower = c(CI_RMSE[1], CI_MAE[1], CI_R2[1], CI_MRE[1], CI_REmax[1]),
  CI_upper = c(CI_RMSE[2], CI_MAE[2], CI_R2[2], CI_MRE[2], CI_REmax[2]),
  IQR = c(iqr$RMSE, iqr$MAE, iqr$R2, iqr$MRE, iqr$REmax)
)
metrics <- metrics %>% mutate(across(where(is.numeric), ~ round(.x, 3)))