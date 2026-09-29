## Results from the reconstruction
not_normalized <- data.frame(
  Model = c("Original_Model", "RS42", "RS96", "RS23", "RS14", "RS97", "RS56", "RS78", "RS12", "RS62", "RS83"),
  RMSE = c(0.97, 0.9590082, 0.8257291, 0.9030714, 0.8766148, 0.8445143, 1.0352365, 1.0101544, 0.99385, 0.8404565, 0.7984785),
  MAE = c(0.785, 0.7968548, 0.7252967, 0.8011899, 0.7749198, 0.7543971, 0.9110405, 0.8466437, 0.8804312, 0.7312622, 0.6958504),
  RSquared = c(0.949, 0.9497634, 0.9627565, 0.9554528, 0.9580247, 0.9610426, 0.9414597, 0.944262, 0.9460468, 0.9614161, 0.9651741),
  MRE = c(11.138, 16.4281668, 13.0585175, 15.4687968, 15.8059001, 15.1752186, 18.1649176, 17.3834385, 17.6132565, 13.8107595, 15.0196954),
  REmax = c(28.687, 47.8953938, 26.3942219, 44.1524028, 46.5398477, 53.8374919, 59.5461606, 55.7939376, 48.0410224, 33.4799174, 56.6200744)
)

actb_normalized <- data.frame(
  Model = c("Original_Model", "RS42", "RS96", "RS23", "RS14", "RS97", "RS56", "RS78", "RS12", "RS62", "RS83"),
  RMSE = c(0.97, 0.9785549, 0.8890049, 0.9925517, 0.906599, 1.0755243, 1.0467913, 1.0534075, 0.9622348, 1.0848004, 0.9912945),
  MAE = c(0.785, 0.8697101, 0.7309728, 0.8863158, 0.7979444, 0.9492748, 0.9461661, 0.9231825, 0.8577676, 0.9761635, 0.90631),
  RSquared = c(0.949, 0.9476946, 0.9568298, 0.9461876, 0.9551041, 0.9368147, 0.9401456, 0.9393866, 0.9494248, 0.9357201, 0.9463239),
  MRE = c(11.138, 16.9375149, 12.1183711, 17.3946305, 14.9243027, 18.5402708, 17.1514452, 18.5633467, 14.5426096, 18.6953352, 16.3513918),
  REmax  = c(28.687, 59.1901526 , 30.7570529 , 57.0417206 , 35.3727026 , 56.0210236 , 54.8476289 , 67.0355538 , 29.7240522 , 54.429246 , 45.374453)
)

actb_normalized <- actb_normalized %>% column_to_rownames("Model")
not_normalized <- not_normalized %>% column_to_rownames("Model")

## Descriptive stats
mean_nn <- not_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) mean(x)))
sd_nn <- not_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) sd(x)))
median_nn <- not_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) median(x)))
iqr_nn <- not_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) IQR(x)))
mean_an <- actb_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) mean(x)))
sd_an <- actb_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) sd(x)))
median_an <- actb_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) median(x)))
iqr_an <- actb_normalized %>% slice(-1) %>%
  summarise(across(where(is.numeric), \(x) IQR(x)))

CI_RMSE_nn <- confidence_interval(not_normalized, "RMSE")
CI_RMSE_an <- confidence_interval(actb_normalized, "RMSE")
CI_MAE_nn <- confidence_interval(not_normalized, "MAE")
CI_MAE_an <- confidence_interval(actb_normalized, "MAE")
CI_R2_nn <- confidence_interval(not_normalized, "RSquared")
CI_R2_an <- confidence_interval(actb_normalized, "RSquared")
CI_MRE_nn <- confidence_interval(not_normalized, "MRE")
CI_MRE_an <- confidence_interval(actb_normalized, "MRE")
CI_REmax_nn <- confidence_interval(not_normalized, "REmax")
CI_REmax_an <- confidence_interval(actb_normalized, "REmax")

t.test(not_normalized$RMSE[2:11], actb_normalized$RMSE[2:11], paired = TRUE)
t.test(not_normalized$MAE[2:11], actb_normalized$MAE[2:11], paired = TRUE)
t.test(not_normalized$RSquared[2:11], actb_normalized$RSquared[2:11], paired = TRUE)
t.test(not_normalized$MRE[2:11], actb_normalized$MRE[2:11], paired = TRUE)
t.test(not_normalized$REmax[2:11], actb_normalized$REmax[2:11], paired = TRUE)

## Seed metrics
metrics_nn <- data.frame(
  Metric = c("RMSE", "MAE", "RSquared", "MRE", "REmax"),
  Original = c(0.97, 0.785, 0.949, 11.138, 28.687),
  Mean = c(mean_nn$RMSE, mean_nn$MAE, mean_nn$RSquared, mean_nn$MRE, mean_nn$REmax),
  SD = c(sd_nn$RMSE, sd_nn$MAE, sd_nn$RSquared, sd_nn$MRE, sd_nn$REmax),
  Median = c(median_nn$RMSE, median_nn$MAE, median_nn$RSquared, median_nn$MRE, median_nn$REmax),
  CI_lower = c(CI_RMSE_nn[1], CI_MAE_nn[1], CI_R2_nn[1], CI_MRE_nn[1], CI_REmax_nn[1]),
  CI_upper = c(CI_RMSE_nn[2], CI_MAE_nn[2], CI_R2_nn[2], CI_MRE_nn[2], CI_REmax_nn[2]),
  IQR = c(iqr_nn$RMSE, iqr_nn$MAE, iqr_nn$RSquared, iqr_nn$MRE, iqr_nn$REmax)
)

metrics_an <- data.frame(
  Metric = c("RMSE", "MAE", "RSquared", "MRE", "REmax"),
  Original = c(0.97, 0.785, 0.949, 11.138, 28.687),
  Mean = c(mean_an$RMSE, mean_an$MAE, mean_an$RSquared, mean_an$MRE, mean_an$REmax),
  SD = c(sd_an$RMSE, sd_an$MAE, sd_an$RSquared, sd_an$MRE, sd_an$REmax),
  Median = c(median_an$RMSE, median_an$MAE, median_an$RSquared, median_an$MRE, median_an$REmax),
  CI_lower = c(CI_RMSE_an[1], CI_MAE_an[1], CI_R2_an[1], CI_MRE_an[1], CI_REmax_an[1]),
  CI_upper = c(CI_RMSE_an[2], CI_MAE_an[2], CI_R2_an[2], CI_MRE_an[2], CI_REmax_an[2]),
  IQR = c(iqr_an$RMSE, iqr_an$MAE, iqr_an$RSquared, iqr_an$MRE, iqr_an$REmax)
)
