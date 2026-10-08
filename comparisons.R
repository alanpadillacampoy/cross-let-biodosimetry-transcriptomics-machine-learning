## this is to do the comparisons between different experimental summaries

original_signature_performance_summary <- read.csv("100_seeds_original_performance_summary.csv")
reduced_signature_performance_summary <- read.csv("100_seeds_reduced_performance_summary.csv")
normalized_signature_performance_summary <- read.csv("100_seeds_normalized_performance_summary.csv")
original_signature_performance <- read.csv("100_seeds_original_performance.csv")
reduced_signature_performance <- read.csv("100_seeds_reduced_performance.csv")
normalized_signature_performance <- read.csv("100_seeds_normalized_performance.csv")
# Paired t-tests ----
# RMSE
rmse_ttest <- t.test(unlist(normal_signature_performance[1,3:102]), 
       unlist(reduced_signature_performance[1,3:102]), paired = TRUE)
# MAE
mae_ttest <- t.test(unlist(normal_signature_performance[2,3:102]), 
                     unlist(reduced_signature_performance[2,3:102]), paired = TRUE)
# R2
r2_ttest <- t.test(unlist(normal_signature_performance[3,3:102]), 
                     unlist(reduced_signature_performance[3,3:102]), paired = TRUE)
# MRE
mre_ttest <- t.test(unlist(normal_signature_performance[4,3:102]), 
                     unlist(reduced_signature_performance[4,3:102]), paired = TRUE)
# REmax
remax_ttest <- t.test(unlist(normal_signature_performance[5,3:102]), 
                     unlist(reduced_signature_performance[5,3:102]), paired = TRUE)

rmse_ttest
mae_ttest
r2_ttest
mre_ttest
remax_ttest

## Paired differences
remax_complete <- unlist(normal_signature_performance[5, 3:102])
remax_reduced   <- unlist(reduced_signature_performance[5, 3:102])

remax_diff <- remax_complete - remax_reduced

mae_complete <- unlist(normal_signature_performance[2, 3:102])
mae_reduced   <- unlist(reduced_signature_performance[2, 3:102])

mae_diff <- mae_complete - mae_reduced

sum(remax_diff > 0)
sum(remax_diff < 0)
sum(remax_diff == 0)

sum(mae_diff > 0)
sum(mae_diff < 0)
sum(mae_diff == 0)
