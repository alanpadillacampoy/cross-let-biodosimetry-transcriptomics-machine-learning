## this is to do the comparisons between different experimental summaries

original_signature_performance_summary <- read.csv("100_seeds_original_performance_summary.csv")
reduced_signature_performance_summary <- read.csv("100_seeds_reduced_performance_summary.csv")
normalized_signature_performance_summary <- read.csv("100_seeds_normalized_performance_summary.csv")
original_signature_performance <- read.csv("100_seeds_original_performance.csv")
reduced_signature_performance <- read.csv("100_seeds_reduced_performance.csv")
normalized_signature_performance <- read.csv("100_seeds_normalized_performance.csv")
## Paired t-tests ----
# Normalized vs non-normalized ----
# RMSE
rmse_ttest_norm <- t.test(unlist(original_signature_performance[1,3:102]), 
       unlist(normalized_signature_performance[1,3:102]), paired = TRUE)
# MAE
mae_ttest_norm <- t.test(unlist(original_signature_performance[2,3:102]), 
                     unlist(normalized_signature_performance[2,3:102]), paired = TRUE)
# R2
r2_ttest_norm <- t.test(unlist(original_signature_performance[3,3:102]), 
                     unlist(normalized_signature_performance[3,3:102]), paired = TRUE)
# MRE
mre_ttest_norm <- t.test(unlist(original_signature_performance[4,3:102]), 
                     unlist(normalized_signature_performance[4,3:102]), paired = TRUE)
# REmax
remax_ttest_norm <- t.test(unlist(original_signature_performance[5,3:102]), 
                     unlist(normalized_signature_performance[5,3:102]), paired = TRUE)

rmse_ttest_norm
mae_ttest_norm
r2_ttest_norm
mre_ttest_norm
remax_ttest_norm

# Paired differences
remax_complete_norm <- unlist(original_signature_performance[5, 3:102])
remax_reduced_norm   <- unlist(normalized_signature_performance[5, 3:102])

remax_diff_norm <- remax_complete_norm - remax_reduced_norm

mae_complete_norm <- unlist(original_signature_performance[2, 3:102])
mae_reduced_norm   <- unlist(normalized_signature_performance[2, 3:102])

mae_diff_norm <- mae_complete_norm - mae_reduced_norm

sum(remax_diff_norm > 0)
sum(remax_diff_norm < 0)
sum(remax_diff_norm == 0)

sum(mae_diff_norm > 0)
sum(mae_diff_norm < 0)
sum(mae_diff_norm == 0)

# Original vs reduced signature ----
# RMSE
rmse_ttest_red <- t.test(unlist(original_signature_performance[1,3:102]), 
                     unlist(reduced_signature_performance[1,3:102]), paired = TRUE)
# MAE
mae_ttest_red <- t.test(unlist(original_signature_performance[2,3:102]), 
                    unlist(reduced_signature_performance[2,3:102]), paired = TRUE)
# R2
r2_ttest_red <- t.test(unlist(original_signature_performance[3,3:102]), 
                   unlist(reduced_signature_performance[3,3:102]), paired = TRUE)
# MRE
mre_ttest_red <- t.test(unlist(original_signature_performance[4,3:102]), 
                    unlist(reduced_signature_performance[4,3:102]), paired = TRUE)
# REmax
remax_ttest_red <- t.test(unlist(original_signature_performance[5,3:102]), 
                      unlist(reduced_signature_performance[5,3:102]), paired = TRUE)

rmse_ttest_red
mae_ttest_red
r2_ttest_red
mre_ttest_red
remax_ttest_red

# Paired differences
remax_complete_red <- unlist(original_signature_performance[5, 3:102])
remax_reduced_red   <- unlist(reduced_signature_performance[5, 3:102])

remax_diff_red <- remax_complete_red - remax_reduced_red

mae_complete_red <- unlist(original_signature_performance[2, 3:102])
mae_reduced_red   <- unlist(reduced_signature_performance[2, 3:102])

mae_diff_red <- mae_complete_red - mae_reduced_red

sum(remax_diff_red > 0)
sum(remax_diff_red < 0)
sum(remax_diff_red == 0)

sum(mae_diff_red > 0)
sum(mae_diff_red < 0)
sum(mae_diff_red == 0)
