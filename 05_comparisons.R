## this is to do the comparisons between different experimental summaries

raw_signature_performance_summary <- read.csv("100_seeds_original_performance_summary.csv")
reduced_signature_performance_summary <- read.csv("100_seeds_reduced_performance_summary.csv")
normalized_signature_performance_summary <- read.csv("100_seeds_normalized_performance_summary.csv")
raw_signature_performance <- read.csv("100_seeds_original_performance.csv")
reduced_signature_performance <- read.csv("100_seeds_reduced_performance.csv")
normalized_signature_performance <- read.csv("100_seeds_normalized_performance.csv")
## Paired t-tests ----
# Normalized vs non-normalized ----
# RMSE
rmse_ttest_norm <- t.test(
  unlist(normalized_signature_performance[1, 3:102]),
  unlist(raw_signature_performance[1, 3:102]),
  paired = TRUE
)

# MAE
mae_ttest_norm <- t.test(
  unlist(normalized_signature_performance[2, 3:102]),
  unlist(raw_signature_performance[2, 3:102]),
  paired = TRUE
)

# R2
r2_ttest_norm <- t.test(
  unlist(normalized_signature_performance[3, 3:102]),
  unlist(raw_signature_performance[3, 3:102]),
  paired = TRUE
)

# MRE
mre_ttest_norm <- t.test(
  unlist(normalized_signature_performance[4, 3:102]),
  unlist(raw_signature_performance[4, 3:102]),
  paired = TRUE
)

# REmax
remax_ttest_norm <- t.test(
  unlist(normalized_signature_performance[5, 3:102]),
  unlist(raw_signature_performance[5, 3:102]),
  paired = TRUE
)

rmse_ttest_norm
mae_ttest_norm
r2_ttest_norm
mre_ttest_norm
remax_ttest_norm

results <- tibble(
  Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
  Norm = c(0.993, 0.87, 0.946, 15.994, 50.56),
  Unnorm = c(0.898, 0.78, 0.956, 14.724, 41.219),
  Difference = c(
    rmse_ttest_norm$estimate[[1]],
    mae_ttest_norm$estimate[[1]],
    r2_ttest_norm$estimate[[1]],
    mre_ttest_norm$estimate[[1]],
    remax_ttest_norm$estimate[[1]]
  ),
  CI_lower = c(
    rmse_ttest_norm$conf.int[[1]],
    mae_ttest_norm$conf.int[[1]],
    r2_ttest_norm$conf.int[[1]],
    mre_ttest_norm$conf.int[[1]],
    remax_ttest_norm$conf.int[[1]]
  ),
  CI_upper = c(
    rmse_ttest_norm$conf.int[[2]],
    mae_ttest_norm$conf.int[[2]],
    r2_ttest_norm$conf.int[[2]],
    mre_ttest_norm$conf.int[[2]],
    remax_ttest_norm$conf.int[[2]]
  ),
  p_value = c(
    rmse_ttest_norm$p.value,
    mae_ttest_norm$p.value,
    r2_ttest_norm$p.value,
    mre_ttest_norm$p.value,
    remax_ttest_norm$p.value
  )
)

view(results)

p_values <- c(rmse_ttest_norm$p.value, mae_ttest_norm$p.value, r2_ttest_norm$p.value,
              mre_ttest_norm$p.value, remax_ttest_norm$p.value)
adjusted_p_values <- p.adjust(p_values, method = "holm")

# Paired differences
remax_complete_norm <- unlist(raw_signature_performance[5, 3:102])
remax_reduced_norm   <- unlist(normalized_signature_performance[5, 3:102])

remax_diff_norm <- remax_complete_norm - remax_reduced_norm

mae_complete_norm <- unlist(raw_signature_performance[2, 3:102])
mae_reduced_norm   <- unlist(normalized_signature_performance[2, 3:102])

mae_diff_norm <- mae_complete_norm - mae_reduced_norm

sum(remax_diff_norm > 0)
sum(remax_diff_norm < 0)
sum(remax_diff_norm == 0)

sum(mae_diff_norm > 0)
sum(mae_diff_norm < 0)
sum(mae_diff_norm == 0)

# raw vs reduced signature ----
# RMSE
rmse_ttest_red <- t.test(unlist(raw_signature_performance[1,3:102]), 
                     unlist(reduced_signature_performance[1,3:102]), paired = TRUE)
# MAE
mae_ttest_red <- t.test(unlist(raw_signature_performance[2,3:102]), 
                    unlist(reduced_signature_performance[2,3:102]), paired = TRUE)
# R2
r2_ttest_red <- t.test(unlist(raw_signature_performance[3,3:102]), 
                   unlist(reduced_signature_performance[3,3:102]), paired = TRUE)
# MRE
mre_ttest_red <- t.test(unlist(raw_signature_performance[4,3:102]), 
                    unlist(reduced_signature_performance[4,3:102]), paired = TRUE)
# REmax
remax_ttest_red <- t.test(unlist(raw_signature_performance[5,3:102]), 
                      unlist(reduced_signature_performance[5,3:102]), paired = TRUE)

rmse_ttest_red
mae_ttest_red
r2_ttest_red
mre_ttest_red
remax_ttest_red

data.frame(
  Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
  Whole = c(0.898, 0.780, 0.956, 14.724, 41.219),
  Reduced = c(0.895, 0.797, 0.956, 14.293, 33.284),
  Difference = c(
    rmse_ttest_red$estimate[[1]],
    mae_ttest_red$estimate[[1]],
    r2_ttest_red$estimate[[1]],
    mre_ttest_red$estimate[[1]],
    remax_ttest_red$estimate[[1]]
  ),
  CI_lower = c(
    rmse_ttest_red$conf.int[[1]],
    mae_ttest_red$conf.int[[1]],
    r2_ttest_red$conf.int[[1]],
    mre_ttest_red$conf.int[[1]],
    remax_ttest_red$conf.int[[1]]
  ),
  CI_upper = c(
    rmse_ttest_red$conf.int[[2]],
    mae_ttest_red$conf.int[[2]],
    r2_ttest_red$conf.int[[2]],
    mre_ttest_red$conf.int[[2]],
    remax_ttest_red$conf.int[[2]]
  ),
  p_value = c(
    rmse_ttest_red$p.value,
    mae_ttest_red$p.value,
    r2_ttest_red$p.value,
    mre_ttest_red$p.value,
    remax_ttest_red$p.value
  )
)

p_values <- c(rmse_ttest_red$p.value, mae_ttest_red$p.value, r2_ttest_red$p.value,
              mre_ttest_red$p.value, remax_ttest_red$p.value)
adjusted_p_values <- p.adjust(p_values, method = "holm")

# Paired differences
remax_complete_red <- unlist(raw_signature_performance[5, 3:102])
remax_reduced_red   <- unlist(reduced_signature_performance[5, 3:102])

remax_diff_red <- remax_complete_red - remax_reduced_red

mae_complete_red <- unlist(raw_signature_performance[2, 3:102])
mae_reduced_red   <- unlist(reduced_signature_performance[2, 3:102])

mae_diff_red <- mae_complete_red - mae_reduced_red

sum(remax_diff_red > 0)
sum(remax_diff_red < 0)
sum(remax_diff_red == 0)

sum(mae_diff_red > 0)
sum(mae_diff_red < 0)
sum(mae_diff_red == 0)

beep(1)



# 
# raw_signature_performance
# shapiro.test(as.numeric(raw_signature_performance[1,3:102]))
# qqnorm(as.numeric(raw_signature_performance[5,3:102]),  main = "Normal Q-Q Plot", col = "blue", pch = 19)
# qqline(as.numeric(raw_signature_performance[5,3:102]),  col = "red", lwd = 2)
# t.test(
#   as.numeric(raw_signature_performance[5, 3:102]),
#   mu = 28.687
# )

