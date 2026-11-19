## this is to do the comparisons between different experimental summaries

# raw_signature_performance_summary <- read.csv("100_seeds_original_performance_summary.csv")
# reduced_signature_performance_summary <- read.csv("100_seeds_reduced_performance_summary.csv")
# normalized_signature_performance_summary <- read.csv("100_seeds_normalized_performance_summary.csv")
# raw_signature_performance <- read.csv("000_no_time.csv")
# reduced_signature_performance <- read.csv("000_time.csv")
# normalized_signature_performance <- read.csv("100_seeds_normalized_performance.csv")
## Paired t-tests ----
# Normalized vs non-normalized ----
# RMSE
# rmse_ttest_norm <- t.test(
#   unlist(normalized_signature_performance[1, 3:102]),
#   unlist(raw_signature_performance[1, 3:102]),
#   paired = TRUE
# )
# 
# # MAE
# mae_ttest_norm <- t.test(
#   unlist(normalized_signature_performance[2, 3:102]),
#   unlist(raw_signature_performance[2, 3:102]),
#   paired = TRUE
# )
# 
# # R2
# r2_ttest_norm <- t.test(
#   unlist(normalized_signature_performance[3, 3:102]),
#   unlist(raw_signature_performance[3, 3:102]),
#   paired = TRUE
# )
# 
# # MRE
# mre_ttest_norm <- t.test(
#   unlist(normalized_signature_performance[4, 3:102]),
#   unlist(raw_signature_performance[4, 3:102]),
#   paired = TRUE
# )
# 
# # REmax
# remax_ttest_norm <- t.test(
#   unlist(normalized_signature_performance[5, 3:102]),
#   unlist(raw_signature_performance[5, 3:102]),
#   paired = TRUE
# )
# 
# rmse_ttest_norm
# mae_ttest_norm
# r2_ttest_norm
# mre_ttest_norm
# remax_ttest_norm
# 
# results <- tibble(
#   Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
#   Norm = c(0.993, 0.87, 0.946, 15.994, 50.56),
#   Unnorm = c(0.898, 0.78, 0.956, 14.724, 41.219),
#   Difference = c(
#     rmse_ttest_norm$estimate[[1]],
#     mae_ttest_norm$estimate[[1]],
#     r2_ttest_norm$estimate[[1]],
#     mre_ttest_norm$estimate[[1]],
#     remax_ttest_norm$estimate[[1]]
#   ),
#   CI_lower = c(
#     rmse_ttest_norm$conf.int[[1]],
#     mae_ttest_norm$conf.int[[1]],
#     r2_ttest_norm$conf.int[[1]],
#     mre_ttest_norm$conf.int[[1]],
#     remax_ttest_norm$conf.int[[1]]
#   ),
#   CI_upper = c(
#     rmse_ttest_norm$conf.int[[2]],
#     mae_ttest_norm$conf.int[[2]],
#     r2_ttest_norm$conf.int[[2]],
#     mre_ttest_norm$conf.int[[2]],
#     remax_ttest_norm$conf.int[[2]]
#   ),
#   p_value = c(
#     rmse_ttest_norm$p.value,
#     mae_ttest_norm$p.value,
#     r2_ttest_norm$p.value,
#     mre_ttest_norm$p.value,
#     remax_ttest_norm$p.value
#   )
# )
# 
# view(results)
# 
# p_values <- c(rmse_ttest_norm$p.value, mae_ttest_norm$p.value, r2_ttest_norm$p.value,
#               mre_ttest_norm$p.value, remax_ttest_norm$p.value)
# adjusted_p_values <- p.adjust(p_values, method = "holm")
# 
# # Paired differences
# remax_complete_norm <- unlist(raw_signature_performance[5, 3:102])
# remax_reduced_norm   <- unlist(normalized_signature_performance[5, 3:102])
# 
# remax_diff_norm <- remax_complete_norm - remax_reduced_norm
# 
# mae_complete_norm <- unlist(raw_signature_performance[2, 3:102])
# mae_reduced_norm   <- unlist(normalized_signature_performance[2, 3:102])
# 
# mae_diff_norm <- mae_complete_norm - mae_reduced_norm
# 
# sum(remax_diff_norm > 0)
# sum(remax_diff_norm < 0)
# sum(remax_diff_norm == 0)
# 
# sum(mae_diff_norm > 0)
# sum(mae_diff_norm < 0)
# sum(mae_diff_norm == 0)
# 
# # raw vs reduced signature ----
# # RMSE
# rmse_ttest_red <- t.test(unlist(raw_signature_performance[1,3:102]), 
#                      unlist(reduced_signature_performance[1,3:102]), paired = TRUE)
# # MAE
# mae_ttest_red <- t.test(unlist(raw_signature_performance[2,3:102]), 
#                     unlist(reduced_signature_performance[2,3:102]), paired = TRUE)
# # R2
# r2_ttest_red <- t.test(unlist(raw_signature_performance[3,3:102]), 
#                    unlist(reduced_signature_performance[3,3:102]), paired = TRUE)
# # MRE
# mre_ttest_red <- t.test(unlist(raw_signature_performance[4,3:102]), 
#                     unlist(reduced_signature_performance[4,3:102]), paired = TRUE)
# # REmax
# remax_ttest_red <- t.test(unlist(raw_signature_performance[5,3:102]), 
#                       unlist(reduced_signature_performance[5,3:102]), paired = TRUE)
# 
# rmse_ttest_red
# mae_ttest_red
# r2_ttest_red
# mre_ttest_red
# remax_ttest_red
# 
# data.frame(
#   Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
#   Whole = c(0.898, 0.780, 0.956, 14.724, 41.219),
#   Reduced = c(0.895, 0.797, 0.956, 14.293, 33.284),
#   Difference = c(
#     rmse_ttest_red$estimate[[1]],
#     mae_ttest_red$estimate[[1]],
#     r2_ttest_red$estimate[[1]],
#     mre_ttest_red$estimate[[1]],
#     remax_ttest_red$estimate[[1]]
#   ),
#   CI_lower = c(
#     rmse_ttest_red$conf.int[[1]],
#     mae_ttest_red$conf.int[[1]],
#     r2_ttest_red$conf.int[[1]],
#     mre_ttest_red$conf.int[[1]],
#     remax_ttest_red$conf.int[[1]]
#   ),
#   CI_upper = c(
#     rmse_ttest_red$conf.int[[2]],
#     mae_ttest_red$conf.int[[2]],
#     r2_ttest_red$conf.int[[2]],
#     mre_ttest_red$conf.int[[2]],
#     remax_ttest_red$conf.int[[2]]
#   ),
#   p_value = c(
#     rmse_ttest_red$p.value,
#     mae_ttest_red$p.value,
#     r2_ttest_red$p.value,
#     mre_ttest_red$p.value,
#     remax_ttest_red$p.value
#   )
# )
# 
# p_values <- c(rmse_ttest_red$p.value, mae_ttest_red$p.value, r2_ttest_red$p.value,
#               mre_ttest_red$p.value, remax_ttest_red$p.value)
# adjusted_p_values <- p.adjust(p_values, method = "holm")
# 
# # Paired differences
# remax_complete_red <- unlist(raw_signature_performance[5, 3:102])
# remax_reduced_red   <- unlist(reduced_signature_performance[5, 3:102])
# 
# remax_diff_red <- remax_complete_red - remax_reduced_red
# 
# mae_complete_red <- unlist(raw_signature_performance[2, 3:102])
# mae_reduced_red   <- unlist(reduced_signature_performance[2, 3:102])
# 
# mae_diff_red <- mae_complete_red - mae_reduced_red
# 
# sum(remax_diff_red > 0)
# sum(remax_diff_red < 0)
# sum(remax_diff_red == 0)
# 
# sum(mae_diff_red > 0)
# sum(mae_diff_red < 0)
# sum(mae_diff_red == 0)
# 
# beep(1)
# 
# 

# 
# raw_signature_performance
# shapiro.test(as.numeric(raw_signature_performance[1,3:102]))
# qqnorm(as.numeric(raw_signature_performance[5,3:102]),  main = "Normal Q-Q Plot", col = "blue", pch = 19)
# qqline(as.numeric(raw_signature_performance[5,3:102]),  col = "red", lwd = 2)
# t.test(
#   as.numeric(raw_signature_performance[5, 3:102]),
#   mu = 28.687
# )

gse <- "63952"
wang <- paste0("0_RT_GSE", gse, "_wang.csv")
wang <- read.csv(wang)
li <- paste0("0_RT_GSE", gse, "_li.csv")
li <- read.csv(li)
zhao <- paste0("0_RT_GSE", gse, "_zhao.csv")
zhao <- read.csv(zhao)
wang_time <- paste0("0_RT_GSE", gse, "_wang_time.csv")
wang_time <- read.csv(wang_time)
li_time <- paste0("0_RT_GSE", gse, "_li_time.csv")
li_time <- read.csv(li_time)
zhao_time <- paste0("0_RT_GSE", gse, "_zhao_time.csv")
zhao_time <- read.csv(zhao_time)

## NO_TIME vs TIME
#RMSE
wang_v_time_RMSE <- t.test(unlist(wang[1, 3:102]), unlist(wang_time[1, 3:102]), paired = TRUE)
li_v_time_RMSE <- t.test(unlist(li[1, 3:102]), unlist(li_time[1, 3:102]), paired = TRUE)
zhao_v_time_RMSE <- t.test(unlist(zhao[1, 3:102]), unlist(zhao_time[1, 3:102]), paired = TRUE)

#MAE
wang_v_time_MAE <- t.test(unlist(wang[2, 3:102]), unlist(wang_time[2, 3:102]), paired = TRUE)
li_v_time_MAE <- t.test(unlist(li[2, 3:102]), unlist(li_time[2, 3:102]), paired = TRUE)
zhao_v_time_MAE <- t.test(unlist(zhao[2, 3:102]), unlist(zhao_time[2, 3:102]), paired = TRUE)

#R2
wang_v_time_R2 <- t.test(unlist(wang[3, 3:102]), unlist(wang_time[3, 3:102]), paired = TRUE)
li_v_time_R2 <- t.test(unlist(li[3, 3:102]), unlist(li_time[3, 3:102]), paired = TRUE)
zhao_v_time_R2 <- t.test(unlist(zhao[3, 3:102]), unlist(zhao_time[3, 3:102]), paired = TRUE)

#MRE
wang_v_time_MRE <- t.test(unlist(wang[4, 3:102]), unlist(wang_time[4, 3:102]), paired = TRUE)
li_v_time_MRE <- t.test(unlist(li[4, 3:102]), unlist(li_time[4, 3:102]), paired = TRUE)
zhao_v_time_MRE <- t.test(unlist(zhao[4, 3:102]), unlist(zhao_time[4, 3:102]), paired = TRUE)

#REMAX
wang_v_time_REMAX <- t.test(unlist(wang[5, 3:102]), unlist(wang_time[5, 3:102]), paired = TRUE)
li_v_time_REMAX <- t.test(unlist(li[5, 3:102]), unlist(li_time[5, 3:102]), paired = TRUE)
zhao_v_time_REMAX <- t.test(unlist(zhao[5, 3:102]), unlist(zhao_time[5, 3:102]), paired = TRUE)

# Wang table
p_values <- c(wang_v_time_RMSE$p.value, wang_v_time_MAE$p.value, wang_v_time_R2$p.value,
              wang_v_time_MRE$p.value, wang_v_time_REMAX$p.value)
adjusted_p_values <- p.adjust(p_values, method = "holm")
adjusted_p_values

wang_vs_time <- tibble(
  Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
  wang = c(mean(unlist(wang[1, 3:102])),
           mean(unlist(wang[2, 3:102])),
           mean(unlist(wang[3, 3:102])),
           mean(unlist(wang[4, 3:102])),
           mean(unlist(wang[5, 3:102]))),
  wang_time = c(mean(unlist(wang_time[1, 3:102])),
                mean(unlist(wang_time[2, 3:102])),
                mean(unlist(wang_time[3, 3:102])),
                mean(unlist(wang_time[4, 3:102])),
                mean(unlist(wang_time[5, 3:102]))),
  Difference = c(
    wang_v_time_RMSE$estimate[[1]],
    wang_v_time_MAE$estimate[[1]],
    wang_v_time_R2$estimate[[1]],
    wang_v_time_MRE$estimate[[1]],
    wang_v_time_REMAX$estimate[[1]]
  ),
  CI_lower = c(
    wang_v_time_RMSE$conf.int[[1]],
    wang_v_time_MAE$conf.int[[1]],
    wang_v_time_R2$conf.int[[1]],
    wang_v_time_MRE$conf.int[[1]],
    wang_v_time_REMAX$conf.int[[1]]
  ),
  CI_upper = c(
    wang_v_time_RMSE$conf.int[[2]],
    wang_v_time_MAE$conf.int[[2]],
    wang_v_time_R2$conf.int[[2]],
    wang_v_time_MRE$conf.int[[2]],
    wang_v_time_REMAX$conf.int[[2]]
  ),
  p_value = adjusted_p_values[1:5]
) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

# li table
p_values <- c(li_v_time_RMSE$p.value, li_v_time_MAE$p.value, li_v_time_R2$p.value,
              li_v_time_MRE$p.value, li_v_time_REMAX$p.value)
adjusted_p_values <- p.adjust(p_values, method = "holm")
adjusted_p_values

li_vs_time <- tibble(
  Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
  li = c(mean(unlist(li[1, 3:102])),
         mean(unlist(li[2, 3:102])),
         mean(unlist(li[3, 3:102])),
         mean(unlist(li[4, 3:102])),
         mean(unlist(li[5, 3:102]))),
  li_time = c(mean(unlist(li_time[1, 3:102])),
              mean(unlist(li_time[2, 3:102])),
              mean(unlist(li_time[3, 3:102])),
              mean(unlist(li_time[4, 3:102])),
              mean(unlist(li_time[5, 3:102]))),
  Difference = c(
    li_v_time_RMSE$estimate[[1]],
    li_v_time_MAE$estimate[[1]],
    li_v_time_R2$estimate[[1]],
    li_v_time_MRE$estimate[[1]],
    li_v_time_REMAX$estimate[[1]]
  ),
  CI_lower = c(
    li_v_time_RMSE$conf.int[[1]],
    li_v_time_MAE$conf.int[[1]],
    li_v_time_R2$conf.int[[1]],
    li_v_time_MRE$conf.int[[1]],
    li_v_time_REMAX$conf.int[[1]]
  ),
  CI_upper = c(
    li_v_time_RMSE$conf.int[[2]],
    li_v_time_MAE$conf.int[[2]],
    li_v_time_R2$conf.int[[2]],
    li_v_time_MRE$conf.int[[2]],
    li_v_time_REMAX$conf.int[[2]]
  ),
  p_value = adjusted_p_values[1:5]
) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

# zhao table
p_values <- c(zhao_v_time_RMSE$p.value, zhao_v_time_MAE$p.value, zhao_v_time_R2$p.value,
              zhao_v_time_MRE$p.value, zhao_v_time_REMAX$p.value)
adjusted_p_values <- p.adjust(p_values, method = "holm")
adjusted_p_values

zhao_vs_time <- tibble(
  Metric = c("RMSE", "MAE", "R2", "MRE", "REmax"),
  zhao = c(mean(unlist(zhao[1, 3:102])),
           mean(unlist(zhao[2, 3:102])),
           mean(unlist(zhao[3, 3:102])),
           mean(unlist(zhao[4, 3:102])),
           mean(unlist(zhao[5, 3:102]))),
  zhao_time = c(mean(unlist(zhao_time[1, 3:102])),
                mean(unlist(zhao_time[2, 3:102])),
                mean(unlist(zhao_time[3, 3:102])),
                mean(unlist(zhao_time[4, 3:102])),
                mean(unlist(zhao_time[5, 3:102]))),
  Difference = c(
    zhao_v_time_RMSE$estimate[[1]],
    zhao_v_time_MAE$estimate[[1]],
    zhao_v_time_R2$estimate[[1]],
    zhao_v_time_MRE$estimate[[1]],
    zhao_v_time_REMAX$estimate[[1]]
  ),
  CI_lower = c(
    zhao_v_time_RMSE$conf.int[[1]],
    zhao_v_time_MAE$conf.int[[1]],
    zhao_v_time_R2$conf.int[[1]],
    zhao_v_time_MRE$conf.int[[1]],
    zhao_v_time_REMAX$conf.int[[1]]
  ),
  CI_upper = c(
    zhao_v_time_RMSE$conf.int[[2]],
    zhao_v_time_MAE$conf.int[[2]],
    zhao_v_time_R2$conf.int[[2]],
    zhao_v_time_MRE$conf.int[[2]],
    zhao_v_time_REMAX$conf.int[[2]]
  ),
  p_value = adjusted_p_values[1:5]
) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

as.data.frame(wang_vs_time)
as.data.frame(li_vs_time)
as.data.frame(zhao_vs_time)

# Function to compare two signatures
compare_signatures <- function(data1, data2, name1, name2) {
  
  metrics <- c("RMSE", "MAE", "R2", "MRE", "REmax")
  
  # Paired t-tests
  tests <- lapply(1:5, function(i) {
    t.test(
      unlist(data1[i, 3:102]),
      unlist(data2[i, 3:102]),
      paired = TRUE
    )
  })
  
  # Holm-adjusted p-values across five metrics
  adjusted_p_values <- p.adjust(
    sapply(tests, function(x) x$p.value),
    method = "holm"
  )
  
  # Create comparison table
  tibble(
    Metric = metrics,
    
    !!name1 := sapply(1:5, function(i) {
      mean(unlist(data1[i, 3:102]))
    }),
    
    !!name2 := sapply(1:5, function(i) {
      mean(unlist(data2[i, 3:102]))
    }),
    
    Difference = sapply(tests, function(x) x$estimate[[1]]),
    
    CI_lower = sapply(tests, function(x) x$conf.int[[1]]),
    
    CI_upper = sapply(tests, function(x) x$conf.int[[2]]),
    
    p_value = adjusted_p_values
    
  ) %>%
    mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))
}


# Wang vs Li
wang_vs_li <- compare_signatures(
  wang, li, "wang", "li"
)

# Wang vs Zhao
wang_vs_zhao <- compare_signatures(
  wang, zhao, "wang", "zhao"
)

# Li vs Zhao
li_vs_zhao <- compare_signatures(
  li, zhao, "li", "zhao"
)


# Display tables
as.data.frame(wang_vs_li)
as.data.frame(wang_vs_zhao)
as.data.frame(li_vs_zhao)
