## this is to do the comparisons between different experimental summaries

normal_signature_performance_summary <- read.csv("normal_signature_performance_summary.csv")
reduced_signature_performance_summary <- read.csv("reduced_signature_performance_summary.csv")
normal_signature_performance <- read.csv("normal_signature_performance.csv")
reduced_signature_performance <- read.csv("reduced_signature_performance.csv")

# Paired t-tests ----
# RMSE
rmse_ttest <- t.test(unlist(normal_signature_performance[1,3:12]), 
       unlist(reduced_signature_performance[1,3:12]), paired = TRUE)
# MAE
mae_ttest <- t.test(unlist(normal_signature_performance[2,3:12]), 
                     unlist(reduced_signature_performance[2,3:12]), paired = TRUE)
# R2
r2_ttest <- t.test(unlist(normal_signature_performance[3,3:12]), 
                     unlist(reduced_signature_performance[3,3:12]), paired = TRUE)
# MRE
mre_ttest <- t.test(unlist(normal_signature_performance[4,3:12]), 
                     unlist(reduced_signature_performance[4,3:12]), paired = TRUE)
# REmax
remax_ttest <- t.test(unlist(normal_signature_performance[5,3:12]), 
                     unlist(reduced_signature_performance[5,3:12]), paired = TRUE)