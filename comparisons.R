## this is to do the comparisons between different experimental summaries

normal_signature_performance_summary <- read.csv("normal_signature_performance_summary.csv")
reduced_signature_performance_summary <- read.csv("reduced_signature_performance_summary.csv")
normal_signature_performance <- read.csv("normal_signature_performance.csv")
reduced_signature_performance <- read.csv("reduced_signature_performance.csv")

# Paired t-tests ----
# RMSE
t.test(unlist(normal_signature_performance[1,3:12]), 
       unlist(reduced_signature_performance[1,3:12]), paired = TRUE)
t.test(normal_signature_performance_summary$MAE[2:11], 
       actb_normalized$MAE[2:11], paired = TRUE)
t.test(normal_signature_performance_summary$RSquared[2:11], 
       actb_normalized$RSquared[2:11], paired = TRUE)
t.test(normal_signature_performance_summary$MRE[2:11], 
       actb_normalized$MRE[2:11], paired = TRUE)
t.test(normal_signature_performance_summary$REmax[2:11], 
       actb_normalized$REmax[2:11], paired = TRUE)
