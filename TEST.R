df <- read.csv("full_metadata.csv")
library(dplyr)
library(tibble)

# table <- df %>%
#   group_by(dataset_ID) %>%
#   summarise(
#     species = paste(unique(organism), collapse = ", "),
#     time = paste(unique(time_hours), collapse = ", "),
#     radiation = paste(unique(radiation_type), collapse = ", "),
#     product = paste(unique(tissue), collapse = ", "),
#     dose_range = max(dose_Gy),
#     exp_set = paste(unique(experimental_setting), collapse = ", ")
#   )

#ortholog maps

list_se <- readRDS("list_se.rds")
final_matrices <- readRDS("final_matrices.rds")
changes <- readRDS("changes_data_cleaning.rds")

BiocManager::install("orthogene")
a <- final_matrices[[1]]
a <- as.data.frame(t(a))
a$ortologue <- rownames(a)

newgenes <- orthogene::convert_orthologs(a, "rownames", "rownames", 
                                           FALSE, "mmusculus", "hsapiens", 
                                           non121_strategy = "drop_both_species")
