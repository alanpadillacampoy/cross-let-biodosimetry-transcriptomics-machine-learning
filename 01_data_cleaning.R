##This script is to check the quality of the genomic data in the SummarizedExps
#libraries ----
library(DoReMiTra)
library(SummarizedExperiment)
library(stringr)
library(dplyr)
library(BiocManager)
library(GEOquery)
library(tibble)
library(WGCNA)

#Functions
source("00_functions.R")

#Analysis ----


#Extracts the platform from the data set, and launches a GEO query to get the
#platform's metadata
unique_platforms <- unique_gpl_platforms(list_se)

geo_metadata <- GEO_query_list(unique_platforms)

#Standardizes the platform present in dataset 17
geo_metadata <- correct_seventeen(list_se, geo_metadata)

se <- get_expression_matrices(list_se)

gene_column <- find_gene_column(geo_metadata)
id_column <- find_id_column(geo_metadata)

annotated_expression_matrices <- 
  annotate_expression_by_rownames(se, geo_metadata, id_column, gene_column)

log_checked_matrices <- check_log2_transform(annotated_expression_matrices)

complete_matrices <- delete_NAs(annotated_expression_matrices)

collapsed_matrices <- collapse_probes(complete_matrices)

final_matrices <- z_score_matrices(collapsed_matrices)

changes <- track_changes(se, complete_matrices, final_matrices)

# This instruction deletes everything but the final matrices, ----
# which will be used for further analysis. To review each of the individual steps 
# comment the function

clean_environment_1()

saveRDS(final_matrices, file = "final_matrices.rds")
saveRDS(list_se, file = "list_se.rds")
saveRDS(changes, file = "changes_data_cleaning.rds")

