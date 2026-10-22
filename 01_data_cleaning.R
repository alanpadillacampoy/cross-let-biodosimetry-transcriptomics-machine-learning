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

# Signatures ----
wang_mouse_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "H4c3", "H4c9", "Ifit1", 
                                "Igfbp4", "LOC118567921", "Lrrc70", "Ms4a1", "Phlda3", 
                                "Ptprn", "Rps20", "Serpine2", "Thy1", "Actb")
#without the "problem" genes    
wang_mouse_reduced_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "Ifit1", "Igfbp4",
                                        "Ms4a1", "Phlda3", "Ptprn", "Rps20", "Serpine2", "Thy1")

wang_human_signature_genes <- c("CCNG1", "DGKA", "FZR1", "H4C3", "H4C9", "IFIT1", 
                                "IGFBP4", "LOC118567921", "LRRC70", "MS4A1", "PHLDA3", 
                                "PTPRN", "RPS20", "SERPINE2", "THY1", "ACTB")

wang_human_reduced_signature_genes <- c("CCNG1", "DGKA", "FZR1", "IFIT1", 
                                        "IGFBP4", "MS4A1", "PHLDA3", "PTPRN", 
                                        "RPS20", "SERPINE2", "THY1")

li_human_signature_genes <- c("CDKN1A", "BAX", "MDM2", "XPC", "PCNA", "FDXR", 
                              "GDF15", "DDB2", "TNFRSF10B", "PHPT1", "ASTN2", 
                              "RPS27L", "BBC3", "TNFSF4", "POLH", "CCNG1", 
                              "PPM1D", "GADD45A", "ZMAT3", "SESN1", "AEN", 
                              "TRIAP1")

li_mouse_signature_genes <- c("Cdkn1a", "Bax", "Mdm2", "Xpc", "Pcna", "Fdxr", 
                              "Gdf15", "Ddb2", "Tnfrsf10b", "Phpt1", "Astn2", 
                              "Rps27l", "Bbc3", "Tnfsf4", "Polh", "Ccng1", 
                              "Ppm1d", "Gadd45a", "Zmat3", "Sesn1", "Aen", 
                              "Triap1")

zhao_human_signature_genes <- c("DDB2", "CD8A", "TALDO1", "PCNA", "EIF4G2", 
                                "LCN2", "CDKN1A", "PRKCH", "ENO1", "PPM1D")

zhao_mouse_signature_genes <- c("Ddb2", "Cd8a", "Taldo1", "Pcna", "Eif4g2", 
                                "Lcn2", "Cdkn1a", "Prkch", "Eno1", "Ppm1d")

FDXR_human <- "FDXR"
FDXR_mouse <- "Fdxr"
DDB2_human <- "DDB2"
DDB2_mouse <- "Ddb2"
GADD45A_human <- "GADD45A"
GADD45A_mouse <- "Gadd45a"

human_signature_list <- list(
  wang_human = wang_human_reduced_signature_genes,
  li_human   = li_human_signature_genes,
  zhao_human = zhao_human_signature_genes,
  FDXR_human = FDXR_human,
  DDB2_human = DDB2_human,
  GADD45A_human = GADD45A_human
)
mouse_signture_list <- list(
  wang_mouse = wang_mouse_reduced_signature_genes,
  li_mouse   = li_mouse_signature_genes,
  zhao_mouse= zhao_mouse_signature_genes,
  FDXR_mouse = FDXR_mouse,
  DDB2_mouse = DDB2_mouse,
  GADD45A_mouse = GADD45A_mouse
)
gene_space_human <- purrr::reduce(human_signature_list, union)
gene_space_mouse <- purrr::reduce(mouse_signture_list, union)
gene_space_total <- union(gene_space_human, gene_space_mouse)

signature <- wang_mouse_signature_genes
size_signature <- length(signature)

#Analysis ----

#Extracts the platform from the data set, and launches a GEO query to get the
#platform's metadata
unique_platforms <- unique_gpl_platforms(list_se)

geo_metadata <- GEO_query_list(unique_platforms)

#Standardizes the platform present in dataset 17
geo_metadata <- correct_GSE64375(geo_metadata)

se <- get_expression_matrices(list_se)

gene_column <- find_gene_column(geo_metadata)
id_column <- find_id_column(geo_metadata)

annotated_expression_matrices <- 
  annotate_expression_by_rownames(se, geo_metadata, id_column, gene_column)

log_checked_matrices <- check_log2_transform(annotated_expression_matrices)
##
# # Temporary QC: identify datasets exceeding the log2 threshold
# for (dataset_name in names(annotated_expression_matrices)) {
#   expr_df <- annotated_expression_matrices[[dataset_name]]
#   maximum <- max(expr_df[, 3:ncol(expr_df)], na.rm = TRUE)
#   
#   if (maximum > 30) {
#     cat(dataset_name, "- maximum expression:", maximum, "\n")
#   }
# }
# ##
complete_matrices <- delete_NAs(log_checked_matrices)

filtered_matrices <- lapply(names(complete_matrices), function(i){
  x <- complete_matrices[[i]]
  x <- x %>% filter(gene %in% gene_space_total)
  return(x)
})

names(filtered_matrices) <- names(complete_matrices)

collapsed_matrices <- collapse_probes(filtered_matrices)

#scaled_matrices <- z_score_matrices(collapsed_matrices)

#final_matrices <- ortholog_correction(scaled_matrices, list_se)
final_matrices <- collapsed_matrices


#changes <- track_changes(se, complete_matrices, scaled_matrices, final_matrices)

# This instruction deletes everything but the final matrices, ----
# which will be used for further analysis. To review each of the individual steps 
# comment the function
# clean_environment_1()
 
#saveRDS(final_matrices, file = "final_matrices.rds")
#saveRDS(list_se, file = "list_se.rds")
#saveRDS(changes, file = "changes_data_cleaning.rds")
beep(3)