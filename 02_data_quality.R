##This script is to check the quality of the genomic data in the SummarizedExps
#libraries ----
library(DoReMiTra)
library(limma)
library(SummarizedExperiment)
library(readr)
library(stringr)
library(dplyr)
library(hgu133plus2.db)
library(BiocManager)
#BiocManager::install("biomaRt")
library(biomaRt)
library(GEOquery)
#
list_se <- DoReMiTra::list_DoReMiTra_datasets()

platform <- strsplit(list_se$Dataset[1], "_")[[1]][6]
gpl <- GEOquery::getGEO(platform)
gpl_table <- GEOquery::Table(gpl)
head(gpl_table)
se1 <- DoReMiTra::get_DoReMiTra_data(list_se$Dataset[1])
expression_matrix <- as.data.frame(SummarizedExperiment::assay(se1, "exprs"))
probes <- as.data.frame(rownames(expression_matrix))
colnames(probes) <- "probes"
clean_probes <- probes[!grepl("\\(", probes$probes), , drop = FALSE]
head(clean_probes)
new_expression_matrix <- merge(clean_probes, gpl_table[, c("ID", "GENE_SYMBOL")], 
                               by.x = "probes",
                               by.y = "ID",
                               all.x = TRUE
                               ) 
expression_matrix <- merge(expression_matrix, new_expression_matrix,
                           by.x = "row.names",
                           by.y = "probes",
                           all.x = TRUE)

head(expression_matrix)
expression_matrix <- expression_matrix %>% dplyr::relocate("GENE_SYMBOL", .before = 2)
# se1 <- get_DoReMiTra_data(list_se$Dataset[1])
# expr1 <- as.data.frame(assay(se1, "exprs"))
# expr <- tibble::rownames_to_column(expr1, "probe_id")
# for (i in 2:35){
#   se2 <- get_DoReMiTra_data(list_se$Dataset[i])
#   expr2 <- as.data.frame(assay(se2, "exprs"))
#   expr2_rows <- tibble::rownames_to_column(expr2, "probe_id")
#   expr <- left_join(expr, expr2_rows, by = "probe_id")
# }
# 
# write.csv(expr, "expr.csv")
