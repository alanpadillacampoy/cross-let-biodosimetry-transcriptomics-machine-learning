##This script is to check the quality of the genomic data in the SummarizedExps
#libraries ----
library(DoReMiTra)
library(limma)
library(SummarizedExperiment)
library(readr)
library(stringr)
library(dplyr)
library(hgu133plus2.db)

#
list_se <- DoReMiTra::list_DoReMiTra_datasets()

se1 <- get_DoReMiTra_data(list_se$Dataset[1])
expr1 <- as.data.frame(assay(se1, "exprs"))
expr <- tibble::rownames_to_column(expr1, "probe_id")
for (i in 2:35){
  se2 <- get_DoReMiTra_data(list_se$Dataset[i])
  expr2 <- as.data.frame(assay(se2, "exprs"))
  expr2_rows <- tibble::rownames_to_column(expr2, "probe_id")
  expr <- left_join(expr, expr2_rows, by = "probe_id")
}

write.csv(expr, "expr.csv")
