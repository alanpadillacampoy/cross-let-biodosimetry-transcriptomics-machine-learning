##This script is to check the quality of the genomic data in the SummarizedExps
#libraries ----
library(DoReMiTra)
library(limma)
library(SummarizedExperiment)
library(readr)
library(stringr)
library(dplyr)
library(BiocManager)
library(biomaRt)
library(GEOquery)

list_se <- DoReMiTra::list_DoReMiTra_datasets()

#GPL6244 = SE 17 skip it

for (i in (1:35)[-17]){
  print(i)
  platform <- strsplit(list_se$Dataset[i], "_")[[1]][6]
  gpl <- GEOquery::getGEO(platform)
  gpl_table <- GEOquery::Table(gpl)
  se1 <- DoReMiTra::get_DoReMiTra_data(list_se$Dataset[i])
}


expression_matrix <- as.data.frame(SummarizedExperiment::assay(se1, "exprs"))
probes <- as.data.frame(rownames(expression_matrix))
colnames(probes) <- "probes"
clean_probes <- probes[!grepl("\\(", probes$probes), , drop = FALSE]
head(clean_probes)
new_expression_matrix <- merge(clean_probes, 
                               gpl_table[, c("ID", "GENE_SYMBOL")], 
                               by.x = "probes",
                               by.y = "ID",
                               all.x = TRUE
                               ) 
expression_matrix_new <- merge(expression_matrix, new_expression_matrix,
                           by.x = "row.names",
                           by.y = "probes",
                           all.x = TRUE)

head(expression_matrix)
expression_matrix_new <- expression_matrix_new %>% 
  dplyr::relocate(c("GENE" , "GENE_SYMBOL"), .before = 2)

expression_matrix_no_NA <- expression_matrix_new[!is.na(expression_matrix_new$GENE),]

View(expression_matrix_no_NA)
View(expression_matrix_new)

# platform <- c()
# for (i in 1:35){
#   platform <- c(platform, strsplit(list_se$Dataset[i], "_")[[1]][6])
# }
# platform <- levels(as.factor(platform))
# for (i in 1:length(platform)){
#   gpl <- GEOquery::getGEO(platform[i])
#   gpl_table <- GEOquery::Table(gpl)
#   print(platform[i])
#   print(colnames(gpl_table))
# }

se17 <- assay(get_DoReMiTra_data("SE_Rouchka_2015_ExVivo_GSE64375_GPL6244"), "exprs")
gpl <- GEOquery::getGEO("GPL6244")
gpl_table <- GEOquery::Table(gpl)
head(gpl_table)
probe_ids <- rownames(se17)

# Match GPL table rows to your expression data's probe order
matched <- gpl_table[match(probe_ids, gpl_table$ID), ]

gene_assignment_col <- matched$gene_assignment

extract_symbol_unique <- function(x) {
  if (is.na(x) || x == "---") return(NA)
  entries <- strsplit(x, " /// ")[[1]]
  symbols <- sapply(entries, function(e) {
    parts <- strsplit(e, " // ")[[1]]
    if (length(parts) >= 2) trimws(parts[2]) else NA
  })
  unique_symbols <- unique(na.omit(symbols))
  if (length(unique_symbols) == 1) unique_symbols else NA
}

extract_symbol_all <- function(x) {
  if (is.na(x) || x == "---") return(NA)
  entries <- strsplit(x, " /// ")[[1]]
  symbols <- sapply(entries, function(e) {
    parts <- strsplit(e, " // ")[[1]]
    if (length(parts) >= 2) trimws(parts[2]) else NA
  })
  unique_symbols <- unique(na.omit(symbols))
  if (length(unique_symbols) == 0) NA else paste(unique_symbols, collapse = ";")
}

gene_symbol_unique <- sapply(gene_assignment_col, extract_symbol_unique)
gene_symbol_all <- sapply(gene_assignment_col, extract_symbol_all)

# Attach to your expression data as a data frame with the two new columns
expr_annot <- data.frame(probe_id = probe_ids,
                         gene_symbol_unique = gene_symbol_unique,
                         gene_symbol_all = gene_symbol_all)



#write.csv(expression_matrix, "original.csv")
#write.csv(expression_matrix_new, "annotated.csv")

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
