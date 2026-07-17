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
list_problems <- c(4, 5, 14, 20, 33)
list_super_problems <- c(16, 18)

#GPL6244 = SE 17 skip it
#Replaces the names of the probes with the corresponding genes being tested.----
#Skips Data set 17 since platform GPL6244 cannot be generalized in the same way.
for (i in (1:35)[-17]){
 print(i)
 #Extracts the platform from the data set, and launches a GEO query to get the
 #platform's metadata
 platform <- strsplit(list_se$Dataset[i], "_")[[1]][6]
 gpl <- GEOquery::getGEO(platform)
 gpl_table <- GEOquery::Table(gpl)

 #Extracts the expression matrix from the data set
 se1 <- DoReMiTra::get_DoReMiTra_data(list_se$Dataset[i])
 expression_matrix <- as.data.frame(SummarizedExperiment::assay(se1, "exprs"))
 probes <- as.data.frame(SummarizedExperiment::rownames(expression_matrix))
 colnames(probes) <- "probes"
 clean_probes <- probes[!grepl("\\(", probes$probes), , drop = FALSE]
 
 #Selects the possible columns where the probes and symbols may be:
 #PROBES may be in the ID column or in the NAME column
 #The symbol is always a combination of GENE SYMBOL so we search for SYMBOL
 gene_column <- colnames(gpl_table %>% dplyr::select(contains("symbol")))
 id_column <- ifelse(gpl_table[colnames(gpl_table)[1]][1,1] == 1, "NAME", "ID")
 
 #Annotates the probes with the corresponding genes
 annotated_probes <- SummarizedExperiment::merge(clean_probes,
                                                      gpl_table[, c(id_column, gene_column)],
                                                      by.x = "probes",
                                                      by.y = id_column,
                                                      all.x = TRUE)
 
 #Merges the annotated probes with the original expression matrix, and relocates
 #the Gene Symbol column to the second space
 expression_matrix_new <- SummarizedExperiment::merge(expression_matrix,
                                                      annotated_probes,
                                                      by.x = "row.names",
                                                      by.y = "probes",
                                                      all.x = TRUE)
 expression_matrix_new <- expression_matrix_new %>%
   dplyr::relocate(all_of(gene_column), .before = 2)

 #Deletes all NAs or empty spaces in the Gene Symbol column to keep only named genes
 rows_to_keep <- !is.na(expression_matrix_new[[gene_column]]) &
   expression_matrix_new[[gene_column]] != "" &
   expression_matrix_new[[gene_column]] != "-"
 expression_matrix_new <- expression_matrix_new[rows_to_keep, ]
 
 #Assigns the subset number to the new expression matrix
 matrix_name <- paste0("expression_matrix_new_", i)
 assign(matrix_name, expression_matrix_new)
}

for (i in (1:35)[-17]){
 print(i)
 a <- get(paste0("expression_matrix_new_", i))
print(head(a[,c(1,2,3)]))
}

##Deals with the SE17 problem ----

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


# # Ascending order (A to Z)
# df_sorted <- df %>% arrange(column_name)
# 
# # Descending order (Z to A)
# df_sorted <- df %>% arrange(desc(column_name))

