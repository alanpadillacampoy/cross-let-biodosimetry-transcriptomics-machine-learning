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
#Replaces the names of the probes with the corresponding genes being tested.----
#Skips Data set 17 since platform GPL6244 cannot be generalized in the same way.
for (i in 1:35){
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
  
  #Deals with the different structure SE 17 has
  if (i == 17) {
    # 1. Split each row by the triple slashes '///' to isolate the sections
    sections_list <- strsplit(gpl_table$gene_assignment, "///", fixed = TRUE)
    
    # 2. Extract targets and apply the uniqueness condition per row
    extracted_targets <- sapply(sections_list, function(row_sections) {
      
      # Trim spaces from sections
      row_sections <- trimws(row_sections)
      
      # Split each individual section by the double slash '//'
      sub_pieces <- strsplit(row_sections, "//", fixed = TRUE)
      
      # Grab the 2nd element of each section
      targets <- sapply(sub_pieces, function(piece) {
        if (length(piece) >= 2) {
          return(trimws(piece[2]))
        } else {
          return(NA)
        }
      })
      
      # Remove any NAs generated during extraction before checking uniqueness
      targets <- na.omit(targets)
      
      # Apply your condition:
      unique_targets <- unique(targets)
      
      if (length(unique_targets) == 1) {
        return(unique_targets)       
      } else {
        return(NA)                   
      }
    })
    
    # 3. Add it straight back to your table as a new column
    gpl_table$GENE_SYMBOL <- extracted_targets
  }
  
  #Selects the possible columns where the probes and symbols may be:
  #PROBES may be in the ID column or in the NAME column
  #The symbol is always a combination of GENE SYMBOL so we search for SYMBOL
  gene_column <- colnames(gpl_table %>% dplyr::select(contains("symbol")))
  id_column <- ifelse(gpl_table[colnames(gpl_table)[1]][1,1] == 1, "SPOT_ID", "ID")
  

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
  
  #Deletes all NAs or empty spaces in the Gene Symbol column to keep only single 
  #named genes
  rows_to_keep <- !is.na(expression_matrix_new[[gene_column]]) &
    expression_matrix_new[[gene_column]] != "" &
    expression_matrix_new[[gene_column]] != "-" &
    expression_matrix_new[[gene_column]] != "---" &
    !grepl("///", expression_matrix_new[[gene_column]])
  expression_matrix_new <- expression_matrix_new[rows_to_keep, ]
  
  expression_matrix_new <- expression_matrix_new %>% dplyr::rename(Gene = all_of(gene_column))
  expression_matrix_new <- expression_matrix_new %>% dplyr::arrange(Gene)
  
  #Assigns the subset number to the new expression matrix
  matrix_name <- paste0("expression_matrix_new_", i)
  assign(matrix_name, expression_matrix_new)
}

for (i in 1:35){
  print(i)
  a <- get(paste0("expression_matrix_new_", i))
  print(length(unique(na.omit(a$Gene))))
  print(head(a[,c(1,2,3)]))
}

for (i in 1:35) {
  # 1. Print header for the current matrix
  cat("\n=========================================\n")
  cat("  ANALYZING: expression_matrix_new_", i, "\n", sep = "")
  cat("=========================================\n")
  
  # 2. Fetch the matrix dynamically
  a <- get(paste0("expression_matrix_new_", i))
  
  # 3. Isolate the Gene column and drop any missing values (NAs)
  gene_vector <- na.omit(a$Gene)
  
  # 4. Calculate the core metrics
  total_genes  <- length(gene_vector)
  unique_genes <- length(unique(gene_vector))
  
  # Prevent division by zero if a matrix happens to have an empty Gene column
  if (total_genes > 0) {
    pct_unique <- (unique_genes / total_genes) * 100
  } else {
    pct_unique <- 0
  }
  
  # 5. Count how many times each unique gene appears
  gene_counts <- table(gene_vector)
  
  # 6. Print the summary report
  cat("Total Genes (Rows):       ", total_genes, "\n")
  cat("Unique Genes:             ", unique_genes, "\n")
  cat("Percentage of Unique Data:", round(pct_unique, 2), "%\n\n")
  
  # 7. Print the most frequent genes
  cat("Top 6 most frequent genes:\n")
  print(head(sort(gene_counts, decreasing = TRUE)))
  
  # 8. Your original preview of the data frame structure
  cat("\nMatrix preview (First 3 columns):\n")
  print(head(a[, c(1, 2, 3)]))
}


for (i in 1:35) {
  print(i)
  a <- get(paste0("expression_matrix_new_", i))
  max <- max(a[,3:ncol(a)], na.rm = TRUE)
  print(max)
}
