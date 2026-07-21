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
library(tibble)
library(WGCNA)

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
  
  #Extracts the row number of the expression matrix
  rows <- nrow(expression_matrix)
  rows_name <- paste0("rows_original_", i)
  assign(rows_name, rows)
  
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
#Checks for LOG2 transform and if it isnt then transforms it
for (i in 1:35) {
  print(i)
  a <- get(paste0("expression_matrix_new_", i))
  maximum <- max(a[, 3:ncol(a)], na.rm = TRUE)
  print(maximum)
  if (maximum > 30){
    b <- cbind(a[, 1:2], log2(a[, 3:ncol(a)] + 1))
    maximum <- max(b[, 3:ncol(b)], na.rm = TRUE)
    print(maximum)
    #Assigns the subset number to the new expression matrix
    matrix_name <- paste0("expression_matrix_new_", i)
    print(matrix_name)
    assign(matrix_name, b)
  }
}
#TESTING
for (i in 1:35){
  print(i)
  a <- get(paste0("expression_matrix_new_", i))
  print(length(unique(na.omit(a$Gene))))
  print(head(a[,c(1,2,3)]))
}

#DELETE PROBES WITH 15% OR MORE NAs ----
for (i in 1:35) {
  matrix <- get(paste0("expression_matrix_new_", i))
  end_col <- ncol(matrix) - 3
  matrix$NAs <- 
    rowSums(is.na(matrix[,3:end_col])) / 
    end_col * 100
  matrix <- matrix %>%
    dplyr::relocate(NAs, .before = 3)
  matrix_cleaned <- subset.data.frame(matrix, NAs <= 15)
  matrix_name <- paste0("expression_matrix_cleaned_", i)
  assign(matrix_name, matrix_cleaned)
}
#TRACK CHANGES ----
probe_tracker <- NULL  # Run this once before the loop starts
for (i in 1:35) {
  orig <- get(paste0("rows_original_", i))
  n_new <- nrow(get(paste0("expression_matrix_new_", i)))
  n_clean <- nrow(get(paste0("expression_matrix_cleaned_", i)))
  
  probe_tracker <- rbind(probe_tracker, as.data.frame(tibble::tibble(
    dataset = i, original_probes = orig, cleaned_na_genes = n_new,
    difference_1 = original_probes - cleaned_na_genes,
    cleaned_15_percent = n_clean, difference_2 = cleaned_na_genes - cleaned_15_percent
  )))
}

#COLLAPSE PROBES THAT MAP THE SAME GENE INTO ONE BASED ON THE HIGHEST MEAN EXPRESSION ----
for (i in 1:35) {
  print(i)
  cleaned_matrix <- get(paste0("expression_matrix_cleaned_", 1))
  collapsed_data <- WGCNA::collapseRows(datET = cleaned_matrix[,4:ncol(cleaned_matrix)],
                                        rowGroup = cleaned_matrix$Gene,
                                        rowID = rownames(cleaned_matrix),
                                        method = "MaxMean")
  expression_matrix_collapsed <- as.data.frame(collapsed_data$datETcollapsed)
  #expression_matrix_collapsed$Gene <- rownames(expression_matrix_collapsed)
  #expression_matrix_collapsed <- expression_matrix_collapsed %>%
  #  dplyr::relocate(Gene, .before = 1)
  matrix_name <- paste0("expression_matrix_collapsed_", i)
  assign(matrix_name, expression_matrix_collapsed)
}

#Normalizes the data to Z score ----
for (i in 1:35) {
  print(i)
  collapsed_matrix <- get(paste0("expression_matrix_collapsed_", i))
  
  #Replaces NAs for that row mean, effectively turning them to 0 when Z scoring
  row_means <- rowMeans(collapsed_matrix, na.rm = TRUE)
  na_indices <- which(is.na(collapsed_matrix), arr.ind = TRUE)
  collapsed_matrix[na_indices] <- row_means[na_indices[, 1]]
  
  #Checks the variance across genes
  gene_sds <- apply(collapsed_matrix, 1, sd)
  mat_filtered <- collapsed_matrix[gene_sds > 0, ]
  
  #Z scoring
  final_data <- scale(t(mat_filtered))
  print(sum(colMeans(final_data)))
  print(summary(apply(final_data, 2, sd)))
  
  matrix_name <- paste0("expression_matrix_final_", i)
  assign(matrix_name, final_data)
}
