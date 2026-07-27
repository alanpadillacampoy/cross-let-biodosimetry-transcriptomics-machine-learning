#This is the master function document. All functions are stored here for 
#reference in the project

# Libraries ----
library(DoReMiTra)
library(SummarizedExperiment)
library(stringr)
library(dplyr)
library(BiocManager)
library(GEOquery)
library(tibble)
library(WGCNA)

# Start data ----
list_se <- DoReMiTra::get_all_DoReMiTra_datasets()
list_se[["SE_Salah_2025_ExVivo"]] <- NULL

# Data cleaning, quality and validation ----

##Extracts the platform from the data set and checks for uniqueness ----
unique_gpl_platforms <- function(list_se){
  platforms <- vector("list", length = 35)
  for (i in 1:35) {
    platforms[i] <- strsplit(names(list_se)[i], "_")[[1]][6]
  }
  unique_platforms <- unique(platforms)
  
  return(unique_platforms)
}

## Launches a GEO query to obtain the platform's metadata ----
GEO_query_list <- function(unique_platforms){
  
  platform_metadata <- vector("list", length = length(unique_platforms))
  names(platform_metadata) <- unique_platforms
  
  for (i in 1:length(unique_platforms)) {
    unique_platforms[i]
    geo_gpl <- GEOquery::getGEO(unique_platforms[i])
    platform_metadata[[i]] <- GEOquery::Table(geo_gpl)
  }
  return(platform_metadata)
}

## Extracts the expression matrix form the SE ----
get_expression_matrices <- function(list_se) {
  lapply(list_se, function(se) {
    as.data.frame(SummarizedExperiment::assay(se))
  })
}

## Corrects the issue with the platform metadata in the dataset 17 ----
correct_seventeen <- function(list_se, geo_metadata) {
  
  platform <- strsplit(names(list_se)[17], "_")[[1]][6]
  
  sections_list <- strsplit(geo_metadata[[platform]]$gene_assignment, "///", fixed = TRUE)
  
  extracted_targets <- sapply(sections_list, function(row_sections) {
    
    row_sections <- trimws(row_sections)
    sub_pieces <- strsplit(row_sections, "//", fixed = TRUE)
    
    targets <- sapply(sub_pieces, function(piece) {
      if (length(piece) >= 2) {
        return(trimws(piece[2]))
      } else {
        return(NA)
      }
    })
    
    targets <- targets[!is.na(targets)]
    unique_targets <- unique(targets)
    
    if (length(unique_targets) == 1) {
      return(unique_targets)       
    } else {
      return(NA_character_)         
    }
  }, USE.NAMES = FALSE)

  geo_metadata[[platform]]$GENE_SYMBOL <- extracted_targets
  return(geo_metadata)
}
######DELETE?????? ----
# join_probes <- function(se){
#   lapply(se, function(df_expression){
#     probes <- rownames(df_expression)
#     data.frame(probes = probes, stringsAsFactors = FALSE)
#   })
# }
## Selects the possible columns where the probes and symbols may be: ----
find_gene_column <- function(geo_metadata){

  #The symbol is always a combination of GENE SYMBOL so we search for SYMBOL
  gene_column <- lapply(geo_metadata, function(gene_column){
    colnames(gene_column %>% dplyr::select(contains("symbol")))
  }) 
  
  return(gene_column)
}  
find_id_column <- function(geo_metadata){
  
  #PROBES may be in the ID column or in the NAME column  
  id_column <- lapply(geo_metadata, function(id_column){
    ifelse(id_column[colnames(id_column)[1]][1,1] == 1, "SPOT_ID", "ID")
  }) 
  
  return(id_column)
} 

## Annotates the expression matrices with the genes, matching the probes ----
#and the genes they map for
annotate_expression_by_rownames <- function(se, geo_metadata, id_column, gene_column) {
  
  res <- lapply(names(se), function(dataset_name) {
    
    # 1. Get expression matrix as data.frame
    expr_df <- as.data.frame(se[[dataset_name]])
    
    # 2. Extract platform ID (e.g., "GPL11202")
    platform <- strsplit(dataset_name, "_")[[1]][6]
    platform_df <- geo_metadata[[platform]]
    
    # 3. Get target ID and Gene column names for this platform/dataset
    current_id_col   <- as.character(id_column[[platform]])[1]
    current_gene_col <- as.character(gene_column[[platform]])[1]
    
    # Fallback to numeric indexing if id/gene_column are indexed by position rather than platform string
    if (is.na(current_id_col)) {
      idx <- match(dataset_name, names(se))
      current_id_col   <- as.character(id_column[[idx]])[1]
      current_gene_col <- as.character(gene_column[[idx]])[1]
    }
    
    # 4. Subset platform metadata to just the key ID and Gene columns
    meta_subset <- platform_df[, c(current_id_col, current_gene_col), drop = FALSE]
    
    # 5. Merge expression matrix (by rownames) directly with platform metadata
    merged_df <- merge(
      x     = expr_df,
      y     = meta_subset,
      by.x  = "row.names",
      by.y  = current_id_col,
      all.x = TRUE
    )
    
    # Clean up row.names column name
    colnames(merged_df)[1] <- "probes"
    
    # 6. Relocate the gene symbol column directly to position #2
    merged_df <- merged_df %>%
      dplyr::relocate(dplyr::all_of(current_gene_col), .before = 2)
    
    #Renames the gene column to Gene
    merged_df <- merged_df %>% 
      dplyr::rename(gene = all_of(current_gene_col))
    
    return(merged_df)
  })
  
  names(res) <- names(se)
  return(res)
}

check_log2_transform <- function(annotated_expression_matrices) {
  
  res <- lapply(annotated_expression_matrices, function(expr_df) {
    
    # Calculate maximum value across sample columns (columns 3 onwards)
    maximum <- max(expr_df[, 3:ncol(expr_df)], na.rm = TRUE)
    
    # Transform if raw intensities (> 30)
    if (maximum > 30) {
      expr_df[, 3:ncol(expr_df)] <- log2(expr_df[, 3:ncol(expr_df)] + 1)
    }
    
    return(expr_df)
  })
  
  return(res)
}

## This function does two things, deletes all probes that don't map to one ----
#single gene and deletes all probes missing in more than 15% of samples 
delete_NAs <- function(annotated_expression_matrices){
  
  res <- lapply(names(annotated_expression_matrices), function(dataset_name) {
    
    expr_df <- as.data.frame(annotated_expression_matrices[[dataset_name]])
    
    #Deletes all NAs or empty spaces in the Gene column to keep only single named genes
    rows_to_keep <- !is.na(expr_df$gene) &
      expr_df$gene != "" &
      expr_df$gene != "-" &
      expr_df$gene != "---" &
      !grepl("///", expr_df$gene)
    
    expr_df$keep <- rows_to_keep
    expr_df <- subset(expr_df, keep == TRUE)
    expr_df$keep <- NULL
    
    #Deletes the random errors present in samples. If more than 15% of samples 
    #are missing a probe then the probe is deleted. 
    
    end_col <- ncol(expr_df) - 3
    expr_df$NAs <- rowSums(is.na(expr_df[, 3:end_col])) / end_col * 100
    expr_df <- subset(expr_df, NAs <= 15)
    expr_df$NAs <- NULL
    
    #Alphabetizes the gene symbols
    
    expr_df <- expr_df %>% dplyr::arrange(gene)
    
    return(expr_df)

  })
  names(res) <- names(annotated_expression_matrices)
  
  return(res)
}

# Collapses probes that map to the same genes based on the highest mean expression ----
collapse_probes <- function(complete_matrices){

  res <- lapply(names(complete_matrices), function(dataset_name){

    expr_df <- as.data.frame(complete_matrices[[dataset_name]])

    collapsed_data <- WGCNA::collapseRows(datET = expr_df[, 3:ncol(expr_df)],
                                          rowGroup = expr_df$gene,
                                          rowID = rownames(expr_df),
                                          method = "MaxMean")
    expr_df <- as.data.frame(collapsed_data$datETcollapsed)

    return(expr_df)
  })
  names(res) <- names(complete_matrices)

  return(res)
}
## Collapses probes that map to the same genes based on the highest mean expression ----
collapse_probes <- function(complete_matrices) {

  dataset_names <- names(complete_matrices)
  total_datasets <- length(dataset_names)

  res <- lapply(seq_along(dataset_names), function(i) {
    dataset_name <- dataset_names[i]

    # Extract dataframe
    expr_df <- as.data.frame(complete_matrices[[dataset_name]])

    # 2. Extract expression matrix (ensure numeric matrix conversion)
    expr_matrix <- as.matrix(expr_df[, 3:ncol(expr_df)])

    # 3. Perform WGCNA probe collapsing
    collapsed_data <- WGCNA::collapseRows(
      datET    = expr_matrix,
      rowGroup = expr_df$gene,
      rowID    = rownames(expr_df),
      method   = "MaxMean"
    )

    # Convert collapsed data back to data.frame
    collapsed_df <- as.data.frame(collapsed_data$datETcollapsed)

    return(collapsed_df)
  })

  # Preserve dataset names on the resulting list
  names(res) <- dataset_names

  return(res)
}
## Scales the data based on the Z score ----
z_score_matrices <- function(collapsed_matrices){
  res <- lapply(names(collapsed_matrices), function(dataset_name){
    
    expr_df <- as.data.frame(collapsed_matrices[[dataset_name]])
    
    #Replaces NAs for that row mean, effectively turning them to 0 when Z scoring
    row_means <- rowMeans(expr_df, na.rm = TRUE)
    na_indices <- which(is.na(expr_df), arr.ind = TRUE)
    expr_df[na_indices] <- row_means[na_indices[, 1]]
    
    # Checks the variance across genes 
    gene_sd <- apply(expr_df, 1, sd)
    filtered_matrix <- expr_df[gene_sd > 0, ]
    
    # Z Scoring
    expr_df <- scale(t(filtered_matrix))
    
    return(expr_df)
  })
  
  return(res)

}










## Tracks changes ----



# Metadata ----