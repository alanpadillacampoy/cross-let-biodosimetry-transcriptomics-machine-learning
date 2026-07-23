#This is the master function document. All functions are stored here for 
#reference in the project

#Libraries ----
library(DoReMiTra)
library(SummarizedExperiment)
library(stringr)
library(dplyr)
library(BiocManager)
library(GEOquery)
library(tibble)
library(WGCNA)

##Data cleaning, quality and validation ----

#Extracts the platform from the data set and checks for uniqueness 
unique_gpl_platforms <- function(list_se){
  platforms <- vector("list", length = 35)
  for (i in 1:35) {
    platforms[i] <- strsplit(names(list_se)[i], "_")[[1]][6]
  }
  unique_platforms <- unique(platforms)
  
  return(unique_platforms)
}

#Launches a GEO query to obtain the platform's metadata
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

#Extracts the expression matrix form the SE
get_expression_matrices <- function(list_se) {
  lapply(list_se, function(se) {
    as.data.frame(SummarizedExperiment::assay(se))
  })
}

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

# join_probes <- function(se){
#   lapply(se, function(df_expression){
#     probes <- rownames(df_expression)
#     data.frame(probes = probes, stringsAsFactors = FALSE)
#   })
# }


#Selects the possible columns where the probes and symbols may be:
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

#Annotates the expression matrices with the genes, matching the probes and the genes they map for
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
