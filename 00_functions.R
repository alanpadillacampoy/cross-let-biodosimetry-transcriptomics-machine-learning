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
  platforms <- vector("character", length = 35)
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
  
  # 1. Extract the platform ID for dataset #17
  platform <- strsplit(names(list_se)[17], "_")[[1]][6]
  
  # 2. Split each row by the triple slashes '///' to isolate the sections
  sections_list <- strsplit(geo_metadata[[platform]]$gene_assignment, "///", fixed = TRUE)
  
  # 3. Extract targets and apply the uniqueness condition per row
  extracted_targets <- sapply(sections_list, function(row_sections) {
    
    # Trim spaces from sections
    row_sections <- trimws(row_sections)
    
    # Split each individual section by the double slash '//'
    sub_pieces <- strsplit(row_sections, "//", fixed = TRUE)
    
    # Grab the 2nd element (gene symbol) of each section
    targets <- sapply(sub_pieces, function(piece) {
      if (length(piece) >= 2) {
        return(trimws(piece[2]))
      } else {
        return(NA)
      }
    })
    
    # Remove any NAs generated during extraction before checking uniqueness
    targets <- targets[!is.na(targets)]
    
    # Check uniqueness
    unique_targets <- unique(targets)
    
    if (length(unique_targets) == 1) {
      return(unique_targets)       
    } else {
      return(NA_character_)         
    }
  }, USE.NAMES = FALSE)
  
  # 4. Add the new column to the targeted platform data frame
  geo_metadata[[platform]]$GENE_SYMBOL <- extracted_targets
  
  # 5. CRITICAL: Return the modified list so changes persist outside the function!
  return(geo_metadata)
}
