## This script is the primary data cleaner for the metadata,
## it is meant to normalize all the different experiments
## for easier comparisons

#libraries ----
library(DoReMiTra)
library(limma)
library(SummarizedExperiment)
library(readr)
library(stringr)
library(dplyr)

#defines the Summarized-experiment ----
list_se <- list_DoReMiTra_datasets()

#these lists match the names of the specific type of tissue 
#if the tissue col in the metadata is tissue.ch1 then its on list_tissue
#if its in cell.type.ch1
list_tissue = c(1, 2, 3, 4, 5, 7, 8, 11, 13, 17, 21, 22, 23, 25, 29, 30, 32,
                33, 34, 35)
list_cell = c(6, 9, 15, 16, 18)

#creates one matrix per experiment, and fills most of the ---- 
#columns with information from the experiment's metadata 
for (i in 1:35){
  sum_exp <- get_DoReMiTra_data(list_se$Dataset[i])
  number_rows <- length(colnames(sum_exp))
  metadata <- matrix(nrow = number_rows, ncol = 14)
  columns <- colData(sum_exp)
  colnames(metadata) <- c("dataset_name", "dataset_ID", "platform_ID", "sample_ID", 
                          "organism", "dose_Gy", "radiation_type", "sex", 
                          "time_hours", "tissue", "experimental_setting", 
                          "dose_rate_photon", "platform", "condition")
  metadata[,1] <- list_se$Dataset[i]
  metadata[,2] <- strsplit(metadata[,1], "_")[[1]][5]
  metadata[,3] <- strsplit(metadata[,1], "_")[[1]][6]
  metadata[,4] <- columns$geo_accession
  metadata[,5] <- str_to_lower(columns$Organism)
  metadata[,6] <- as.numeric(ifelse(str_detect(columns$Dose, "Control"), 0, 
    str_remove_all(columns$Dose, "[a-zA-Z/_\\s-]")))
  metadata[,7] <- ifelse(metadata[,6] == 0, "control",
                         str_to_lower(columns$Radiation_type))
  metadata[,8] <- str_to_lower(columns$Sex)
  metadata[,9] <- ifelse(str_detect(columns$Time_point, "(?i)pre-irradiation"), 
                         0, 
                         ifelse(
                           str_detect(columns$Time_point, "(?i)day"), 
                           as.numeric(str_remove_all(columns$Time_point, 
                                                     "[a-zA-Z/_-]")) * 24 * 60, 
                           ifelse(
                             str_detect(columns$Time_point, "(?i)hr|hour|h"), 
                             as.numeric(str_remove_all(columns$Time_point, 
                                                       "[a-zA-Z/_-]")) * 60,  
                             as.numeric(str_remove_all(columns$Time_point, 
                                                       "[a-zA-Z/_-]")))))
  metadata[,10] <- if (i %in% list_tissue) {
    metadata[, 10] <- str_to_lower(columns$tissue.ch1)
  } else if (i %in% list_cell) {
    metadata[, 10] <- str_to_lower(columns$cell.type.ch1)
  } else {
    metadata[, 10] <- ""
  }
  metadata[,11] <- str_to_lower(columns$Exp_setting)
  metadata[,13] <- str_to_lower(columns$Platform)
  metadata[,14] <- ifelse(metadata[,6] == 0, "control", "irradiated")
  matrix_name <- paste0("metadata_", i)
  assign(matrix_name, metadata)
}

#Corrects the Dose error in SE_Broustas_2017_ExVivo_GSE90909_GPL13497 ----
se21 <- DoReMiTra::get_DoReMiTra_data("SE_Broustas_2017_ExVivo_GSE90909_GPL13497")

se21$Dose<- ifelse(
  str_detect(se21$title, regex("control", ignore_case = TRUE)),
  "0Gy",
  str_replace_all(
    str_extract(se21$title, "\\d+(\\.\\d+)?\\s*Gy"),
    "\\s+",
    ""
  )
)
number_rows <- length(colnames(se21))
metadata <- matrix(nrow = number_rows, ncol = 14)
columns <- colData(se21)
colnames(metadata) <- c("dataset_name", "dataset_ID", "platform_ID", "sample_ID", 
                        "organism", "dose_Gy", "radiation_type", "sex", 
                        "time_hours", "tissue", "experimental_setting", 
                        "dose_rate_photon", "platform", "condition")
metadata[,1] <- list_se$Dataset[21]
metadata[,2] <- strsplit(metadata[,1], "_")[[1]][5]
metadata[,3] <- strsplit(metadata[,1], "_")[[1]][6]
metadata[,4] <- columns$geo_accession
metadata[,5] <- str_to_lower(columns$Organism)
metadata[,6] <- ifelse(str_detect(columns$Dose, "Control"), 0, 
                       str_remove_all(columns$Dose, "[a-zA-Z/_-]"))
metadata[,7] <- ifelse(metadata[,6] == "0", "control", 
                       str_to_lower(columns$Radiation_type))
metadata[,8] <- str_to_lower(columns$Sex)
metadata[,9] <- ifelse(str_detect(columns$Time_point, "(?i)pre-irradiation"), 
                       0, 
                       ifelse(
                         str_detect(columns$Time_point, "(?i)day"), 
                         as.numeric(str_remove_all(columns$Time_point, 
                                                   "[a-zA-Z/_-]")) * 24 * 60, 
                         ifelse(
                           str_detect(columns$Time_point, "(?i)hr|hour|h"), 
                           as.numeric(str_remove_all(columns$Time_point, 
                                                     "[a-zA-Z/_-]")) * 60,  
                           as.numeric(str_remove_all(columns$Time_point, 
                                                     "[a-zA-Z/_-]")))))
metadata[,10] <- columns$tissue.ch1
metadata[,11] <- str_to_lower(columns$Exp_setting)
metadata[,13] <- str_to_lower(columns$Platform)
metadata[,14] <- ifelse(metadata[,6] == 0, "control", "irradiated")
metadata_21 <- metadata

#adds the photon dose rate information extracted from the literature ----
metadata_1[,12] <- 1.05
metadata_2[,12] <- 1.05
metadata_3[,12] <- 1.05
metadata_4[,12] <- 0.8
metadata_5[,12] <- 0.8
metadata_6[,12] <- 0.67
metadata_7[,12] <- 0.82
metadata_8[,12] <- 0.82
metadata_9[,12] <- 0.05
metadata_10[,12] <- 2.8
metadata_11[,12] <- 0.82
metadata_12[,12] <- 0.45
metadata_13[,12] <- ifelse(metadata_13[, 4] == "0.15 Gy", 6,
                           ifelse(metadata_13[, 4] == "0.30 Gy", 6, 5.6))
metadata_14[,12] <- 0.82
metadata_15[,12] <- 1.16
metadata_16[,12] <- NA
metadata_17[,12] <- ifelse(metadata_17[, 4] == "0.15 Gy", 6,
                            ifelse(metadata_17[, 4] == "0.3 Gy", 6, 5.6))
metadata_18[,12] <- NA
metadata_19[,12] <- ifelse(metadata_19[,7] == "neutron",  0.0258, 0.006667)
metadata_20[,12] <- ifelse(metadata_20[,7] == "neutron",  0.0258, 0.006667)
metadata_21[,12] <- ifelse(metadata_21[,7] == "neutron",  0.0258, 1.23)
metadata_22[,12] <- NA
metadata_23[,12] <- colData(get_DoReMiTra_data(list_se$Dataset[23]))$Dose_rate
metadata_23[,12] <- ifelse(metadata_23[, 12] == "Acute", 1.03,
                            ifelse(metadata_23[, 12] == "Low", 0.0031, 0))
metadata_24[,12] <- NA
metadata_25[,12] <- 0.1
metadata_26[,12] <- 0.86
metadata_27[,12] <- 1.67
metadata_28[,12] <- 1
metadata_29[,12] <- 1
metadata_30[,12] <- NA
metadata_31[,12] <- 1.67
metadata_32[,12] <- 1.45
metadata_33[,12] <- 0.85
metadata_34[,12] <- 6
metadata_35[,12] <- NA

#this was extracted from the other parts of the metadata and are manually input
metadata_10[,10] <- "Human Peripheral blood lymphocytes (PBL)"
metadata_12[,10] <- colData(get_DoReMiTra_data(list_se$Dataset[12]))$Cell_type
metadata_14[,10] <- "whole blood"
metadata_19[,10] <- "whole blood"
metadata_20[,10] <- "whole blood"
metadata_24[,10] <- "peripheral blood lymphocytes"
metadata_26[,10] <- "whole blood"
metadata_27[,10] <- "whole blood"
metadata_28[,10] <- "whole blood"
metadata_31[,10] <- "whole blood"

#compiles all the metadata matrices into one ----
full_metadata <- metadata_1
for (i in 2:35) {
  count <- paste0("metadata_", i)
  full_metadata <- rbind(full_metadata, get(count))
}

#syntax corrections ----
#corrects dose rates
full_metadata[,12] <- ifelse(full_metadata[,6] == 0, 0, full_metadata[,12])

#corrects tissues types, and categorizes them into 4 groups
full_metadata[,10] <- case_when(
  str_detect(full_metadata[,10],
             "(?i)(h1f|hepm|hescs)") ~ "non blood cell lines",
  str_detect(full_metadata[,10],
             "(?i)(cd4|cd8|nk|lymphocyte tcd4)") ~ "isolated lymphocytes",
  str_detect(full_metadata[,10],
             "(?i)(pbmc|pbl|lymphocyte|leukocyte)") ~ "PBMCs PBLs",
  str_detect(full_metadata[,10], 
             "(?i)(whole blood|peripheral blood|^blood$)") ~ "whole blood",
)

write.csv(full_metadata, 
          file = "full_metadata.csv", row.names = FALSE)