df <- read.csv("full_metadata.csv")
library(dplyr)
library(tibble)

# table <- df %>%
#   group_by(dataset_ID) %>%
#   summarise(
#     species = paste(unique(organism), collapse = ", "),
#     time = paste(unique(time_hours), collapse = ", "),
#     radiation = paste(unique(radiation_type), collapse = ", "),
#     product = paste(unique(tissue), collapse = ", "),
#     dose_range = max(dose_Gy),
#     exp_set = paste(unique(experimental_setting), collapse = ", ")
#   )

#ortholog maps

list_se <- readRDS("list_se.rds")
final_matrices <- readRDS("final_matrices.rds")
changes <- readRDS("changes_data_cleaning.rds")
wang_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "H4c3", "H4c9", "Ifit1", "Igfbp4", "LOC118567921",
                "Lrrc70", "Ms4a1", "Phlda3", "Ptprn", "Rps20", "Serpine2", "Thy1")
li_signature_genes <- c("ACTA2", "AEN", "ASTN2", "ATF3", "ATM", "BAX", "BBC3", "BTG2",
                        "CCNG1", "CDKN1A", "DDB2", "EGR1", "EGR4", "FAS", "FDXR", "GADD45A",
                        "GDF15", "JUN", "KU80", "MDM2", "MYC", "PCNA", "PHLDA3", "PHPT1",
                        "PLK3", "POLH", "PPM1D", "RPS27L", "SESN1", "TNFRSR10B", "TNFSF4",
                        "TNFSF9", "TP53I3", "TRIAP1", "TRIM22", "XPC", "ZMAT3")

wang_hs_genes <- orthogene::convert_orthologs(wang_mm_genes, "rownames", "rownames", 
                                                          FALSE, "mmusculus", "hsapiens", 
                                                          non121_strategy = "drop_both_species")

which(rownames(final_matrices[[1]]) == "H4C9")








final_matrices_photon <- final_matrices
final_matrices_photon$SE_Broustas_2018_InVivo_GSE113509_GPL11202 <- NULL
final_matrices_photon$SE_Broustas_2017_InVivo_GSE85323_GPL10333 <- NULL
final_matrices_photon$SE_Broustas_2017_ExVivo_GSE90909_GPL13497 <- NULL
#these are photon but with bad genes
final_matrices_photon$SE_Flores_2009_ExVivo_GSE15341_GPL8332 <- NULL
final_matrices_photon$SE_Ghandhi_2015_ExVivo_GSE65292_GPL13497 <- NULL
final_matrices_photon$SE_Ghandhi_2018_InVivo_GSE84898_GPL13497 <- NULL







count_gene_occurrences <- function(matrices_list) {
  
  # 1. Extract all gene names from all datasets into a single long vector
  all_genes <- unlist(lapply(matrices_list, rownames), use.names = FALSE)
  
  # 2. Count frequency and format into a tidy tibble
  gene_counts <- table(all_genes) %>% 
    as_tibble() %>% 
    dplyr::rename(gene = all_genes, dataset_count = n) %>% 
    arrange(desc(dataset_count))
  
  return(gene_counts)
}

# Usage:
gene_freq <- count_gene_occurrences(final_matrices_photon)

x = 20

# Calculate cumulative counts for datasets >= 18
cumulative_counts <- gene_freq %>%
  filter(dataset_count >= x) %>% 
  group_by(dataset_count) %>%
  summarise(n_genes = n()) %>%
  arrange(desc(dataset_count)) %>%
  mutate(cumulative_genes = cumsum(n_genes))

ggplot(cumulative_counts, aes(x = dataset_count, y = cumulative_genes)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "darkred", size = 2.5) +
  geom_text(aes(label = cumulative_genes), vjust = -1, size = 3) +
  scale_x_continuous(breaks = x:32) +
  labs(
    title = "Cumulative Gene Pool (At Least Half of Datasets)",
    x = "Minimum Dataset Presence Threshold",
    y = "Total Retained Genes"
  ) +
  coord_cartesian(xlim = c(x, 32)) +
  theme_minimal()
# 1. Calculate exact slopes and rate of change of slope
slope_table <- cumulative_counts %>%
  arrange(dataset_count) %>%
  mutate(
    # 1st Difference: Change in retained genes between step n and step n-1
    slope = cumulative_genes - lag(cumulative_genes),
    
    # 2nd Difference: Rate of change of the slope itself (acceleration/curvature)
    slope_change = slope - lag(slope)
  )

# 2. Inspect the slope table around your cutoff (n >= 25)
print(slope_table %>% filter(dataset_count >= 20))

library(dplyr)

# 1. Extract the core gene names (present in >= 31 datasets)
core_genes_31 <- gene_freq %>% 
  filter(dataset_count >= 28) %>% 
  pull(gene)

message(sprintf("Retained %d core genes present in >= 28 datasets.", length(core_genes_31)))

# 2. Subset your final matrices to only keep these core 2,282 genes
final_matrices_filtered <- lapply(final_matrices_photon, function(mat) {
  keep_rows <- intersect(rownames(mat), core_genes_31)
  return(mat[keep_rows, , drop = FALSE])
})

