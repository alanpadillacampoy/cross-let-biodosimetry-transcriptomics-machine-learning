list_se <- readRDS("list_se.rds")
final_matrices <- readRDS("final_matrices.rds")
changes <- readRDS("changes_data_cleaning.rds")
wang_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "H4c3", "H4c9", "Ifit1", "Igfbp4", "LOC118567921",
                          "Lrrc70", "Ms4a1", "Phlda3", "Ptprn", "Rps20", "Serpine2", "Thy1")
#Wang data



#Linear Regression Model
linear_model <- lm()