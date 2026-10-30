salah_dataset <- get_DoReMiTra_data("SE_Salah_2025_ExVivo")

mouse_rna_datasets <- list(read.csv("GSE296450.csv"),read.csv("GSE306805.csv"),
                           read.csv("GSE105266_gamma.csv"))
names(mouse_rna_datasets) <- c("GSE296450", "GSE306805", "GSE105266")

mouse_rna_datasets <- lapply(mouse_rna_datasets, function(x){
  x <- x %>% filter(Sample %in% c("Dose", gene_space_mouse))
  rownames(x) <- x[,1]
  x[,1] <- NULL
  return(x)
})

mouse_rna_datasets[[3]]
