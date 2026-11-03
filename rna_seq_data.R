library(org.Hs.eg.db)
library(AnnotationDbi)

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


## Proton datasets
neutron_set <- read.csv("GSE308422.csv")
neutron_set$symbol <- mapIds(org.Hs.eg.db,
                                          keys    = neutron_set$Sample,
                                          keytype = "ENSEMBL",
                                          column  = "SYMBOL",
                                          multiVals = "first")
neutron_set <- neutron_set %>% relocate(symbol, .before = 2) 
neutron_set[1,2] <- "Dose"
neutron_set$Sample <- NULL
neutron_set <- neutron_set %>% rename( Sample = symbol) %>%
  filter(Sample %in% c("Dose", gene_space_human))
rownames(neutron_set) <- neutron_set[,1]
neutron_set[,1] <- NULL
