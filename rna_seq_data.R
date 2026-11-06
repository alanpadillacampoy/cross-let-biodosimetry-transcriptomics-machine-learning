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


## Neutron whole normalized ----
neutron_set_1<- read.csv("GSE308422.csv")
neutron_set_1$symbol <- mapIds(org.Hs.eg.db,
                                          keys    = neutron_set_1$Sample,
                                          keytype = "ENSEMBL",
                                          column  = "SYMBOL",
                                          multiVals = "first")
neutron_set_1 <- neutron_set_1 %>% relocate(symbol, .before = 2) 
neutron_set_1[1,2] <- "Dose"
neutron_set_1$Sample <- NULL
neutron_set_1 <- neutron_set_1 %>% rename( Sample = symbol) %>%
  filter(Sample %in% c("Dose", gene_space_human))
rownames(neutron_set_1) <- neutron_set_1[,1]
neutron_set_1[,1] <- NULL

## Neutron fast normalized ----
neutron_set_2<- read.csv("GSE308422_RARAF_fast.csv")
neutron_set_2$symbol <- mapIds(org.Hs.eg.db,
                               keys    = neutron_set_2$Sample,
                               keytype = "ENSEMBL",
                               column  = "SYMBOL",
                               multiVals = "first")
neutron_set_2 <- neutron_set_2 %>% relocate(symbol, .before = 2) 
neutron_set_2[1,2] <- "Dose"
neutron_set_2$Sample <- NULL
neutron_set_2 <- neutron_set_2 %>% rename( Sample = symbol) %>%
  filter(Sample %in% c("Dose", gene_space_human))
rownames(neutron_set_2) <- neutron_set_2[,1]
neutron_set_2[,1] <- NULL

## Neutron thermal normalized ----
neutron_set_3<- read.csv("GSE308422_RINSC_thermal.csv")
neutron_set_3$symbol <- mapIds(org.Hs.eg.db,
                               keys    = neutron_set_3$Sample,
                               keytype = "ENSEMBL",
                               column  = "SYMBOL",
                               multiVals = "first")
neutron_set_3 <- neutron_set_3 %>% relocate(symbol, .before = 2) 
neutron_set_3[1,2] <- "Dose"
neutron_set_3$Sample <- NULL
neutron_set_3 <- neutron_set_3 %>% rename( Sample = symbol) %>%
  filter(Sample %in% c("Dose", gene_space_human))
rownames(neutron_set_3) <- neutron_set_3[,1]
neutron_set_3[,1] <- NULL

## Neutron whole RAW ----
neutron_set_4<- read.csv("GSE308422_raw.csv")
neutron_set_4$symbol <- mapIds(org.Hs.eg.db,
                               keys    = neutron_set_4$Sample,
                               keytype = "ENSEMBL",
                               column  = "SYMBOL",
                               multiVals = "first")
neutron_set_4 <- neutron_set_4 %>% relocate(symbol, .before = 2) 
neutron_set_4[1,2] <- "Dose"
neutron_set_4$Sample <- NULL
neutron_set_4 <- neutron_set_4 %>% rename( Sample = symbol) %>%
  filter(Sample %in% c("Dose", gene_space_human))
rownames(neutron_set_4) <- neutron_set_4[,1]
neutron_set_4[,1] <- NULL

## Neutron fast RAW ----
neutron_set_5<- read.csv("GSE308422_RARAF_fast_raw.csv")
neutron_set_5$symbol <- mapIds(org.Hs.eg.db,
                               keys    = neutron_set_5$Sample,
                               keytype = "ENSEMBL",
                               column  = "SYMBOL",
                               multiVals = "first")
neutron_set_5 <- neutron_set_5 %>% relocate(symbol, .before = 2) 
neutron_set_5[1,2] <- "Dose"
neutron_set_5$Sample <- NULL
neutron_set_5 <- neutron_set_5 %>% rename( Sample = symbol) %>%
  filter(Sample %in% c("Dose", gene_space_human))
rownames(neutron_set_5) <- neutron_set_5[,1]
neutron_set_5[,1] <- NULL

## Neutron thermal RAW ----
neutron_set_6<- read.csv("GSE308422_RINSC_thermal_raw.csv")
neutron_set_6$symbol <- mapIds(org.Hs.eg.db,
                               keys    = neutron_set_6$Sample,
                               keytype = "ENSEMBL",
                               column  = "SYMBOL",
                               multiVals = "first")
neutron_set_6 <- neutron_set_6 %>% relocate(symbol, .before = 2) 
neutron_set_6[1,2] <- "Dose"
neutron_set_6$Sample <- NULL
neutron_set_6 <- neutron_set_6 %>% rename( Sample = symbol) %>%
  filter(Sample %in% c("Dose", gene_space_human))
rownames(neutron_set_6) <- neutron_set_6[,1]
neutron_set_6[,1] <- NULL