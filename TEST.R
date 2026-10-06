full_metadata <- read.csv("full_metadata.csv")

full_metadata %>% group_by(organism) %>%
  count(dataset_name) %>%
  print(n = Inf)

homo_sapiens_datasets <- c(
  "SE_Amundson_2008_ExVivo_GSE8917_GPL1708",
  "SE_Amundson_2011_ExVivo_GSE23515_GPL6480",
  "SE_Amundson_2011_InVivo_GSE20162_GPL6480",
  "SE_Ankermit_2015_ExVivo_GSE55953_GPL14550",
  "SE_Flores_2009_ExVivo_GSE15341_GPL8332",
  "SE_Ghandhi_2015_ExVivo_GSE65292_GPL13497",
  "SE_Girardi_2012_ExVivo_GSE20173_GPL6480",
  "SE_Gruel_2008_ExVivo_GSE6978_GPL4803",
  "SE_Lee_2013_ExVivo_GSE44245_GPL570",
  "SE_Manikandan_2014_ExVivo_GSE36355_GPL6883",
  "SE_Nosel_2013_ExVivo_GSE43151_GPL13497",
  "SE_Park_2017_ExVivo_GSE102971_GPL10332_HomoSapiens",
  "SE_Paul_2010_InVivo_GSE23393_GPL6480",
  "SE_Paul_2013_ExVivo_GSE44201_GPL6480",
  "SE_Paul_2013_ExVivo_GSE44201_GPL6848",
  "SE_Rouchka_2015_ExVivo_GSE64375_GPL6244",
  "SE_Rouchka_2019_ExVivo_GSE63952_GPL15207",
  "SE_Vasilyev_2017_ExVivo_GSE97000_GPL17077"
)

macaca_mulatta_datasets <- c(
  "SE_Ghandhi_2018_InVivo_GSE84898_GPL13497",
  "SE_Park_2017_ExVivo_GSE102971_GPL10332_MacacaMulatta"
)

mus_musculus_datasets <- c(
  "SE_Amundson_2018_InVivo_GSE101402_GPL11202",
  "SE_Amundson_2018_InVivo_GSE99176_GPL11202",
  "SE_Amundson_2019_InVivo_GSE124612_GPL11202",
  "SE_Aryankalayil_2018_InVivo_GSE104121_GPL10787",
  "SE_Aryankalayil_2018_InVivo_GSE104121_GPL21163",
  "SE_Broustas_2021_InVivo_GSE132559_GPL11202",
  "SE_Broustas_2021_InVivo_GSE133451_GPL11202",
  "SE_Broustas_2022_InVivo_GSE184361_GPL11202",
  "SE_Broustas_2023_InVivo_GSE196400_GPL11202",
  "SE_Mukherjee_2019_InVivo_GSE114142_GPL11202",
  "SE_Paul_2015_InVivo_GSE62623_GPL10333",
  "SE_Yamaguchi_2020_InVivo_GSE137192_GPL1261"
)

neutron_datasets <- c("SE_Broustas_2017_ExVivo_GSE90909_GPL13497",
                      "SE_Broustas_2017_InVivo_GSE85323_GPL10333",
                      "SE_Broustas_2018_InVivo_GSE113509_GPL11202")


for (i in 1:length(mus_musculus_datasets)){
  print(unique(colData(list_se[[mus_musculus_datasets[1]]])$Organism))
}

lapply(mus_musculus_datasets, function(dataset){
  a <- final_matrices[[dataset]]
  x <- as.data.frame(t(a))
  genes_in_x <- colnames(select(x, any_of(wang_murine_signature_genes)))
  genes_not_in_signature <- symdiff(genes_in_x, wang_murine_signature_genes)
})

