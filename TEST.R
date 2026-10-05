full_metadata <- read.csv("full_metadata.csv")

full_metadata %>% group_by(dataset_ID) %>%
  count() %>%
  print(n = Inf)

