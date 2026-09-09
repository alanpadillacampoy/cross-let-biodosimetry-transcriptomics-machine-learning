library(dplyr)
library(randomForest)
library(tibble)
library(ranger)
library(glmnet)

list_se <- readRDS("list_se.rds")
final_matrices <- readRDS("final_matrices.rds")
changes <- readRDS("changes_data_cleaning.rds")

wang_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "H4c3", "H4c9", "Ifit1", "Igfbp4", "LOC118567921",
                          "Lrrc70", "Ms4a1", "Phlda3", "Ptprn", "Rps20", "Serpine2", "Thy1")


# Read Wang Data
wang_training_data <- read.csv("wang_training_dataset.csv")
wang_training_data <- as.data.frame(t(wang_training_data))
colnames(wang_training_data) <- wang_training_data[1,]
wang_training_data <- wang_training_data[-1,]
wang_training_data <- type.convert(wang_training_data, as.is = TRUE)

# Gene Selection
wang_training_data <- wang_training_data %>% dplyr::select(c("Dose", all_of(wang_signature_genes)))

# Linear Regression Model
linear_model <- lm(Dose ~ ., data = wang_training_data)

# Random Forest Regression (version with ranger)
random_forest_model <- ranger::ranger(
  Dose ~ .,
  data = wang_training_data,
  num.trees = 60,
  seed = 42,
  max.depth = 2,
  splitrule = "variance"
)

# Polynomial Linear Model
wang_poly <- wang_training_data %>% 
  dplyr::mutate(across(2:last_col(), ~ .x^2))
polynomial_model <- lm(Dose ~ ., data = wang_poly)

# Square Root Linear Model
wang_root <- wang_training_data %>% 
  dplyr::mutate(across(2:last_col(), ~ sqrt(.x)))
square_root_model <- lm(Dose ~., data = wang_root)

# Overall

predictions <- tibble::tibble(Dose = wang_training_data$Dose, 
                              Linear = predict(linear_model),
                              RandomForest = random_forest_model$predictions,
                              Polynomial = predict(polynomial_model),
                              SquareRoot = predict(square_root_model)
                              )
predictions
