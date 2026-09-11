library(dplyr)
library(randomForest)
library(tibble)
library(ranger)
library(glmnet)
library(catboost)

list_se <- readRDS("list_se.rds")
final_matrices <- readRDS("final_matrices.rds")
changes <- readRDS("changes_data_cleaning.rds")

wang_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "H4c3", "H4c9", "Ifit1", "Igfbp4", "LOC118567921",
                          "Lrrc70", "Ms4a1", "Phlda3", "Ptprn", "Rps20", "Serpine2", "Thy1")
set.seed(42)

# Read Wang Training Data
wang_training_data <- read.csv("wang_training_dataset.csv")
wang_training_data <- as.data.frame(t(wang_training_data))
colnames(wang_training_data) <- wang_training_data[1,]
wang_training_data <- wang_training_data[-1,]
wang_training_data <- type.convert(wang_training_data, as.is = TRUE)

# Read Wang Test Data
wang_testing_data <- read.csv("wang_testing_dataset.csv")
wang_testing_data <- as.data.frame(t(wang_testing_data))
colnames(wang_testing_data) <- wang_testing_data[1,]
wang_testing_data <- wang_testing_data[-1,]
wang_testing_data <- type.convert(wang_testing_data, as.is = TRUE)

# Gene Selection
wang_training_data <- wang_training_data %>% dplyr::select(c("Dose", all_of(wang_signature_genes)))
wang_testing_data <- wang_testing_data %>% dplyr::select(c("Dose", all_of(wang_signature_genes)))

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

# ElasticNet Regression Model
elastic_model <- glmnet::glmnet(
  x = as.matrix(wang_training_data[, 2:ncol(wang_training_data)]),
  y = wang_training_data$Dose,
  alpha = 0.5,
  lambda = 1.0,
  standardize = FALSE,
  thresh = 0.0002
)

# Data augmentation

wang_training_data_augmented <- wang_training_data %>%
  mutate(Target = case_when(
    Dose == 0  ~ 0,
    Dose < 4  ~ 2,
    Dose < 6  ~ 4,
    Dose < 9  ~ 6,
    Dose < 12  ~ 9,
    TRUE  ~ 12,
  )) %>% 
  relocate(Target, 1) %>%
  dplyr::group_by(Target) %>%
  dplyr::mutate(across(-Dose, ~ .x + runif(n = 1, min = -0.1, max = 0.1)*sd(.x))) %>%
  dplyr::ungroup()

augmented_set <- rbind(wang_training_data_augmented, wang_training_data)
augmented_set$Target <- NULL



























# Overall

predictions <- tibble::tibble(Dose = wang_testing_data$Dose, 
                              Linear = predict(linear_model, newdata = wang_testing_data),
                              RandomForest = predict(random_forest_model, 
                                                     data = wang_testing_data)$predictions,
                              Polynomial = predict(polynomial_model, newdata = wang_testing_data),
                              SquareRoot = predict(square_root_model, newdata = wang_testing_data),
                              ElasticNet = 
                                as.numeric(glmnet::predict.glmnet(elastic_model,
                                                       newx = as.matrix(wang_testing_data[, 2:ncol(wang_testing_data)])))
                              )
predictions

# Augmentation





# Meta-model: CatBoost Regression

#dese_pool <- catboost::catboost.load_pool()
#meta_model <- catboost::catboost.train(predictions, )