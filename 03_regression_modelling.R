library(dplyr)
library(randomForest)
library(tibble)
library(ranger)
library(glmnet)
library(catboost)
library(rsample)
library(purrr)

list_se <- readRDS("list_se.rds")
final_matrices <- readRDS("final_matrices.rds")
changes <- readRDS("changes_data_cleaning.rds")

wang_signature_genes <- c("Ccng1", "Dgka", "Fzr1", "H4c3", "H4c9", "Ifit1", "Igfbp4", "LOC118567921",
                          "Lrrc70", "Ms4a1", "Phlda3", "Ptprn", "Rps20", "Serpine2", "Thy1")
set.seed(42)

# Read Wang Training Data
wang_training_data <- read.csv("wang_training_dataset.csv", header = FALSE)
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

wang_training_data_augmented$Target <- NULL

augmented_set <- rbind(wang_training_data_augmented, wang_training_data)
rownames(augmented_set) <- NULL
augmented_set$Sample_id <- seq_len(nrow(augmented_set))
augmented_set <- augmented_set %>% relocate(Sample_id, 1)

# Cross Validation

cross_validation <- rsample::vfold_cv(augmented_set, v = 5, repeats = 10)

predictions <- lapply(cross_validation$splits, function(splits){
  testing_folds <- rsample::assessment(splits)
  training_folds <- rsample::training(splits)
  
  # Simple Linear model
  linear_model <- lm(Dose ~ . - Sample_id, data = training_folds)

  # Polynomial Linear Model
  wang_poly <- training_folds %>% 
    dplyr::mutate(across(3:last_col(), ~ .x^2))
  polynomial_model <- lm(Dose ~ . - Sample_id, data = wang_poly)
  
  # Square Root Linear Model
  wang_root <- training_folds %>% 
    dplyr::mutate(across(3:last_col(), ~ sqrt(.x)))
  square_root_model <- lm(Dose ~. - Sample_id, data = wang_root)
  
  # Random Forest Regression
  random_forest_model <- ranger::ranger(
    Dose ~ . - Sample_id,
    data = training_folds,
    num.trees = 60,
    seed = 42,
    max.depth = 2,
    splitrule = "variance"
  )
  
  # ElasticNet Regression Model
  elastic_model <- glmnet::glmnet(
    x = as.matrix(training_folds[, 3:ncol(training_folds)]),
    y = training_folds$Dose,
    alpha = 0.5,
    lambda = 1.0,
    standardize = FALSE,
    thresh = 0.0002,
    seed = 42
  )
  
  # Create the data frame MOVE THIS TO THE END
  matrix_general <- data.frame(
    Sample_id = testing_folds$Sample_id,
    Dose = testing_folds$Dose,
    Linear = predict(linear_model, newdata = testing_folds),
    Polynomial = predict(polynomial_model, newdata = testing_folds),
    SquareRoot = predict(square_root_model, newdata = testing_folds),
    RandomForest = predict(random_forest_model, 
                           data = testing_folds)$predictions,
    ElasticNet = 
      as.numeric(glmnet::predict.glmnet(elastic_model,
                                        newx = as.matrix(testing_folds[, 3:ncol(testing_folds)])))
    )
  
  catboost_matrix <- data.frame(
    Sample_id = training_folds$Sample_id,
    Dose = training_folds$Dose,
    Linear = predict(linear_model, newdata = training_folds),
    Polynomial = predict(polynomial_model, newdata = training_folds),
    SquareRoot = predict(square_root_model, newdata = training_folds),
    RandomForest = predict(random_forest_model, 
                           data = training_folds)$predictions,
    ElasticNet = 
      as.numeric(glmnet::predict.glmnet(elastic_model,
                                        newx = as.matrix(training_folds[, 3:ncol(training_folds)])))
  )
  
  #Meta-model: CatBoost Regression Model
  catboost_model <- catboost::catboost.train(
    learn_pool = catboost::catboost.load_pool(
      data = as.matrix(catboost_matrix[, 3:7]),
      label = catboost_matrix$Dose
    ),
    params = list(
      loss_function = "RMSE",
      depth = 2,
      iterations = 50,
      learning_rate = 0.08,
      l2_leaf_reg = 4.5,
      subsample = 0.75,
      random_seed = 42
    )
  )
  
  #predict testing_folds with the base trained on training_folds
  #predict testing_folds with the catboost trained on training_folds

})


################################################################################
################################################################################

final_predictions_matrix <- predictions %>% purrr::list_rbind() %>%
  dplyr::arrange(Sample_id) %>%
  tibble::remove_rownames() %>%
  dplyr::group_by(Sample_id) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  dplyr::ungroup()



# Meta-model: CatBoost Regression

catboost_model <- catboost::catboost.train(
  learn_pool = catboost::catboost.load_pool(
    data = as.matrix(final_predictions_matrix[, 3:7]),
    label = final_predictions_matrix$Dose
  ),
  params = list(
    loss_function = "RMSE",
    depth = 2,
    iterations = 50,
    learning_rate = 0.08,
    l2_leaf_reg = 4.5,
    subsample = 0.75,
    random_seed = 42
  )
)
catboost_predictions <- catboost.predict(
  catboost_model,
  catboost.load_pool(
    data = as.matrix(final_predictions_matrix[, 3:7])
  )
)
final_predictions_matrix$CatBoost <- catboost_predictions
