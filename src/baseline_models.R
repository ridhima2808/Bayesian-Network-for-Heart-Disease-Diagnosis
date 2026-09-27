# Cardio-AI baseline models
#
# Prerequisite:
#   Rscript src/preprocessing.R
#
# Run from the project root with:
#   Rscript src/baseline_models.R

required_packages <- c("rpart")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before training baseline models: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

random_seed <- 509L
classification_threshold <- 0.50

find_project_root <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  candidates <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[[1L]])
    script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
    candidates <- unique(c(candidates, script_dir, dirname(script_dir)))
  }

  for (candidate in candidates) {
    current <- candidate
    repeat {
      if (file.exists(file.path(current, "README.md")) &&
          dir.exists(file.path(current, "src"))) {
        return(current)
      }

      parent <- dirname(current)
      if (identical(parent, current)) {
        break
      }
      current <- parent
    }
  }

  stop("Could not locate the Cardio-AI project root.", call. = FALSE)
}

validate_model_data <- function(train, test) {
  if (!identical(names(train), names(test))) {
    stop("Training and test schemas do not match.", call. = FALSE)
  }
  if (!"target" %in% names(train)) {
    stop("The target column is missing.", call. = FALSE)
  }
  if (anyNA(train) || anyNA(test)) {
    stop("Training or test data contains missing values.", call. = FALSE)
  }
  if (!is.factor(train$target) || !identical(levels(train$target), c("0", "1"))) {
    stop("Training target must be a factor with levels 0 and 1.", call. = FALSE)
  }
  if (!is.factor(test$target) || !identical(levels(test$target), c("0", "1"))) {
    stop("Test target must be a factor with levels 0 and 1.", call. = FALSE)
  }

  factor_columns <- names(train)[vapply(train, is.factor, logical(1))]
  predictor_factor_columns <- setdiff(factor_columns, "target")
  for (column_name in predictor_factor_columns) {
    observed_train_levels <- unique(as.character(train[[column_name]]))
    unseen_test_levels <- setdiff(
      unique(as.character(test[[column_name]])),
      observed_train_levels
    )
    if (length(unseen_test_levels) > 0L) {
      stop(
        "Test column '", column_name,
        "' contains level(s) absent from training: ",
        paste(unseen_test_levels, collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

align_factor_levels <- function(train, test) {
  predictor_factor_columns <- setdiff(
    names(train)[vapply(train, is.factor, logical(1))],
    "target"
  )
  for (column_name in predictor_factor_columns) {
    observed_levels <- levels(droplevels(train[[column_name]]))
    train[[column_name]] <- factor(
      as.character(train[[column_name]]),
      levels = observed_levels
    )
    test[[column_name]] <- factor(
      as.character(test[[column_name]]),
      levels = observed_levels
    )
  }

  list(train = train, test = test)
}

project_root <- find_project_root()
training_path <- file.path(project_root, "data", "processed", "heart_train.rds")
test_path <- file.path(project_root, "data", "processed", "heart_test.rds")
model_dir <- file.path(project_root, "results", "models")

if (!file.exists(training_path) || !file.exists(test_path)) {
  stop(
    "Processed train/test data not found. Run src/preprocessing.R first.",
    call. = FALSE
  )
}

heart_train <- readRDS(training_path)
heart_test <- readRDS(test_path)
validate_model_data(heart_train, heart_test)
aligned_data <- align_factor_levels(heart_train, heart_test)
heart_train <- aligned_data$train
heart_test <- aligned_data$test

set.seed(random_seed)

# Logistic regression uses the original numeric measurements and categorical
# factors. The explicit binary response makes the positive class unambiguous.
logistic_train <- heart_train
logistic_train$target_binary <- as.integer(logistic_train$target == "1")
logistic_train$target <- NULL

logistic_model <- stats::glm(
  target_binary ~ .,
  data = logistic_train,
  family = stats::binomial(link = "logit")
)
logistic_probability <- as.numeric(stats::predict(
  logistic_model,
  newdata = heart_test[setdiff(names(heart_test), "target")],
  type = "response"
))
if (anyNA(logistic_probability) || any(!is.finite(logistic_probability))) {
  stop("Logistic regression produced invalid test probabilities.", call. = FALSE)
}
logistic_class <- factor(
  ifelse(logistic_probability >= classification_threshold, "1", "0"),
  levels = c("0", "1")
)

# The decision tree uses the exact same predictors and training/test rows.
decision_tree_model <- rpart::rpart(
  target ~ .,
  data = heart_train,
  method = "class",
  control = rpart::rpart.control(
    cp = 0.01,
    minsplit = 10,
    xval = 10
  )
)
tree_probabilities <- stats::predict(
  decision_tree_model,
  newdata = heart_test[setdiff(names(heart_test), "target")],
  type = "prob"
)
if (!"1" %in% colnames(tree_probabilities)) {
  stop("Decision tree did not return probabilities for class 1.", call. = FALSE)
}
tree_probability <- as.numeric(tree_probabilities[, "1"])
if (anyNA(tree_probability) || any(!is.finite(tree_probability))) {
  stop("Decision tree produced invalid test probabilities.", call. = FALSE)
}
tree_class <- factor(
  ifelse(tree_probability >= classification_threshold, "1", "0"),
  levels = c("0", "1")
)

predictions <- data.frame(
  test_row = seq_len(nrow(heart_test)),
  actual = as.character(heart_test$target),
  logistic_probability = logistic_probability,
  logistic_class = as.character(logistic_class),
  decision_tree_probability = tree_probability,
  decision_tree_class = as.character(tree_class),
  stringsAsFactors = FALSE
)

logistic_coefficients <- data.frame(
  term = rownames(summary(logistic_model)$coefficients),
  summary(logistic_model)$coefficients,
  row.names = NULL,
  check.names = FALSE
)

tree_importance <- decision_tree_model$variable.importance
if (is.null(tree_importance)) {
  tree_importance_table <- data.frame(
    variable = character(0),
    importance = numeric(0)
  )
} else {
  tree_importance_table <- data.frame(
    variable = names(tree_importance),
    importance = as.numeric(tree_importance),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  tree_importance_table <- tree_importance_table[
    order(tree_importance_table$importance, decreasing = TRUE),
    ,
    drop = FALSE
  ]
}

baseline_bundle <- list(
  logistic_model = logistic_model,
  decision_tree_model = decision_tree_model,
  predictions = predictions,
  classification_threshold = classification_threshold,
  random_seed = random_seed,
  training_rows = nrow(heart_train),
  test_rows = nrow(heart_test),
  predictor_columns = setdiff(names(heart_train), "target"),
  target_levels = levels(heart_train$target)
)

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(baseline_bundle, file.path(model_dir, "baseline_models.rds"))
write.csv(
  predictions,
  file.path(model_dir, "baseline_predictions.csv"),
  row.names = FALSE
)
write.csv(
  logistic_coefficients,
  file.path(model_dir, "logistic_coefficients.csv"),
  row.names = FALSE
)
write.csv(
  tree_importance_table,
  file.path(model_dir, "decision_tree_importance.csv"),
  row.names = FALSE
)

message("Baseline model training complete.")
message("Logistic regression and decision tree saved to: ", model_dir)
