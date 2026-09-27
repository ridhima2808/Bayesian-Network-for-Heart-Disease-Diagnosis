# Cardio-AI model evaluation
#
# Prerequisites, in order:
#   Rscript src/preprocessing.R
#   Rscript src/bayesian_network.R
#   Rscript src/baseline_models.R
#
# Run from the project root with:
#   Rscript src/evaluate.R

required_packages <- c("bnlearn", "ggplot2", "pROC")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before evaluation: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

classification_threshold <- 0.50
inference_particles <- 20000L
random_seed <- 509L

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

safe_ratio <- function(numerator, denominator) {
  if (denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

calculate_classification_metrics <- function(actual, predicted, probability) {
  actual <- factor(actual, levels = c("0", "1"))
  predicted <- factor(predicted, levels = c("0", "1"))
  confusion <- table(actual = actual, predicted = predicted)

  true_negative <- unname(confusion["0", "0"])
  false_positive <- unname(confusion["0", "1"])
  false_negative <- unname(confusion["1", "0"])
  true_positive <- unname(confusion["1", "1"])

  roc_object <- pROC::roc(
    response = actual,
    predictor = probability,
    levels = c("0", "1"),
    direction = "<",
    quiet = TRUE
  )

  precision <- safe_ratio(true_positive, true_positive + false_positive)
  recall <- safe_ratio(true_positive, true_positive + false_negative)

  list(
    metrics = c(
      accuracy = safe_ratio(
        true_positive + true_negative,
        sum(confusion)
      ),
      auc = as.numeric(pROC::auc(roc_object)),
      sensitivity = recall,
      specificity = safe_ratio(
        true_negative,
        true_negative + false_positive
      ),
      precision = precision,
      f1 = if (is.na(precision) || is.na(recall) || precision + recall == 0) {
        NA_real_
      } else {
        2 * precision * recall / (precision + recall)
      }
    ),
    confusion = confusion,
    roc = roc_object
  )
}

confusion_to_data_frame <- function(confusion, model_name) {
  result <- as.data.frame(confusion, stringsAsFactors = FALSE)
  names(result) <- c("actual", "predicted", "count")
  result$model <- model_name
  result[c("model", "actual", "predicted", "count")]
}

roc_to_data_frame <- function(roc_object, model_name) {
  result <- data.frame(
    false_positive_rate = 1 - roc_object$specificities,
    true_positive_rate = roc_object$sensitivities,
    model = model_name,
    stringsAsFactors = FALSE
  )
  result[order(result$false_positive_rate, result$true_positive_rate), ]
}

project_root <- find_project_root()
inference_script <- file.path(project_root, "src", "inference.R")
test_path <- file.path(project_root, "data", "processed", "heart_test.rds")
bayesian_path <- file.path(project_root, "results", "models", "bayesian_network.rds")
baseline_path <- file.path(project_root, "results", "models", "baseline_models.rds")
output_dir <- file.path(project_root, "results", "evaluation")

required_files <- c(inference_script, test_path, bayesian_path, baseline_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required evaluation input(s) are missing: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

# Defines load_bayesian_model() and diagnose_patient(). Worked examples do not
# run when inference.R is sourced.
source(inference_script, local = environment())

heart_test <- readRDS(test_path)
bayesian_bundle <- load_bayesian_model(bayesian_path)
baseline_bundle <- readRDS(baseline_path)

if (anyNA(heart_test) || !"target" %in% names(heart_test)) {
  stop("Test data is missing the target or contains missing values.", call. = FALSE)
}
if (!"predictions" %in% names(baseline_bundle)) {
  stop("Baseline model bundle does not contain test predictions.", call. = FALSE)
}

baseline_predictions <- baseline_bundle$predictions
if (nrow(baseline_predictions) != nrow(heart_test)) {
  stop("Baseline predictions do not align with the saved test split.", call. = FALSE)
}
actual <- as.character(heart_test$target)
if (!identical(actual, as.character(baseline_predictions$actual))) {
  stop("Baseline prediction labels do not match the saved test targets.", call. = FALSE)
}

# Use every non-target field as evidence for each test patient. A row-specific
# seed makes the Monte Carlo predictions reproducible across runs.
predictor_columns <- setdiff(names(heart_test), "target")
bayesian_probability <- vapply(seq_len(nrow(heart_test)), function(row_index) {
  patient <- as.list(heart_test[row_index, predictor_columns, drop = FALSE])
  diagnosis <- diagnose_patient(
    patient,
    model_bundle = bayesian_bundle,
    n = inference_particles,
    seed = random_seed + row_index - 1L
  )
  diagnosis$probability_disease
}, numeric(1))

if (anyNA(bayesian_probability) || any(!is.finite(bayesian_probability))) {
  stop("Bayesian network produced invalid test probabilities.", call. = FALSE)
}

prediction_table <- data.frame(
  test_row = seq_len(nrow(heart_test)),
  actual = actual,
  bayesian_network_probability = bayesian_probability,
  bayesian_network_class = ifelse(
    bayesian_probability >= classification_threshold,
    "1",
    "0"
  ),
  logistic_regression_probability = baseline_predictions$logistic_probability,
  logistic_regression_class = baseline_predictions$logistic_class,
  decision_tree_probability = baseline_predictions$decision_tree_probability,
  decision_tree_class = baseline_predictions$decision_tree_class,
  stringsAsFactors = FALSE
)

model_inputs <- list(
  `Bayesian network` = list(
    predicted = prediction_table$bayesian_network_class,
    probability = prediction_table$bayesian_network_probability
  ),
  `Logistic regression` = list(
    predicted = prediction_table$logistic_regression_class,
    probability = prediction_table$logistic_regression_probability
  ),
  `Decision tree` = list(
    predicted = prediction_table$decision_tree_class,
    probability = prediction_table$decision_tree_probability
  )
)

evaluation <- lapply(model_inputs, function(model_input) {
  calculate_classification_metrics(
    actual = actual,
    predicted = model_input$predicted,
    probability = model_input$probability
  )
})

metrics_table <- do.call(rbind, lapply(names(evaluation), function(model_name) {
  data.frame(
    model = model_name,
    as.list(evaluation[[model_name]]$metrics),
    row.names = NULL,
    check.names = FALSE
  )
}))

confusion_table <- do.call(rbind, lapply(names(evaluation), function(model_name) {
  confusion_to_data_frame(evaluation[[model_name]]$confusion, model_name)
}))

roc_table <- do.call(rbind, lapply(names(evaluation), function(model_name) {
  roc_to_data_frame(evaluation[[model_name]]$roc, model_name)
}))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  metrics_table,
  file.path(output_dir, "model_metrics.csv"),
  row.names = FALSE
)
write.csv(
  confusion_table,
  file.path(output_dir, "confusion_matrices.csv"),
  row.names = FALSE
)
write.csv(
  prediction_table,
  file.path(output_dir, "test_predictions.csv"),
  row.names = FALSE
)
write.csv(
  roc_table,
  file.path(output_dir, "roc_curve_points.csv"),
  row.names = FALSE
)

# Accuracy and AUC comparison ----------------------------------------------

metric_plot_data <- rbind(
  data.frame(
    model = metrics_table$model,
    metric = "Accuracy",
    value = metrics_table$accuracy
  ),
  data.frame(
    model = metrics_table$model,
    metric = "AUC",
    value = metrics_table$auc
  )
)

metric_plot <- ggplot2::ggplot(
  metric_plot_data,
  ggplot2::aes(x = model, y = value, fill = metric)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", value)),
    position = ggplot2::position_dodge(width = 0.8),
    vjust = -0.35,
    size = 3.5
  ) +
  ggplot2::scale_fill_manual(values = c("Accuracy" = "#4C78A8", "AUC" = "#F58518")) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, by = 0.2),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Test-set model performance",
    x = NULL,
    y = "Score",
    fill = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 15, hjust = 1),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(output_dir, "model_metric_comparison.png"),
  metric_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# Confusion matrices -------------------------------------------------------

confusion_plot <- ggplot2::ggplot(
  confusion_table,
  ggplot2::aes(x = predicted, y = actual, fill = count)
) +
  ggplot2::geom_tile(color = "white", linewidth = 1) +
  ggplot2::geom_text(ggplot2::aes(label = count), size = 5) +
  ggplot2::facet_wrap(~model, nrow = 1) +
  ggplot2::scale_fill_gradient(low = "#E8F1F8", high = "#2C699A") +
  ggplot2::scale_x_discrete(labels = c("0" = "No disease", "1" = "Disease")) +
  ggplot2::scale_y_discrete(labels = c("0" = "No disease", "1" = "Disease")) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "Test-set confusion matrices",
    x = "Predicted class",
    y = "Actual class",
    fill = "Patients"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

ggplot2::ggsave(
  file.path(output_dir, "confusion_matrices.png"),
  confusion_plot,
  width = 12,
  height = 4.8,
  dpi = 300
)

# ROC curves ---------------------------------------------------------------

auc_labels <- setNames(
  sprintf("%s (AUC %.3f)", metrics_table$model, metrics_table$auc),
  metrics_table$model
)
roc_table$model_label <- unname(auc_labels[roc_table$model])

roc_plot <- ggplot2::ggplot(
  roc_table,
  ggplot2::aes(
    x = false_positive_rate,
    y = true_positive_rate,
    color = model_label
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "#888888"
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::coord_equal() +
  ggplot2::scale_x_continuous(limits = c(0, 1)) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(
    title = "Test-set ROC curves",
    x = "False positive rate",
    y = "True positive rate",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  file.path(output_dir, "roc_curves.png"),
  roc_plot,
  width = 8,
  height = 7,
  dpi = 300
)

evaluation_bundle <- list(
  metrics = metrics_table,
  confusion_matrices = lapply(evaluation, `[[`, "confusion"),
  roc_objects = lapply(evaluation, `[[`, "roc"),
  predictions = prediction_table,
  classification_threshold = classification_threshold,
  inference_particles = inference_particles,
  random_seed = random_seed
)
saveRDS(evaluation_bundle, file.path(output_dir, "evaluation_results.rds"))

message("Evaluation complete. Results saved to: ", output_dir)
print(metrics_table, row.names = FALSE)
