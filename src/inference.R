# Cardio-AI Bayesian network inference
#
# Prerequisites:
#   Rscript src/preprocessing.R
#   Rscript src/bayesian_network.R
#
# Run the worked examples from the project root with:
#   Rscript src/inference.R
#
# This file can also be sourced by the Shiny application. Its primary public
# functions are load_bayesian_model(), diagnose_patient(), and
# sample_target_distribution().

required_packages <- c("bnlearn")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before running inference: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

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

load_bayesian_model <- function(model_path = NULL) {
  if (is.null(model_path)) {
    model_path <- file.path(
      find_project_root(),
      "results",
      "models",
      "bayesian_network.rds"
    )
  }
  if (!file.exists(model_path)) {
    stop(
      "Bayesian network model not found at ", model_path,
      ". Run src/bayesian_network.R first.",
      call. = FALSE
    )
  }

  bundle <- readRDS(model_path)
  required_fields <- c(
    "selected_name", "selected_fit", "discretization", "training_columns",
    "node_levels"
  )
  if (!all(required_fields %in% names(bundle))) {
    stop("The model bundle is missing required inference metadata.", call. = FALSE)
  }
  bundle
}

as_patient_list <- function(patient) {
  if (is.data.frame(patient)) {
    if (nrow(patient) != 1L) {
      stop("Patient data must contain exactly one row.", call. = FALSE)
    }
    patient <- as.list(patient[1L, , drop = FALSE])
  }
  if (!is.list(patient) || is.null(names(patient)) ||
      any(names(patient) == "")) {
    stop("Patient evidence must be a named list or one-row data frame.", call. = FALSE)
  }
  if (anyDuplicated(names(patient))) {
    stop("Patient evidence contains duplicate field names.", call. = FALSE)
  }
  patient
}

prepare_patient_evidence <- function(patient, model_bundle) {
  patient <- as_patient_list(patient)
  allowed_nodes <- setdiff(model_bundle$training_columns, "target")
  unexpected <- setdiff(names(patient), allowed_nodes)
  if (length(unexpected) > 0L) {
    stop(
      "Unknown or disallowed evidence field(s): ",
      paste(unexpected, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(patient) == 0L) {
    stop("Provide at least one item of patient evidence.", call. = FALSE)
  }

  evidence <- lapply(names(patient), function(column_name) {
    value <- patient[[column_name]]
    if (length(value) != 1L || is.na(value)) {
      stop(
        "Evidence field '", column_name, "' must contain one non-missing value.",
        call. = FALSE
      )
    }

    if (column_name %in% names(model_bundle$discretization)) {
      numeric_value <- suppressWarnings(as.numeric(as.character(value)))
      if (is.na(numeric_value) || !is.finite(numeric_value)) {
        stop(
          "Evidence field '", column_name, "' must be a finite number.",
          call. = FALSE
        )
      }
      specification <- model_bundle$discretization[[column_name]]
      return(as.character(cut(
        numeric_value,
        breaks = specification$breaks,
        labels = specification$labels,
        include.lowest = TRUE,
        ordered_result = FALSE
      )))
    }

    categorical_value <- as.character(value)
    allowed_levels <- model_bundle$node_levels[[column_name]]
    if (!categorical_value %in% allowed_levels) {
      stop(
        "Evidence field '", column_name, "' must be one of: ",
        paste(allowed_levels, collapse = ", "),
        call. = FALSE
      )
    }
    categorical_value
  })
  names(evidence) <- names(patient)
  evidence
}

diagnose_patient <- function(
    patient,
    model_bundle = load_bayesian_model(),
    n = 50000L,
    seed = 509L) {
  if (length(n) != 1L || is.na(n) || n < 1000) {
    stop("n must be a single value of at least 1000.", call. = FALSE)
  }
  evidence <- prepare_patient_evidence(patient, model_bundle)

  set.seed(seed)
  disease_probability <- bnlearn::cpquery(
    model_bundle$selected_fit,
    event = (target == "1"),
    evidence = evidence,
    method = "lw",
    n = as.integer(n)
  )
  if (length(disease_probability) != 1L ||
      !is.finite(disease_probability)) {
    stop("Inference did not return a finite probability.", call. = FALSE)
  }

  structure(
    list(
      probability_disease = unname(disease_probability),
      probability_no_disease = unname(1 - disease_probability),
      evidence = evidence,
      method = "likelihood weighting via bnlearn::cpquery",
      particles = as.integer(n),
      seed = as.integer(seed),
      selected_structure = model_bundle$selected_name
    ),
    class = "cardio_diagnosis"
  )
}

sample_target_distribution <- function(
    patient,
    model_bundle = load_bayesian_model(),
    n = 50000L,
    seed = 509L) {
  if (length(n) != 1L || is.na(n) || n < 1000) {
    stop("n must be a single value of at least 1000.", call. = FALSE)
  }
  evidence <- prepare_patient_evidence(patient, model_bundle)

  set.seed(seed)
  samples <- bnlearn::cpdist(
    model_bundle$selected_fit,
    nodes = "target",
    evidence = evidence,
    method = "lw",
    n = as.integer(n)
  )
  weights <- attr(samples, "weights")
  if (is.null(weights) || length(weights) != nrow(samples) || sum(weights) <= 0) {
    stop("Conditional sampling did not return usable likelihood weights.", call. = FALSE)
  }

  target_levels <- model_bundle$node_levels$target
  weighted_probabilities <- vapply(target_levels, function(target_level) {
    sum(weights[as.character(samples$target) == target_level]) / sum(weights)
  }, numeric(1))

  data.frame(
    target = target_levels,
    probability = unname(weighted_probabilities),
    stringsAsFactors = FALSE
  )
}

print.cardio_diagnosis <- function(x, ...) {
  cat("Bayesian network diagnosis\n")
  cat(sprintf("P(heart disease | evidence): %.3f\n", x$probability_disease))
  cat(sprintf("P(no heart disease | evidence): %.3f\n", x$probability_no_disease))
  cat("Discretized evidence:\n")
  print(unlist(x$evidence))
  invisible(x)
}

run_worked_examples <- function() {
  model_bundle <- load_bayesian_model()

  example_a <- list(
    age = 45,
    sex = 0,
    cp = 1,
    trestbps = 120,
    chol = 210,
    fbs = 0,
    restecg = 0,
    thalach = 165,
    exang = 0,
    oldpeak = 0.2,
    slope = 0,
    ca = 0,
    thal = 2
  )
  example_b <- list(
    age = 63,
    sex = 1,
    cp = 3,
    trestbps = 145,
    chol = 280,
    fbs = 1,
    restecg = 1,
    thalach = 115,
    exang = 1,
    oldpeak = 2.8,
    slope = 2,
    ca = 2,
    thal = 3
  )

  diagnosis_a <- diagnose_patient(example_a, model_bundle, seed = 509L)
  diagnosis_b <- diagnose_patient(example_b, model_bundle, seed = 510L)

  cat("\nWorked example A\n")
  print(diagnosis_a)
  cat("\nWorked example B\n")
  print(diagnosis_b)

  example_results <- data.frame(
    example = c("A", "B"),
    probability_disease = c(
      diagnosis_a$probability_disease,
      diagnosis_b$probability_disease
    ),
    probability_no_disease = c(
      diagnosis_a$probability_no_disease,
      diagnosis_b$probability_no_disease
    )
  )
  output_dir <- file.path(find_project_root(), "results", "inference")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    example_results,
    file.path(output_dir, "worked_examples.csv"),
    row.names = FALSE
  )

  invisible(example_results)
}

if (sys.nframe() == 0L) {
  run_worked_examples()
}
