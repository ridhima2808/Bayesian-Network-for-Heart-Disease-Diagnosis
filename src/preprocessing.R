# Cardio-AI data preprocessing
#
# Run from the project root with:
#   Rscript src/preprocessing.R
#
# Outputs are written to data/processed/. RDS is used so factor types and levels
# are preserved for the Bayesian network and baseline models.

expected_columns <- c(
  "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
  "thalach", "exang", "oldpeak", "slope", "ca", "thal", "target"
)

numeric_columns <- c("age", "trestbps", "chol", "thalach", "oldpeak")

categorical_levels <- list(
  sex = c("0", "1"),
  cp = c("0", "1", "2", "3"),
  fbs = c("0", "1"),
  restecg = c("0", "1", "2"),
  exang = c("0", "1"),
  slope = c("0", "1", "2"),
  ca = c("0", "1", "2", "3", "4"),
  thal = c("0", "1", "2", "3"),
  target = c("0", "1")
)

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
      if (file.exists(file.path(current, "data", "heart.csv"))) {
        return(current)
      }

      parent <- dirname(current)
      if (identical(parent, current)) {
        break
      }
      current <- parent
    }
  }

  stop(
    "Could not locate the project root containing data/heart.csv.",
    call. = FALSE
  )
}

validate_raw_data <- function(data) {
  if (nrow(data) != 303L) {
    stop("Expected 303 rows; found ", nrow(data), ".", call. = FALSE)
  }
  if (ncol(data) != 14L) {
    stop("Expected 14 columns; found ", ncol(data), ".", call. = FALSE)
  }
  if (!identical(names(data), expected_columns)) {
    stop(
      "Column names or order do not match the expected dataset schema.",
      call. = FALSE
    )
  }
  if (anyNA(data)) {
    stop("The raw dataset contains missing values.", call. = FALSE)
  }

  target_values <- sort(unique(as.character(data$target)))
  if (!identical(target_values, c("0", "1"))) {
    stop("Target must contain only 0 and 1.", call. = FALSE)
  }

  invisible(TRUE)
}

coerce_numeric_column <- function(values, column_name) {
  converted <- suppressWarnings(as.numeric(as.character(values)))
  if (anyNA(converted) || any(!is.finite(converted))) {
    stop(
      "Column '", column_name, "' contains a non-numeric or non-finite value.",
      call. = FALSE
    )
  }
  converted
}

coerce_factor_column <- function(values, column_name, allowed_levels) {
  character_values <- trimws(as.character(values))
  unexpected <- setdiff(unique(character_values), allowed_levels)

  if (length(unexpected) > 0L) {
    stop(
      "Column '", column_name, "' contains unexpected value(s): ",
      paste(sort(unexpected), collapse = ", "),
      call. = FALSE
    )
  }

  factor(character_values, levels = allowed_levels)
}

stratified_split <- function(target, train_fraction = 0.80, seed = 509L) {
  if (train_fraction <= 0 || train_fraction >= 1) {
    stop("train_fraction must be between 0 and 1.", call. = FALSE)
  }

  set.seed(seed)
  indices_by_class <- split(seq_along(target), target, drop = TRUE)
  train_indices <- unlist(
    lapply(indices_by_class, function(indices) {
      sample(indices, size = floor(length(indices) * train_fraction))
    }),
    use.names = FALSE
  )

  sort(train_indices)
}

project_root <- find_project_root()
input_path <- file.path(project_root, "data", "heart.csv")
output_dir <- file.path(project_root, "data", "processed")

heart_raw <- read.csv(
  input_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A", "?"),
  strip.white = TRUE,
  fileEncoding = "UTF-8-BOM"
)

validate_raw_data(heart_raw)

heart_clean <- heart_raw
for (column_name in numeric_columns) {
  heart_clean[[column_name]] <- coerce_numeric_column(
    heart_clean[[column_name]],
    column_name
  )
}
for (column_name in names(categorical_levels)) {
  heart_clean[[column_name]] <- coerce_factor_column(
    heart_clean[[column_name]],
    column_name,
    categorical_levels[[column_name]]
  )
}

if (anyNA(heart_clean)) {
  stop("Preprocessing introduced missing values.", call. = FALSE)
}

train_indices <- stratified_split(heart_clean$target)
test_indices <- setdiff(seq_len(nrow(heart_clean)), train_indices)

heart_train <- heart_clean[train_indices, , drop = FALSE]
heart_test <- heart_clean[test_indices, , drop = FALSE]

# Preserve all predefined factor levels, even when a level is absent in a split.
stopifnot(
  identical(lapply(heart_train[names(categorical_levels)], levels), categorical_levels),
  identical(lapply(heart_test[names(categorical_levels)], levels), categorical_levels)
)

split_manifest <- data.frame(
  source_row = seq_len(nrow(heart_clean)),
  split = ifelse(seq_len(nrow(heart_clean)) %in% train_indices, "train", "test")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(heart_clean, file.path(output_dir, "heart_clean.rds"))
saveRDS(heart_train, file.path(output_dir, "heart_train.rds"))
saveRDS(heart_test, file.path(output_dir, "heart_test.rds"))
write.csv(
  split_manifest,
  file.path(output_dir, "split_manifest.csv"),
  row.names = FALSE
)

message("Preprocessing complete.")
message("Clean rows: ", nrow(heart_clean))
message("Training rows: ", nrow(heart_train))
message("Test rows: ", nrow(heart_test))
message("Outputs: ", output_dir)
