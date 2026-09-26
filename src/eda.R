# Cardio-AI exploratory data analysis
#
# Prerequisite:
#   Rscript src/preprocessing.R
#
# Run from the project root with:
#   Rscript src/eda.R

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before running EDA: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

numeric_columns <- c("age", "trestbps", "chol", "thalach", "oldpeak")
categorical_columns <- c(
  "sex", "cp", "fbs", "restecg", "exang", "slope", "ca", "thal"
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

numeric_summary <- function(data, columns) {
  rows <- lapply(columns, function(column_name) {
    values <- data[[column_name]]
    data.frame(
      variable = column_name,
      n = length(values),
      mean = mean(values),
      sd = stats::sd(values),
      min = min(values),
      q1 = unname(stats::quantile(values, 0.25)),
      median = stats::median(values),
      q3 = unname(stats::quantile(values, 0.75)),
      max = max(values),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

categorical_summary <- function(data, columns) {
  rows <- lapply(columns, function(column_name) {
    counts <- as.data.frame(
      table(value = data[[column_name]], target = data$target),
      stringsAsFactors = FALSE
    )
    counts$variable <- column_name
    counts$within_variable_percent <- ave(
      counts$Freq,
      counts$variable,
      FUN = function(values) 100 * values / sum(values)
    )
    names(counts)[names(counts) == "Freq"] <- "count"
    counts[c("variable", "value", "target", "count", "within_variable_percent")]
  })

  do.call(rbind, rows)
}

project_root <- find_project_root()
input_path <- file.path(project_root, "data", "processed", "heart_clean.rds")
output_dir <- file.path(project_root, "results", "eda")

if (!file.exists(input_path)) {
  stop(
    "Cleaned data not found at ", input_path,
    ". Run src/preprocessing.R first.",
    call. = FALSE
  )
}

heart <- readRDS(input_path)
required_columns <- c(numeric_columns, categorical_columns, "target")
if (!all(required_columns %in% names(heart))) {
  stop("Cleaned data does not contain all required columns.", call. = FALSE)
}
if (anyNA(heart)) {
  stop("Cleaned data contains missing values.", call. = FALSE)
}
if (!all(vapply(heart[c(categorical_columns, "target")], is.factor, logical(1)))) {
  stop("Expected categorical columns and target to be factors.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Summary tables -----------------------------------------------------------

numeric_stats <- numeric_summary(heart, numeric_columns)
write.csv(
  numeric_stats,
  file.path(output_dir, "numeric_summary.csv"),
  row.names = FALSE
)

categorical_stats <- categorical_summary(heart, categorical_columns)
write.csv(
  categorical_stats,
  file.path(output_dir, "categorical_summary_by_target.csv"),
  row.names = FALSE
)

class_balance <- as.data.frame(table(target = heart$target))
names(class_balance)[names(class_balance) == "Freq"] <- "count"
class_balance$percent <- 100 * class_balance$count / sum(class_balance$count)
class_balance$diagnosis <- ifelse(
  class_balance$target == "1",
  "Heart disease",
  "No heart disease"
)
write.csv(
  class_balance[c("target", "diagnosis", "count", "percent")],
  file.path(output_dir, "class_balance.csv"),
  row.names = FALSE
)

correlation_matrix <- stats::cor(heart[numeric_columns], use = "complete.obs")
write.csv(
  correlation_matrix,
  file.path(output_dir, "numeric_correlations.csv"),
  row.names = TRUE
)

# Class-balance plot -------------------------------------------------------

class_balance$label <- sprintf(
  "%d (%.1f%%)",
  class_balance$count,
  class_balance$percent
)

class_balance_plot <- ggplot2::ggplot(
  class_balance,
  ggplot2::aes(x = diagnosis, y = count, fill = target)
) +
  ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    vjust = -0.4,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = c("0" = "#4C78A8", "1" = "#E45756")) +
  ggplot2::scale_y_continuous(
    limits = c(0, max(class_balance$count) * 1.15),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Heart disease class balance",
    x = NULL,
    y = "Number of patients"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(output_dir, "class_balance.png"),
  class_balance_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# Numeric distributions and target comparison -----------------------------

numeric_long <- utils::stack(heart[numeric_columns])
names(numeric_long) <- c("value", "variable")

numeric_distribution_plot <- ggplot2::ggplot(
  numeric_long,
  ggplot2::aes(x = value)
) +
  ggplot2::geom_histogram(
    bins = 25,
    fill = "#4C78A8",
    color = "white",
    linewidth = 0.2
  ) +
  ggplot2::facet_wrap(~variable, scales = "free", ncol = 2) +
  ggplot2::labs(
    title = "Distributions of numeric patient measurements",
    x = "Value",
    y = "Number of patients"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  file.path(output_dir, "numeric_distributions.png"),
  numeric_distribution_plot,
  width = 10,
  height = 9,
  dpi = 300
)

numeric_long$target <- rep(heart$target, times = length(numeric_columns))

numeric_by_target_plot <- ggplot2::ggplot(
  numeric_long,
  ggplot2::aes(x = target, y = value, fill = target)
) +
  ggplot2::geom_boxplot(outlier.alpha = 0.35, show.legend = FALSE) +
  ggplot2::facet_wrap(~variable, scales = "free_y", ncol = 2) +
  ggplot2::scale_fill_manual(values = c("0" = "#4C78A8", "1" = "#E45756")) +
  ggplot2::scale_x_discrete(labels = c("0" = "No disease", "1" = "Disease")) +
  ggplot2::labs(
    title = "Numeric measurements by diagnosis",
    x = NULL,
    y = "Value"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  file.path(output_dir, "numeric_by_target.png"),
  numeric_by_target_plot,
  width = 10,
  height = 9,
  dpi = 300
)

# Categorical distributions -----------------------------------------------

categorical_long <- do.call(
  rbind,
  lapply(categorical_columns, function(column_name) {
    data.frame(
      variable = column_name,
      value = as.character(heart[[column_name]]),
      target = heart$target,
      stringsAsFactors = FALSE
    )
  })
)

categorical_distribution_plot <- ggplot2::ggplot(
  categorical_long,
  ggplot2::aes(x = value, fill = target)
) +
  ggplot2::geom_bar(position = "dodge") +
  ggplot2::facet_wrap(~variable, scales = "free_x", ncol = 3) +
  ggplot2::scale_fill_manual(
    name = "Diagnosis",
    values = c("0" = "#4C78A8", "1" = "#E45756"),
    labels = c("0" = "No disease", "1" = "Disease")
  ) +
  ggplot2::labs(
    title = "Categorical feature distributions by diagnosis",
    x = "Category code",
    y = "Number of patients"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  file.path(output_dir, "categorical_distributions.png"),
  categorical_distribution_plot,
  width = 12,
  height = 9,
  dpi = 300
)

# Numeric correlation heatmap ---------------------------------------------

correlation_long <- as.data.frame(as.table(correlation_matrix))
names(correlation_long) <- c("variable_x", "variable_y", "correlation")

correlation_plot <- ggplot2::ggplot(
  correlation_long,
  ggplot2::aes(x = variable_x, y = variable_y, fill = correlation)
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", correlation)),
    size = 3.5
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "Correlation among numeric variables",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

ggplot2::ggsave(
  file.path(output_dir, "numeric_correlation_heatmap.png"),
  correlation_plot,
  width = 8,
  height = 7,
  dpi = 300
)

message("EDA complete. Tables and plots saved to: ", output_dir)
