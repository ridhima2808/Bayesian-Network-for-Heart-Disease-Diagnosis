# Cardio-AI Bayesian network structure and parameter learning
#
# Prerequisite:
#   Rscript src/preprocessing.R
#
# Run from the project root with:
#   Rscript src/bayesian_network.R

required_packages <- c("bnlearn")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before training the Bayesian network: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

numeric_columns <- c("age", "trestbps", "chol", "thalach", "oldpeak")
random_seed <- 509L
equivalent_sample_size <- 10

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

learn_discretization <- function(data, columns, bins = 3L) {
  if (bins != 3L) {
    stop("This project expects exactly three ordered bins.", call. = FALSE)
  }

  specifications <- lapply(columns, function(column_name) {
    values <- data[[column_name]]
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
      stop(
        "Numeric column '", column_name,
        "' is missing, non-numeric, or contains invalid values.",
        call. = FALSE
      )
    }

    internal_breaks <- unname(stats::quantile(
      values,
      probs = seq(0, 1, length.out = bins + 1L)[-c(1L, bins + 1L)],
      type = 7
    ))
    if (anyDuplicated(internal_breaks)) {
      stop(
        "Cannot form distinct tertile boundaries for '", column_name, "'.",
        call. = FALSE
      )
    }

    list(
      breaks = c(-Inf, internal_breaks, Inf),
      labels = c("low", "medium", "high"),
      training_range = range(values)
    )
  })
  names(specifications) <- columns
  specifications
}

apply_discretization <- function(data, specifications) {
  discretized <- data
  for (column_name in names(specifications)) {
    specification <- specifications[[column_name]]
    discretized[[column_name]] <- cut(
      data[[column_name]],
      breaks = specification$breaks,
      labels = specification$labels,
      include.lowest = TRUE,
      ordered_result = FALSE
    )
  }

  if (anyNA(discretized)) {
    stop("Discretization introduced missing values.", call. = FALSE)
  }
  if (!all(vapply(discretized, is.factor, logical(1)))) {
    stop("Bayesian network input must contain factors only.", call. = FALSE)
  }

  discretized
}

build_expert_network <- function(node_names) {
  # Upstream risk factors point toward the diagnosis. The diagnosis then points
  # to symptoms and test findings, allowing evidence on those observations to
  # update posterior disease probability without an excessively large CPT.
  expert_arcs <- data.frame(
    from = c(
      "age", "age", "age", "sex", "sex", "fbs", "trestbps", "chol",
      "target", "target", "target", "target", "target", "target", "target",
      "cp", "exang", "oldpeak"
    ),
    to = c(
      "trestbps", "chol", "target", "chol", "target", "target", "target",
      "target", "cp", "restecg", "thalach", "exang", "oldpeak", "ca",
      "thal", "exang", "oldpeak", "slope"
    ),
    stringsAsFactors = FALSE
  )

  unknown_nodes <- setdiff(unique(c(expert_arcs$from, expert_arcs$to)), node_names)
  if (length(unknown_nodes) > 0L) {
    stop(
      "Expert structure references unknown node(s): ",
      paste(unknown_nodes, collapse = ", "),
      call. = FALSE
    )
  }

  network <- bnlearn::empty.graph(node_names)
  for (row_index in seq_len(nrow(expert_arcs))) {
    network <- bnlearn::set.arc(
      network,
      from = expert_arcs$from[[row_index]],
      to = expert_arcs$to[[row_index]],
      check.cycles = TRUE
    )
  }
  network
}

plot_network_base <- function(network, title) {
  node_names <- bnlearn::nodes(network)
  node_count <- length(node_names)
  angles <- seq(pi / 2, pi / 2 - 2 * pi, length.out = node_count + 1L)
  angles <- angles[-length(angles)]
  coordinates <- cbind(x = cos(angles), y = sin(angles))
  rownames(coordinates) <- node_names
  network_arcs <- bnlearn::arcs(network)

  graphics::plot.new()
  graphics::plot.window(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), asp = 1)
  graphics::title(main = title, line = 0.5, cex.main = 1.05)

  if (nrow(network_arcs) > 0L) {
    for (arc_index in seq_len(nrow(network_arcs))) {
      from <- coordinates[network_arcs[arc_index, "from"], ]
      to <- coordinates[network_arcs[arc_index, "to"], ]
      direction <- to - from
      distance <- sqrt(sum(direction^2))
      unit_direction <- direction / distance
      start <- from + 0.13 * unit_direction
      end <- to - 0.16 * unit_direction
      graphics::arrows(
        start[[1L]], start[[2L]], end[[1L]], end[[2L]],
        length = 0.07,
        angle = 22,
        col = "#6B7280",
        lwd = 1.1
      )
    }
  }

  node_colors <- ifelse(node_names == "target", "#E45756", "#DCEAF7")
  graphics::symbols(
    coordinates[, "x"],
    coordinates[, "y"],
    circles = rep(0.13, node_count),
    inches = FALSE,
    add = TRUE,
    bg = node_colors,
    fg = "#263238"
  )
  graphics::text(
    coordinates[, "x"],
    coordinates[, "y"],
    labels = node_names,
    cex = 0.68,
    font = ifelse(node_names == "target", 2, 1)
  )
}

save_network_panel <- function(networks, titles, output_path) {
  grDevices::png(output_path, width = 2100, height = 760, res = 180)
  tryCatch(
    {
      graphics::par(mfrow = c(1, length(networks)), mar = c(1, 1, 3, 1))
      for (index in seq_along(networks)) {
        plot_network_base(networks[[index]], titles[[index]])
      }
    },
    finally = grDevices::dev.off()
  )
}

save_single_network <- function(network, title, output_path) {
  grDevices::png(output_path, width = 1200, height = 1100, res = 180)
  tryCatch(
    plot_network_base(network, title),
    finally = grDevices::dev.off()
  )
}

project_root <- find_project_root()
training_path <- file.path(project_root, "data", "processed", "heart_train.rds")
model_dir <- file.path(project_root, "results", "models")
plot_dir <- file.path(project_root, "results", "bayesian_network")

if (!file.exists(training_path)) {
  stop(
    "Training data not found at ", training_path,
    ". Run src/preprocessing.R first.",
    call. = FALSE
  )
}

heart_train <- readRDS(training_path)
if (anyNA(heart_train)) {
  stop("Training data contains missing values.", call. = FALSE)
}
if (!all(numeric_columns %in% names(heart_train)) ||
    !"target" %in% names(heart_train)) {
  stop("Training data does not match the expected schema.", call. = FALSE)
}

# Numeric variables are discretized from training data only. This avoids test
# leakage and produces a fully discrete network with inspectable CPTs.
discretization <- learn_discretization(heart_train, numeric_columns)
bn_train <- apply_discretization(heart_train, discretization)

node_names <- names(bn_train)
target_blacklist <- data.frame(
  from = rep("target", length(node_names) - 1L),
  to = setdiff(node_names, "target"),
  stringsAsFactors = FALSE
)

set.seed(random_seed)
hc_network <- bnlearn::hc(
  bn_train,
  score = "bde",
  iss = equivalent_sample_size,
  blacklist = target_blacklist,
  restart = 20,
  perturb = 5
)

set.seed(random_seed)
tabu_network <- bnlearn::tabu(
  bn_train,
  score = "bde",
  iss = equivalent_sample_size,
  blacklist = target_blacklist,
  tabu = 10,
  max.tabu = 100
)

expert_network <- build_expert_network(node_names)

networks <- list(
  hill_climbing = hc_network,
  tabu = tabu_network,
  expert = expert_network
)

network_scores <- vapply(
  networks,
  function(network) {
    bnlearn::score(
      network,
      data = bn_train,
      type = "bde",
      iss = equivalent_sample_size
    )
  },
  numeric(1)
)

learned_names <- c("hill_climbing", "tabu")
selected_name <- learned_names[[which.max(network_scores[learned_names])]]
selected_network <- networks[[selected_name]]

fitted_networks <- lapply(networks, function(network) {
  bnlearn::bn.fit(
    network,
    data = bn_train,
    method = "bayes",
    iss = equivalent_sample_size
  )
})

comparison <- data.frame(
  model = names(networks),
  structure = c("data-driven", "data-driven", "expert-defined"),
  bde_score = unname(network_scores),
  arc_count = vapply(
    networks,
    function(network) as.integer(bnlearn::narcs(network)),
    integer(1)
  ),
  selected_for_inference = names(networks) == selected_name,
  stringsAsFactors = FALSE
)

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  comparison,
  file.path(model_dir, "bn_structure_comparison.csv"),
  row.names = FALSE
)
write.csv(
  bnlearn::arcs(hc_network),
  file.path(model_dir, "hill_climbing_arcs.csv"),
  row.names = FALSE
)
write.csv(
  bnlearn::arcs(tabu_network),
  file.path(model_dir, "tabu_arcs.csv"),
  row.names = FALSE
)
write.csv(
  bnlearn::arcs(expert_network),
  file.path(model_dir, "expert_arcs.csv"),
  row.names = FALSE
)

model_bundle <- list(
  selected_name = selected_name,
  selected_network = selected_network,
  selected_fit = fitted_networks[[selected_name]],
  networks = networks,
  fitted_networks = fitted_networks,
  comparison = comparison,
  hc_tabu_structural_hamming_distance = bnlearn::shd(hc_network, tabu_network),
  discretization = discretization,
  training_columns = node_names,
  node_levels = lapply(bn_train, levels),
  target_levels = levels(bn_train$target),
  random_seed = random_seed,
  equivalent_sample_size = equivalent_sample_size
)
saveRDS(model_bundle, file.path(model_dir, "bayesian_network.rds"))

save_network_panel(
  networks,
  c(
    "Hill-climbing structure",
    "Tabu-search structure",
    "Expert-defined structure"
  ),
  file.path(plot_dir, "network_structure_comparison.png")
)
save_single_network(
  selected_network,
  paste("Selected structure:", gsub("_", " ", selected_name)),
  file.path(plot_dir, "selected_network.png")
)

message("Bayesian network training complete.")
message("Selected learned structure: ", selected_name)
message("Model bundle: ", file.path(model_dir, "bayesian_network.rds"))
message("Comparison plots: ", plot_dir)
