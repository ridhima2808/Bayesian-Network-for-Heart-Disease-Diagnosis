# Cardio-AI Shiny application
#
# Prerequisites:
#   Rscript src/preprocessing.R
#   Rscript src/bayesian_network.R
#
# Launch from the project root with:
#   shiny::runApp("app")

required_packages <- c("bnlearn", "ggplot2", "shiny")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required package(s) before launching the app: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

find_app_project_root <- function() {
  candidates <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[[1L]])
    script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
    candidates <- unique(c(candidates, script_dir, dirname(script_dir)))
  }

  for (candidate in candidates) {
    current <- candidate
    repeat {
      if (file.exists(file.path(current, "src", "inference.R")) &&
          file.exists(file.path(current, "README.md"))) {
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

project_root <- find_app_project_root()
inference_script <- file.path(project_root, "src", "inference.R")
model_path <- file.path(project_root, "results", "models", "bayesian_network.rds")

# Import the shared validation and inference functions without running the
# command-line worked examples.
source(inference_script, local = TRUE)
model_bundle <- load_bayesian_model(model_path)

inference_particles <- 15000L
inference_seed <- 509L
numeric_input_ranges <- lapply(
  model_bundle$discretization,
  function(specification) specification$training_range
)

factor_labels <- c(
  age = "Age",
  sex = "Sex",
  cp = "Chest-pain type",
  trestbps = "Resting blood pressure",
  chol = "Serum cholesterol",
  fbs = "Fasting blood sugar",
  restecg = "Resting ECG",
  thalach = "Maximum heart rate",
  exang = "Exercise-induced angina",
  oldpeak = "Exercise ST depression",
  slope = "Peak ST-segment slope",
  ca = "Major vessels observed",
  thal = "Thallium stress-test result"
)

patient_from_input <- function(input) {
  list(
    age = input$age,
    sex = input$sex,
    cp = input$cp,
    trestbps = input$trestbps,
    chol = input$chol,
    fbs = input$fbs,
    restecg = input$restecg,
    thalach = input$thalach,
    exang = input$exang,
    oldpeak = input$oldpeak,
    slope = input$slope,
    ca = input$ca,
    thal = input$thal
  )
}

validate_patient_ranges <- function(patient) {
  for (column_name in names(numeric_input_ranges)) {
    value <- patient[[column_name]]
    allowed_range <- numeric_input_ranges[[column_name]]
    if (length(value) != 1L || is.na(value) || !is.finite(value) ||
        value < allowed_range[[1L]] || value > allowed_range[[2L]]) {
      stop(
        factor_labels[[column_name]], " must be between ",
        allowed_range[[1L]], " and ", allowed_range[[2L]], ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

calculate_factor_contributions <- function(
    patient,
    full_diagnosis,
    model_bundle,
    n = inference_particles,
    seed = inference_seed) {
  feature_names <- names(patient)
  probability_without_feature <- vapply(seq_along(feature_names), function(index) {
    reduced_patient <- patient[-index]
    reduced_diagnosis <- diagnose_patient(
      reduced_patient,
      model_bundle = model_bundle,
      n = n,
      seed = seed
    )
    reduced_diagnosis$probability_disease
  }, numeric(1))

  contributions <- data.frame(
    feature = feature_names,
    factor = unname(factor_labels[feature_names]),
    entered_value = vapply(patient, as.character, character(1)),
    model_value = unname(unlist(full_diagnosis$evidence[feature_names])),
    probability_without = probability_without_feature,
    contribution = full_diagnosis$probability_disease - probability_without_feature,
    stringsAsFactors = FALSE
  )
  contributions$absolute_contribution <- abs(contributions$contribution)
  contributions$direction <- ifelse(
    contributions$contribution >= 0,
    "Raised estimate",
    "Lowered estimate"
  )
  contributions <- contributions[
    order(contributions$absolute_contribution, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  rownames(contributions) <- NULL
  contributions
}

risk_band <- function(probability) {
  if (probability < 0.33) {
    return(list(label = "Lower model estimate", class = "risk-low"))
  }
  if (probability < 0.66) {
    return(list(label = "Intermediate model estimate", class = "risk-medium"))
  }
  list(label = "Higher model estimate", class = "risk-high")
}

app_css <- paste(
  "body { background: #f4f7fa; color: #24313d; }",
  ".app-header { margin: 20px 0 18px; }",
  ".app-header h1 { margin-bottom: 4px; font-weight: 700; }",
  ".subtitle { color: #5e6b76; font-size: 16px; }",
  ".well { background: #ffffff; border: 1px solid #dce3e8; box-shadow: none; }",
  ".diagnose-button { width: 100%; margin-top: 12px; font-weight: 600; }",
  paste0(
    ".result-card { background: #ffffff; border: 1px solid #dce3e8; ",
    "border-radius: 8px; padding: 24px; margin-bottom: 18px; text-align: center; }"
  ),
  ".probability { font-size: 52px; line-height: 1.1; font-weight: 700; }",
  paste0(
    ".risk-label { display: inline-block; border-radius: 14px; padding: 5px 12px; ",
    "margin-top: 8px; font-weight: 600; }"
  ),
  ".risk-low { color: #176b42; background: #dff3e8; }",
  ".risk-medium { color: #825d09; background: #fff0c7; }",
  ".risk-high { color: #9b2c2c; background: #fde2e2; }",
  ".method-note { color: #64717d; margin-top: 12px; }",
  paste0(
    ".disclaimer { border-left: 4px solid #d97706; background: #fff8e8; ",
    "padding: 12px 14px; margin: 18px 0; }"
  ),
  paste0(
    ".section-card { background: #ffffff; border: 1px solid #dce3e8; ",
    "border-radius: 8px; padding: 18px; margin-bottom: 18px; }"
  ),
  sep = "\n"
)

ui <- shiny::fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML(app_css))
  ),
  shiny::div(
    class = "app-header",
    shiny::h1("Cardio-AI"),
    shiny::div(
      class = "subtitle",
      "Bayesian network demonstration for heart disease diagnosis"
    )
  ),
  shiny::div(
    class = "disclaimer",
    shiny::strong("Educational use only. "),
    "This course project is not a medical device and must not replace evaluation by a qualified clinician."
  ),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 4,
      shiny::h3("Patient information"),
      shiny::numericInput(
        "age", "Age (years)", value = 54,
        min = numeric_input_ranges$age[[1L]],
        max = numeric_input_ranges$age[[2L]],
        step = 1
      ),
      shiny::selectInput(
        "sex", "Sex",
        choices = c("Female" = "0", "Male" = "1"),
        selected = "1"
      ),
      shiny::selectInput(
        "cp", "Chest-pain type",
        choices = c(
          "Typical angina" = "0",
          "Atypical angina" = "1",
          "Non-anginal pain" = "2",
          "Asymptomatic" = "3"
        ),
        selected = "0"
      ),
      shiny::numericInput(
        "trestbps", "Resting blood pressure (mm Hg)",
        value = 130,
        min = numeric_input_ranges$trestbps[[1L]],
        max = numeric_input_ranges$trestbps[[2L]],
        step = 1
      ),
      shiny::numericInput(
        "chol", "Serum cholesterol (mg/dL)",
        value = 246,
        min = numeric_input_ranges$chol[[1L]],
        max = numeric_input_ranges$chol[[2L]],
        step = 1
      ),
      shiny::selectInput(
        "fbs", "Fasting blood sugar",
        choices = c("120 mg/dL or below" = "0", "Above 120 mg/dL" = "1")
      ),
      shiny::selectInput(
        "restecg", "Resting ECG",
        choices = c(
          "Normal" = "0",
          "ST-T wave abnormality" = "1",
          "Left-ventricular hypertrophy" = "2"
        )
      ),
      shiny::numericInput(
        "thalach", "Maximum heart rate achieved",
        value = 150,
        min = numeric_input_ranges$thalach[[1L]],
        max = numeric_input_ranges$thalach[[2L]],
        step = 1
      ),
      shiny::selectInput(
        "exang", "Exercise-induced angina",
        choices = c("No" = "0", "Yes" = "1")
      ),
      shiny::numericInput(
        "oldpeak", "Exercise-induced ST depression",
        value = 1,
        min = numeric_input_ranges$oldpeak[[1L]],
        max = numeric_input_ranges$oldpeak[[2L]],
        step = 0.1
      ),
      shiny::selectInput(
        "slope", "Peak ST-segment slope",
        choices = c("Upsloping" = "0", "Flat" = "1", "Downsloping" = "2")
      ),
      shiny::selectInput(
        "ca", "Major vessels observed by fluoroscopy",
        choices = setNames(as.character(0:4), as.character(0:4))
      ),
      shiny::selectInput(
        "thal", "Thallium stress-test result",
        choices = c(
          "Unspecified source code" = "0",
          "Fixed defect" = "1",
          "Normal" = "2",
          "Reversible defect" = "3"
        ),
        selected = "2"
      ),
      shiny::actionButton(
        "diagnose",
        "Estimate probability",
        class = "btn-primary diagnose-button"
      )
    ),
    shiny::mainPanel(
      width = 8,
      shiny::uiOutput("probability_card"),
      shiny::div(
        class = "section-card",
        shiny::h3("Most influential entered factors"),
        shiny::p(
          "Influence is estimated by removing one item of evidence at a time. ",
          "Positive values raised the model estimate; negative values lowered it."
        ),
        shiny::plotOutput("contribution_plot", height = "430px")
      ),
      shiny::div(
        class = "section-card",
        shiny::h3("Factor details"),
        shiny::tableOutput("contribution_table")
      )
    )
  )
)

server <- function(input, output, session) {
  diagnosis_result <- shiny::eventReactive(input$diagnose, {
    patient <- patient_from_input(input)

    tryCatch(
      {
        validate_patient_ranges(patient)
        shiny::withProgress(message = "Running Bayesian inference", value = 0, {
          full_diagnosis <- diagnose_patient(
            patient,
            model_bundle = model_bundle,
            n = inference_particles,
            seed = inference_seed
          )
          shiny::incProgress(1 / (length(patient) + 1), detail = "Estimating factor influence")

          contributions <- calculate_factor_contributions(
            patient,
            full_diagnosis,
            model_bundle,
            n = inference_particles,
            seed = inference_seed
          )
          shiny::incProgress(length(patient) / (length(patient) + 1))

          list(
            diagnosis = full_diagnosis,
            contributions = contributions
          )
        })
      },
      error = function(error) {
        shiny::showNotification(
          paste("Diagnosis could not be completed:", conditionMessage(error)),
          type = "error",
          duration = 10
        )
        NULL
      }
    )
  }, ignoreNULL = TRUE)

  output$probability_card <- shiny::renderUI({
    result <- diagnosis_result()
    shiny::req(result)

    probability <- result$diagnosis$probability_disease
    band <- risk_band(probability)
    shiny::div(
      class = "result-card",
      shiny::div("Estimated probability of heart disease"),
      shiny::div(class = "probability", sprintf("%.1f%%", 100 * probability)),
      shiny::div(class = paste("risk-label", band$class), band$label),
      shiny::div(
        class = "method-note",
        paste0(
          "Approximate likelihood-weighted inference using ",
          format(inference_particles, big.mark = ","),
          " particles."
        )
      )
    )
  })

  output$contribution_plot <- shiny::renderPlot({
    result <- diagnosis_result()
    shiny::req(result)

    plot_data <- utils::head(result$contributions, 8L)
    plot_data$factor <- factor(
      plot_data$factor,
      levels = rev(plot_data$factor)
    )
    plot_data$percentage_points <- 100 * plot_data$contribution

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = factor, y = percentage_points, fill = direction)
    ) +
      ggplot2::geom_col(width = 0.68) +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(yintercept = 0, color = "#4b5563", linewidth = 0.4) +
      ggplot2::scale_fill_manual(
        values = c("Raised estimate" = "#E45756", "Lowered estimate" = "#4C78A8")
      ) +
      ggplot2::labs(
        x = NULL,
        y = "Change in disease probability (percentage points)",
        fill = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.major.y = ggplot2::element_blank()
      )
  })

  output$contribution_table <- shiny::renderTable({
    result <- diagnosis_result()
    shiny::req(result)

    table_data <- utils::head(result$contributions, 8L)
    data.frame(
      Factor = table_data$factor,
      `Entered value` = table_data$entered_value,
      `BN category` = table_data$model_value,
      Direction = table_data$direction,
      `Probability change` = sprintf(
        "%+.2f percentage points",
        100 * table_data$contribution
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
}

shiny::shinyApp(ui = ui, server = server)
