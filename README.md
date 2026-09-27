# Cardio-AI: Bayesian Network for Heart Disease Diagnosis

Cardio-AI is a college machine learning and data mining project that uses the Cleveland subset of the UCI Heart Disease dataset to develop and evaluate a Bayesian network for heart disease diagnosis. The repository contains reproducible R scripts for preprocessing, exploratory analysis, structure and parameter learning, probabilistic inference, baseline comparison, evaluation, and an interactive Shiny application.

> **Educational use only:** This project is not a medical device and its output must not replace evaluation by a qualified healthcare professional.

## Dataset

The project uses exactly the redistributed Cleveland CSV specified for the course project:

<https://raw.githubusercontent.com/halekpetigo/BIOF509/main/heart%202.csv>

The checked-in file at `data/heart.csv` contains 303 observations and 14 columns, has no missing values, and has 138 records with `target = 0` and 165 with `target = 1`. See [`data/README.md`](data/README.md) for the checksum, verification results, encodings, and complete data dictionary.

## Requirements

- R 4.4.0 or newer
- An internet connection for the initial package installation
- A browser for the Shiny application

Install the packages used by the project from CRAN:

```r
install.packages(c(
  "bnlearn",
  "ggplot2",
  "pROC",
  "rpart",
  "shiny"
))
```

`rpart` is normally included with the recommended R packages, but listing it explicitly makes the setup requirement clear. The project uses base R for its stratified split and metric calculations, so `caret` and `caTools` are not required.

## Reproduce the project

Clone the repository, open a terminal in its root directory, and run the scripts in this order.

### 1. Preprocess and split the data

```bash
Rscript src/preprocessing.R
```

This script revalidates the source CSV, converts categorical fields to factors, and creates a reproducible stratified 80/20 training/test split using seed `509`. Factor-preserving RDS files and a split manifest are written under `data/processed/`.

### 2. Run exploratory data analysis

```bash
Rscript src/eda.R
```

Summary tables, class-balance checks, distributions, target comparisons, and the numeric correlation heatmap are written to `results/eda/`.

### 3. Train the Bayesian networks

```bash
Rscript src/bayesian_network.R
```

Numeric features are discretized into tertiles learned from training data only. The script compares hill-climbing and tabu-search structures using the BDe score, constructs an expert-defined structure, fits smoothed conditional probability tables, and saves the selected learned network to `results/models/bayesian_network.rds`. Network plots are written to `results/bayesian_network/`.

### 4. Run the inference examples

```bash
Rscript src/inference.R
```

The two worked patients demonstrate `P(target = 1 | evidence)` using likelihood-weighted Monte Carlo inference. Their probabilities are written to `results/inference/worked_examples.csv`. The same file exposes reusable `load_bayesian_model()`, `diagnose_patient()`, and `sample_target_distribution()` functions.

### 5. Train baseline models

```bash
Rscript src/baseline_models.R
```

Logistic regression and a decision tree are trained on the same split used by the Bayesian network. Models, test probabilities, logistic coefficients, and tree importance values are written to `results/models/`.

### 6. Evaluate all models

```bash
Rscript src/evaluate.R
```

The evaluation compares the Bayesian network, logistic regression, and decision tree on the held-out test set. It writes accuracy, AUC, sensitivity, specificity, precision, F1, confusion matrices, test predictions, ROC coordinates, and comparison plots to `results/evaluation/`.

Bayesian predictions use Monte Carlo inference, so the script fixes row-specific random seeds for reproducibility. Small numerical differences may still occur across R or package versions.

## Run the Shiny application

Complete preprocessing and Bayesian-network training first, then launch the app from the repository root:

```r
shiny::runApp("app")
```

Alternatively, run it from a shell:

```bash
R -e 'shiny::runApp("app")'
```

Enter the patient measurements and click **Estimate probability**. The app displays the inferred probability of heart disease and the eight most influential entered factors. Factor influence is computed by removing one item of evidence at a time and measuring the change in model probability; it is an associative model explanation, not a causal or clinical conclusion.

## Modeling approach

The project keeps continuous measurements numeric for EDA and the two baseline models. For the Bayesian network, those measurements are converted to low, medium, and high tertiles using cut points derived only from the training set. This produces a fully discrete network with inspectable conditional probability tables while avoiding test-set leakage.

Two data-driven structures are learned with hill-climbing and tabu search. The better of these by BDe score is selected for inference. A compact expert-defined network is fitted and retained for structural comparison. Bayesian parameter estimation uses an equivalent sample size of 10 to smooth sparse conditional probability tables.

All model comparisons use the same held-out test rows and a classification threshold of `0.50`.

## Project structure

```text
.
├── README.md
├── .gitignore
├── data/
│   ├── README.md
│   ├── heart.csv
│   └── processed/                  # Generated cleaned data and split files
├── src/
│   ├── preprocessing.R
│   ├── eda.R
│   ├── bayesian_network.R
│   ├── inference.R
│   ├── baseline_models.R
│   └── evaluate.R
├── app/
│   └── app.R
├── results/
│   ├── eda/                        # Generated EDA tables and figures
│   ├── bayesian_network/           # Generated network diagrams
│   ├── models/                     # Generated fitted models and metadata
│   ├── inference/                  # Generated worked-example probabilities
│   └── evaluation/                 # Generated metrics and comparison plots
└── report/
    └── notes.md                    # Short project findings and limitations
```

Generated directories appear after their corresponding scripts run.

## Troubleshooting

- **`Rscript: command not found`:** Install R 4.4.0 or newer and ensure its executable directory is on `PATH`.
- **Missing package error:** Run the `install.packages()` command above in the same R library used to execute the scripts.
- **Missing processed data:** Run `src/preprocessing.R` before any analysis or modeling script.
- **Missing Bayesian model:** Run `src/bayesian_network.R` before inference, evaluation, or the Shiny app.
- **App starts slowly after clicking Diagnose:** The probability and factor contributions use repeated Monte Carlo queries. This is expected on modest hardware.

To record the exact environment used for a report or reproducibility check, run:

```r
sessionInfo()
```
