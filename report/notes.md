# Cardio-AI findings and notes

## Reproducibility snapshot

The full workflow was executed successfully on 27 September 2026 with R 4.5.3, `bnlearn` 5.2.1, and random seed 509. The stratified split contains 242 training patients and 61 test patients. The complete dataset contains 138 patients without heart disease (45.5%) and 165 with heart disease (54.5%).

## Exploratory findings

The largest absolute numeric correlation is between age and maximum achieved heart rate (`r = -0.399`). Maximum heart rate and exercise-induced ST depression are also moderately negatively correlated (`r = -0.344`). The other numeric correlations are comparatively weak, so no numeric pair is close to being redundant.

## Bayesian-network structure

Hill-climbing and tabu search each learned 13 arcs and produced effectively identical BDe scores (`-2968.73`). The implementation selected the tabu structure using the full-precision scores. The 18-arc expert network had a lower BDe score (`-3013.32`), suggesting that the learned sparse structure fits these training data better under this scoring criterion.

## Held-out test results

| Model | Accuracy | AUC | Sensitivity | Specificity | F1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Bayesian network | 0.738 | 0.809 | 0.939 | 0.500 | 0.795 |
| Logistic regression | 0.803 | 0.861 | 0.758 | 0.857 | 0.806 |
| Decision tree | 0.803 | 0.744 | 0.788 | 0.821 | 0.813 |

Logistic regression has the strongest test AUC and ties the decision tree for the highest accuracy. The Bayesian network detects the largest share of positive cases, but its lower specificity produces more false positives. This sensitivity-specificity tradeoff is more informative than accuracy alone for this dataset.

## Limitations

- The dataset is small, and results come from one held-out split rather than repeated cross-validation.
- Discretizing numeric variables into tertiles makes CPTs tractable and interpretable but discards within-bin information.
- Bayesian-network probabilities use approximate Monte Carlo inference and may vary slightly across software versions.
- The learned relationships are probabilistic associations and should not be interpreted as causal medical effects.
- No probability calibration study or external clinical validation was performed.
- The Shiny app is an educational demonstration, not a medical device or clinical decision-support system.
