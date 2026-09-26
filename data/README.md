# Dataset

`heart.csv` is the Cleveland subset of the UCI Heart Disease dataset, using the redistributed CSV specified for this project:

<https://raw.githubusercontent.com/halekpetigo/BIOF509/main/heart%202.csv>

The downloaded file is retained unchanged. Its SHA-256 checksum is `7c3014365675306819510a49ff289efbec1d1a6a666a2dc7652f1547b383d859`.

## Verification

Verification performed after download confirmed:

- 303 data rows and one header row
- 14 columns in the expected order
- no empty, `NA`, `N/A`, or `?` values
- target values limited to `0` and `1`
- target distribution: 138 records with `target = 0` and 165 records with `target = 1`

The source file includes a UTF-8 byte-order mark and CRLF line endings. These encoding markers do not affect the data and were left unchanged.

## Data dictionary

All fields are numeric in the raw CSV. Variables marked categorical will be converted to factors during preprocessing.

| Column | Raw type | Modeling type | Meaning and encoding |
| --- | --- | --- | --- |
| `age` | Integer | Numeric | Age in years. |
| `sex` | Integer | Categorical | Sex: `0` = female, `1` = male. |
| `cp` | Integer | Categorical | Chest-pain type: `0` = typical angina, `1` = atypical angina, `2` = non-anginal pain, `3` = asymptomatic. |
| `trestbps` | Integer | Numeric | Resting blood pressure in mm Hg on admission. |
| `chol` | Integer | Numeric | Serum cholesterol in mg/dL. |
| `fbs` | Integer | Categorical | Fasting blood sugar above 120 mg/dL: `0` = no, `1` = yes. |
| `restecg` | Integer | Categorical | Resting ECG: `0` = normal, `1` = ST-T wave abnormality, `2` = probable or definite left-ventricular hypertrophy. |
| `thalach` | Integer | Numeric | Maximum heart rate achieved. |
| `exang` | Integer | Categorical | Exercise-induced angina: `0` = no, `1` = yes. |
| `oldpeak` | Decimal | Numeric | ST depression induced by exercise relative to rest. |
| `slope` | Integer | Categorical | Slope of the peak exercise ST segment: `0` = upsloping, `1` = flat, `2` = downsloping. |
| `ca` | Integer | Categorical | Number of major vessels colored by fluoroscopy, encoded `0`–`4` in this file. |
| `thal` | Integer | Categorical | Thallium stress-test result: `0` = unspecified source code, `1` = fixed defect, `2` = normal, `3` = reversible defect. |
| `target` | Integer | Categorical | Diagnosis: `0` = no heart disease, `1` = heart disease. |
