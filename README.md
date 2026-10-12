# Heart Disease Regression Analysis

**Reproducible generalized linear model analysis of a clinical benchmark**, based on a
2023 university group project by **David Yang, Arman Jahangiri, and Mojtaba Kanani Sarcheshmeh**.

This repository preserves the original collaborative analysis and adds a
self-contained, reproducible reanalysis: exploratory plots, full logistic GLM,
stepwise logistic regression, LASSO feature selection, and held-out model comparison.

[**Read the rendered analysis report**](results/report.html) ·
[**Results and statistics**](results/RESULTS.md) ·
[**Reproduction workflow**](.github/workflows/reproduce.yml)

![Held-out ROC curves](results/figures/test_roc_curves.png)

![LASSO CV plot](results/figures/lasso_cv.png)

## Reproduced results

**1,025** source rows · **302** unique records · **723** exact duplicate rows removed ·
**241** training observations · **61** held-out test observations.

| Model | Holdout ROC AUC | Accuracy | F1 |
| --- | ---: | ---: | ---: |
| Full logistic GLM | **0.936** | **85.2%** | **0.870** |
| Stepwise logistic GLM | 0.922 | 82.0% | 0.836 |
| LASSO logistic GLM | 0.930 | 83.6% | 0.861 |

These are **reproduced measurements**, using a single, stratified 80/20 split
(random seed 2026) after removing exact duplicates. They are **not**
independent external validation or estimates of clinical performance.

See [the full metric table](results/model_metrics.csv),
[data audit](results/data_audit.csv), and
[analysis report](results/report.html).

![Comparison of holdout model performance](results/figures/model_comparison.png)

## Key methodological improvement: avoid duplicate leakage

The source file contains **1,025 rows but numerous identical records**.
The original team notebooks combine extensive exploratory modelling with
cross-validation on a prepared dataset; their original data dependencies include
`heart.xlsx` and `heart_processed.csv`, which were not in the source repository.

This reanalysis reads the included `data/heart.csv` directly, drops exact
duplicates *before* an 80/20 class-stratified split, uses five-fold
cross-validation on training data only for LASSO, and reports performance on
the untouched holdout. It makes no promise that a small historical dataset
represents prospective patients.

## What is included

| Location | Description |
| --- | --- |
| [scripts/run_analysis.R](scripts/run_analysis.R) | Full executable R analysis and figure generation |
| [report.Rmd](report.Rmd) | Reproducible analytical report |
| [results/](results/) | Generated metrics, model coefficients, predictions, PNG plots, HTML report |
| [data/heart.csv](data/heart.csv) | Original 1,025-row source dataset |
| [legacy/](legacy/) | Historical R Markdown sources, unchanged from the group repository |
| [.github/workflows/reproduce.yml](.github/workflows/reproduce.yml) | Automated execution, report rendering and results publishing |

## Reproduce

Requires **R**, the **glmnet** package, and optionally **rmarkdown** for HTML.

```bash
Rscript -e 'install.packages(c("glmnet", "rmarkdown"), repos="https://cloud.r-project.org")'
Rscript scripts/run_analysis.R
Rscript -e "rmarkdown::render('report.Rmd', output_file='report.html', output_dir='results')"
```

Open `results/report.html` locally for the compiled analysis. GitHub Actions
runs the same code on each relevant update and commits the generated results.
A fixed seed (2026) makes the single held-out split reproducible.

## Evaluation design

- **Predictors:** 13 clinical features; categorical integer codes modelled as factors.
- **Models:** full logistic GLM, AIC-reduced GLM, and L1-regularized logistic regression.
- **Data partition:** remove exact duplicates, then class-stratified holdout.
- **Selection:** stepwise feature selection and LASSO tuning restricted to training.
- **Metrics:** ROC AUC, accuracy, precision, recall, specificity, F1, Brier score.
- **Output:** numerical audit, held-out predictions, odds/log-odds diagnostics, figures.

## Provenance and collaborators

This is a collaborative academic project and **not individual-only work**.
The original 2023 team members were **David Yang, Arman Jahangiri, and
Mojtaba Kanani Sarcheshmeh**. Original files are preserved with attribution
from [davidhy8/heart-regression](https://github.com/davidhy8/heart-regression).
The modernized reproduction workflow is a subsequent implementation and does
not overwrite the original project methodology.

Data source: [Kaggle Heart Disease Dataset](https://www.kaggle.com/datasets/johnsmith88/heart-disease-dataset),
a repackaged Cleveland-derived clinical benchmark.

**Research and education only.** These estimates are not clinically validated,
cannot establish causation, and must not be used for diagnosis or treatment.
