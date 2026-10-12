#!/usr/bin/env Rscript
# Reproducible heart-disease GLM analysis, adapted from a 2023 collaborative project.
# All metrics are computed using a duplicate-free stratified holdout.
options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026)

cmd <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(cmd)) stop("Run with Rscript scripts/run_analysis.R from this repository.")
script <- normalizePath(sub("^--file=", "", cmd[1]))
root <- normalizePath(file.path(dirname(script), ".."))
dir.create(file.path(root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(root, "data", "heart.csv")
raw <- read.csv(source_file, check.names = FALSE)
cols <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
          "thalach", "exang", "oldpeak", "slope", "ca", "thal", "target")
if (!identical(names(raw), cols)) stop("Unexpected dataset schema.")
if (anyNA(raw)) stop("Dataset unexpectedly contains missing values.")
if (!all(raw$target %in% c(0, 1))) stop("Target is not binary.")
if (any(!vapply(raw, is.numeric, logical(1)))) stop("All source columns must be numeric.")

n_original <- nrow(raw)
data <- unique(raw)
n_unique <- nrow(data)
duplicates <- n_original - n_unique
if (n_unique < 150) stop("Insufficient distinct records for the demo.")

# The source uses codes for categorical clinical attributes. Treat these as
# factors rather than assuming that the codes represent continuous scales.
factor_columns <- c("sex", "cp", "fbs", "restecg", "exang", "slope", "ca", "thal")
for (column in factor_columns) {
  data[[column]] <- factor(data[[column]], levels = sort(unique(data[[column]])))
}

# Split AFTER deduplicating so exact duplicate rows never appear across
# train and test. Use a fixed random seed and class-stratified sampling.
train_idx <- sort(unlist(lapply(split(seq_len(nrow(data)), data$target), function(i) {
  sample(i, size = floor(length(i) * 0.8), replace = FALSE)
})))
test_idx <- setdiff(seq_len(nrow(data)), train_idx)
train <- data[train_idx, , drop = FALSE]
test <- data[test_idx, , drop = FALSE]
if (length(intersect(train_idx, test_idx)) > 0) stop("Overlapping train/test indices.")
if (any(duplicated(raw)) && duplicates == 0) stop("Duplicate check failed.")
if (anyDuplicated(rbind(train, test)) != 0) stop("Duplicated records leaked across folds.")

# Model selection and penalization are fitted entirely inside train.
full <- suppressWarnings(glm(target ~ ., data = train, family = binomial(),
                             control = glm.control(maxit = 100)))
selected <- suppressWarnings(step(full, direction = "backward", trace = 0))
if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("The glmnet package is required (see README and CI workflow).")
}
x <- model.matrix(target ~ ., data = data)[, -1, drop = FALSE]
x_train <- x[train_idx, , drop = FALSE]
x_test <- x[test_idx, , drop = FALSE]
y_train <- train$target
y_test <- test$target

# Class-stratified inner CV, with no information from the test partition.
set.seed(2026)
fold_id <- integer(nrow(train))
for (class in sort(unique(y_train))) {
  positions <- which(y_train == class)
  fold_id[positions] <- sample(rep(seq_len(5), length.out = length(positions)))
}
cv_lasso <- glmnet::cv.glmnet(x_train, y_train, family = "binomial", alpha = 1,
                             foldid = fold_id, type.measure = "auc",
                             standardize = TRUE)

probabilities <- list(
  "Full logistic GLM" = as.numeric(predict(full, newdata = test, type = "response")),
  "Stepwise logistic GLM" = as.numeric(predict(selected, newdata = test, type = "response")),
  "LASSO logistic GLM" = as.numeric(predict(cv_lasso, newx = x_test,
                                            s = "lambda.1se", type = "response"))
)

auc_score <- function(y, p) {
  positives <- sum(y == 1)
  negatives <- sum(y == 0)
  if (positives == 0 || negatives == 0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) -
     positives * (positives + 1) / 2) / (positives * negatives)
}
safe_divide <- function(a, b) if (b == 0) NA_real_ else a / b
metrics_one <- function(name, p) {
  predicted <- as.integer(p >= 0.5)
  tp <- sum(predicted == 1 & y_test == 1)
  fp <- sum(predicted == 1 & y_test == 0)
  tn <- sum(predicted == 0 & y_test == 0)
  fn <- sum(predicted == 0 & y_test == 1)
  precision <- safe_divide(tp, tp + fp)
  recall <- safe_divide(tp, tp + fn)
  data.frame(
    model = name,
    auc = auc_score(y_test, p),
    accuracy = (tp + tn) / length(y_test),
    precision = precision,
    recall = recall,
    specificity = safe_divide(tn, tn + fp),
    f1 = safe_divide(2 * precision * recall, precision + recall),
    brier = mean((p - y_test)^2),
    test_n = length(y_test),
    threshold = 0.5,
    stringsAsFactors = FALSE
  )
}
metrics <- do.call(rbind, Map(metrics_one, names(probabilities), probabilities))
write.csv(metrics, file.path(root, "results", "model_metrics.csv"),
          row.names = FALSE, na = "")
audit <- data.frame(
  original_rows = n_original,
  unique_rows = n_unique,
  exact_duplicates_removed = duplicates,
  train_rows = nrow(train),
  test_rows = nrow(test),
  original_prevalence = mean(raw$target),
  deduplicated_prevalence = mean(data$target),
  test_prevalence = mean(y_test)
)
write.csv(audit, file.path(root, "results", "data_audit.csv"), row.names = FALSE)
write.csv(data.frame(
  model = rep(names(probabilities), each = length(y_test)),
  observation = rep(seq_along(y_test), times = length(probabilities)),
  actual = rep(y_test, times = length(probabilities)),
  predicted_probability = unlist(probabilities, use.names = FALSE)
), file.path(root, "results", "test_predictions.csv"), row.names = FALSE)

fig_dir <- file.path(root, "results", "figures")
save_png <- function(name, draw) {
  png(file.path(fig_dir, name), width = 1200, height = 780, res = 130)
  tryCatch(draw(), finally = dev.off())
}
save_png("outcome_balance.png", function() {
  counts <- table(factor(data$target, levels = c(0, 1),
                         labels = c("No disease", "Disease")))
  par(mar = c(5, 5, 4, 2))
  bars <- barplot(counts, col = c("#657b92", "#c77c64"), ylim = c(0, max(counts) * 1.2),
                  ylab = "Distinct records", main = "Target distribution after deduplication")
  text(bars, as.numeric(counts) + 5, labels = as.integer(counts), cex = 1.2)
})
save_png("age_distribution.png", function() {
  par(mar = c(5, 5, 4, 2))
  boxplot(age ~ factor(target, labels = c("No disease", "Disease")), data = data,
          col = c("#d7e1ea", "#f1d1bd"), ylab = "Age (years)",
          xlab = "Dataset outcome", main = "Age and heart disease label")
  stripchart(age ~ factor(target, labels = c("No disease", "Disease")),
             data = data, add = TRUE, method = "jitter", jitter = 0.13,
             pch = 16, cex = 0.35, col = rgb(0.1, 0.2, 0.3, alpha = 0.24))
})
save_png("numeric_correlations.png", function() {
  numeric_fields <- c("age", "trestbps", "chol", "thalach", "oldpeak")
  mat <- cor(data[, numeric_fields])
  par(mar = c(9, 9, 4, 4))
  image(seq_along(numeric_fields), seq_along(numeric_fields), mat,
        xlab = "", ylab = "", axes = FALSE, col = colorRampPalette(
          c("#355c7d", "#f7f7f7", "#ce6a58"))(100), zlim = c(-1, 1),
        main = "Pearson correlations (continuous attributes)")
  axis(1, at = seq_along(numeric_fields), labels = numeric_fields, las = 2)
  axis(2, at = seq_along(numeric_fields), labels = numeric_fields, las = 2)
  for (i in seq_along(numeric_fields)) {
    for (j in seq_along(numeric_fields)) {
      text(i, j, sprintf("%.2f", mat[i, j]), cex = 0.75)
    }
  }
})
roc_points <- function(y, p) {
  ord <- order(p, decreasing = TRUE)
  p <- p[ord]
  y <- y[ord]
  cbind(fpr = c(0, cumsum(y == 0) / sum(y == 0)),
        tpr = c(0, cumsum(y == 1) / sum(y == 1)))
}
save_png("test_roc_curves.png", function() {
  par(mar = c(5, 5, 4, 2))
  plot(c(0, 1), c(0, 1), type = "n", xlab = "False positive rate",
       ylab = "True positive rate", main = "ROC curves on untouched holdout test set",
       xaxs = "i", yaxs = "i")
  abline(0, 1, lty = 3, col = "gray60")
  palette <- c("#416e90", "#bb704b", "#4c997d")
  for (i in seq_along(probabilities)) {
    roc <- roc_points(y_test, probabilities[[i]])
    lines(roc[, "fpr"], roc[, "tpr"], col = palette[i], lwd = 3)
  }
  legend("bottomright", legend = sprintf("%s (AUC %.3f)", names(probabilities), metrics$auc),
         col = palette, lwd = 3, cex = 0.8, bg = "white")
})
save_png("lasso_cv.png", function() {
  par(mar = c(5, 5, 4, 2))
  plot(cv_lasso, main = "LASSO: five-fold CV on training set only")
})
coefs <- as.matrix(coef(cv_lasso, s = "lambda.1se"))
coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
coefs <- coefs[abs(coefs[, 1]) > 1e-8, , drop = FALSE]
coefficients <- data.frame(feature = rownames(coefs), coefficient = coefs[, 1],
                           row.names = NULL)
write.csv(coefficients, file.path(root, "results", "lasso_coefficients.csv"),
          row.names = FALSE)
save_png("lasso_coefficients.png", function() {
  par(mar = c(5, 13, 4, 2))
  if (!nrow(coefficients)) {
    plot.new()
    title(main = "No retained LASSO coefficients at lambda.1se")
  } else {
    order_idx <- order(coefficients$coefficient)
    vals <- coefficients$coefficient[order_idx]
    labs <- coefficients$feature[order_idx]
    barplot(vals, horiz = TRUE, names.arg = labs, las = 1, cex.names = 0.8,
            col = ifelse(vals >= 0, "#4c997d", "#bb704b"),
            main = "LASSO log-odds coefficients (lambda.1se)",
            xlab = "Coefficient on model-matrix scale")
    abline(v = 0, col = "gray30")
  }
})
save_png("model_comparison.png", function() {
  par(mar = c(8, 5, 4, 2))
  barplot(rbind(metrics$auc, metrics$accuracy), beside = TRUE,
          names.arg = c("Full logistic", "Stepwise logistic", "LASSO logistic"),
          las = 2, col = c("#416e90", "#c77c64"), ylim = c(0, 1),
          ylab = "Holdout score", main = "Held-out predictive performance")
  legend("topright", legend = c("ROC AUC", "Accuracy"),
         fill = c("#416e90", "#c77c64"), cex = 0.85)
})
fmt <- function(x) sprintf("%.3f", x)
lines <- c(
  "# Reproduced results",
  "",
  sprintf("**Source dataset:** %s rows; **distinct records:** %s; **exact duplicates removed:** %s.",
          n_original, n_unique, duplicates),
  sprintf("**Stratified split:** %s training records, %s held-out testing records (seed 2026).",
          nrow(train), nrow(test)),
  "",
  "All model selection and cross-validation occur on the training set. The test set is used",
  "only for final evaluation. No clinical deployment or external validation is claimed.",
  "",
  "| Model | Holdout ROC AUC | Accuracy | Precision | Recall | F1 | Brier |",
  "|---|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(metrics))) {
  m <- metrics[i, ]
  lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                            m$model, fmt(m$auc), fmt(m$accuracy),
                            fmt(m$precision), fmt(m$recall), fmt(m$f1), fmt(m$brier)))
}
lines <- c(lines,
  "",
  "## Interpretation and limitations",
  "",
  "- AUC evaluates ranking; accuracy at a fixed 0.50 threshold can conceal calibration problems.",
  "- Confidence intervals are not supplied: this is a small, single held-out split.",
  "- The 1,025-row dataset has many exact repeats. Splitting those repeats across folds",
  "  can overstate generalization; only distinct records are analysed here.",
  "- The source is a repackaged Cleveland-derived dataset, not a current clinical cohort.",
  "- Predictors are observational and do not support causal risk-factor claims.",
  "- All models are educational; predictions must not be used for diagnosis.",
  "",
  "## Figures",
  "",
  "![Held-out ROC curves](figures/test_roc_curves.png)",
  "![Model comparison](figures/model_comparison.png)",
  "![LASSO cross-validation](figures/lasso_cv.png)",
  "![LASSO coefficients](figures/lasso_coefficients.png)",
  "![Age distribution](figures/age_distribution.png)",
  "![Outcome balance](figures/outcome_balance.png)",
  "![Numeric correlations](figures/numeric_correlations.png)"
)
writeLines(lines, file.path(root, "results", "RESULTS.md"))
cat("Analysis complete: ", n_original, " original; ", n_unique,
    " unique; ", nrow(test), " held-out test records.\n", sep = "")
print(metrics, row.names = FALSE)
