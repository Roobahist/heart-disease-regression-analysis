#!/usr/bin/env Rscript
# Publication-style plots computed from the data and fitted results, never
# from hardcoded performance numbers. Adapted for the portfolio reproduction.
render_figures <- function(root) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install the ggplot2 package to render figures.")
  }
  library(ggplot2)
  out <- file.path(root, "results", "figures")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)

  # Accessible, restrained palette and one shared chart system.
  ink <- "#1A2B40"
  muted <- "#627286"
  teal <- "#136F82"
  orange <- "#C47B50"
  purple <- "#7773A5"
  grid <- "#E5EBF1"
  palette <- c("Full logistic GLM" = teal, "Stepwise logistic GLM" = orange,
               "LASSO logistic GLM" = purple)
  models <- names(palette)
  theme_report <- function() {
    theme_minimal(base_size = 14, base_family = "sans") +
      theme(
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        text = element_text(colour = ink),
        plot.title = element_text(size = 21, face = "bold", colour = ink,
                                  margin = margin(b = 8)),
        plot.subtitle = element_text(size = 11.7, colour = muted,
                                     lineheight = 1.2, margin = margin(b = 16)),
        plot.caption = element_text(size = 9.8, colour = muted,
                                    hjust = 0, margin = margin(t = 16)),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        panel.grid.major = element_line(colour = grid, linewidth = 0.45),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11, colour = muted),
        axis.text = element_text(size = 11, colour = ink),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        plot.margin = margin(22, 28, 23, 25)
      )
  }
  save_figure <- function(name, p, w = 11, h = 6.7) {
    ggsave(file.path(out, paste0(name, ".png")), plot = p, width = w,
           height = h, units = "in", dpi = 220, bg = "white")
    ggsave(file.path(out, paste0(name, ".svg")), plot = p, width = w,
           height = h, units = "in", bg = "white")
  }
  source_data <- read.csv(file.path(root, "data", "heart.csv"))
  dat <- unique(source_data)
  metrics <- read.csv(file.path(root, "results", "model_metrics.csv"))
  preds <- read.csv(file.path(root, "results", "test_predictions.csv"))
  coef_df <- read.csv(file.path(root, "results", "lasso_coefficients.csv"))
  cv <- read.csv(file.path(root, "results", "lasso_cv_curve.csv"))
  audit <- read.csv(file.path(root, "results", "data_audit.csv"))
  stopifnot(nrow(dat) == audit$unique_rows[1],
            setequal(metrics$model, models), setequal(unique(preds$model), models),
            nrow(preds) == sum(metrics$test_n),
            all(preds$predicted_probability >= 0 & preds$predicted_probability <= 1))
  sample_note <- sprintf("%d unique records · %d duplicates removed",
                         audit$unique_rows[1], audit$exact_duplicates_removed[1])
  test_note <- sprintf("%d held-out test observations · one split", metrics$test_n[1])

  # Outcome counts: absolute length, and counts/percentages labeled directly.
  count <- data.frame(outcome = c("No disease", "Disease"),
                      n = c(sum(dat$target == 0), sum(dat$target == 1)))
  count$outcome <- factor(count$outcome, levels = c("Disease", "No disease"))
  count$label <- sprintf("%d records  ·  %.1f%%", count$n, 100 * count$n / nrow(dat))
  p <- ggplot(count, aes(n, outcome, fill = outcome)) +
    geom_col(width = 0.52, show.legend = FALSE) +
    geom_text(aes(label = label), hjust = -0.12, size = 4.2,
              colour = ink, fontface = "bold") +
    scale_fill_manual(values = c("No disease" = teal, "Disease" = orange)) +
    scale_x_continuous(limits = c(0, max(count$n) * 1.72),
                       breaks = seq(0, 220, 40), expand = c(0, 0)) +
    labs(title = "Outcome composition", subtitle = sample_note,
         x = "Number of distinct records", y = NULL,
         caption = "Observed labels in the benchmark dataset; not population disease prevalence.") +
    theme_report() + theme(panel.grid.major.y = element_blank())
  save_figure("outcome_balance", p, h = 5)

  # Observed age distributions, with transparent raw points and IQR.
  dat$outcome <- factor(dat$target, levels = c(0, 1),
                        labels = c("No disease", "Disease"))
  group_labels <- c(sprintf("No disease · n=%d", sum(dat$target == 0)),
                    sprintf("Disease · n=%d", sum(dat$target == 1)))
  levels(dat$outcome) <- group_labels
  color_groups <- setNames(c(teal, orange), group_labels)
  set.seed(2026)
  p <- ggplot(dat, aes(outcome, age, fill = outcome)) +
    geom_violin(width = 0.77, alpha = 0.17, colour = NA,
                trim = FALSE, show.legend = FALSE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white",
                 colour = ink, linewidth = 0.75, show.legend = FALSE) +
    geom_jitter(aes(colour = outcome), width = 0.12, height = 0,
                alpha = 0.38, size = 1.6, show.legend = FALSE) +
    scale_fill_manual(values = color_groups) +
    scale_colour_manual(values = color_groups) +
    scale_y_continuous(breaks = seq(20, 80, 10)) +
    labs(title = "Age distribution across labels",
         subtitle = "Each point is a distinct record; violin, median and interquartile range",
         x = NULL, y = "Age in years",
         caption = "Descriptive comparison, unadjusted for other variables.") +
    theme_report() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(face = "bold", margin = margin(t = 9)))
  save_figure("age_distribution", p, h = 6.3)

  # Lower-triangle correlation matrix, with real coefficients and full scale.
  variables <- c("age", "trestbps", "chol", "thalach", "oldpeak")
  labels <- c("Age", "Resting BP", "Cholesterol", "Max heart rate", "ST depression")
  correlation <- cor(dat[, variables, drop = FALSE])
  tiles <- expand.grid(column = seq_along(variables), row = seq_along(variables))
  tiles <- tiles[tiles$row >= tiles$column, , drop = FALSE]
  tiles$r <- correlation[cbind(tiles$row, tiles$column)]
  p <- ggplot(tiles, aes(column, row, fill = r)) +
    geom_tile(colour = "white", linewidth = 2) +
    geom_text(aes(label = sprintf("%+.2f", r)),
              colour = ink, size = 4.6, fontface = "bold") +
    scale_fill_gradient2(low = teal, mid = "#F4F7FA", high = orange,
                         midpoint = 0, limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         name = "Pearson r") +
    scale_x_continuous(breaks = seq_along(variables), labels = labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(variables), labels = labels, expand = c(0, 0)) +
    coord_equal() +
    labs(title = "Relationships among continuous features",
         subtitle = paste("Pearson correlation ·", sample_note),
         x = NULL, y = NULL,
         caption = "Only numerical predictors are shown. Correlation does not establish causality.") +
    theme_report() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "right",
          axis.text.x = element_text(angle = 35, hjust = 1),
          legend.title = element_text(size = 11))
  save_figure("numeric_correlations", p, w = 11.7, h = 7.8)

  # ROC data collapse tied prediction thresholds into the correct step function.
  roc_one <- function(d, name) {
    o <- order(d$predicted_probability, decreasing = TRUE)
    y <- d$actual[o]
    prob <- d$predicted_probability[o]
    last_tie <- !duplicated(prob, fromLast = TRUE)
    data.frame(fpr = c(0, cumsum(y == 0)[last_tie] / sum(y == 0)),
               tpr = c(0, cumsum(y == 1)[last_tie] / sum(y == 1)),
               model = name)
  }
  roc <- do.call(rbind, lapply(models, function(name)
    roc_one(preds[preds$model == name, ], name)))
  roc$model <- factor(roc$model, levels = models)
  auc_labels <- setNames(paste0(models, "   |   AUC ",
                               sprintf("%.3f", metrics$auc[match(models, metrics$model)])),
                         models)
  p <- ggplot(roc, aes(fpr, tpr, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.8, colour = "#B5C0CC") +
    geom_step(linewidth = 1.12, direction = "hv") +
    scale_colour_manual(values = palette, breaks = models, labels = auc_labels) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    labs(title = "How well do the models separate the two labels?",
         subtitle = paste("ROC curves ·", test_note),
         x = "False positive rate", y = "True positive rate",
         caption = "AUC calculated from held-out predictions. Diagonal: chance-level discrimination.") +
    theme_report()
  save_figure("test_roc_curves", p, w = 11.6, h = 7.2)

  # Small multiple dot plots share a zero baseline; no misleading truncation.
  comparison <- do.call(rbind, lapply(c("auc", "accuracy", "f1"), function(key) {
    data.frame(model = metrics$model, metric = key, value = metrics[[key]])
  }))
  comparison$model <- factor(comparison$model, levels = rev(models))
  comparison$metric <- factor(comparison$metric,
                              levels = c("auc", "accuracy", "f1"),
                              labels = c("ROC AUC", "Accuracy", "F1 score"))
  p <- ggplot(comparison, aes(value, model, colour = model)) +
    geom_segment(aes(x = 0, xend = value, yend = model),
                 linewidth = 1.4, alpha = 0.28) +
    geom_point(size = 4.3) +
    geom_text(aes(x = value + 0.022, label = sprintf("%.3f", value)),
              hjust = 0, colour = ink, fontface = "bold", size = 3.8) +
    facet_wrap(~metric, nrow = 1) +
    scale_colour_manual(values = palette, guide = "none") +
    scale_x_continuous(limits = c(0, 1.17), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", ".25", ".50", ".75", "1.0"), expand = c(0, 0)) +
    labs(title = "Held-out performance, side by side",
         subtitle = test_note, x = "Score (0–1)", y = NULL,
         caption = "Accuracy and F1 use a 0.50 probability threshold. No statistical difference test was performed.") +
    theme_report() +
    theme(panel.grid.major.y = element_blank(), panel.spacing.x = unit(0.6, "lines"),
          strip.text = element_text(face = "bold", size = 12.3, colour = ink))
  save_figure("model_comparison", p, w = 14.8, h = 5.6)

  # Cross-validated LASSO AUC with real fold uncertainty and both lambdas.
  stopifnot(length(unique(cv$lambda_min)) == 1, length(unique(cv$lambda_1se)) == 1)
  cv$log_lambda <- log10(cv$lambda)
  cv <- cv[order(cv$log_lambda), ]
  selected_lambda <- unique(cv$lambda_1se)
  minimum_lambda <- unique(cv$lambda_min)
  p <- ggplot(cv, aes(log_lambda, cvm)) +
    geom_ribbon(aes(ymin = cvm - cvsd, ymax = cvm + cvsd),
                fill = teal, alpha = 0.15) +
    geom_line(colour = teal, linewidth = 1.4) +
    geom_vline(xintercept = log10(minimum_lambda), colour = orange,
               linetype = "dotted", linewidth = 1.0) +
    geom_vline(xintercept = log10(selected_lambda), colour = purple,
               linetype = "longdash", linewidth = 1.0) +
    labs(title = "Selecting the LASSO penalty",
         subtitle = "Five-fold training CV · mean ROC AUC with ±1 standard error",
         x = expression(log[10](lambda)), y = "Cross-validated ROC AUC",
         caption = sprintf(
           "Orange dotted: lambda.min = %.4f. Purple dashed: selected lambda.1se = %.4f. Test set excluded.",
           minimum_lambda, selected_lambda)) +
    theme_report()
  save_figure("lasso_cv", p, w = 11.6, h = 6.5)

  # Signed LASSO coefficients, labeled by the dataset's coded categories.
  term_labels <- c(
    sex1 = "Sex: code 1", cp1 = "Chest pain: code 1", cp2 = "Chest pain: code 2",
    cp3 = "Chest pain: code 3", thalach = "Maximum heart rate",
    exang1 = "Exercise angina: code 1", oldpeak = "ST depression",
    slope1 = "ST slope: code 1", slope2 = "ST slope: code 2",
    ca1 = "Major vessels: code 1", ca2 = "Major vessels: code 2",
    ca3 = "Major vessels: code 3", ca4 = "Major vessels: code 4",
    thal1 = "Thal status: code 1", thal2 = "Thal status: code 2",
    thal3 = "Thal status: code 3")
  coef_df$label <- ifelse(coef_df$feature %in% names(term_labels),
                          unname(term_labels[coef_df$feature]), coef_df$feature)
  coef_df$label <- factor(coef_df$label,
                          levels = coef_df$label[order(coef_df$coefficient)])
  coef_df$direction <- ifelse(coef_df$coefficient >= 0, "Positive", "Negative")
  extent <- max(1, max(abs(coef_df$coefficient)))
  p <- ggplot(coef_df, aes(coefficient, label, fill = direction)) +
    geom_col(width = 0.63, show.legend = FALSE) +
    geom_vline(xintercept = 0, colour = ink, linewidth = 0.85) +
    geom_text(aes(label = sprintf("%+.2f", coefficient),
                  hjust = ifelse(coefficient >= 0, -0.16, 1.17)),
              colour = ink, size = 3.7, fontface = "bold") +
    scale_fill_manual(values = c(Positive = teal, Negative = orange)) +
    scale_x_continuous(limits = c(-1.45 * extent, 1.45 * extent)) +
    labs(title = "LASSO coefficients retained by regularization",
         subtitle = sprintf("%d non-zero terms · lambda.1se selection", nrow(coef_df)),
         x = "Log-odds coefficient", y = NULL,
         caption = "Coded categorical terms compare with a reference level. Coefficients are not causal effects.") +
    theme_report() +
    theme(panel.grid.major.y = element_blank())
  save_figure("lasso_coefficients", p, w = 12.5,
              h = max(6.7, 3.3 + nrow(coef_df) * 0.36))

  # New evidence chart: every held-out prediction shown, not a smoothed
  # probability curve that obscures the limited test sample size.
  full <- preds[preds$model == "Full logistic GLM", ]
  full$outcome <- factor(full$actual, levels = c(0, 1),
                         labels = c("No disease", "Disease"))
  set.seed(2026)
  p <- ggplot(full, aes(outcome, predicted_probability, colour = outcome)) +
    geom_hline(yintercept = 0.50, linetype = "dashed",
               colour = "#8A99A9", linewidth = 0.75) +
    geom_boxplot(width = 0.28, colour = ink, fill = "white",
                 outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(width = 0.12, height = 0, size = 2.4, alpha = 0.78) +
    scale_colour_manual(values = c("No disease" = teal, Disease = orange),
                        guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.20),
                       labels = function(v) sprintf("%.0f%%", v * 100)) +
    labs(title = "Predicted probabilities on held-out records",
         subtitle = sprintf("Full logistic model · %s", test_note),
         x = "Observed outcome", y = "Predicted probability of disease",
         caption = "Each point is one test record; dashed line is the 0.50 classification threshold.") +
    theme_report() + theme(panel.grid.major.x = element_blank())
  save_figure("test_probability_distribution", p, h = 6.7)
  message("Regenerated 8 charts as high-resolution PNG and vector SVG.")
}

if (sys.nframe() == 0) {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) stop("Run Rscript scripts/plot_results.R")
  this_script <- normalizePath(sub("^--file=", "", arg[1]))
  render_figures(normalizePath(file.path(dirname(this_script), "..")))
}
