# ---------------------------------------------------------------------------
# make_figures.R
# Reproduces every figure and every number used in README.Rmd.
# Run from the package root:  Rscript data-raw/make_figures.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(randomForest)
  library(ggplot2)
  library(scales)
})

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
set.seed(215)

## ---------------------------------------------------------------------------
## 1. Load data, out-of-bag predictions and the fitted model
## ---------------------------------------------------------------------------
dat <- read.csv("data-raw/m6A_imbalanced.csv", stringsAsFactors = FALSE)
oob <- readRDS("data-raw/oob_predictions.rds")
fit <- readRDS("inst/extdata/rf_fit.rds")

dat$m6A_status <- factor(dat$m6A_status, levels = c("Negative", "Positive"))
oob$truth      <- factor(oob$truth,      levels = c("Negative", "Positive"))

pal_status <- c(Negative = "#2C7FB8", Positive = "#D7301F")

## ---------------------------------------------------------------------------
## 2. Summary statistics written to a text file
## ---------------------------------------------------------------------------
sink("data-raw/analysis_summary.txt")

cat("=== DATA ===\n")
cat("rows:", nrow(dat), " cols:", ncol(dat), "\n")
print(table(dat$m6A_status))
cat("prevalence:", round(mean(dat$m6A_status == "Positive") * 100, 1), "%\n\n")

cat("=== RNA_type ===\n")
print(table(dat$RNA_type))
cat("\npositivity rate by RNA_type (%):\n")
print(round(100 * prop.table(table(dat$RNA_type, dat$m6A_status), 1)[, "Positive"], 1))

cat("\n=== RNA_region ===\n")
print(table(dat$RNA_region))
cat("\npositivity rate by RNA_region (%):\n")
print(round(100 * prop.table(table(dat$RNA_region, dat$m6A_status), 1)[, "Positive"], 1))

cat("\n=== numeric summaries ===\n")
num_cols <- c("m6A_level", "gc_content", "exon_length", "distance_to_junction",
              "evolutionary_conservation")
print(round(sapply(dat[num_cols], function(x) c(min = min(x), q1 = quantile(x, .25),
                                                median = median(x), mean = mean(x),
                                                q3 = quantile(x, .75), max = max(x))), 3))

cat("\n=== DNA 5-mer composition (top 10) ===\n")
print(head(sort(table(dat$DNA_5mer), decreasing = TRUE), 10))
cat("distinct 5-mers:", length(unique(dat$DNA_5mer)), "\n")
cat("fraction with A at centre (position 3):",
    round(mean(substr(dat$DNA_5mer, 3, 3) == "A") * 100, 1), "%\n")

cat("\n=== FEATURE IMPORTANCE (MeanDecreaseGini) ===\n")
imp <- fit$importance   # stored MeanDecreaseGini matrix
imp_df <- data.frame(feature = rownames(imp),
                     gini = as.numeric(imp[, "MeanDecreaseGini"]),
                     row.names = NULL)
imp_df <- imp_df[order(imp_df$gini), ]
print(imp_df, row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. ROC, PRC, threshold sweep, calibration, confusion matrix
## ---------------------------------------------------------------------------
truth_pos <- oob$truth == "Positive"
p <- oob$prob_positive
n_pos <- sum(truth_pos)
n_neg <- sum(!truth_pos)

ord <- order(p, decreasing = TRUE)
tp <- cumsum(truth_pos[ord])
fp <- cumsum(!truth_pos[ord])
tpr <- tp / n_pos
fpr <- fp / n_neg
prec <- tp / (tp + fp)
rec <- tpr

# AUC by trapezoid over the ROC points
auc <- sum(diff(c(0, fpr)) * (c(0, tpr)[-1] + c(0, tpr)[-length(c(0, tpr))]) / 2)
# Average precision: step-wise integral of the precision-recall curve
ap <- sum(diff(c(0, rec)) * prec)

cat("\n=== OOB PERFORMANCE ===\n")
cat("AUC:", round(auc, 4), "\n")
cat("Average precision:", round(ap, 4), "\n")
cat("baseline prevalence:", round(n_pos / (n_pos + n_neg), 4), "\n")

thresholds <- seq(0.05, 0.95, by = 0.05)
sweep <- do.call(rbind, lapply(thresholds, function(th) {
  pred_pos <- p > th
  TP <- sum(pred_pos & truth_pos);  FP <- sum(pred_pos & !truth_pos)
  TN <- sum(!pred_pos & !truth_pos); FN <- sum(!pred_pos & truth_pos)
  data.frame(threshold = th,
             TP = TP, FP = FP, TN = TN, FN = FN,
             sensitivity = TP / (TP + FN),
             specificity = TN / (TN + FP),
             precision   = ifelse(TP + FP == 0, NA, TP / (TP + FP)),
             accuracy    = (TP + TN) / length(p),
             f1          = ifelse(TP + FP == 0, NA,
                                  2 * TP / (2 * TP + FP + FN)),
             predicted_positive = sum(pred_pos))
}))

cat("\n=== THRESHOLD SWEEP ===\n")
print(sweep, row.names = FALSE)

cm <- sweep[sweep$threshold == 0.5, ]
cat("\n=== CONFUSION MATRIX AT 0.5 ===\n")
print(matrix(c(cm$TN, cm$FP, cm$FN, cm$TP), nrow = 2, byrow = TRUE,
             dimnames = list(Predicted = c("Negative", "Positive"),
                             Truth     = c("Negative", "Positive"))))
cat("accuracy:", round(cm$accuracy, 4),
    " sensitivity:", round(cm$sensitivity, 4),
    " specificity:", round(cm$specificity, 4),
    " precision:", round(cm$precision, 4), "\n")

# calibration: 12 equal-width probability bins
cal <- data.frame(prob = p, truth = truth_pos)
cal$bin <- cut(cal$prob, breaks = seq(0, 1, by = 1 / 12),
               include.lowest = TRUE, right = FALSE)
cal_df <- do.call(rbind, lapply(split(cal, cal$bin), function(d) {
  data.frame(mid = mean(range(d$prob)),
             observed = mean(d$truth),
             n = nrow(d))
}))
cal_df <- cal_df[order(cal_df$mid), ]
cat("\n=== CALIBRATION (12 bins) ===\n")
print(cal_df, row.names = FALSE)

# recommended threshold maximising F1
best_f1 <- sweep[which.max(sweep$f1), ]
cat("\nbest F1 threshold:", best_f1$threshold, " F1:", round(best_f1$f1, 4),
    " sensitivity:", round(best_f1$sensitivity, 4),
    " specificity:", round(best_f1$specificity, 4), "\n")

sink()

## ---------------------------------------------------------------------------
## 4. Figures
## ---------------------------------------------------------------------------
theme_report <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

save_png <- function(plot, file, width, height) {
  ggsave(file, plot, width = width, height = height, dpi = 160, bg = "white")
}

# dependency-free 2 x 2 panel built with grid (patchwork is not installed)
save_2x2 <- function(p1, p2, p3, p4, file, width, height) {
  grDevices::png(file, width = width * 160, height = height * 160, res = 160)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p3, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p4, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  grDevices::dev.off()
}

## --- Figure 1: dataset overview -------------------------------------------
p1a <- ggplot(dat, aes(m6A_status, fill = m6A_status)) +
  geom_bar(width = 0.6) +
  geom_text(stat = "count", aes(label = comma(after_stat(count))),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = pal_status, guide = "none") +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .12))) +
  labs(title = "A. Class balance", x = NULL, y = "Number of sites")

p1b <- ggplot(dat, aes(RNA_type, fill = m6A_status)) +
  geom_bar(position = "dodge", width = 0.7) +
  scale_fill_manual(values = pal_status) +
  scale_y_continuous(labels = comma) +
  labs(title = "B. Transcript biotype", x = NULL, y = "Number of sites", fill = NULL)

p1c <- ggplot(dat, aes(RNA_region, fill = m6A_status)) +
  geom_bar(position = "dodge", width = 0.7) +
  scale_fill_manual(values = pal_status) +
  scale_y_continuous(labels = comma) +
  labs(title = "C. Transcript region", x = NULL, y = "Number of sites", fill = NULL)

p1d <- ggplot(dat, aes(m6A_level, fill = m6A_status)) +
  geom_histogram(bins = 40, colour = NA, alpha = 0.85) +
  scale_fill_manual(values = pal_status) +
  labs(title = "D. m6A modification level", x = "m6A_level", y = "Number of sites", fill = NULL)

p1a <- p1a + theme_report; p1b <- p1b + theme_report
p1c <- p1c + theme_report; p1d <- p1d + theme_report
save_2x2(p1a, p1b, p1c, p1d, "man/figures/fig01_dataset_overview.png", 10, 8)

## --- Figure 2: predictor distributions ------------------------------------
p2a <- ggplot(dat, aes(gc_content, colour = m6A_status, fill = m6A_status)) +
  geom_density(alpha = 0.18, linewidth = 0.7) +
  scale_colour_manual(values = pal_status) +
  scale_fill_manual(values = pal_status) +
  labs(title = "A. GC content", x = "gc_content", y = "Density", colour = NULL, fill = NULL)

p2b <- ggplot(dat, aes(exon_length, colour = m6A_status, fill = m6A_status)) +
  geom_density(alpha = 0.18, linewidth = 0.7) +
  scale_colour_manual(values = pal_status) +
  scale_fill_manual(values = pal_status) +
  labs(title = "B. Exon length (log scale)", x = "exon_length", y = "Density",
       colour = NULL, fill = NULL)

p2c <- ggplot(dat, aes(distance_to_junction, colour = m6A_status, fill = m6A_status)) +
  geom_density(alpha = 0.18, linewidth = 0.7) +
  scale_colour_manual(values = pal_status) +
  scale_fill_manual(values = pal_status) +
  labs(title = "C. Distance to junction (log scale)",
       x = "distance_to_junction", y = "Density", colour = NULL, fill = NULL)

p2d <- ggplot(dat, aes(evolutionary_conservation, colour = m6A_status, fill = m6A_status)) +
  geom_density(alpha = 0.18, linewidth = 0.7) +
  scale_colour_manual(values = pal_status) +
  scale_fill_manual(values = pal_status) +
  labs(title = "D. Evolutionary conservation", x = "evolutionary_conservation",
       y = "Density", colour = NULL, fill = NULL)

p2a <- p2a + theme_report; p2b <- p2b + theme_report
p2c <- p2c + theme_report; p2d <- p2d + theme_report
save_2x2(p2a, p2b, p2c, p2d, "man/figures/fig02_feature_distributions.png", 10, 8)

## --- Figure 3: feature importance -----------------------------------------
imp_df$feature <- factor(imp_df$feature, levels = imp_df$feature)
fig3 <- ggplot(imp_df, aes(gini, feature)) +
  geom_col(fill = "#3182BD", width = 0.7) +
  geom_text(aes(label = round(gini, 0)), hjust = -0.15, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Random forest variable importance",
       subtitle = "Mean decrease in Gini impurity across 500 trees",
       x = "MeanDecreaseGini", y = NULL) +
  theme_report
save_png(fig3, "man/figures/fig03_feature_importance.png", 8, 5)

## --- Figure 4: ROC ---------------------------------------------------------
roc_df <- data.frame(fpr = c(0, fpr, 1), tpr = c(0, tpr, 1))
fig4 <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_line(colour = "#D7301F", linewidth = 1) +
  annotate("text", x = 0.62, y = 0.18, hjust = 0,
           label = paste0("AUC = ", sprintf("%.3f", auc)), size = 4.2) +
  coord_equal() +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(title = "Out-of-bag ROC curve",
       subtitle = paste0(nrow(oob), " sites; ", comma(n_pos), " positive / ",
                         comma(n_neg), " negative"),
       x = "False positive rate", y = "True positive rate") +
  theme_report
save_png(fig4, "man/figures/fig04_roc_curve.png", 6.4, 6)

## --- Figure 5: precision-recall -------------------------------------------
pr_df <- data.frame(recall = c(0, rec), precision = c(1, prec))
fig5 <- ggplot(pr_df, aes(recall, precision)) +
  geom_hline(yintercept = mean(truth_pos), linetype = "dashed", colour = "grey60") +
  geom_line(colour = "#238B45", linewidth = 1) +
  annotate("text", x = 0.05, y = mean(truth_pos) + 0.03, hjust = 0,
           label = paste0("baseline = ", sprintf("%.3f", mean(truth_pos))),
           size = 3.4, colour = "grey35") +
  annotate("text", x = 0.62, y = 0.40, hjust = 0,
           label = paste0("AP = ", sprintf("%.3f", ap)), size = 4.2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(title = "Out-of-bag precision-recall curve",
       subtitle = "Classifier scores are class probabilities",
       x = "Recall (sensitivity)", y = "Precision") +
  theme_report
save_png(fig5, "man/figures/fig05_prc_curve.png", 6.4, 6)

## --- Figure 6: threshold trade-off ----------------------------------------
sweep_long <- rbind(
  data.frame(threshold = sweep$threshold, value = sweep$sensitivity, metric = "Sensitivity"),
  data.frame(threshold = sweep$threshold, value = sweep$specificity, metric = "Specificity"),
  data.frame(threshold = sweep$threshold, value = sweep$precision,   metric = "Precision"),
  data.frame(threshold = sweep$threshold, value = sweep$f1,          metric = "F1")
)
fig6 <- ggplot(sweep_long, aes(threshold, value, colour = metric)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey30") +
  annotate("text", x = 0.51, y = 0.06, hjust = 0, size = 3.2, colour = "grey30",
           label = "default 0.5") +
  scale_colour_manual(values = c(Sensitivity = "#D7301F", Specificity = "#2C7FB8",
                                 Precision = "#238B45", F1 = "#6A51A3")) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(title = "Decision-threshold trade-off",
       subtitle = "How the operating point changes with the positive-class cut-off",
       x = "Probability threshold", y = "Metric value", colour = NULL) +
  theme_report
save_png(fig6, "man/figures/fig06_threshold_tradeoff.png", 8, 5)

## --- Figure 7: calibration ------------------------------------------------
fig7 <- ggplot(cal_df, aes(mid, observed)) +
  geom_abline(colour = "grey60", linetype = "dashed") +
  geom_point(aes(size = n), colour = "#D7301F", alpha = 0.85) +
  geom_line(colour = "#D7301F", linewidth = 0.7) +
  scale_size_continuous(labels = comma, range = c(1.5, 6), name = "Sites") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(title = "Calibration of predicted probabilities",
       subtitle = "Observed positive rate versus mean predicted probability (12 bins)",
       x = "Mean predicted probability", y = "Observed positive rate") +
  theme_report
save_png(fig7, "man/figures/fig07_calibration.png", 6.4, 5.4)

## --- Figure 8: confusion matrix -------------------------------------------
cm_df <- data.frame(
  predicted = factor(rep(c("Negative", "Positive"), each = 2),
                     levels = c("Negative", "Positive")),
  truth = factor(rep(c("Negative", "Positive"), times = 2),
                 levels = c("Negative", "Positive")),
  n = c(cm$TN, cm$FN, cm$FP, cm$TP)
)
cm_max <- max(cm_df$n)
fig8 <- ggplot(cm_df, aes(predicted, truth, fill = n)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = comma(n), colour = n > 0.55 * cm_max),
            size = 6, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_fill_gradient(low = "#9ECAE1", high = "#08306B",
                      breaks = c(0, 10000, 20000, 30000), labels = comma) +
  labs(title = "Confusion matrix at threshold 0.5",
       subtitle = paste0("OOB accuracy = ", sprintf("%.1f%%", 100 * cm$accuracy)),
       x = "Predicted class", y = "Observed class", fill = "Sites") +
  guides(fill = guide_colourbar(barwidth = unit(4, "cm"),
                                barheight = unit(0.35, "cm"))) +
  theme_report +
  theme(panel.grid = element_blank(),
        legend.title.position = "top")
save_png(fig8, "man/figures/fig08_confusion_matrix.png", 6, 4.6)

cat("\nAll figures written to man/figures/.\n")
