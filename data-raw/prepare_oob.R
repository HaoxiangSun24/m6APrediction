# Extract out-of-bag predictions from the original random forest fitted in
# Practical 6. The full model (with OOB votes) lives outside the package; only
# the compact prediction table is kept here for the analysis report.
suppressMessages(library(randomForest))

fit <- readRDS("C:/Users/Lenovo/Tianxi/tmp_m6a_app/rf_fit.rds")

oob <- data.frame(
  truth = fit$y,
  prob_positive = fit$votes[, "Positive"],
  stringsAsFactors = FALSE
)

saveRDS(oob, "data-raw/oob_predictions.rds", compress = "xz")

cat("rows:", nrow(oob), "\n")
print(table(oob$truth))
cat("size MB:", round(file.info("data-raw/oob_predictions.rds")$size / 1024^2, 2), "\n")
