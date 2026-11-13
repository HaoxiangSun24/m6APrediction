test_that("dna_encoding returns five nucleotide factors", {
  encoded <- dna_encoding(c("GGACA", "GGAUA"))

  expect_s3_class(encoded, "data.frame")
  expect_equal(dim(encoded), c(2L, 5L))
  expect_equal(names(encoded), paste0("nt_pos", 1:5))
  expect_true(all(vapply(encoded, is.factor, logical(1L))))
  expect_equal(levels(encoded$nt_pos1), c("A", "T", "C", "G"))
  # U is converted to T
  expect_equal(as.character(encoded$nt_pos5[2]), "A")
  expect_equal(as.character(encoded$nt_pos4[2]), "T")
})

test_that("single-site prediction returns the expected named vector", {
  result <- prediction_single(
    gc_content = 0.65,
    RNA_type = "mRNA",
    RNA_region = "CDS",
    exon_length = 10.63,
    distance_to_junction = 7.66,
    evolutionary_conservation = 0.18,
    DNA_5mer = "GGACA"
  )

  expect_type(result, "character")
  expect_named(result, c("predicted_m6A_prob", "predicted_m6A_status"))
  expect_gte(as.numeric(result["predicted_m6A_prob"]), 0)
  expect_lte(as.numeric(result["predicted_m6A_prob"]), 1)
  expect_true(result["predicted_m6A_status"] %in% c("Positive", "Negative"))
})

test_that("batch prediction matches single-site prediction", {
  sites <- data.frame(
    gc_content = c(0.65, 0.52),
    RNA_type = c("mRNA", "lncRNA"),
    RNA_region = c("CDS", "3'UTR"),
    exon_length = c(10.63, 9.13),
    distance_to_junction = c(7.66, 7.59),
    evolutionary_conservation = c(0.18, 0.62),
    DNA_5mer = c("GGACA", "GAACA"),
    stringsAsFactors = FALSE
  )

  batch <- prediction_multiple(feature_df = sites)
  single <- prediction_single(
    gc_content = sites$gc_content[2],
    RNA_type = sites$RNA_type[2],
    RNA_region = sites$RNA_region[2],
    exon_length = sites$exon_length[2],
    distance_to_junction = sites$distance_to_junction[2],
    evolutionary_conservation = sites$evolutionary_conservation[2],
    DNA_5mer = sites$DNA_5mer[2]
  )

  expect_equal(nrow(batch), 2L)
  expect_true(all(
    c("predicted_m6A_prob", "predicted_m6A_status") %in% names(batch)
  ))
  expect_equal(
    batch$predicted_m6A_prob[2],
    as.numeric(single["predicted_m6A_prob"])
  )
})

test_that("invalid inputs produce informative errors", {
  expect_error(
    prediction_single(
      gc_content = 0.65, RNA_type = "mRNA", RNA_region = "CDS",
      exon_length = 10.63, distance_to_junction = 7.66,
      evolutionary_conservation = 0.18, DNA_5mer = "GGAA"
    ),
    "five characters"
  )

  expect_error(
    prediction_single(
      gc_content = 0.65, RNA_type = "miRNA", RNA_region = "CDS",
      exon_length = 10.63, distance_to_junction = 7.66,
      evolutionary_conservation = 0.18, DNA_5mer = "GGACA"
    ),
    "Invalid value"
  )

  expect_error(
    prediction_single(
      gc_content = 0.65, RNA_type = "mRNA", RNA_region = "CDS",
      exon_length = 10.63, distance_to_junction = 7.66,
      evolutionary_conservation = 0.18, DNA_5mer = "GTGCA"
    ),
    "third nucleotide"
  )

  expect_error(
    prediction_multiple(feature_df = data.frame(gc_content = 0.65)),
    "Missing required column"
  )
})

test_that("the bundled model can be loaded explicitly", {
  model_path <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
  expect_true(nzchar(model_path))
  expect_true(file.exists(model_path))

  model <- readRDS(model_path)
  expect_s3_class(model, "randomForest")

  result <- prediction_single(
    ml_fit = model,
    gc_content = 0.65, RNA_type = "mRNA", RNA_region = "CDS",
    exon_length = 10.63, distance_to_junction = 7.66,
    evolutionary_conservation = 0.18, DNA_5mer = "GGACA"
  )
  expect_equal(names(result), c("predicted_m6A_prob", "predicted_m6A_status"))
})
