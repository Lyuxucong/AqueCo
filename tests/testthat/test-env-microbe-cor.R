test_that("environment-microbe correlations match samples and return associations", {
  otu <- matrix(
    c(10, 0, 5, 0,
      8, 2, 0, 0,
      0, 7, 3, 0,
      0, 5, 4, 1,
      2, 0, 0, 8,
      0, 0, 2, 9),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("S", 1:6), paste0("ASV", 1:4))
  )
  tax <- data.frame(
    Phylum = c("P1", "P1", "P2", "P2"),
    Genus = c("G1", "G2", "G3", "G4"),
    row.names = paste0("ASV", 1:4)
  )
  env <- data.frame(
    pH = c(6.1, 6.3, 6.8, 7.0, 7.3, 7.5),
    TN = c(10, 11, 12, 15, 16, 17),
    row.names = paste0("S", 1:6)
  )
  md <- data.frame(Group = rep(c("A", "B"), each = 3),
                   row.names = paste0("S", 1:6))
  dat <- aq_data(otu, taxonomy = tax, metadata = md, env = env)

  res <- aq_env_microbe_cor(
    dat, level = "Genus", method = "spearman",
    prevalence = 0, group = NULL
  )
  expect_s3_class(res, "aq_env_microbe_cor")
  expect_equal(nrow(res$results), 8)
  expect_true(all(c("Taxon", "Environment", "Correlation",
                    "p_value", "p_adjusted", "N") %in% names(res$results)))

  by_group <- aq_env_microbe_cor(
    dat, level = "Genus", method = "spearman",
    prevalence = 0, group = "Group"
  )
  expect_setequal(unique(by_group$results$Group), c("A", "B"))

  plt <- aq_env_microbe_cor_plot(res)
  expect_s3_class(plt, "aq_plot")
  expect_false(is.null(plt$plot))
})


test_that("group-specific prevalence filtering binds abundance tables safely", {
  otu <- matrix(
    c(10, 0, 0,
      9, 1, 0,
      8, 2, 0,
      0, 0, 10,
      0, 1, 9,
      0, 2, 8),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("S", 1:6), paste0("ASV", 1:3))
  )
  tax <- data.frame(
    Genus = c("G1", "G2", "G3"),
    row.names = paste0("ASV", 1:3)
  )
  env <- data.frame(
    pH = c(6.1, 6.2, 6.3, 7.1, 7.2, 7.3),
    row.names = paste0("S", 1:6)
  )
  md <- data.frame(
    Group = rep(c("A", "B"), each = 3),
    row.names = paste0("S", 1:6)
  )
  dat <- aq_data(otu, taxonomy = tax, metadata = md, env = env)
  res <- aq_env_microbe_cor(
    dat, level = "Genus", group = "Group",
    prevalence = 0.5, method = "spearman"
  )
  expect_s3_class(res, "aq_env_microbe_cor")
  expect_true(all(c("SampleID", "Group", "G1", "G2", "G3") %in%
                    names(res$abundance)))
})
