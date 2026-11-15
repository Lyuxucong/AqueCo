test_that("environment analysis aligns env and metadata by sample IDs", {
  otu <- matrix(
    c(10, 4, 5, 8, 3, 7),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), c("ASV1", "ASV2"))
  )
  md <- data.frame(
    Group = factor(c("KO", "OE", "WT")),
    row.names = c("S3", "S1", "S2")
  )
  env <- data.frame(
    pH = c(7.0, 6.5, 7.5),
    TN = c(10, 20, 30),
    Note = c("a", "b", "c"),
    row.names = c("S2", "S3", "S1")
  )
  dat <- aq_data(otu, metadata = md, env = env)
  res <- aq_env(dat, group = "Group", posthoc = "both")

  expect_s3_class(res, "aq_env")
  expect_setequal(res$parameters$variables, c("pH", "TN"))
  expect_equal(as.character(res$data$Group[res$data$SampleID == "S1"][1]), "OE")
  expect_true(all(c("dunn_p_adjusted", "wilcox_p_adjusted") %in%
                    names(res$pairwise)))
})

test_that("environment plot supports boxplot and mean plus SE barplot", {
  skip_if_not_installed("ggplot2")
  otu <- matrix(
    c(10, 4, 5, 8, 3, 7, 9, 2, 6, 4, 8, 5),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("S", 1:6), c("ASV1", "ASV2"))
  )
  md <- data.frame(
    Group = factor(rep(c("KO", "OE", "WT"), each = 2)),
    row.names = paste0("S", 1:6)
  )
  env <- data.frame(
    pH = c(6.5, 6.7, 7.1, 7.2, 7.5, 7.6),
    TN = c(10, 11, 13, 14, 16, 17),
    row.names = paste0("S", 1:6)
  )
  dat <- aq_data(otu, metadata = md, env = env)
  res <- aq_env(dat, group = "Group", posthoc = "both")

  p1 <- aq_env_plot(res, type = "boxplot", facet = TRUE,
                    pairwise_label = TRUE)
  p2 <- aq_env_plot(res, type = "bar", bar_error = "SE", facet = TRUE)

  expect_s3_class(p1, "aq_plot")
  expect_s3_class(p2, "aq_plot")
  expect_false(is.null(p1$plot))
  expect_false(is.null(p2$plot))
  expect_equal(p2$parameters$bar_error, "SE")
})
