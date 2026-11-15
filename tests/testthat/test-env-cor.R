test_that("environment correlation returns matrices and adjusted pairwise results", {
  otu <- matrix(
    c(10, 4, 5, 8, 3, 7, 9, 2, 6, 4, 8, 5),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("S", 1:6), c("ASV1", "ASV2"))
  )
  env <- data.frame(
    pH = c(6, 6.2, 6.4, 7, 7.2, 7.4),
    TN = c(10, 11, 12, 15, 16, 17),
    TP = c(3, 2, 4, 7, 6, 8),
    row.names = paste0("S", 1:6)
  )
  md <- data.frame(Group = rep(c("A", "B"), each = 3),
                   row.names = paste0("S", 1:6))
  dat <- aq_data(otu, metadata = md, env = env)

  res <- aq_env_cor(dat, method = "spearman")
  expect_s3_class(res, "aq_env_cor")
  expect_equal(dim(res$correlation$All), c(3, 3))
  expect_equal(nrow(res$results), 3)

  by_group <- aq_env_cor(dat, method = "pearson", group = "Group")
  expect_equal(length(by_group$correlation), 2)
  expect_true(all(c("p_value", "p_adjusted", "N") %in% names(by_group$results)))

  plt <- aq_env_cor_plot(res)
  expect_s3_class(plt, "aq_plot")
  expect_false(is.null(plt$plot))
})
