test_that("RDA and CCA functions are exported and validate AqueCo inputs", {
  expect_true("aq_rda_cca" %in% getNamespaceExports("aqueco"))
  expect_true("aq_rda_cca_plot" %in% getNamespaceExports("aqueco"))
  otu <- matrix(
    c(10, 2, 1, 0, 8, 3, 1, 0, 7, 4, 2, 0,
      0, 2, 10, 1, 0, 3, 8, 2, 1, 4, 7, 3,
      9, 1, 0, 2, 8, 2, 1, 1, 1, 5, 4, 2,
      0, 1, 9, 3, 1, 2, 8, 2, 2, 4, 5, 3),
    nrow = 12, byrow = TRUE,
    dimnames = list(paste0("S", 1:12), paste0("ASV", 1:4))
  )
  env <- data.frame(
    pH = c(6.1, 6.2, 6.3, 6.4, 7.1, 7.2, 7.3, 7.4, 6.0, 6.5, 7.0, 7.5),
    TN = c(10, 11, 12, 13, 20, 21, 22, 23, 9, 14, 19, 24),
    TP = c(3, 3.2, 3.4, 3.6, 8, 8.2, 8.4, 8.6, 2.8, 4, 7.5, 9),
    row.names = paste0("S", 1:12)
  )
  md <- data.frame(
    Group = factor(rep(c("KO", "OE", "WT"), each = 4)),
    row.names = paste0("S", 1:12)
  )
  dat <- aq_data(otu, metadata = md, env = env)
  expect_error(
    aq_rda_cca(dat, method = "RDA", env_variables = "missing",
               group = "Group", permutations = 19),
    "Unknown environmental"
  )
  expect_error(
    aq_rda_cca(dat, method = "CCA", transform = "hellinger",
               group = "Group", permutations = 19),
    "Hellinger transformation"
  )
})
