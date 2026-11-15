test_that("VPA supports environmental groups and validates inputs", {
  otu <- matrix(
    c(10, 2, 1, 0, 8, 3, 1, 0, 7, 4, 2, 0,
      0, 2, 10, 1, 0, 3, 8, 2, 1, 4, 7, 3,
      9, 1, 0, 2, 8, 2, 1, 1, 1, 5, 4, 2,
      0, 1, 9, 3, 1, 2, 8, 2, 2, 4, 5, 3,
      7, 2, 2, 1, 6, 3, 2, 1, 5, 4, 3, 2),
    nrow = 15, byrow = TRUE,
    dimnames = list(paste0("S", 1:15), paste0("ASV", 1:4))
  )
  env <- data.frame(
    Temp = seq(10, 24, length.out = 15),
    TN = seq(1, 4, length.out = 15),
    TP = seq(2, 8, length.out = 15),
    row.names = paste0("S", 1:15)
  )
  dat <- aq_data(otu, env = env)
  expect_error(
    aq_vpa(dat, groups = list(Physical = "Temp", Chemical = "missing")),
    "not found"
  )
  expect_error(
    aq_vpa(dat, groups = list(Physical = "Temp")),
    "2 to 4 explanatory sets"
  )
})

test_that("VPA returns adjusted-R2 fractions for two environmental groups", {
  set.seed(42)
  n <- 18
  otu <- matrix(
    stats::rpois(n * 5, lambda = 8),
    nrow = n, ncol = 5,
    dimnames = list(paste0("S", seq_len(n)), paste0("ASV", 1:5))
  )
  env <- data.frame(
    Temp = rnorm(n), Depth = rnorm(n),
    TN = rnorm(n), TP = rnorm(n),
    row.names = paste0("S", seq_len(n))
  )
  dat <- aq_data(otu, env = env)
  res <- aq_vpa(
    dat,
    groups = list(Physical = c("Temp", "Depth"),
                  Chemical = c("TN", "TP"))
  )
  expect_s3_class(res, "aq_vpa")
  expect_equal(length(res$group_names), 2)
  expect_true(is.data.frame(res$fractions))
  expect_true("Adj.R.square" %in% names(res$fractions))
  expect_true(is.data.frame(res$table_fractions))
  expect_equal(res$varpart$tables, c("Physical", "Chemical"))
})
