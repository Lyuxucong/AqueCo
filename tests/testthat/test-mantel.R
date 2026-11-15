test_that("group-wise Mantel analysis returns results and a combined plot", {
  otu <- matrix(
    c(10, 2, 0, 1,
      8, 3, 0, 1,
      6, 4, 1, 0,
      0, 2, 10, 1,
      0, 3, 8, 2,
      1, 4, 6, 3,
      9, 1, 0, 2,
      0, 1, 9, 3,
      2, 5, 4, 1),
    nrow = 9, byrow = TRUE,
    dimnames = list(paste0("S", 1:9), paste0("ASV", 1:4))
  )
  env <- data.frame(
    pH = c(6.1, 6.3, 6.5, 7.1, 7.3, 7.5, 6.2, 7.4, 7.6),
    TN = c(10, 12, 14, 20, 22, 24, 11, 23, 25),
    TP = c(3, 4, 5, 8, 9, 10, 3.5, 9.5, 10.5),
    row.names = paste0("S", 1:9)
  )
  md <- data.frame(
    Group = factor(c("KO", "KO", "KO", "OE", "OE", "OE", "WT", "WT", "WT")),
    row.names = paste0("S", 1:9)
  )
  dat <- aq_data(otu, metadata = md, env = env)
  res <- aq_mantel(
    dat, group = "Group", permutations = 19,
    prevalence = 0, community_method = "bray"
  )
  expect_s3_class(res, "aq_mantel")
  expect_true(is.data.frame(res$env_data))
  expect_true(all(c("Group", "Environment", "Mantel_r", "Mantel_p",
                    "Mantel_p_adjusted") %in% names(res$results)))
  expect_equal(nrow(res$results), 9)
  skip_if_not_installed("linkET")
  p <- aq_mantel_plot(res)
  expect_s3_class(p, "aq_plot")
  expect_false(is.null(p$plot))
  expect_equal(p$plot$theme$axis.text.x$angle, 90)
  p45 <- aq_mantel_plot(res, x_label_angle = 45, x_label_size = 8)
  expect_equal(p45$plot$theme$axis.text.x$angle, 45)
  expect_equal(p45$plot$theme$axis.text.x$size, 8)
  expect_error(aq_mantel_plot(res, x_label_angle = -5), "x_label_angle")
})


test_that("group-wise Mantel plot retains one link per environment-group pair", {
  otu <- matrix(
    c(10, 2, 0, 1, 8, 3, 0, 1, 6, 4, 1, 0,
      0, 2, 10, 1, 0, 3, 8, 2, 1, 4, 6, 3,
      9, 1, 0, 2, 0, 1, 9, 3, 2, 5, 4, 1),
    nrow = 9, byrow = TRUE,
    dimnames = list(paste0("S", 1:9), paste0("ASV", 1:4))
  )
  env <- data.frame(
    pH = c(6.1, 6.3, 6.5, 7.1, 7.3, 7.5, 6.2, 7.4, 7.6),
    TN = c(10, 12, 14, 20, 22, 24, 11, 23, 25),
    TP = c(3, 4, 5, 8, 9, 10, 3.5, 9.5, 10.5),
    row.names = paste0("S", 1:9)
  )
  md <- data.frame(
    Group = factor(rep(c("KO", "OE", "WT"), each = 3)),
    row.names = paste0("S", 1:9)
  )
  dat <- aq_data(otu, metadata = md, env = env)
  res <- aq_mantel(dat, group = "Group", permutations = 19, prevalence = 0)
  skip_if_not_installed("linkET")
  plt <- aq_mantel_plot(res)
  expect_s3_class(plt, "aq_plot")
  expect_equal(nrow(plt$results), length(unique(plt$results$Environment)) * 3)
})


test_that("Mantel plot exposes group labels and significance-only correlation marks", {
  fml <- names(formals(aq_mantel_plot))
  expect_true("show_group_labels" %in% fml)
  expect_true("group_label_size" %in% fml)
})
