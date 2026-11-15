test_that("beta decomposition calculates pairwise components within every group", {
  otu <- matrix(
    c(
      10, 0, 5, 0,
      8, 2, 0, 0,
      0, 7, 3, 0,
      0, 5, 4, 1,
      2, 0, 0, 8,
      0, 0, 2, 9
    ),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("S", 1:6), paste0("ASV", 1:4))
  )
  md <- data.frame(
    Group = factor(rep(c("KO", "OE", "WT"), each = 2)),
    row.names = paste0("S", 1:6)
  )
  dat <- aq_data(otu, metadata = md)

  res <- aq_beta_decomp(dat, group = "Group", method = "sorensen")

  expect_s3_class(res, "aq_beta_decomp")
  expect_equal(nrow(res$data), 3)
  expect_setequal(as.character(res$data$Group), c("KO", "OE", "WT"))
  expect_equal(nrow(res$means), 3)
  expect_true(all(abs(
    res$data$Similarity_pct +
      res$data$Turnover_pct +
      res$data$Nestedness_pct - 100
  ) < 1e-6))
  expect_true(all(res$data$Turnover_pct >= 0))
  expect_true(all(res$data$Nestedness_pct >= 0))
})

test_that("ternary plotting returns a combined plot and one plot per group", {
  skip_if_not_installed("ggtern")
  otu <- matrix(
    c(
      10, 0, 5, 0,
      8, 2, 0, 0,
      0, 7, 3, 0,
      0, 5, 4, 1,
      2, 0, 0, 8,
      0, 0, 2, 9
    ),
    nrow = 6, byrow = TRUE,
    dimnames = list(paste0("S", 1:6), paste0("ASV", 1:4))
  )
  md <- data.frame(
    Group = factor(rep(c("KO", "OE", "WT"), each = 2)),
    row.names = paste0("S", 1:6)
  )
  dat <- aq_data(otu, metadata = md)
  res <- aq_beta_decomp(dat, group = "Group")
  plt <- aq_beta_decomp_plot(res)

  expect_s3_class(plt, "aq_plot")
  expect_false(is.null(plt$plot))
  expect_true(inherits(plt$plot, "ggplot"))
  expect_true(all(vapply(plt$group_plots, inherits, logical(1), "ggplot")))
  expect_equal(length(plt$group_plots), 3)
  expect_setequal(names(plt$group_plots), c("KO", "OE", "WT"))
})


test_that("ternary mean point layer works without non-exported ggtern geom", {
  skip_if_not_installed("ggtern")
  otu <- matrix(c(10,0,5,0, 8,2,0,0, 0,7,3,0,
                  0,5,4,1, 2,0,0,8, 0,0,2,9),
                nrow=6,byrow=TRUE,
                dimnames=list(paste0("S",1:6),paste0("ASV",1:4)))
  md <- data.frame(Group=rep(c("KO","OE","WT"),each=2),
                   row.names=paste0("S",1:6))
  res <- aq_beta_decomp(aq_data(otu,metadata=md),group="Group")
  plotted <- aq_beta_decomp_plot(res)
  expect_true(inherits(plotted$plot,"ggplot"))
  expect_true(length(plotted$group_plots)==3L)
})
