test_that("aq_alpha_plot runs statistical tests without Value lookup errors", {
  otu <- matrix(
    c(
      10, 5, 3, 8, 4, 6,
      8, 7, 3, 6, 5, 7,
      3, 8, 7, 2, 9, 6,
      4, 7, 8, 3, 8, 7,
      9, 4, 5, 7, 3, 8,
      2, 9, 7, 3, 8, 6
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("S", 1:6), paste0("ASV", 1:6))
  )
  taxonomy <- data.frame(
    Phylum = rep(c("Firmicutes", "Proteobacteria"), each = 3),
    row.names = colnames(otu)
  )
  metadata <- data.frame(
    Group = factor(rep(c("A", "B", "C"), each = 2)),
    Time = factor(rep(c("T1", "T2"), 3)),
    row.names = rownames(otu)
  )
  dat <- aq_data(otu, taxonomy = taxonomy, metadata = metadata)
  alpha <- aq_alpha(dat, index = "all", group = c("Group", "Time"))

  p_kw <- aq_alpha_plot(
    alpha, index = "Shannon", group = "Group",
    test = "kruskal", show_p = TRUE
  )
  p_wilcox <- aq_alpha_plot(
    alpha, index = "Shannon", group = "Group",
    test = "wilcox", show_p = FALSE
  )
  p_anova <- aq_alpha_plot(
    alpha, index = "Shannon", group = "Group",
    test = "anova", show_p = FALSE
  )
  p_multi <- aq_alpha_plot(
    alpha, index = c("Shannon", "Simpson"), group = "Group",
    facet = "Index", test = "kruskal", show_p = TRUE
  )

  expect_s3_class(p_kw, "aq_plot")
  expect_s3_class(p_wilcox, "aq_plot")
  expect_s3_class(p_anova, "aq_plot")
  expect_s3_class(p_multi, "aq_plot")
  expect_s3_class(p_kw$plot, "ggplot")
  expect_true(is.data.frame(p_kw$statistics))
  expect_true(all(c("p_value", "p_adjusted") %in% names(p_kw$statistics)))
})
