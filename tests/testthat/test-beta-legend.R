test_that("ordination legends use metadata names, not conditional expressions", {
  scores <- data.frame(
    Axis1 = c(-0.2, 0.1, 0.3, -0.1, 0.2, 0.0),
    Axis2 = c(0.1, -0.1, 0.2, -0.2, 0.0, 0.15),
    SampleID = paste0("S", 1:6),
    Group = factor(rep(c("KO", "OE", "WT"), each=2))
  )
  beta <- structure(
    list(scores=scores, parameters=list(ordination="NMDS"),
         stress=0.117, variance=NULL, permanova=NULL),
    class="aq_beta"
  )

  out <- aq_beta_plot(
    beta, method="NMDS", color="Group", shape="Group",
    centroid=TRUE, ellipse=TRUE, test="none"
  )

  expect_identical(out$plot$labels$colour, "Group")
  expect_identical(out$plot$labels$shape, "Group")
  expect_false(grepl("^if\\s*\\(", out$plot$labels$colour))
  expect_false(grepl("^if\\s*\\(", out$plot$labels$shape))
})
