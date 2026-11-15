test_that("marginal ordination boxplots use aplot panels", {
  scores <- data.frame(
    Axis1 = c(-0.25, -0.18, 0.08, 0.12, 0.21, 0.27),
    Axis2 = c(0.12, -0.08, 0.18, -0.12, 0.05, -0.02),
    SampleID = paste0("S", 1:6),
    Group = factor(rep(c("KO", "OE", "WT"), each = 2),
                   levels = c("KO", "OE", "WT"))
  )
  beta <- structure(
    list(scores=scores, parameters=list(ordination="NMDS"),
         stress=0.117, variance=NULL, permanova=NULL),
    class="aq_beta"
  )

  if (requireNamespace("aplot", quietly = TRUE)) {
    out <- aq_beta_plot(
      beta, method="NMDS", color="Group", shape="Group",
      marginal=TRUE, test="none"
    )
    expect_s3_class(out, "aq_plot")
    expect_false(is.null(out$plot))
    expect_null(out$plot$labels$title)
  } else {
    expect_error(
      aq_beta_plot(beta, method="NMDS", color="Group", marginal=TRUE),
      "requires the 'aplot' package"
    )
  }
})
