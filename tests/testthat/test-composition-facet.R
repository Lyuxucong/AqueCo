test_that("aq_comp supports free x scales and readable sample labels", {
  otu <- matrix(
    c(10, 5, 3, 8, 4, 6, 7, 2, 9),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("S1", "S2", "S3"), paste0("ASV", 1:3))
  )
  taxonomy <- data.frame(
    Phylum = c("Firmicutes", "Proteobacteria", "Firmicutes"),
    row.names = colnames(otu)
  )
  metadata <- data.frame(
    Group = factor(c("A", "A", "B")),
    row.names = rownames(otu)
  )
  dat <- aq_data(otu = otu, taxonomy = taxonomy, metadata = metadata)

  p <- aq_comp(
    dat, level = "Phylum", summary = "sample",
    facet = "Group", facet_scales = "free_x",
    x_angle = 45, x_size = 7
  )
  expect_s3_class(p, "ggplot")
})
