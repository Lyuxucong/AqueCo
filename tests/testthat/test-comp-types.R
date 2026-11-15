test_that("aq_comp supports bar, heatmap, and bubble types", {
  otu <- matrix(
    c(10, 5, 3, 8, 4, 6, 7, 2, 9, 6, 8, 3),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("S", 1:4), paste0("ASV", 1:3))
  )
  taxonomy <- data.frame(
    Phylum = c("Firmicutes", "Proteobacteria", "Firmicutes"),
    Genus = c("GenusA", "GenusB", "GenusC"),
    row.names = colnames(otu)
  )
  metadata <- data.frame(
    Group = factor(c("A", "A", "B", "B")),
    row.names = rownames(otu)
  )
  dat <- aq_data(otu = otu, taxonomy = taxonomy, metadata = metadata)

  p_bar <- aq_comp(dat, type = "bar", level = "Phylum",
                   facet = "Group", facet_scales = "free_x")
  p_heat <- aq_comp(dat, type = "heatmap", level = "Genus",
                    facet = "Group", heatmap_scale = "row_z")
  p_heat_mean <- aq_comp(dat, type = "heatmap", level = "Genus",
                         summary = "mean", group = "Group",
                         heatmap_scale = "row_z")
  p_bubble <- aq_comp(dat, type = "bubble", level = "Genus",
                      facet = "Group", bubble_color = "taxon")
  p_bubble_mean <- aq_comp(dat, type = "bubble", level = "Genus",
                           summary = "mean", group = "Group",
                           bubble_color = "group")

  expect_s3_class(p_bar, "ggplot")
  expect_s3_class(p_heat, "ggplot")
  expect_s3_class(p_heat_mean, "ggplot")
  expect_s3_class(p_bubble, "ggplot")
  expect_s3_class(p_bubble_mean, "ggplot")
})
