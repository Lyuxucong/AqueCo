test_that("aq_import accepts OTU data already read with row.names = 1", {
  otu <- data.frame(
    S1 = c(10, 4),
    S2 = c(5, 8),
    row.names = c("ASV1", "ASV2"),
    check.names = FALSE
  )
  tax <- data.frame(Genus = c("G1", "G2"),
                    row.names = c("ASV1", "ASV2"))
  meta <- data.frame(Group = c("A", "B"),
                     row.names = c("S1", "S2"))

  x <- aq_import(otu, taxonomy = tax, metadata = meta)
  expect_equal(dim(x$otu), c(2L, 2L))
  expect_equal(rownames(x$otu), c("S1", "S2"))
  expect_equal(colnames(x$otu), c("ASV1", "ASV2"))
  expect_equal(unname(x$otu["S1", "ASV1"]), 10)
})

test_that("aq_import does not remove the first sample column", {
  otu <- data.frame(
    S1 = c(10, 4),
    S2 = c(5, 8),
    row.names = c("ASV1", "ASV2"),
    check.names = FALSE
  )
  x <- aq_import(otu)
  expect_equal(ncol(x$otu), 2L)
  expect_true(all(c("S1", "S2") %in% rownames(x$otu)))
})

test_that("aq_import rejects OTU data without taxon row names", {
  otu <- data.frame(S1 = c(10, 4), S2 = c(5, 8))
  expect_error(aq_import(otu), "row.names = 1")
})

test_that("aq_theme_set sets a consistent theme_bw-based default", {
  old <- ggplot2::theme_get()
  on.exit(ggplot2::theme_set(old), add = TRUE)
  aq_theme_set()
  expect_s3_class(ggplot2::theme_get(), "theme")
})

test_that("aq_comp venn type dispatches to aq_venn", {
  skip_if_not_installed("ggVennDiagram")
  otu <- matrix(c(5, 2, 4, 8, 3, 7), nrow = 3,
                dimnames = list(c("S1", "S2", "S3"), c("ASV1", "ASV2")))
  tax <- data.frame(Genus = c("G1", "G2"), row.names = colnames(otu))
  meta <- data.frame(Group = c("A", "A", "B"), row.names = rownames(otu))
  dat <- aq_data(otu, taxonomy = tax, metadata = meta)
  v <- aq_comp(dat, type = "venn", level = "Genus", group = "Group")
  expect_s3_class(v, "aq_venn")
  expect_s3_class(v$plot, "ggplot")
})


test_that("aq_import uses fixed row.names = 1 input format", {
  expect_false("otu_orientation" %in% names(formals(aq_import)))

  otu_tab <- data.frame(
    S1 = c(10, 4),
    S2 = c(5, 8),
    row.names = c("ASV1", "ASV2"),
    check.names = FALSE
  )
  tax <- data.frame(Genus = c("G1", "G2"),
                    row.names = c("ASV1", "ASV2"))
  meta <- data.frame(Group = c("A", "B"),
                     row.names = c("S1", "S2"))

  x <- aq_import(otu_tab, taxonomy = tax, metadata = meta)
  expect_equal(dim(x$otu), c(2L, 2L))
  expect_equal(rownames(x$otu), c("S1", "S2"))
  expect_equal(colnames(x$otu), c("ASV1", "ASV2"))
  expect_equal(unname(x$otu["S1", "ASV1"]), 10)
})

test_that("aq_import transposes matrix input from taxa-rows format", {
  otu <- matrix(
    c(10, 5, 4, 8),
    nrow = 2,
    dimnames = list(c("ASV1", "ASV2"), c("S1", "S2"))
  )
  tax <- data.frame(
    Phylum = c("Firmicutes", "Proteobacteria"),
    row.names = c("ASV1", "ASV2")
  )
  meta <- data.frame(
    Group = c("A", "B"),
    row.names = c("S1", "S2")
  )

  x <- aq_import(otu, taxonomy = tax, metadata = meta)
  expect_equal(dim(x$otu), c(2L, 2L))
  expect_equal(rownames(x$otu), c("S1", "S2"))
  expect_equal(colnames(x$otu), c("ASV1", "ASV2"))
  expect_equal(unname(x$otu["S1", "ASV1"]), 10)
})


test_that("aq_import accepts data.frames with taxon IDs already in row names", {
  otu <- data.frame(
    S1 = c(10, 4),
    S2 = c(5, 8),
    row.names = c("ASV1", "ASV2"),
    check.names = FALSE
  )
  tax <- data.frame(Genus = c("G1", "G2"),
                    row.names = c("ASV1", "ASV2"))
  meta <- data.frame(Group = c("A", "B"),
                     row.names = c("S1", "S2"))

  x <- aq_import(otu, taxonomy = tax, metadata = meta)
  expect_equal(dim(x$otu), c(2L, 2L))
  expect_equal(rownames(x$otu), c("S1", "S2"))
  expect_equal(colnames(x$otu), c("ASV1", "ASV2"))
})

test_that("aq_import gives an informative error for ambiguous data.frames", {
  otu <- data.frame(S1 = c(10, 4), S2 = c(5, 8))
  expect_error(aq_import(otu), "row.names = 1")
})
