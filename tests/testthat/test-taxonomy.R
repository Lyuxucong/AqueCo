test_that("aq_taxa recognizes canonical taxonomic rank names", {
  otu <- matrix(
    c(10, 5, 3, 8, 4, 6),
    nrow = 2,
    dimnames = list(c("S1", "S2"), paste0("ASV", 1:3))
  )

  taxonomy <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Bacteria"),
    Phylum = c("Firmicutes", "Proteobacteria", "Firmicutes"),
    Genus = c("GenusA", "GenusB", "GenusC"),
    row.names = colnames(otu)
  )

  dat <- aq_data(otu = otu, taxonomy = taxonomy)

  expect_s3_class(aq_taxa(dat, level = "Phylum"), "aq_taxa")
  expect_s3_class(aq_taxa(dat, level = "phylum"), "aq_taxa")
  expect_s3_class(aq_taxa(dat, level = "PHYLUM"), "aq_taxa")
  expect_s3_class(aq_taxa(dat, level = "Genus"), "aq_taxa")
})

test_that("aq_taxa recognizes common rank abbreviations", {
  otu <- matrix(
    c(10, 5, 3, 8, 4, 6),
    nrow = 2,
    dimnames = list(c("S1", "S2"), paste0("ASV", 1:3))
  )

  taxonomy <- data.frame(
    p = c("Firmicutes", "Proteobacteria", "Firmicutes"),
    g = c("GenusA", "GenusB", "GenusC"),
    row.names = colnames(otu)
  )

  dat <- aq_data(otu = otu, taxonomy = taxonomy)

  expect_s3_class(aq_taxa(dat, level = "p"), "aq_taxa")
  expect_s3_class(aq_taxa(dat, level = "g"), "aq_taxa")
})

test_that("aq_taxa aggregation preserves sample rows and taxon columns", {
  otu <- matrix(
    c(10, 5, 3, 8, 4, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("S1", "S2"), paste0("ASV", 1:3))
  )

  taxonomy <- data.frame(
    Phylum = c("Firmicutes", "Proteobacteria", "Firmicutes"),
    Genus = c("GenusA", "GenusB", "GenusC"),
    row.names = colnames(otu)
  )

  dat <- aq_data(otu = otu, taxonomy = taxonomy)
  res <- aq_taxa(dat, level = "Phylum", abundance = "count")

  expect_equal(rownames(res$count), rownames(otu))
  expect_equal(ncol(res$count), 2L)
  expect_equal(unname(res$count["S1", "Firmicutes"]), 13)
  expect_equal(unname(res$count["S1", "Proteobacteria"]), 5)
  expect_equal(rowSums(res$relative), c(S1 = 1, S2 = 1))
})
