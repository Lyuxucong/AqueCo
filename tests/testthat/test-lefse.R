
test_that("LEfSe wrapper checks data and optional microeco dependency", {
  otu <- matrix(c(10,20,30,40,12,22,32,42),nrow=4,
    dimnames=list(paste0("S",1:4),c("A","B")))
  dat <- aq_data(otu)
  expect_error(aq_lefse(dat,group="Group"),"microeco|taxonomy")
})
test_that("LEfSe wrapper exports data and validates plots", {
  x <- structure(list(res_diff=data.frame(Taxon="A",LDA=3),
    group="Group"),class="aq_lefse")
  f <- tempfile(fileext=".csv")
  aq_lefse_export(x,f)
  expect_true(file.exists(f))
  expect_error(aq_lefse_plot(x,type="other"))
})

test_that("microeco receives a data.frame OTU table", {
  skip_if_not_installed("microeco")
  otu <- matrix(c(12, 5, 9, 20, 18, 6, 11, 15),
                nrow=4, dimnames=list(paste0("S",1:4),c("ASV1","ASV2")))
  otu_df <- as.data.frame(t(otu), check.names=FALSE)
  expect_s3_class(otu_df,"data.frame")
  expect_equal(rownames(otu_df),colnames(otu))
  expect_equal(colnames(otu_df),rownames(otu))
  expect_equal(as.matrix(otu_df),unname(t(otu)),ignore_attr=TRUE)
})
