test_that("aq_data validates and stores samples", {
  x <- matrix(1:12, nrow=3, dimnames=list(paste0("S",1:3),paste0("T",1:4)))
  d <- aq_data(x)
  expect_s3_class(d,"aq_data")
  expect_equal(nrow(d$otu),3)
})

test_that("alpha and beta analyses return result classes", {
  skip_if_not_installed("vegan")
  set.seed(1)
  x <- matrix(rpois(60,10),nrow=6,dimnames=list(paste0("S",1:6),paste0("T",1:10)))
  md <- data.frame(Group=rep(c("A","B"),each=3),row.names=rownames(x))
  d <- aq_data(x,metadata=md)
  a <- aq_alpha(d,index=c("Observed","Shannon"))
  expect_s3_class(a,"aq_alpha")
  b <- aq_beta(d,ordination="PCoA",group="Group")
  expect_s3_class(b,"aq_beta")
})
