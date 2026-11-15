
test_that("fast RDA skips selection and optional permutation tests", {
  skip_if_not_installed("vegan")
  set.seed(7)
  otu <- matrix(sample(1:20,60,TRUE),10,6,
    dimnames=list(paste0("S",1:10),paste0("T",1:6)))
  env <- data.frame(pH=seq_len(10),row.names=rownames(otu))
  dat <- aq_data(otu=otu,env=env)
  fit <- aq_rda_cca(dat,method="RDA",select="none",mode="fast",
                    tests="global",permutations=9)
  expect_s3_class(fit,"aq_ordination")
  expect_null(fit$term_test)
  expect_null(fit$axis_test)
})
test_that("phylogenetic null model accepts cached distance and workers", {
  skip_if_not_installed("ape")
  skip_if_not_installed("picante")
  set.seed(2)
  otu <- matrix(sample(1:15,30,TRUE),6,5,
    dimnames=list(paste0("S",1:6),paste0("T",1:5)))
  tree <- ape::rtree(5,tip.label=colnames(otu))
  dat <- aq_data(otu=otu,tree=tree)
  x <- aq_null(dat,group=NULL,filter_thres=0,cache=TRUE)
  expect_true(is.matrix(x$phylo_dist))
  res <- aq_beta_nti(x,runs=4,workers=1,seed=2)
  expect_equal(dim(res$matrix),c(6L,6L))
})
