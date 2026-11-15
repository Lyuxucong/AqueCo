
test_that("null model preparation and assembly classifications work", {
  otu <- matrix(c(10,0,5,2, 0,8,1,3, 4,2,7,0, 1,6,0,9),
                nrow=4,byrow=TRUE,
                dimnames=list(paste0("S",1:4),paste0("T",1:4)))
  md <- data.frame(Group=c("KO","KO","WT","WT"),
                   row.names=rownames(otu))
  dat <- aq_data(otu=otu,metadata=md)
  x <- aq_null(dat,filter_thres=0,seed=123)
  expect_s3_class(x,"aq_null")
  expect_equal(nrow(x$otu),4L)
  m <- matrix(NA_real_,4,4,dimnames=list(rownames(otu),rownames(otu)))
  m[1,2] <- m[2,1] <- 2.5
  m[3,4] <- m[4,3] <- -1
  r <- m
  r[3,4] <- r[4,3] <- 0.99
  nti <- structure(list(matrix=m,metadata=md),class="aq_null_result")
  rc <- structure(list(matrix=r,metadata=md),class="aq_null_result")
  a <- aq_assembly(nti,rc)
  expect_equal(nrow(a$pairs),2L)
  expect_true("Heterogeneous Selection" %in% a$pairs$process)
  expect_true("Dispersal Limitation" %in% a$pairs$process)
  expect_equal(as.numeric(tapply(a$summary$percentage,
               a$summary$Group,sum)),c(100,100))
})
test_that("RCbray is symmetric and bounded for raw counts", {
  skip_if_not_installed("vegan")
  otu <- matrix(c(10,3,2, 4,8,1, 3,2,9, 7,4,5),4,3,byrow=TRUE,
                dimnames=list(paste0("S",1:4),paste0("T",1:3)))
  md <- data.frame(Group=c("KO","KO","WT","WT"),row.names=rownames(otu))
  x <- aq_null(aq_data(otu=otu,metadata=md),filter_thres=0)
  a <- aq_rcbray(x,runs=20,seed=7)
  expect_equal(a$matrix,t(a$matrix))
  expect_true(all(abs(a$matrix[upper.tri(a$matrix)])<=1))
})
