
test_that("distance-decay supports overall, group and selected groups", {
  set.seed(42)
  ids <- paste0("S",seq_len(12))
  otu <- matrix(sample(1:100,12*15,TRUE),12,15,
                dimnames=list(ids,paste0("T",seq_len(15))))
  md <- data.frame(Group=rep(c("IW","CW","TW"),each=4),
                   Longitude=seq(110,111,length.out=12),
                   Latitude=seq(30,31,length.out=12),row.names=ids)
  env <- data.frame(pH=seq(6,8,length.out=12),
                    TN=seq(1,3,length.out=12)+rnorm(12,0,.2),
                    row.names=ids)
  dat <- aq_data(otu,metadata=md,env=env)
  res <- aq_ddr(dat,group="Group",scope="both",permutations=9)
  expect_s3_class(res,"aq_ddr")
  expect_equal(sort(unique(res$summary$Group)),
               sort(c("Overall","IW","CW","TW")))
  expect_equal(nrow(res$summary),8L)
  expect_true(all(is.finite(res$summary$R2)))
  expect_true(all(res$summary$Mantel_p>=0 & res$summary$Mantel_p<=1))
  expect_s3_class(aq_ddr_plot(res),"ggplot")
  expect_s3_class(aq_ddr_plot(res,groups="IW",distance="geographic"),"ggplot")
  expect_s3_class(aq_ddr_plot(res,compare=FALSE),"ggplot")
  single <- aq_ddr(dat,group="Group",groups="IW",scope="group",
                   distance="environmental",permutations=9)
  expect_equal(unique(single$summary$Group),"IW")
  expect_equal(single$summary$Pairs,6)
  f <- tempfile()
  files <- aq_ddr_export(single,f)
  expect_true(all(file.exists(files)))
})

test_that("geographic distance uses metadata only and preserves matrix dimensions", {
  set.seed(10)
  ids <- paste0("S",1:6)
  otu <- matrix(sample(1:20,36,TRUE),6,6,
                dimnames=list(ids,paste0("ASV",1:6)))
  md <- data.frame(Group=rep(c("IW","CW"),each=3),
                   Longitude=c(110,110.1,110.2,111,111.1,111.2),
                   Latitude=c(30,30.1,30.2,31,31.1,31.2),
                   row.names=ids)
  env <- data.frame(pH=seq(6,7,length.out=6),
                    TN=c(1,3,2,4,5,6),row.names=ids)
  dat <- aq_data(otu,metadata=md,env=env)
  result <- aq_ddr(dat,distance=c("geographic","environmental"),
                   group="Group",scope="both",permutations=9)
  expect_equal(nrow(result$summary),6L)
  expect_equal(nrow(result$pairs[result$pairs$Group=="Overall" &
                   result$pairs$Distance=="Geographic",]),15L)
  expect_true(all(result$pairs$Predictor>=0))
  expect_s3_class(aq_ddr_plot(result),"ggplot")
})
