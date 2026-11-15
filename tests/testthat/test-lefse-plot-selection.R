
test_that("LDA plotting applies total and per-group caps", {
  x <- structure(list(res_diff=data.frame(
    Taxa=paste0("g__Taxon",1:90),
    Group=rep(c("IW","CW","TW"),each=30),
    LDA=seq(4.5,2.1,length.out=90))),class="aq_lefse")
  p <- aq_lefse_plot(x,top_n_per_group=5)
  expect_s3_class(p,"ggplot")
  expect_equal(nrow(p$data),15L)
  expect_equal(as.integer(table(p$data$Group)),c(5L,5L,5L))
  expect_equal(nrow(aq_lefse_plot(x,top_n=12)$data),12L)
  expect_equal(nrow(aq_lefse_plot(x,use_number=9)$data),9L)
  expect_equal(nrow(aq_lefse_plot(x,tax_level="Genus",top_n=8)$data),8L)
  expect_error(aq_lefse_plot(x,tax_level="Family"),"No LEfSe features")
})
