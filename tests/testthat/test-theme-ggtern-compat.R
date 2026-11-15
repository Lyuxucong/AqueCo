test_that("AqueCo themes keep ggtern tick lengths compatible", {
  th <- theme_aq("publication")
  expect_s3_class(th$tern.axis.ticks.length.major, "rel")
  expect_s3_class(th$tern.axis.ticks.length.minor, "rel")
})
