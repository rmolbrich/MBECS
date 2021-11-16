test_that("capitalization works", {
  expect_identical(mbecUpperCase("muh"), "Muh")
  expect_identical(mbecUpperCase("MUH"), "MUH")
})

test_that("linear modeling works", {
  expect_identical(length(mbecLM(datadummy, "lm", c("group","batch"))), dim(datadummy$cnts)[2])
  expect_identical(class(mbecLM(datadummy, "lm", c("group","batch"))), "numeric")
  expect_identical(typeof(mbecLM(datadummy, "lm", c("group","batch"))), "double")
})

test_that("linear mixed modeling works", {
  expect_identical(length(mbecLM(datadummy, "lmm", c("group","batch"))), dim(datadummy$cnts)[2])
  expect_identical(class(mbecLM(datadummy, "lmm", c("group","batch"))), "numeric")
  expect_identical(typeof(mbecLM(datadummy, "lmm", c("group","batch"))), "double")
})


test_that("Percentile normalization works", {
  expect_identical(dim(percentileNorm(cnts=datadummy$cnts,meta=datadummy$meta[,c("group","batch")])), dim(datadummy$cnts))



})




#mtx.pn_counts <- percentileNorm(cnts=datadummy$cnts,meta=datadummy$meta[,c("group","batch")])











