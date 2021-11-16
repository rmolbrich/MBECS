test_that("capitalization works", {
  expect_identical(mbecUpperCase("muh"), "Muh")
  expect_identical(mbecUpperCase("MUH"), "MUH")
})

test_that("linear modeling works", {
  lm.res <- evaluate_promise(mbecLM(datadummy, "lmm", c("group","batch")))
  expect_identical(length(lm.res$result), dim(datadummy$cnts)[2])
  expect_identical(class(lm.res$result), "numeric")
  expect_identical(typeof(lm.res$result), "double")
})

test_that("linear mixed modeling works", {
  lmm.res <- evaluate_promise(mbecLM(datadummy, "lmm", c("group","batch")))
  expect_identical(length(lmm.res$result), dim(datadummy$cnts)[2])
  expect_identical(class(lmm.res$result), "numeric")
  expect_identical(typeof(lmm.res$result), "double")
})


test_that("Percentile normalization works", {
  # just use 'evaluate_promise()' to get all the relevant events for testing
  pn.res <- evaluate_promise(percentileNorm(cnts=datadummy$cnts,meta=datadummy$meta[,c("group","batch")]))
  expect_identical(dim(pn.res$result), dim(datadummy$cnts))
  expect_identical(pn.res$messages[1], "Group 0-0.5 is considered control group, i.e., reference for normalization procedure. To change reference please 'relevel()' grouping factor accordingly.\n")
  # test deterministic result!?
  expect_identical(percentileNorm(cnts=datadummy$cnts,meta=datadummy$meta[,c("group","batch")]) ,pn.res$result)
})

test_that("percentile of score works", {
  expect_identical(poscore(c(1:50), 42, type="rank"), 84)
  expect_identical(poscore(c(1:50), 42, type="weak"), 84)
  expect_identical(poscore(c(1:50), 42, type="strict"), 82)
  expect_identical(poscore(c(1:50), 42, type="mean"), 83)
})






