test_that("class methods work", {
  constructor.res <- evaluate_promise(MbecData(input.obj=datadummy$cnts, meta.obj = datadummy$meta))
  # correct class?
  expect_identical(class(constructor.res$result)[1], "MbecData")
  # correct package - which seems redundant but then again.. why not
  expect_identical(attr(class(constructor.res$result), "package"), "MBECS")
  # all the slots present?
  expect_identical(names(attributes(constructor.res$result)), c("type","log","transformations","otu_table","tax_table","sam_data","phy_tree","refseq","class"))
  # worked without a hitch?
  expect_identical(constructor.res$warnings, character(0))
  expect_identical(constructor.res$messages, character(0))

  ## ToDo: check if function fails correctly.


})

