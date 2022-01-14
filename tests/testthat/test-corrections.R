# TEST ASSESSMENT FUNCTIONS -----------------------------------------------

test_that("mbecSVA works", {
  model.vars <- c("batch","group")
  # test 'sva' modelling
  sva.test <- evaluate_promise(mbecSVA(input.obj=dummy.mbec, model.vars=model.vars))
  # expect message about surrogate variables
  expect_true(grepl("Number of significant surrogate", sva.test$output, fixed=TRUE))
  # expect a p-value for each feature
  expect_equal(length(sva.test$result), ncol(dummy.mbec@otu_table))
  # expect values <= 1
  expect_true(all(sva.test$result <= 1))

})



# TEST RUV FUNCTIONS ------------------------------------------------------

test_that("mbecRUV2 works", {
  model.vars <- c("batch","group")
  ruv2.test <- evaluate_promise(mbecRUV2(input.obj=dummy.mbec, model.vars, type="otu", nc.features=NULL))
  # expect no warning
  expect_warning(ruv2.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying Remove Unwanted Variantion v2 (RUV-II) for batch-correction.\n", ruv2.test$messages, fixed=TRUE)))
  expect_true(any(grepl("No negative control features provided. Using pseudo-negative controls.\n", ruv2.test$messages, fixed=TRUE)))

  # And with pre-selected features
  otu.select <- c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
  ruv2.test <- evaluate_promise(mbecRUV2(input.obj=dummy.mbec, model.vars, type="otu", nc.features=otu.select))
  # expect no warning
  expect_warning(ruv2.test, NA)
  # expect feedback message
  expect_false(any(grepl("No negative control features provided. Using pseudo-negative controls.\n", ruv2.test$messages, fixed=TRUE)))
})


test_that("mbecRUV3 works", {

  model.vars <- c("batch","group")
  ruv3.test <- evaluate_promise(mbecRUV3(input.obj=dummy.mbec, model.vars, type="clr", nc.features=NULL))
  # expect no warning
  expect_warning(ruv3.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying Remove Unwanted Variantion v3 (RUV-III) for batch-correction.\n", ruv3.test$messages, fixed=TRUE)))
  expect_true(any(grepl("No negative control features provided. Using pseudo-negative controls.\n", ruv3.test$messages, fixed=TRUE)))

  # And with pre-selected features
  otu.select <- c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
  ruv3.test <- evaluate_promise(mbecRUV2(input.obj=dummy.mbec, model.vars, type="otu", nc.features=otu.select))
  # expect no warning
  expect_warning(ruv3.test, NA)
  # expect feedback message
  expect_false(any(grepl("No negative control features provided. Using pseudo-negative controls.\n", ruv3.test$messages, fixed=TRUE)))

  # expect errors without replicate column - after I removed the column..
  dummy.noRep <- dummy.mbec
  dummy.noRep@sam_data <- dummy.mbec@sam_data[,-3]
  expect_warning(mbecRUV3(input.obj=dummy.noRep, model.vars, type="clr", nc.features=NULL), "No technical replicates found. RUV-3 is not available for this data-set!")
})


test_that("mbecRUV4 works", {
  model.vars <- c("batch","group")
  ruv4.test <- evaluate_promise(mbecRUV4(input.obj=dummy.mbec, model.vars, type="otu", nc.features=NULL))
  # expect no warning
  expect_warning(ruv4.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying Remove Unwanted Variantion v4 (RUV-IV) for batch-correction.\n" , ruv4.test$messages, fixed=TRUE)))
  expect_true(any(grepl("No negative control features provided. Using pseudo-negative controls.\n", ruv4.test$messages, fixed=TRUE)))
  # expect p-values for every feature
  expect_equal(length(ruv4.test$result), ncol(dummy.mbec@otu_table))

  # And with pre-selected features
  otu.select <- c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
  ruv4.test <- evaluate_promise(mbecRUV4(input.obj=dummy.mbec, model.vars, type="otu", nc.features=otu.select))
  # expect no warning
  expect_warning(ruv4.test, NA)
  # expect feedback message
  expect_false(any(grepl("No negative control features provided. Using pseudo-negative controls.\n", ruv4.test$messages, fixed=TRUE)))
  # expect p-values for every feature
  expect_equal(length(ruv4.test$result), ncol(dummy.mbec@otu_table))
})


# TEST CORRECTION FUNCTIONS -----------------------------------------------

test_that("mbecBMC works", {
  cnts <- matrix(1:24, nrow=4, ncol=6,
                 dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4","F5","F6")))
  meta <- data.frame("sID"=c("A","B","C","D"),
                     "group"=factor(c("X","X","Y","Y")),
                     "batch"=factor(c(1,2,1,2)), row.names = "sID")


  test.frame <- mbecTransform(list(cnts,meta), method="clr", required.col = c("batch","group"))
  bmc.test <- evaluate_promise(mbecBMC(input.obj=test.frame, "batch", type="otu"))
  # expect no warning
  expect_warning(bmc.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying Batch Mean-Centering", bmc.test$messages, fixed=TRUE)))
  # for mockup values expect identical result
  expect_identical(bmc.test$result, as.data.frame(matrix(rep(c(-1,-1,1,1), times=6), nrow=4, ncol=6,
                                           dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4","F5","F6")))))
})


test_that("mbecBat works", {
  model.vars <- c("batch","group")

  bat.test <- evaluate_promise(mbecBat(input.obj=dummy.mbec, model.vars, type="otu"))
  # expect no warning
  expect_warning(bat.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying ComBat (sva) for batch-correction.\n", bat.test$messages, fixed=TRUE)))

  # test for single covariate
  bat.single.test <- evaluate_promise(mbecBat(input.obj=dummy.mbec, model.vars[1], type="clr"))
  # expect no warning
  expect_warning(bat.single.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying ComBat (sva) for batch-correction.\n", bat.single.test$messages, fixed=TRUE)))
  # expect that values differ
  expect_true(!all(bat.test$result == bat.single.test$result))

})



test_that("mbecRBE works", {
  model.vars <- c("batch","group")

  rbe.test <- evaluate_promise(mbecRBE(input.obj=dummy.mbec, model.vars, type="otu"))
  # expect no warning
  expect_warning(rbe.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying 'removeBatchEffect' (limma) for batch-correction.\n", rbe.test$messages, fixed=TRUE)))

  # test for single covariate
  rbe.clr.test <- evaluate_promise(mbecRBE(input.obj=dummy.mbec, model.vars[1], type="clr"))
  # expect no warning
  expect_warning(rbe.clr.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying 'removeBatchEffect' (limma) for batch-correction.\n", rbe.clr.test$messages, fixed=TRUE)))
  # expect that values differ
  expect_true(!all(rbe.test$result == rbe.clr.test$result))

})


test_that("mbecPN works", {
  model.vars <- c("batch","group")
  cnts <- matrix(1:24, nrow=4, ncol=6,
                 dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4","F5","F6")))
  meta <- data.frame("sID"=c("A","B","C","D"),
                     "group"=factor(c("X","X","Y","Y")),
                     "batch"=factor(c(1,2,1,2)), row.names = "sID")
  test.frame <- mbecTransform(list(cnts,meta), method="tss", required.col = c("batch","group"))

  pn.test <- evaluate_promise(mbecRBE(input.obj=test.frame, model.vars, type="tss"))
  # expect no warning
  expect_warning(pn.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying 'removeBatchEffect' (limma) for batch-correction.\n", pn.test$messages, fixed=TRUE)))
  # expect result for mockup data
  expect_equal(pn.test$result, matrix(c(0.02059746,0.07902514,0.13745282,0.19588051,0.25430819,0.31273588,
                                                           0.02233183,0.08006577,0.13779970,0.19553363,0.25326757,0.31100150,
                                                           0.04390748,0.09301116,0.14211483,0.19121850,0.24032218,0.28942585,
                                                           0.04217310,0.09197053,0.14176795,0.19156538,0.24136280,0.29116023), nrow=6, ncol=4,
                                                         dimnames=list(c("F1","F2","F3","F4","F5","F6"),c("A","B","C","D"))), tolerance=0.0000001)

})


test_that("mbecSVD works", {

  model.vars <- c("batch","group")
  cnts <- matrix(1:24, nrow=4, ncol=6,
                 dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4","F5","F6")))
  meta <- data.frame("sID"=c("A","B","C","D"),
                     "group"=factor(c("X","X","Y","Y")),
                     "batch"=factor(c(1,2,1,2)), row.names = "sID")
  test.frame <- mbecTransform(list(cnts,meta), method="tss", required.col = c("batch","group"))

  svd.test <- evaluate_promise(mbecSVD(input.obj=test.frame, model.vars, type="otu"))
  # expect no warning
  expect_warning(svd.test, NA)
  # expect feedback message
  expect_true(any(grepl("Applying Singular Value Decomposition (SVD) for batch-correction.\n", svd.test$messages, fixed=TRUE)))
  # expect result for mockup data
  expect_equal(svd.test$result, as.data.frame(matrix(c(rep(2.5,4), rep(6.5, 4), rep(10.5,4), rep(14.5,4),rep(18.5,4),rep(22.5,4)), nrow=4, ncol=6,
                                      dimnames=list(c("A","B","C","D"),c("F1","F2","F3","F4","F5","F6")))))


})


# TEST CORRECTION WRAPPER -------------------------------------------------


test_that("mbecCorrection works", {

  # ToDo:

})






test_that("mbecRunCorrections wrapper works", {

  # ToDo:

})









