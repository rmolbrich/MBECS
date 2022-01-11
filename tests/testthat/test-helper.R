# TEST HELPER SECTION -----------------------------------------------------

test_that("mbecTestModel works", {
  # Works with list, phyloseq and MbecDdata input due to 'mbecProcessInput'
  # model is estimable and return value is NULL
  expect_identical(mbecTestModel(input.obj=datadummy,
                                 model.vars=c("group","batch")), NULL)
  # 'model.form' is class formula and is estimable
  expect_identical(mbecTestModel(input.obj=datadummy,
                                 model.form=stats::as.formula("y ~ group + batch")), NULL)
  # problem with estimability and return value is a character vector
  expect_vector(mbecTestModel(input.obj=datadummy,
                              model.vars=c("sample","group","batch")),
                character())
  # covariates and model-formula are missing
  expect_error(mbecTestModel(input.obj=datadummy),
               "Please supply covariates and/or model-formula.")

  # will construct generic model-formula if input is not class 'formula'
  form.res <- evaluate_promise(mbecTestModel(input.obj=datadummy,
                                             model.vars=c("group","batch"),
                                             model.form="y ~ group + batch"))
  expect_true(any(grepl("lm-formula", form.res$messages)))
})


test_that("mbecPCTest works", {
  # ToDo: meh
})


test_that("capitalization works", {
  expect_identical(mbecUpperCase("muh"), "Muh")
  expect_identical(mbecUpperCase("MUH"), "MUH")
})


test_that("linear modeling works", {
  lm.res <- evaluate_promise(mbecLM(datadummy, "lmm", c("group","batch"), type="otu"))
  expect_identical(length(lm.res$result), dim(datadummy$cnts)[2])
  expect_identical(class(lm.res$result), "numeric")
  expect_identical(typeof(lm.res$result), "double")
})


test_that("linear mixed modeling works", {
  lmm.res <- evaluate_promise(mbecLM(datadummy, "lmm", c("group","batch"), type="otu"))
  expect_identical(length(lmm.res$result), dim(datadummy$cnts)[2])
  expect_identical(class(lmm.res$result), "numeric")
  expect_identical(typeof(lmm.res$result), "double")
})



# TEST TRANSFORMATION SECTION ---------------------------------------------

test_that("Transformation wrapper works", {
  # creat small dataset for testing
  cnts <- matrix(1:16, nrow=4, ncol=4,
                 dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4")))
  meta <- data.frame("sID"=c("A","B","C","D"),
                     "group"=factor(c("X","X","Y","Y")),
                     "batch"=factor(c(1,2,1,2)), row.names = "sID")

  tss.res <- evaluate_promise(mbecTransform(input.obj=list(cnts,meta), method = "tss"))
  expect_identical(tss.res$result@tss, apply(cnts, 1, function(x){x/sum(x)}))
  expect_identical(colnames(tss.res$result@tss),c("A","B","C","D"))
  expect_identical(rownames(tss.res$result@tss),c("F1","F2","F3","F4"))

  clr.res <- evaluate_promise(mbecTransform(input.obj=list(cnts,meta), method = "clr"))
  expect_identical(clr.res$result@clr, t(mbecCLR(cnts)))
  expect_identical(colnames(clr.res$result@clr),c("A","B","C","D"))
  expect_identical(rownames(clr.res$result@clr),c("F1","F2","F3","F4"))

})


test_that("Percentile normalization works", {
  cnts <- matrix(1:16, nrow=4, ncol=4,
                 dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4")))
  cnts.norm <- matrix(rep.int(c(50,100), times=8), nrow=4, ncol=4,
                      dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4")))
  meta <- data.frame("sID"=c("A","B","C","D"),
                     "group"=factor(c("X","X","Y","Y")),
                     "batch"=factor(c(1,2,1,2)), row.names = "sID")
  # just use 'evaluate_promise()' to get all the relevant events for testing
  pn.res <- evaluate_promise(percentileNorm(cnts,meta=meta[,c("group","batch")]))
  expect_identical(dim(pn.res$result), dim(cnts))
  expect_true(grepl("Group.*",pn.res$messages))
  # test deterministic result!?
  expect_identical(cnts.norm ,pn.res$result)
})


test_that("percentile of score works", {
  expect_identical(poscore(c(1:50), 42, type="rank"), 84)
  expect_identical(poscore(c(1:50), 42, type="weak"), 84)
  expect_identical(poscore(c(1:50), 42, type="strict"), 82)
  expect_identical(poscore(c(1:50), 42, type="mean"), 83)
})


test_that("CLR transformation works", {
  cnts <- matrix(c(1,2,5,6), nrow=2, ncol=2,
                 dimnames=list(c("A","B"), c("F1","F2")))

  cnts.norm <- matrix(c(-0.80471896, -0.54930614, 0.80471896, 0.54930614), nrow=2, ncol=2,
                      dimnames=list(c("A","B"), c("F1","F2")))
  # test deterministic result!?
  expect_equal(mbecCLR(cnts) ,cnts.norm)
})

