# TEST ANALYSIS FUNCTIONS -------------------------------------------------


test_that("mbecRLE works", {
  model.vars <- c("batch","group")
  # check out the returned data
  rle.test <- evaluate_promise(mbecRLE(dummy.mbec, model.vars = eval(model.vars), type="otu", return.data = TRUE))
  # expect class data.frame
  expect_s3_class(rle.test$result, "data.frame")
  # expect #n.sample rows and #n.sample + #covariate columns
  expect_equal(dim(rle.test$result), c(nrow(dummy.mbec@sam_data) * ncol(dummy.mbec@otu_table), (dim(dummy.mbec@sam_data)[2] + 3)))
  # expect #n.sample rows and 3 columns
  expect_identical(colnames(rle.test$result), c("specimen", "values", colnames(dummy.mbec@sam_data), "plot.order"))
  # expect no warnings
  expect_warning(rle.test, NA)
  # now put in plotting function and test that this produces a 'gtable'
  expect_s3_class(mbecRLEPlot(rle.df=rle.test$result,
                              model.vars=eval(model.vars)),
                  "ggplot")
})


test_that("mbecPCA works", {
  # check out the returned data
  pca.otu.test <- evaluate_promise(mbecPCA(dummy.mbec, model.vars = c("batch","group"), pca.axes = c(1,2), type="otu", return.data = TRUE))
  # expect #n.sample rows and #n.sample + #covariate columns
  expect_equal(dim(pca.otu.test$result[[1]]), c(nrow(dummy.mbec@sam_data), nrow(dummy.mbec@sam_data)+ncol(dummy.mbec@sam_data)))
  # expect #n.sample rows and 3 columns
  expect_equal(dim(pca.otu.test$result[[2]]), c(nrow(dummy.mbec@sam_data), 3))
  # expect colnames to be c("var.explained","axis.min","axis.max")
  expect_identical(colnames(pca.otu.test$result[[2]]), c("var.explained", "axis.min", "axis.max"))
  # expect that used axes are listed correctly
  expect_equal(pca.otu.test$result[[3]], c(1,2))
  # expect no warnings
  expect_warning(pca.otu.test, NA)
  # now put in plotting function and test that this produces a 'gtable'
  expect_s3_class(mbecPCAPlot(plot.df=pca.otu.test$result[[1]],
              metric.df=pca.otu.test$result[[2]],
              model.vars=c("batch","group"), pca.axes=pca.otu.test$result[[3]]),
              "gtable")

  # crepeat with different data, axes and only single covariate
  pca.clr.test <- evaluate_promise(mbecPCA(dummy.mbec, model.vars = c("batch"), pca.axes = c(3,4), type="clr", return.data = TRUE))
  # expect #n.sample rows and #n.sample + #covariate columns - but it will use all of them
  expect_equal(dim(pca.clr.test$result[[1]]), c(nrow(dummy.mbec@sam_data), nrow(dummy.mbec@sam_data)+ncol(dummy.mbec@sam_data)))
  # expect #n.sample rows and 3 columns
  expect_equal(dim(pca.clr.test$result[[2]]), c(nrow(dummy.mbec@sam_data), 3))
  # expect colnames to be c("var.explained","axis.min","axis.max")
  expect_identical(colnames(pca.clr.test$result[[2]]), c("var.explained", "axis.min", "axis.max"))
  # expect that used axes are listed correctly
  expect_equal(pca.clr.test$result[[3]], c(3,4))
  # expect no warnings
  expect_warning(pca.clr.test, NA)
  # now put in plotting function and test that this produces a 'gtable'
  expect_s3_class(mbecPCAPlot(plot.df=pca.clr.test$result[[1]],
                              metric.df=pca.clr.test$result[[2]],
                              model.vars=c("batch"), pca.axes=pca.clr.test$result[[3]]),
                  "gtable")
})


test_that("mbecBox works", {
  # check out the returned data
  box.test <- evaluate_promise(mbecBox(dummy.mbec, method = "TOP", n = 10, model.var = "batch", type="otu", label=character(),
                                           return.data = TRUE))
  # expect #n.sample rows and #n + 2
  expect_equal(dim(box.test$result[[1]]), c(nrow(dummy.mbec@sam_data), 12))
  # expect #n.sample rows and 3 columns
  expect_equal(length(box.test$result[[2]]), 10)
  # expect colnames
  expect_true(all(box.test$result[[2]] %in% colnames(box.test$result[[1]])))
  # expect no warnings
  expect_warning(box.test, NA)
  # now put in plotting function for further testing
  box.plot.test <- evaluate_promise(mbecBoxPlot(box.test$result[[1]], box.test$result[[2]], "batch"))
  # expect list of 10 grobs
  expect_equal(length(box.plot.test$result), 10)
  # expect no warnings
  expect_warning(box.plot.test, NA)

  # repeat for method 'ALL'
  # check out the returned data
  box.test <- evaluate_promise(mbecBox(dummy.mbec, method = "ALL", model.var = "batch", type="otu", label=character(),
                                       return.data = TRUE))
  # expect #n.sample rows and #n + 2
  expect_equal(dim(box.test$result[[1]]), c(nrow(dummy.mbec@sam_data), ncol(dummy.mbec@otu_table) + ncol(dummy.mbec@sam_data) + 1))
  # expect #n.sample rows and 3 columns
  expect_equal(length(box.test$result[[2]]), ncol(dummy.mbec@otu_table))
  # expect colnames
  expect_true(all(box.test$result[[2]] %in% colnames(box.test$result[[1]])))
  # expect no warnings
  expect_warning(box.test, NA)
  # now put in plotting function for further testing
  box.plot.test <- evaluate_promise(mbecBoxPlot(box.test$result[[1]], box.test$result[[2]], "batch"))
  # expect list of 500 grobs
  expect_equal(length(box.plot.test$result), ncol(dummy.mbec@otu_table))
  # expect no warnings
  expect_warning(box.plot.test, NA)

  # repeat for method selection of OTUs
  # check out the returned data
  otu.select <- c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
  box.test <- evaluate_promise(mbecBox(dummy.mbec, method = otu.select, model.var = "batch", type="otu", label=character(),
                                       return.data = TRUE))
  # expect #n.sample rows and #n + 2
  expect_equal(dim(box.test$result[[1]]), c(nrow(dummy.mbec@sam_data), 12))
  # expect #n.sample rows and 3 columns
  expect_equal(length(box.test$result[[2]]), 10)
  # expect colnames
  expect_true(all(box.test$result[[2]] %in% otu.select))
  # expect no warnings
  expect_warning(box.test, NA)
  # now put in plotting function for further testing
  box.plot.test <- evaluate_promise(mbecBoxPlot(box.test$result[[1]], box.test$result[[2]], "batch"))
  # expect list of 10 grobs
  expect_equal(length(box.plot.test$result), 10)
  # expect no warnings
  expect_warning(box.plot.test, NA)
})


test_that("mbecHeat works", {
  # check out the returned data
  heat.test <- evaluate_promise(mbecHeat(dummy.mbec, model.vars = c("batch", "group"),
                                         method = "TOP", n = 10, type="otu", label=character(),
                                       return.data = TRUE))
  # expect #n.sample rows and #n + 2
  expect_equal(dim(heat.test$result[[1]]), c(10, nrow(dummy.mbec@sam_data)))
  # expect #n.sample rows and 3 columns
  expect_equal(dim(heat.test$result[[2]]), dim(dummy.mbec@sam_data))
  # expect colnames
  expect_true(all(heat.test$result$sID %in% colnames(heat.test$result[[1]])))
  # expect no warnings
  expect_warning(heat.test, NA)
  # now put in plotting function for further testing
  heat.plot.test <- evaluate_promise(mbecHeatPlot(center=T, scale=T, heat.test$result[[1]], heat.test$result[[2]], c("batch", "group")))
  # expect list of 10 grobs
  expect_s3_class(heat.plot.test$result, "pheatmap")
  # expect no warnings
  expect_warning(heat.plot.test, NA)

  # check out the returned data for method 'ALL'
  heat.test <- evaluate_promise(mbecHeat(dummy.mbec, model.vars = c("batch", "group"),
                                         method = "ALL", type="otu", label=character(),
                                         return.data = TRUE))
  # expect #n.sample rows and #n + 2
  expect_equal(dim(heat.test$result[[1]]), c(ncol(dummy.mbec@otu_table), nrow(dummy.mbec@sam_data)))
  # expect #n.sample rows and 3 columns
  expect_equal(dim(heat.test$result[[2]]), dim(dummy.mbec@sam_data))
  # expect colnames
  expect_true(all(heat.test$result$sID %in% colnames(heat.test$result[[1]])))
  # expect no warnings
  expect_warning(heat.test, NA)
  # now put in plotting function for further testing
  heat.plot.test <- evaluate_promise(mbecHeatPlot(center=T, scale=T, heat.test$result[[1]], heat.test$result[[2]], c("batch", "group")))
  # expect list of 10 grobs
  expect_s3_class(heat.plot.test$result, "pheatmap")
  # expect no warnings
  expect_warning(heat.plot.test, NA)

  # check out the returned data for selected OTUs
  otu.select <- c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
  heat.test <- evaluate_promise(mbecHeat(dummy.mbec, model.vars = c("batch", "group"),
                                         method = otu.select, type="otu", label=character(),
                                         return.data = TRUE))
  # expect #n.sample rows and #n + 2
  expect_equal(dim(heat.test$result[[1]]), c(length(otu.select), nrow(dummy.mbec@sam_data)))
  # expect #n.sample rows and 3 columns
  expect_equal(dim(heat.test$result[[2]]), dim(dummy.mbec@sam_data))
  # expect colnames
  expect_true(all(heat.test$result$sID %in% colnames(heat.test$result[[1]])))
  # expect no warnings
  expect_warning(heat.test, NA)
  # now put in plotting function for further testing
  heat.plot.test <- evaluate_promise(mbecHeatPlot(center=T, scale=T, heat.test$result[[1]], heat.test$result[[2]], c("batch", "group")))
  # expect list of 10 grobs
  expect_s3_class(heat.plot.test$result, "pheatmap")
  # expect no warnings
  expect_warning(heat.plot.test, NA)
})


test_that("mbecMosaic works", {
  mosaic.test <- evaluate_promise(mbecMosaic(dummy.mbec,
                                             model.vars = c("batch", "group"),
                                             return.data = TRUE))
  expect_true(all(colnames(mosaic.test$result) %in% c("Var1", "Var2", "Freq", "Freq.scaled")))
  # expect no warnings
  expect_warning(mosaic.test, NA)
  # now put in plotting function for further testing
  mosaic.plot.test <- evaluate_promise(mbecMosaicPlot(mosaic.test$result, c("batch", "group")))
  # expect list of 10 grobs
  expect_s3_class(mosaic.plot.test$result, "gtable")
  # expect no warnings
  expect_warning(mosaic.plot.test, NA)
})



# TEST VARIANCE CALCULATIONS ----------------------------------------------

test_that("mbecModelVariance LM works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'lm' modelling
  mvar.test <- evaluate_promise(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "lm",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL))
  # expect dimension
  expect_equal(dim(mvar.test$result), c(ncol(dummy.mbec@otu_table), 4))
  # expect colnames
  expect_true(all(colnames(mvar.test$result) %in% c("batch","group","Residuals","type")))
  # expect no warnings
  expect_warning(mvar.test, NA)

  # build directly from 'mbecModelVarianceLM(model.form, model.vars, tmp.cnts, tmp.meta, type)'
  lm.test <- evaluate_promise(mbecModelVarianceLM(model.form=NULL, model.vars = c("batch", "group"), tmp.cnts=tmp[[1]], tmp.meta=tmp[[2]], type="clr"))

  expect_identical(mvar.test$result, lm.test$result)
  # expect no warnings
  expect_warning(lm.test, NA)
})


test_that("mbecModelVariance LMM works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'lmm' modelling
  mvar.test <- evaluate_promise(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "lmm",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL))
  # expect dimension
  expect_equal(dim(mvar.test$result), c(ncol(dummy.mbec@otu_table), 4))
  # expect colnames
  expect_true(all(colnames(mvar.test$result) %in% c("batch","groupB","Residuals","type")))
  # expect no warnings
  expect_warning(mvar.test, NA)

  # build directly from 'mbecModelVarianceLM(model.form, model.vars, tmp.cnts, tmp.meta, type)'
  lmm.test <- evaluate_promise(mbecModelVarianceLMM(model.form=NULL, model.vars = c("batch", "group"), tmp.cnts=tmp[[1]], tmp.meta=tmp[[2]], type="clr"))

  expect_identical(mvar.test$result, lmm.test$result)
  # expect no warnings
  expect_warning(lmm.test, NA)

})


test_that("mbecModelVariance RDA works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'rda' modelling
  mvar.test <- evaluate_promise(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "rda",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL))
  # expect dimension
  expect_equal(dim(mvar.test$result), c(1, 3))
  # expect colnames
  expect_true(all(colnames(mvar.test$result) %in% c("batch","group","type")))
  # expect no warnings
  expect_warning(mvar.test, NA)

  # build directly from 'mbecModelVarianceLM(model.form, model.vars, tmp.cnts, tmp.meta, type)'
  rda.test <- evaluate_promise(mbecModelVarianceRDA(model.vars = c("batch", "group"), tmp.cnts=tmp[[1]], tmp.meta=tmp[[2]], type="clr"))

  expect_identical(mvar.test$result, rda.test$result)
  # expect no warnings
  expect_warning(rda.test, NA)
})


test_that("mbecModelVariance PVCA works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'pvca' modelling
  mvar.test <- evaluate_promise(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "pvca",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL))
  # expect dimension
  expect_equal(dim(mvar.test$result), c(1, 5))
  # expect colnames
  expect_true(all(colnames(mvar.test$result) %in% c("batch.group","group","batch","Residual","type")))
  # expect no warnings
  expect_warning(mvar.test, NA)

  # build directly from 'mbecModelVarianceLM(model.form, model.vars, tmp.cnts, tmp.meta, type)'
  pvca.test <- evaluate_promise(mbecModelVariancePVCA(model.vars = c("batch", "group"), tmp.cnts=tmp[[1]], tmp.meta=tmp[[2]], type="clr",pct_threshold = 0.5876))

  expect_identical(mvar.test$result, pvca.test$result)
  # expect no warnings
  expect_warning(pvca.test, NA)
})


test_that("mbecModelVariance S.COEF works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 's.coef' modelling
  mvar.test <- evaluate_promise(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "s.coef",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL))
  # expect dimension
  expect_equal(dim(mvar.test$result), c(4, 4))
  # expect colnames
  expect_true(all(colnames(mvar.test$result) %in% c("variable","cluster","sil.coefficient","type")))
  # expect no warnings
  expect_warning(mvar.test, NA)

  # build directly from 'mbecModelVarianceSCOEF(model.form, model.vars, tmp.cnts, tmp.meta, type)'
  scoef.test <- evaluate_promise(mbecModelVarianceSCOEF(model.vars = c("batch", "group"), tmp.cnts=tmp[[1]], tmp.meta=tmp[[2]], type="clr"))

  expect_identical(mvar.test$result, scoef.test$result)
  # expect no warnings
  expect_warning(scoef.test, NA)
})



# TEST VARIANCE STATISTICS ------------------------------------------------

test_that("mbecVarianceStats LM works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'lm' modelling
  mvar.test <- evaluate_promise(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "lm",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL))

  tmp.formula = stats::as.formula(paste("y", " ~ ", paste(c("batch", "group"), collapse = " + ")))
  model.variances <- NULL
  y <- tmp[[1]][[1]]
  model.variances <- mbecVarianceStats(stats::lm(tmp.formula, data = tmp[[2]]))
  # expect same result
  expect_equal(mvar.test$result[1,1:3], data.frame(model.variances))

})


test_that("mbecVarianceStats LMM works", {
  model.vars <- c("batch", "group")
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = model.vars, type="clr", label=character())
  # test 'lm' modelling
  mvar.test <- evaluate_promise(mbecModelVariance(dummy.mbec, model.vars = model.vars,
                                                  method = "lmm",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL))

  control = lme4::lmerControl(calc.derivs = TRUE, check.rankX = "stop.deficient")
  f.terms <- paste("(1|", model.vars, ")", sep = "")
  tmp.formula <- stats::as.formula(paste(paste("y", paste(model.vars[-1], collapse = " + "), sep = " ~ "),
                                           paste(f.terms[1], collapse = " + "), sep = " + "))

    y <- tmp[[1]][[1]]
    model.variances <- evaluate_promise(mbecVarianceStats(lme4::lmer(tmp.formula, data = tmp[[2]], control = control)))
    expect_equal(mvar.test$result[1,1:3], data.frame(model.variances$result))

})


test_that("mbecMixedvariance works", {
  model.vars <- c("batch", "group")
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = model.vars, type="clr", label=character())

  control = lme4::lmerControl(calc.derivs = TRUE, check.rankX = "stop.deficient")
  f.terms <- paste("(1|", model.vars, ")", sep = "")
  tmp.formula <- stats::as.formula(paste(paste("y", paste(model.vars[-1], collapse = " + "), sep = " ~ "),
                                         paste(f.terms[1], collapse = " + "), sep = " + "))

  y <- tmp[[1]][[1]]

  mvar.test <- evaluate_promise(mbecMixedVariance(lme4::lmer(tmp.formula, data = tmp[[2]], control = control)))

  expect_equal(length(mvar.test$result), 3)
  # expect no warnings
  expect_warning(mvar.test, NA)
})


test_that("mbecValidateModel works", {
  model.vars <- c("batch", "group")
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = model.vars, type="clr", label=character())

  # expect warning
  mval.test <- evaluate_promise(mbecValidateModel(lme4::lmer(tmp[[1]][,1] ~ group + (1|batch),data=tmp[[2]]), colinearityThreshold = 0.5))
  expect_true(grepl("Some covariates are strongly *", mval.test$result))

  # expect no problems
  mval.test <- evaluate_promise(mbecValidateModel(lme4::lmer(tmp[[1]][,1] ~ group + (1|batch),data=tmp[[2]]), colinearityThreshold = 0.999))
  expect_identical(mval.test$result, NULL)

})


test_that("colinScore works", {

  cnts <- matrix(1:16, nrow=4, ncol=4,
                 dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4")))
  meta <- data.frame("sID"=c("A","B","C","D"),
                     "group"=factor(c("X","X","Y","Y")),
                     "batch"=factor(c(1,2,1,2)), row.names = "sID")
  # expect warning
  mval.test <- evaluate_promise(colinScore(lme4::lmer(cnts[,1] ~ group + (1|batch),data=meta)))
  expect_equal(mval.test$result[1], 6.3220273e-08)
  expect_identical(names(attributes(mval.test$result)), "vcor")
  expect_warning(mval.test, NA)
})





# VARIANCE CALCULATION ----------------------------------------------------

test_that("SOMETHING works", {
  # Works with list, phyloseq and MbecDdata input due to 'mbecProcessInput'
  # model is estimable and return value is NULL
  expect_identical(mbecTestModel(input.obj=dummy.mbec,
                                 model.vars=c("group","batch")), NULL)
  # 'model.form' is class formula and is estimable
  expect_identical(mbecTestModel(input.obj=dummy.mbec,
                                 model.form=stats::as.formula("y ~ group + batch")), NULL)
  # problem with estimability and return value is a character vector
  expect_vector(mbecTestModel(input.obj=dummy.mbec,
                              model.vars=c("sID","group","batch")),
                character())
  # covariates and model-formula are missing
  expect_error(mbecTestModel(input.obj=dummy.mbec),
               "Please supply covariates and/or model-formula.")

  # will construct generic model-formula if input is not class 'formula'
  form.res <- evaluate_promise(mbecTestModel(input.obj=dummy.mbec,
                                             model.vars=c("group","batch"),
                                             model.form="y ~ group + batch"))
  expect_true(any(grepl("lm-formula", form.res$messages)))
})

