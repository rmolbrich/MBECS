# TEST PRELIMINARY PLOTS --------------------------------------------------

## already part of 'test-analyses.R'


# TEST VARIANCE PLOTTATION ------------------------------------------------

test_that("mbecVarianceStatsPlot LM works", {
  # test for lm variances
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'lm' modelling
  mvar.plot.test <- evaluate_promise(mbecVarianceStatsPlot(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "lm",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL)))
  # expect ggplot object
  expect_s3_class(mvar.plot.test$result, "ggplot")
  # expect no warnings, messages or errors
  expect_warning(mvar.plot.test, NA)
  expect_message(mvar.plot.test, NA)
  expect_error(mvar.plot.test, NA)

})


test_that("mbecVarianceStatsPlot LMM works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test for lmm variances
  mvar.plot.test <- evaluate_promise(mbecVarianceStatsPlot(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                                             method = "lmm",
                                                                             model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                                             na.action = NULL)))
  # expect ggplot object
  expect_s3_class(mvar.plot.test$result, "ggplot")
  # expect no warnings, messages or errors
  expect_warning(mvar.plot.test, NA)
  expect_message(mvar.plot.test, NA)
  expect_error(mvar.plot.test, NA)
})


test_that("mbecRDAStatsPlot works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'rda' modelling
  rda.plot.test <- evaluate_promise(mbecRDAStatsPlot(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "rda",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL)))
  expect_s3_class(rda.plot.test$result, "ggplot")
  # expect no warnings, messages or errors
  expect_warning(rda.plot.test, NA)
  expect_message(rda.plot.test, NA)
  expect_error(rda.plot.test, NA)
})


test_that("mbecPVCAStatsPlot works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 'pvca' modelling
  pvca.plot.test <- evaluate_promise(mbecPVCAStatsPlot(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "pvca",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL)))
  expect_s3_class(pvca.plot.test$result, "ggplot")
  # expect no warnings, messages or errors
  expect_warning(pvca.plot.test, NA)
  expect_message(pvca.plot.test, NA)
  expect_error(pvca.plot.test, NA)

})


test_that("mbecSCOEFStatsPlot works", {
  tmp <- mbecGetData(input.obj = dummy.mbec, orientation = "sxf",
                     required.col = c("batch", "group"), type="clr", label=character())
  # test 's.coef' modelling
  scoef.plot.test <- evaluate_promise(mbecSCOEFStatsPlot(mbecModelVariance(dummy.mbec, model.vars = c("batch", "group"),
                                                  method = "s.coef",
                                                  model.form = NULL, type = "clr", label=character(), no.warning = TRUE,
                                                  na.action = NULL)))

  expect_s3_class(scoef.plot.test$result, "ggplot")
  # expect no warnings, messages or errors
  expect_warning(scoef.plot.test, NA)
  expect_message(scoef.plot.test, NA)
  expect_error(scoef.plot.test, NA)
})









