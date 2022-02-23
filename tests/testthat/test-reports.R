# REPORTING WRAPPER -------------------------------------------------------

test_that("mbecReportPrelim works", {
    data(dummy.list)

    dummy.test <- mbecTransform( list(dummy.list$cnts[,seq(20)],
                                      dummy.list$meta), method="tss")
    dummy.test <- mbecTransform(dummy.test, method="clr")

    # have the data returned to do some testing
    PreRep.test <- evaluate_promise(
        mbecReportPrelim(input.obj=dummy.test, model.vars=c("batch","group"),
                         type="clr", file.name=NULL, file.dir=NULL,
                         return.data = TRUE))
    # all analyses present
    expect_true(
        all(names(PreRep.test$result) %in%
                c("mosaic","pca","rle","heat","box",
                  "linmod","linmixmod","rda","pvca","scoef")))

    # mosaic works
    ds.test <- mbecMosaic(
        input.obj=dummy.test, model.vars=c("batch","group"), return.data=TRUE)

    expect_equal(PreRep.test$result$mosaic, ds.test)
    # pca works
    ds.test <- mbecPCA(
        input.obj=dummy.test, model.vars=c("batch","group"), return.data=TRUE)

    expect_equal(PreRep.test$result$pca, ds.test)
    # rle works
    ds.test <- mbecRLE(
        input.obj=dummy.test, model.vars=c("batch","group"), return.data=TRUE)

    expect_equal(PreRep.test$result$rle, ds.test)
    # heatmap works
    ds.test <- mbecHeat(
        input.obj=dummy.test, model.vars=c("batch","group"), return.data=TRUE)

    expect_equal(PreRep.test$result$heat, ds.test)
    # box works
    ds.test <- mbecBox(
        input.obj=dummy.test, method="TOP", n=5, model.var="batch",
        return.data=TRUE)

    expect_equal(PreRep.test$result$box, ds.test)
    # linmod works
    ds.test <- mbecModelVariance(
        input.obj=dummy.test, model.vars=c("batch","group"), method="lm",
        type="clr")

    expect_equal(PreRep.test$result$linmod, ds.test)
    # linmixmod works
    ds.test <- mbecModelVariance(
        input.obj=dummy.test, model.vars=c("batch","group"), method="lmm",
        type="clr")

    expect_equal(PreRep.test$result$linmixmod, ds.test)
    # rda works
    ds.test <- mbecModelVariance(
        input.obj=dummy.test, model.vars=c("batch","group"), method="rda",
        type="clr")

    expect_equal(PreRep.test$result$rda, ds.test)
    # pvca works
    ds.test <- mbecModelVariance(
        input.obj=dummy.test, model.vars=c("batch","group"), method="pvca",
        type="clr")

    expect_equal(PreRep.test$result$pvca, ds.test)
    # scoef works
    ds.test <- mbecModelVariance(
        input.obj=dummy.test, model.vars=c("batch","group"), method="s.coef",
        type="clr")

    expect_equal(PreRep.test$result$scoef, ds.test)

    # check if it will create clr transformed values from list input
    PreRep.test <- evaluate_promise(
        mbecReportPrelim(input.obj=list(dummy.list$cnts[,seq(20)],
                                        dummy.list$meta),
                         model.vars=c("batch","group"),
                         type="clr", file.name=NULL, file.dir=NULL,
                         return.data=TRUE))

    # all analyses present
    expect_true(
        all(names(PreRep.test$result) %in%
                c("mosaic","pca","rle","heat","box",
                  "linmod","linmixmod","rda","pvca","scoef")))

    expect_true(all(unique(PreRep.test$result$linmod$type) %in% "clr"))

})


test_that("mbecReportPost works", {
    data(dummy.list)

    dummy.test <- mbecTransform( list(dummy.list$cnts[,seq(20)],
                                      dummy.list$meta), method="tss")
    dummy.test <- mbecTransform(dummy.test, method="clr")

    # expect an error if no corrections have been performed
    expect_error(
        mbecReportPost(input.obj=dummy.test, model.vars=c("batch","group"),
                       type="clr", file.name=NULL, file.dir=NULL,
                       return.data = TRUE), "No corrections available.")
    # Apply ComBat to the dummy
    dummy.corrected <- mbecCorrection(input.obj=dummy.test,
                                      model.vars=c("batch","group"),
                                      method="bat", type="clr" )
    # Get Post correction report data object
    PosRep.test <- evaluate_promise(
        mbecReportPost(input.obj=dummy.corrected, model.vars=c("batch","group"),
                       type="clr", file.name=NULL, file.dir=NULL,
                       return.data = TRUE))
    # all analyses present
    expect_true(
        all(names(PosRep.test$result) %in%
                c("mosaic","pca","rle","heat","box",
                  "linmod","linmixmod","rda","pvca","scoef")))
    # outputs contain clr and bat values
    expect_true(all(unique(PosRep.test$result$linmod$type) %in% c("clr","bat")))
    # mosaic stay the same
    expect_true(
        all(names(PosRep.test$result$mosaic) %in%
                c("Var1","Var2","Freq","Freq.scaled")))
    # rest contains data for clr and bat values
    expect_true(all(names(PosRep.test$result$pca) %in% c("pre","bat")))
    expect_true(all(names(PosRep.test$result$rle) %in% c("pre","bat")))
    expect_true(all(names(PosRep.test$result$heat) %in% c("pre","bat")))
    expect_true(all(names(PosRep.test$result$box) %in% c("pre","bat")))

})






