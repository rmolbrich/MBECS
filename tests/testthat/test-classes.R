# CONSTRUCTOR -------------------------------------------------------------

test_that("class constructor works", {
    data(dummy.list)
    constructor.res <-
        evaluate_promise(MbecData(cnt_table=dummy.list$cnts,
                                  meta_data=dummy.list$meta[,-4]))
    # correct class?
    expect_identical(class(constructor.res$result)[1], "MbecData")
    # correct package - which seems redundant but then again.. why not
    expect_identical(attr(class(constructor.res$result), "package"), "MBECS")
    # all the slots present?
    expect_identical(names(attributes(constructor.res$result)),
                     c("assessments","corrections","tss","clr","otu_table",
                       "tax_table","sam_data","phy_tree","refseq","class"))
    # worked without a hitch?
    expect_identical(constructor.res$warnings, character(0))
    expect_identical(constructor.res$messages,
                     "No 'sID' column present, creating from rownames now.\n")
    # are 'assessments' and 'corrections' initialized as empty lists?
    expect_identical(constructor.res$result@assessments, list())
    expect_identical(constructor.res$result@corrections, list())

    # Test correct failure
    # 1. test for non-numerical counts failure
    cnts.char <- dummy.list$cnts
    cnts.char[1,1] <- as.character(cnts.char[1,1])
    expect_error(MbecData(cnt_table=cnts.char,
                          meta_data=dummy.list$meta),
                 "All the values in your count-table need to be numeric!")
    # 2. no metadata supplied
    expect_error(MbecData(cnt_table=dummy.list$cnts, meta_data=NULL),
                 "You need to supply a meta-frame.")
})



# SETTER ------------------------------------------------------------------

test_that("tss setter works", {
    data(dummy.list)
    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta)
    # test tss setter
    set.res <- evaluate_promise(mbecSetData(input.obj=testdummy,
                                            new.cnts=dummy.list$cnts,
                                            type="tss", label=NULL))
    expect_identical(as.data.frame(t(set.res$result@tss)), dummy.list$cnts)
    # still works after rotation?
    set.res <- evaluate_promise(mbecSetData(input.obj=testdummy,
                                            new.cnts=t(dummy.list$cnts),
                                            type="tss", label=NULL))
    expect_identical(as.data.frame(t(set.res$result@tss)), dummy.list$cnts)
})


test_that("clr setter works", {

    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta)
    # test clr setter
    set.res <- evaluate_promise(mbecSetData(input.obj=testdummy,
                                            new.cnts=dummy.list$cnts,
                                            type="clr", label=NULL))
    expect_identical(as.data.frame(t(set.res$result@clr)), dummy.list$cnts)
    # still works after rotation?
    set.res <- evaluate_promise(mbecSetData(input.obj=testdummy,
                                            new.cnts=t(dummy.list$cnts),
                                            type="clr", label=NULL))
    expect_identical(as.data.frame(t(set.res$result@clr)), dummy.list$cnts)
})


test_that("assessments setter works", {

    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta)
    # test assessments setter
    set.res <- evaluate_promise(mbecSetData(input.obj=testdummy,
                                            new.cnts=dummy.list$cnts[1,],
                                            type="ass", label="asstest"))
    expect_identical(set.res$result@assessments[["asstest"]],
                     dummy.list$cnts[1,])
})


test_that("corrections setter works", {

    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta)
    # test corrections setter
    set.res <- evaluate_promise(mbecSetData(input.obj=testdummy,
                                            new.cnts=dummy.list$cnts,
                                            type="cor", label="cortest"))
    expect_identical(as.data.frame(t(set.res$result@corrections[["cortest"]])),
                     dummy.list$cnts)
    # still works after rotation?
    set.res <- evaluate_promise(mbecSetData(input.obj=testdummy,
                                            new.cnts=t(dummy.list$cnts),
                                            type="cor", label="cortest"))
    expect_identical(as.data.frame(t(set.res$result@corrections[["cortest"]])),
                     dummy.list$cnts)
})



# GETTER ------------------------------------------------------------------

test_that("otu getter works", {

    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta)
    # test otu getter
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy,
                    orientation="fxs",
                    required.col=c("group", "batch", "replicate"),
                    type="otu", label=NULL))

    expect_identical(
        get.res$result[[1]],
        data.frame(as(t(phyloseq::otu_table(testdummy)),"matrix"),
                   check.names=FALSE))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
    # still works after rotation?
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="sxf",
                    required.col=c("group","batch", "replicate"),
                    type="otu", label=NULL))

    expect_identical(
        get.res$result[[1]],
        data.frame(as(phyloseq::otu_table(testdummy),"matrix"),
                   check.names=FALSE))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
})


test_that("tss getter works", {

    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta,
                          tss=dummy.list$cnts)
    # test tss getter
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="fxs",
                    required.col=c("group","batch", "replicate"),
                    type="tss", label=NULL))

    expect_identical(get.res$result[[1]], as.data.frame(t(dummy.list$cnts)))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
    # still works after rotation?
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="sxf",
                    required.col=c("group","batch", "replicate"),
                    type="tss", label=NULL))

    expect_identical(get.res$result[[1]], as.data.frame(dummy.list$cnts))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
})


test_that("clr getter works", {

    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta,
                          clr=dummy.list$cnts)
    # test tss getter
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="fxs",
                    required.col=c("group","batch", "replicate"),
                    type="clr", label=NULL))

    expect_identical(get.res$result[[1]], as.data.frame(t(dummy.list$cnts)))

    expect_identical(names(
        get.res$result[[2]]), c("group","batch","replicate","sID"))
    # still works after rotation?
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="sxf",
                    required.col=c("group","batch", "replicate"),
                    type="clr", label=NULL))

    expect_identical(get.res$result[[1]], as.data.frame(dummy.list$cnts))
    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
})


test_that("assessment getter works", {

    testdummy <- mbecSetData(
        input.obj=MbecData(cnt_table=dummy.list$cnts,
                           meta_data=dummy.list$meta),
        new.cnts=dummy.list$cnts, type="ass", label="asstest")
    # test tss getter
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="fxs",
                    required.col=c("group","batch", "replicate"),
                    type="ass", label="asstest"))

    expect_identical(get.res$result[[1]], as.data.frame(t(dummy.list$cnts)))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
    # still works after rotation?
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="sxf",
                    required.col=c("group","batch", "replicate"),
                    type="ass", label="asstest"))

    expect_identical(get.res$result[[1]], as.data.frame(dummy.list$cnts))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
})


test_that("correction getter works", {

    testdummy <- mbecSetData(
        input.obj=MbecData(cnt_table=dummy.list$cnts,
                           meta_data=dummy.list$meta), new.cnts=dummy.list$cnts,
        type="cor", label="cortest")
    # test tss getter
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="fxs",
                    required.col=c("group","batch", "replicate"),
                    type="cor", label="cortest"))

    expect_identical(get.res$result[[1]], as.data.frame(t(dummy.list$cnts)))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
    # still works after rotation?
    get.res <- evaluate_promise(
        mbecGetData(input.obj=testdummy, orientation="sxf",
                    required.col=c("group","batch", "replicate"),
                    type="cor", label="cortest"))

    expect_identical(get.res$result[[1]], as.data.frame(dummy.list$cnts))

    expect_identical(
        names(get.res$result[[2]]), c("group","batch","replicate","sID"))
})



# PROCESSING --------------------------------------------------------------

test_that("MbecData process input works", {

    # MbecData input, i.e., check required.col test
    testdummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta)
    # no error
    expect_error(
        mbecProcessInput(testdummy,
                         required.col=c("group","batch", "replicate")), NA)
    # error
    expect_error(
        mbecProcessInput(testdummy,
                         required.col=c("missing")),
        "You need to supply a meta-frame that contains the columns: missing")
})


test_that("list process input works", {
    # 2. check list input
    pi.res <- mbecProcessInput(dummy.list,
                               required.col=c("group","batch", "replicate"))
    # is the same as the testdummy
    expect_identical(
        pi.res, MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta))
    # error with column verification
    expect_error(
        mbecProcessInput(dummy.list, required.col=c("missing")))
    # error wrong input
    expect_error(
        mbecProcessInput(dummy.list[1]))
})


test_that("phyloseq process input works", {
    data("dummy.ps")
    data("dummy.mbec")
    # 3. check phyloseq input
    mbec.dummy <- MbecData(cnt_table=dummy.list$cnts, meta_data=dummy.list$meta)

    # no error
    expect_identical(
        mbecProcessInput(dummy.ps,
                         required.col=c("group", "batch", "replicate")),
        mbec.dummy)
    # error missing column
    expect_error(mbecProcessInput(dummy.ps, required.col=c("missing")))
})

