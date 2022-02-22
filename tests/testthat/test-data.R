# TEST DATA GENERATION ----------------------------------------------------

test_that("mbecDummy works", {

    mockup <- evaluate_promise(mbecDummy(n.otus=500, n.samples=40))
    # check output object
    expect_type(mockup$result, "list")
    expect_true(length(mockup$result) == 2)
    expect_true(all(names(mockup$result) %in% c("cnts","meta")))
    # check abundance table dimensions
    expect_true(all(dim(mockup$result$cnts) == c(40,500)))
    # check meta-data
    expect_true(all(names(mockup$result$meta) %in% c("group","batch",
                                                     "replicate")))
    expect_true(all(dim(mockup$result$meta) == c(40,3)))
})


# TEST DATA-SETS ----------------------------------------------------------

test_that("dummy.list works", {

    data("dummy.list")
    expect_type(dummy.list, "list")
    expect_true(length(dummy.list) == 2)
    expect_true(all(names(dummy.list) %in% c("cnts","meta")))
    # check abundance table dimensions
    expect_true(all(dim(dummy.list$cnts) == c(40,500)))
    # check meta-data
    expect_true(all(names(dummy.list$meta) %in% c("group","batch",
                                                     "replicate", "sID")))
    expect_true(all(dim(dummy.list$meta) == c(40,4)))
})


test_that("dummy.ps works", {

    data("dummy.ps")
    expect_s4_class(dummy.ps, "phyloseq")
    # check abundance table dimensions
    expect_true(all(dim(otu_table(dummy.ps)) == c(40,500)))
    # check meta-data
    expect_true(all(names(sample_data(dummy.ps)) %in% c("group","batch",
                                                  "replicate", "sID")))
    expect_true(all(dim(sample_data(dummy.ps)) == c(40,4)))
})


test_that("dummy.mbec works", {

    data("dummy.mbec")
    expect_s4_class(dummy.mbec, "MbecData")
    # check abundance table dimensions
    expect_true(all(dim(otu_table(dummy.mbec)) == c(40,500)))
    # check meta-data
    expect_true(all(names(sample_data(dummy.mbec)) %in% c("group","batch",
                                                        "replicate", "sID")))
    expect_true(all(dim(sample_data(dummy.mbec)) == c(40,4)))

    expect_true(all(dim(dummy.mbec@clr) == c(500,40)))
    expect_true(all(dim(dummy.mbec@tss) == c(500,40)))
})


