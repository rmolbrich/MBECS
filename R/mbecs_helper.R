# HELPER FUNCTIONS --------------------------------------------------------

#' Check If Model Is Estimable
#'
#' Applies Limma's 'nonEstimable()' to a given model and returns NULL if
#' everything works out, or a warning and a vector of problematic covariates in
#' case there is a problem.
#'
#' The usefull part is that you can just put in all the covariates of interest
#' as model.vars and the function will build a simple linear model and its
#' model.matrix for testing. You can also provide more complex linear models
#' and the function will do the rest.
#'
#' @keywords limma nonEstimable Wrapper
#' @param input.obj MbecData, phyloseq or list (counts, meta-data).
#' @param model.vars Names of covariates to construct formula from.
#' @param model.form Formula for a linear model to test.
#' @return Either NULL if everything is fine or a vector of strings that denote
#' covariates and their respective problematic levels.
#'
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return NULL because it is estimable.
#' data(dummy.mbec)
#' eval.obj <- mbecTestModel(input.obj=dummy.mbec,
#' model.vars=c("group","batch"))
mbecTestModel <- function(input.obj, model.vars=NULL, model.form=NULL) {
    if( is.null(model.vars) && is.null(model.form) )
        stop("Please supply covariates and/or model-formula.")
    input.obj <- mbecProcessInput(input.obj, required.col=model.vars)

    # 1. extract covariate information
    tmp.meta <- mbecGetData(input.obj)[[2]]

    # 2. check if model-formula was supplied
    if( is.null(model.form) || !is(model.form, "formula") ) {
        # construct linear model from covariates
        message("Construct lm-formula from covariates.")

        model.form <- stats::as.formula(paste("y", " ~ ",
                                              paste(model.vars,
                                                    collapse=" + ")))
    }
    # if model.form is complete --> remove LHS
    if( length(model.form) == 3 ) model.form <- model.form[-2]

    # create model matrix from RHS
    model.mtx <- stats::model.matrix(model.form, tmp.meta)

    res.est <- limma::nonEstimable(model.mtx)

    if( !is.null(res.est) ) {
        message("There is a problem with the estimatibility of your model.
            Check out covariate: ", paste("'",res.est, "'", sep="",
                                          collapse=", "))
    }
    return(res.est)
}





#' Capitalize Word Beginning
#'
#' Change the first letter of the input to uppercase. Used in plotting functions
#' to make covariates, i.e., axis-labels look nicer.
#'
#' @keywords uppercase
#' @param input Any string whose first letter should be capitalized.
#' @return Input with first letter capitalized
mbecUpperCase <- function(input=character()) {
    return(gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",input,perl = TRUE))
}


#' Linear (Mixed) Model Feature to Batch Fit
#'
#' Helper function that fits lm/lmm with covariates 'treatment' and 'batch' to
#' every feature in the data-set. Returns the fdr corrected significance value
#' for the "treatment" variable. The method 'lm' will fit the linear model
#' \code{y ~ model.vars[1] + model.vars[2]} and the linear mixed model will
#' consider the second term as random effect, i.e.,
#' \code{y ~ model.vars[1] + (1|model.vars[2])}.
#'
#' The function returns either a plot-frame or the finished ggplot object.
#' Input for th data-set can be an MbecData-object, a phyloseq-object or a list
#' that contains counts and covariate data. The covariate table requires an
#' 'sID' column that contains sample IDs equal to the sample naming in the
#' counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Significance Linear Mixed Model Batch
#' @param input.obj MbecData object
#' @param model.vars Covariates of interest, first relates to batch and second
#' to treatment.
#' @param method Either 'lm' or 'lmm' for linear models and linear mixed models.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr, cor'.
#' DEFAULT is clr' and the use of 'cor' requires the parameter label to be
#' set as well.
#' @param label Which corrected abundance matrix to use for analysis in case
#' 'cor' was selected as type.
#' @return A vector of fdr corrected p-values that show significance of
#' treatment for every feature
#'
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return p-value for the linear model fit of every feature.
#' data(dummy.mbec)
#' val.score <- mbecLM(input.obj=dummy.mbec, model.vars=c("batch","group"),
#' method="lm")
mbecLM <- function(input.obj, method=c("lm","lmm"),
                   model.vars=c("batch","group"),
                   type=c("clr","otu","tss","cor"),
                   label=character()) {

    type <- match.arg(type)
    ## check and prepare inputs
    input.obj <- mbecProcessInput(input.obj, required.col=model.vars)
    tmp <- mbecGetData(input.obj = input.obj, orientation = "sxf",
                       required.col = eval(model.vars), type=eval(type),
                       label=eval(label))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    if( method == "auto" ) {
        # ToDo: sth.
        message("This feature is supposed to detect unbalanced designs and
                select lmm instead of lm. Alas, it has not been implemented.")

    } else if( method == "lm" ) {
        # fit lm to every feature with treatment and batch as as model
        # parameters then extract p-value for treatment
        tmp.group.p <- apply(tmp.cnts, 2, FUN = function(y) {
            # y ~ group + batch
            nc.lm <- stats::lm(y ~ get(model.vars[2]) + get(model.vars[1]),
                               data=tmp.meta)
            nc.lm.summary <- summary(nc.lm)
            # extract p-value of group (treatment)
            p <- nc.lm.summary$coefficients[2,4]
        })

    } else if( method == "lmm" ) {
        # overkill, but just keep it for now
        f.terms <- paste("(1|",model.vars[1],")", sep="")

        tmp.group.p <- apply(tmp.cnts, 2, FUN = function(y) {
            tmp.formula <- stats::as.formula(paste(
                paste("y", model.vars[-1], sep=" ~ "), paste(f.terms[],
                                                             collapse=" + "),
                sep=" + "))
            nc.lmm <- eval(bquote(lmerTest::lmer(.(tmp.formula),
                                                 data=tmp.meta)))
            nc.lmm.summary <- summary(nc.lmm)
            p <- nc.lmm.summary$coefficients[2,5]
        })
    }

    # correct for multiple testing
    tmp.group.p <- stats::p.adjust(tmp.group.p, method = 'fdr')

    return(tmp.group.p)
}


# TRANSFORMATION FUNCTIONS ------------------------------------------------


#' Normalizing Transformations
#'
#' Wrapper to help perform cumulative log-ratio  and total sum-scaling
#' transformations ,adapted from packages 'mixOmics' and robCompositions' to
#' work on matrices and Phyloseq objects alike.
#'
#' The function returns an MbecData object with transformed counts and covariate
#' information. Input for the data-set can be of type MbecData, phyloseq or a
#' list that contains counts and covariate data. Correct orientation of counts
#' will be handled internally, as long as both abundance table contain sample
#' names.
#'
#' @keywords CLR TSS Transformation
#' @param input.obj MbecData, phyloseq, list(counts, meta-data)
#' @param method one of 'CLR' or 'TSS'
#' @param offset (OPTIONAL) Offset in case of sparse matrix, for DEFAULT (0) an
#' offset will be calculated if required.
#' @param required.col (OPTIONAL) A vector of column names in the meta-data that
#' need to be present. Sanity check for subsequent steps.
#' @return MbecData with transformed counts in 'clr' and 'tss' attributes
#' respectively.
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return the cumulative log-ratio transformed counts in an
#' # MbecData object.
#' data(dummy.mbec)
#' mbec.CLR <- mbecTransform(input.obj=dummy.mbec, method="clr", offset=0,
#' required.col=c("batch","group"))
#'
#' # This will return total sum-scaled counts in an MbecData object.
#' mbec.CLR <- mbecTransform(input.obj=dummy.mbec, method="tss", offset=0,
#' required.col=c("batch","group"))
mbecTransform <- function(input.obj, method=c("clr","tss"),
                          offset=0, required.col=NULL) {
    ## 00. Check if 'method' was chosen correctly and get optional arguments
    method <- match.arg(method)
    ## VALIDATE input and change to 'MbecData' if needed
    input.obj <- mbecProcessInput(input.obj, required.col=eval(required.col))

    ## needs sxf orientation
    tmp <- mbecGetData(input.obj, orientation="sxf", type="otu")

    if( method == "clr" ) {

        tmp.cnts <- mbecCLR(tmp[[1]], offset = offset)
        # rebuild sample AND feature names for reassembly
        colnames(tmp.cnts) <- colnames(tmp[[1]])
        rownames(tmp.cnts) <- rownames(tmp[[1]])

        input.obj <- mbecSetData(input.obj, new.cnts=tmp.cnts, type="clr")
    } else if( method == "tss" ) {
        tmp.cnts <- t(apply(tmp[[1]], 1, function(x){x/sum(x)}))
        # rebuild sample AND feature names for reassembly
        colnames(tmp.cnts) <- colnames(tmp[[1]])
        rownames(tmp.cnts) <- rownames(tmp[[1]])

        input.obj <- mbecSetData(input.obj, new.cnts=tmp.cnts, type="tss")
    }

    return(input.obj)
}


#' Percentile Normalization
#'
#' Wrapper to help perform percentile normalization on a matrix of counts.
#' Takes counts and a data-frame of grouping variables and returns a matrix of
#' transformed counts. This is designed (by the Developers of the procedure) to
#' work with case/control experiments by taking the untreated group as reference
#' and adjusting the other groupings of TRT x Batch to it.
#'
#' The function returns a matrix of normalized abundances.
#'
#' @keywords Percentile Normalization
#' @param cnts A numeric matrix  of abundances (samples x features).
#' @param meta Data-frame of covariate columns, first column contains batches,
#' second column contains grouping.
#' @return Numeric matrix of corrected/normalized counts.
#'
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalized counts, according to the covariate
#' # information in meta
#' data(dummy.list)
#' mtx.pn_counts <- percentileNorm(cnts=dummy.list$cnts,
#' meta=dummy.list$meta[,c("batch","group")])
percentileNorm <- function(cnts, meta) {

    ref.group <- levels(meta[,2])[1]
    message("Group ",ref.group, " is considered control group, i.e., reference
          for normalization procedure. To change reference please 'relevel()'
          grouping factor accordingly.")

    norm.cnts <- cnts; norm.cnts[,] <- NA

    # for every batch
    for( b.idx in levels(meta[,1]) ) {
        # for every feature
        for( f.idx in seq_len(ncol(cnts)) ) {
            # which are the control-group values
            ctrl.group.vec <- cnts[which((meta[,2] %in% ref.group) &
                                             (meta[,1] %in% b.idx)), f.idx]
            # for every sample in the batch
            for( s.idx in which(meta[,1] %in% b.idx) ) {
                # call 'poscore' and get normalized value
                norm.cnts[s.idx, f.idx] <- poscore(ctrl.group.vec,
                                                   cnts[s.idx, f.idx],
                                                   "mean")
            }
        }
    }
    return(norm.cnts)
}


#' Percentile of Score
#'
#' Helper function that calculates percentiles of scores for batch-correction
#' method 'pn' (percentile normalization). R-implementation of Claire Duvallet's
#' 'percentileofscore()' for python.
#'
#' Calculates the number of values that bigger than reference (left) and the
#' number of values that are smaller than the reference (right). Percentiles of
#' scores are given in the interval \eqn{I:[0,100]}. Depending on type of
#' calculation, the score will be computed as follows:
#'
#' \code{rank = (right + left + ifelse(right > left, 1, 0)) * 50/n}
#'
#' \code{weak = right / n*100}
#'
#' \code{strict = left / n*100}
#'
#' \code{mean = (right + left) * 50/n)}
#'
#' @keywords Percentile Score
#' @param cnt.vec A vector of counts that acts as reference for score
#' calculation.
#' @param cnt A numeric value to calculate percentile-score for.
#' @param type One of 'rank', 'weak', 'strict' or 'mean' to determine how the
#' score is calculated.
#' @return A score for given count in relation to reference counts.
#'
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a score for the supplied vector with default evaluation
#' # (strict).
#' val.score <- poscore(cnt.vec=runif(100, min=0, max=100), cnt=42,
#' type="strict")
poscore <- function( cnt.vec, cnt, type=c("rank","weak","strict","mean") ) {
    # check argument
    type <- match.arg(type)
    # get number of cnts
    n <- length(cnt.vec)

    # if nothing to compare to, then return max value
    if( n == 0 ) {
        return(100)

    } else {
        left <- sum(cnt.vec < cnt)
        right <- sum(cnt.vec <= cnt)

        pos <- switch(type,
                      rank = (right + left + ifelse(right > left, 1, 0)) * 50/n,
                      weak = right / n*100,
                      strict = left / n*100,
                      mean = (right + left) * 50/n)
    }
    return(pos)
}



#' Centered Log-Ratio Transformation
#'
#' Internal function that performs CLR-transformation on input-matrix.
#' Formula is: clr(mtx) = ln( mtx / geometric_mean(mtx_samples))
#'
#' @keywords Log Ratio Transformation
#' @param input.mtx A matrix of counts (samples x features).
#' @param offset An (OPTIONAL) offset in case of sparse matrix. Function will
#' add an offset of 1/#features if matrix is sparse and offset not provided.
#' @return A matrix of transformed counts of same size and orientation as the
#' input.
mbecCLR <- function(input.mtx, offset = 0) {
    if( dim(input.mtx)[2] < 2 ) {
        message("No basis for transformation. Matrix contains less than 2
                features, returning unchanged.")
        return(input.mtx)
    }
    # 1. stop for negative values and NAs
    if( any(input.mtx < 0 | is.na(input.mtx)) ) {
        stop("Examine your data for NAs and negative values, CLR transformation
         requires complete positive values.\n")
    }
    if( any(input.mtx==0) & offset==0 ) {
        message("Found zeros, function will add a small pseudo-count
                (1/#features) for log-ratio transformation.")
        offset <- 1/ncol(input.mtx)
    }
    # add the offset
    input.mtx <- input.mtx + offset
    # 1: geometric-mean is n-th root of the product of all values for a sample
    gman <- apply(input.mtx, 1, function(s.row) exp(mean(log(s.row))))
    # 2: clr(mtx) <- [ ln(x_sample.i / geometric-mean(x_sample.i))]
    input.mtx <- log((input.mtx / gman))

    return(input.mtx)
}







