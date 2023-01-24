# CORRECTIONS WRAPPER -----------------------------------------------------

#' Run Correction Pipeline
#'
#' Run all correction algorithms selected by method and add corrected counts
#' as matrices to the data-set.
#'
#' @keywords Batch-Effect Correction Pipeline
#' @param input.obj Phyloseq object or a list that contains numeric matrix and
#' meta-data table. Requires sample names as row/col-names to handle correct
#' orientation.
#' @param model.vars Two covariates of interest to select by first variable
#' selects panels and second one determines coloring.
#' @param type One of 'otu', 'tss' or 'clr' to determine the abundance matrix
#' to use for evaluation.
#' @param method algorithms to use
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-3. If not supplied, the algorithm will use an 'lm'
#' to find pseudo-negative controls
#' @return an object of class MbecDataSet
#'
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This call will use 'ComBat' for batch effect correction and store the new
#' # counts in a list-obj in the output.
#' data(dummy.mbec)
#' study.obj <- mbecRunCorrections(input.obj=dummy.mbec,
#' model.vars=c("batch","group"), method=c("bat","bmc"))
#'
#' # This call will use 'Percentile Normalization' for batch effect correction
#' # and replace the old count matrix.
#' study.obj <- mbecRunCorrections(dummy.mbec, model.vars=c("batch","group"),
#' method=c("pn"))
mbecRunCorrections <- function(input.obj,
                               model.vars=c("batch","group"), type="clr",
                               method=c("ruv3","bmc","bat","rbe",
                                        "pn","svd","pls"),
                               nc.features=NULL) {

    input.obj <- mbecProcessInput(input.obj, required.col=eval(model.vars))

    if( all(dim(input.obj@clr) == 1) )
        input.obj <- mbecTransform(input.obj, method="clr")

    if( all(dim(input.obj@tss) == 1) )
        input.obj <- mbecTransform(input.obj, method="tss")

    for( m.idx in method ) {
        if( m.idx == "pn" ) {
            # because percentile norm is supposed to run on tss data
            input.obj <- mbecCorrection(input.obj = input.obj,
                                        model.vars = model.vars,
                                        method = "pn", type="tss")
        } else {
            input.obj <- mbecCorrection(input.obj = input.obj, type=eval(type),
                                        model.vars = model.vars,
                                        method=eval(m.idx),
                                        nc.features=eval(nc.features))
        }
    }

    return(input.obj)

}


#' Batch Effect Correction Wrapper
#'
#' Either corrects or accounts for (known) batch effects with one of several
#' algorithms.
#'
#' ASSESSMENT METHODS
#' The assessment methods 'lm, lmm, sva, ruv-2 and ruv-4" estimate the
#' significance of the batch effect and update the attribute 'assessments' with
#' vectors of p-values.
#'
#' Linear (Mixed) Models: A simple linear mixed model with covariates
#' 'treatment' and 'batch', or respective variables in your particular data-set,
#' will be fitted to each feature and the significance for the treatment
#' variable extracted.
#'
#' Surrogate variable Analysis (SVA): Surrogate Variable Analysis (SVA): Two
#' step approach that (1.) identify the number of latent factors to be estimated
#' by fitting a full-model with effect of interest and a null-model with no
#' effects. The function num.sv then calculates the number of latent factors.
#' In the next (2.) step, the sva function will estimate the surrogate
#' variables. And adjust for them in full/null-model . Subsequent F-test gives
#' significance values for each feature - these P-values and Q-values are
#' accounting for surrogate variables (estimated BEs).
#'
#' Remove unwanted Variation 2 (RUV-2): Estimates unknown BEs by using negative
#' control variables that, in principle, are unaffected by treatment/biological
#' effect, i.e., aka the effect of interest in an experiment. These variables
#' are generally determined prior to the experiment. An approach to RUV-2
#' without the presence of negative control variables is the estimation of
#' pseudo-negative controls. To that end, an lm or lmm (depending on whether or
#' not the study design is balanced) with treatment is fitted to each feature
#' and the significance calculated. The features that are not significantly
#' affected by treatment are considered as pseudo-negative control variables.
#' Subsequently, the actual RUV-2 function is applied to the data and returns
#' the p-values for treatment, considering unwanted BEs (whatever that means).
#'
#' Remove Unwanted Variation 4 (RUV-4): The updated version of RUV-2 also
#' incorporates the residual matrix (w/o treatment effect) to estimate the
#' unknown BEs. To that end it follows the same procedure in case there are no
#' negative control variables and computes pseudo-controls from the data via
#' l(m)m. As RUV-2, this algorithm also uses the parameter 'k' for the number of
#' latent factors. RUV-4 brings the function 'getK()' that estimates this factor
#' from the data itself. The calculated values are however not always reliable.
#' A value of k=0 fo example can occur and should be set to 1 instead. The
#' output is the same as with RUV-2.
#'
#' CORRECTION METHODS
#' The correction methods 'ruv3, bmc, bat, rbe, pn, svd' attempt to mitigate
#' the batch effect and update the attribute 'corrections' with the resulting
#' abundance matrices of corrected counts.
#'
#' Remove Unwanted Variation 3 (RUV-3): This algorithm requires negative
#' control-features, i.e., OTUs that are known to be unaffected by the batch
#' effect, as well as technical replicates. The algorithm will check for the
#' existence of a replicate column in the covariate data. If the column is not
#' present, the execution stops and a warning message will be displayed.
#'
#' Batch Mean Centering (BMC): For known BEs, this method takes the batches,
#' i.e., subgroup of samples within a particular batch, and centers them to
#' their mean.
#'
#' Combat Batch Effects (ComBat): This method uses an non-/parametric empirical
#' Bayes framework to correct for BEs. Described by Johnson et al. 2007 this
#' method was initially conceived to work with gene expression data and is part
#' of the sva-package in R.
#'
#' Remove Batch Effects (RBE): As part of the limma-package this method was
#' designed to remove BEs from Microarray Data. The algorithm fits the full-
#' model to the data, i.e., all relevant covariates whose effect should not be
#' removed, and a model that only contains the known BEs. The difference between
#' these models produces a residual matrix that (should) contain only the full-
#' model-effect, e.g., treatment. As of now the mbecs-correction only uses the
#' first input for batch-effect grouping. ToDo: think about implementing a
#' version for more complex models.
#'
#' Percentile Normalization (PN): This method was actually developed
#' specifically to facilitate the integration of microbiome data from different
#' studies/experimental set-ups. This problem is similar to the mitigation of
#' BEs, i.e., when collectively analyzing two or more data-sets, every study is
#' effectively a batch on its own (not withstanding the probable BEs within
#' studies). The algorithm iterates over the unique batches and converts the
#' relative abundance of control samples into their percentiles. The relative
#' abundance of case-samples within the respective batches is then transformed
#' into percentiles of the associated control-distribution. Basically, the
#' procedure assumes that the control-group is unaffected by any effect of
#' interest, e.g., treatment or sickness, but both groups within a batch are
#' affected by that BE. The switch to percentiles (kinda) flattens the effective
#' difference in count values due to batch - as compared to the other batches.
#' This also introduces the two limiting aspects in percentile normalization. It
#' can only be applied to case/control designs because it requires a reference
#' group. In addition, the transformation into percentiles removes information
#' from the data.
#'
#' Singular Value Decomposition (SVD): Basically perform matrix factorization
#' and compute singular eigenvectors (SEV). Assume that the first SEV captures
#' the batch-effect and remove this effect from the data. The interesting thing
#' is that this works pretty well (with the test-data anyway) But since the SEVs
#' are latent factors that are (most likely) confounded with other effects it is
#' not obvious that this is the optimal approach to solve this issue.
#'
#' Principal Least Squares Discriminant Analysis (PLSDA)
#' This function estimates latent dimensions from the explanatory matrix
#' \code{X}. The latent dimensions are maximally associated with the outcome
#' matrix \code{Y}. It is a built-in function of \code{PLSDA_batch} and has been
#' adjusted to work in the MBECS-package. To that end, the function
#' \code{mixOmics::explained_variance} was replaced with a computation based on
#' \code{vegan::cca} since this is already used in the MBECS package.
#' Additionally, the matrix deflation function was replaced with own code. The
#' credit for algorithm and implementation goes to
#' 'https://github.com/EvaYiwenWang/PLSDAbatch' and the associated publication
#' that is referenced in the documentation and vignette.
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be as input, but assessments or corrections-lists will
#' contain the result of the respective chosen method.
#'
#' @keywords Batch-Effect Correction and Assessment
#' @param input.obj An MbecData object with 'tss' and 'clr' matrices.
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param method Denotes the algorithms to use. One of 'lm, lmm, sva, ruv2,
#' ruv4' for assessment methods or one of 'ruv3, bmc, bat, rbe, pn, svd', 'cqr'
#' for correction algorithms.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr' but percentile normalization is supposed to work on tss-abundances.
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-2/3/4. If not supplied, the algorithm will use a
#' linear model to find pseudo-negative controls
#' @return An updated object of class MbecData.
#'
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This call will use 'ComBat' for batch effect correction on CLR-transformed
#' # abundances and store the new counts in the 'corrections' attribute.
#' data(dummy.mbec)
#' study.obj <- mbecCorrection(input.obj=dummy.mbec,
#' model.vars=c("batch","group"), method="bat", type="clr")
#'
#' # This call will use 'Percentile Normalization' for batch effect correction
#' # on TSS-transformed counts and store the new counts in the 'corrections'
#' # attribute.
#' study.obj <- mbecCorrection(dummy.mbec, model.vars=c("batch","group"),
#' method="pn", type="tss")
mbecCorrection <- function(input.obj, model.vars=c("batch","group"),
                           method=c("lm","lmm","sva","ruv2","ruv4","ruv3","bmc",
                                    "bat","rbe","pn","svd","pls"),
                           type=c("clr","otu","tss"), nc.features=NULL) {

    input.obj <- mbecProcessInput(input.obj, required.col=eval(model.vars))

    ##  Check if 'method' was chosen correctly.
    method <- match.arg(method)
    type <- match.arg(type)
    ## START WITH ASSESSMENT METHODS
    if( method == "lm" ) {
        message("Applying Linear Model (LM).")
        res.assess <- mbecLM(input.obj = input.obj, method = "lm",
                             model.vars=model.vars,
                             type=eval(type))
    } else if( method == "lmm" ) {
        message("Applying Linear Mixed Model (LMM)")
        res.assess <- mbecLM(input.obj = input.obj, method = "lmm",
                             model.vars=model.vars,
                             type=eval(type))
    } else if( method == "sva" ) {
        res.assess <- mbecSVA(input.obj, model.vars[-1],
                              type=eval(type))
    } else if( method == "ruv2" ) {
        res.assess <- mbecRUV2(input.obj, model.vars,
                               type=eval(type), nc.features)
    } else if( method == "ruv4" ) {
        res.assess <- mbecRUV4(input.obj, model.vars,
                               type=eval(type), nc.features)
        ## CORRECTION METHODS FROM HERE
    } else if( method == "ruv3" ) {
        res.correct <- mbecRUV3(input.obj, model.vars,
                                type=eval(type), nc.features)
    } else if( method == "bmc" ) {
        res.correct <- mbecBMC(input.obj, model.vars,
                               type=eval(type))
    } else if( method == "bat" ) {
        res.correct <- mbecBat(input.obj, model.vars,
                               type=eval(type))
    } else if( method == "rbe" ) {
        res.correct <- mbecRBE(input.obj, model.vars,
                               type=eval(type))
    } else if( method == "pn" ) {
        # !This is supposed to work on TSS counts
        res.correct <- mbecPN(input.obj, model.vars,
                              type=eval(type))
    } else if( method == "svd" ) {
        res.correct <- mbecSVD(input.obj, model.vars,
                               type=eval(type))
    } else if( method == "pls" ) {
        res.correct <- mbecPLSDA(input.obj, model.vars,
                                 type=eval(type))
    }

    ## figure out where to put the result
    if( method %in% c("lm","lmm","sva","ruv2","ruv4") ) {
        return.obj <- mbecSetData(input.obj = input.obj, new.cnts = res.assess,
                                  type = "ass", label = eval(method))
    } else if( method %in% c("ruv3","bmc","bat","rbe","pn","svd","pls") ) {
        return.obj <- mbecSetData(input.obj = input.obj, new.cnts = res.correct,
                                  type = "cor", label = eval(method))
    }

    return(return.obj)
    ## Point of no return. *badum tsss*
}


# ASSESSMENT FUNCTIONS ----------------------------------------------------

#' Surrogate variable Analysis (SVA)
#'
#' Two step approach that (1.) identify the number of latent factors to be
#' estimated by fitting a full-model with effect of interest and a null-model
#' with no effects. The function 'num.sv()' then calculates the number of latent
#' factors. In the next (2.) step, the sva function will estimate the surrogate
#' variables. And adjust for them in full/null-model . Subsequent F-test gives
#' significance values for each feature - these P-values and Q-values are
#' accounting for surrogate variables (estimated BEs).
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a vector of p-values.
#'
#' @keywords Batch-Effect Assessment SVA
#' @param input.obj MbecData object
#' @param model.vars Vector of covariate names. First element relates to
#' variable of interest.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @return A vector of p-values that indicate significance of the batch-effect
#' for the features.
#'
#' @include mbecs_classes.R
mbecSVA <- function(input.obj, model.vars, type=c("clr","otu","tss")) {

    message("Applying Surrogate Variable Analysis (SVA).")
    # this should only match the uncorrected abundance matrices
    type <- match.arg(type)
    tmp <- mbecGetData(input.obj=input.obj, orientation="fxs", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    tmp.mod <- stats::model.matrix( ~ tmp.meta[[model.vars[1]]])
    tmp.mod0 <- stats::model.matrix( ~ 1, data=tmp.meta[[model.vars[1]]])
    tmp.sva.n <- sva::num.sv(dat = as.matrix(tmp.cnts), mod = tmp.mod)
    tmp.sva <- sva::sva(as.matrix(tmp.cnts), tmp.mod, tmp.mod0, n.sv=tmp.sva.n)

    message("Number of significant surrogate variables is: ", tmp.sva$n.sv)
    # do sth. that I have to understand yet
    tmp.mod.ba <- cbind(tmp.mod, tmp.sva$sv)
    tmp.mod0.ba <- cbind(tmp.mod0, tmp.sva$sv)
    tmp.sva.trt_p <- sva::f.pvalue(as.matrix(tmp.cnts), tmp.mod.ba, tmp.mod0.ba)
    tmp.sva.trt_adjp <- stats::p.adjust(tmp.sva.trt_p, method = 'fdr')

    return(tmp.sva.trt_adjp)
}


# RUV FUNCTIONS -----------------------------------------------------------

#' Remove unwanted Variation 2 (RUV-2)
#'
#' Estimates unknown BEs by using negative control variables that, in principle,
#' are unaffected by treatment/study/biological effect (aka the effect of
#' interest in an experiment). These variables are generally determined prior
#' to the experiment. An approach to RUV-2 without the presence of negative
#' control variables is the estimation of pseudo-negative controls. To that end
#' an lm or lmm (depending on whether or not the study design is balanced) with
#' treatment is fitted to each feature and the significance calculated. The
#' features that are not significantly affected by treatment are considered as
#' pseudo-negative control variables. Subsequently, the actual RUV-2 function
#' is applied to the data and returns the p-values for treatment, considering
#' unwanted BEs (whatever that means).
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a vector of p-values.
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-3. If not supplied, the algorithm will use an 'lm'
#' to find pseudo-negative controls
#' @return A vector of p-values that indicate significance of the batch-effect
#' for the features.
#'
#' @include mbecs_classes.R
mbecRUV2 <- function(input.obj, model.vars, type=c("clr","otu","tss"),
                     nc.features=NULL) {

    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="sxf", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    ## check for neg. controls or generate pseudo negatives
    if( !is.null(nc.features) ) {
        tmp.nc <- colnames(tmp.cnts) %in% nc.features
        names(tmp.nc) <- colnames(tmp.cnts)
    } else {
        message("No negative control features provided.
                Using pseudo-negative controls.")
        tmp.group.p <- mbecLM(input.obj=input.obj, method="lm",
                              model.vars=model.vars, type=eval(type))
        tmp.nc <- tmp.group.p > 0.05
    }

    message("Applying Remove Unwanted Variantion v2 (RUV-II).")
    tmp.ruv2 <- ruv::RUV2(Y = tmp.cnts, X = tmp.meta[[model.vars[1]]],
                          ctl = tmp.nc, k = 3) # k is subjective

    tmp.ruv2.trt_p <- tmp.ruv2$p
    tmp.ruv2.trt_adjp <- stats::p.adjust(tmp.ruv2.trt_p, method="fdr")
    names(tmp.ruv2.trt_adjp) <- phyloseq::taxa_names(input.obj)
    return(tmp.ruv2.trt_adjp)
}


#' Remove Unwanted Variation 4 (RUV-4)
#'
#' The updated version of RUV-2 also incorporates the residual matrix
#' (w/o treatment effect) to estimate the unknown BEs. To that end it follows
#' the same procedure in case there are no negative control variables and
#' computes pseudo-controls from the data via l(m)m. As RUV-2, this algorithm
#' also uses the parameter 'k' for the number of latent factors. RUV-4 brings
#' the function 'getK()' that estimates this factor from the data itself. The
#' calculated values are however not always reliable. A value of k=0 fo example
#' can occur and should be set to 1 instead.
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a vector of p-values.
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-3. If not supplied, the algorithm will use an 'lm'
#' to find pseudo-negative controls
#' @return A vector of p-values that indicate significance of the batch-effect
#' for the features.
#'
#' @include mbecs_classes.R
mbecRUV4 <- function(input.obj, model.vars, type=c("clr","otu","tss"),
                     nc.features=NULL) {
    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="sxf", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    ## check for neg. controls or generate pseudo negatives
    if( !is.null(nc.features) ) {
        tmp.nc <- colnames(tmp.cnts) %in% nc.features
        names(tmp.nc) <- colnames(tmp.cnts)
    } else {
        message("No negative control features provided.
                Using pseudo-negative controls.")
        tmp.group.p <- mbecLM(input.obj=input.obj, method="lm",
                              model.vars=model.vars, type=eval(type))
        tmp.nc <- tmp.group.p > 0.05
    }

    message("Applying Remove Unwanted Variantion v4 (RUV-IV).")
    tmp.k <- ruv::getK(Y=tmp.cnts, X=tmp.meta$group, ctl=tmp.nc)
    tmp.k <- ifelse(tmp.k$k !=0, tmp.k$k, 1)

    tmp.ruv4 <- ruv::RUV4(Y=tmp.cnts, X=tmp.meta$group, ctl=tmp.nc, k=tmp.k)
    tmp.ruv4.trt_p <- tmp.ruv4$p
    tmp.ruv4.trt_adjp <- stats::p.adjust(tmp.ruv4.trt_p, method="fdr")

    names(tmp.ruv4.trt_adjp) <- phyloseq::taxa_names(input.obj)


    return(tmp.ruv4.trt_adjp)
}


#' Remove Unwanted Variation 3 (RUV-3)
#'
#' This algorithm requires negative control-features, i.e., OTUs that are known
#' to be unaffected by the batch effect, as well as technical replicates. The
#' algorithm will check for the existence of a replicate column in the covariate
#' data. If the column is not present, the execution stops and a warning message
#' will be displayed.
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a matrix of corrected abundances.
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-3. If not supplied, the algorithm will use an 'lm'
#' to find pseudo-negative controls
#' @return A matrix of batch-effect corrected counts
#'
#' @include mbecs_classes.R
mbecRUV3 <- function(input.obj, model.vars, type=c("clr","otu","tss"),
                     nc.features=NULL) {

    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="sxf", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    ## check for neg. controls or generate pseudo negatives
    if( !is.null(nc.features) ) {
        tmp.nc <- colnames(tmp.cnts) %in% nc.features
        names(tmp.nc) <- colnames(tmp.cnts)
    } else {
        message("No negative control features provided.
                Using pseudo-negative controls.")
        tmp.group.p <- mbecLM(input.obj=input.obj, method="lm",
                              model.vars=model.vars, type=eval(type))
        tmp.nc <- tmp.group.p > 0.05
    }

    message("Applying Remove Unwanted Variantion v3 (RUV-III).")
    tmp.replicate <- tmp.meta$replicate
    if( length(tmp.replicate) == length(unique(tmp.replicate)) ) {
        warning("No technical replicates found. RUV-3 is not available!")
        return(NULL)
    } else {
        tmp.replicate.matrix <- ruv::replicate.matrix(tmp.replicate)
    }

    corrected.cnts <- ruv::RUVIII(Y=tmp.cnts, M=tmp.replicate.matrix,
                                  ctl=tmp.nc)
    rownames(corrected.cnts) <- rownames(tmp.cnts)

    return(corrected.cnts)
}


# CORRECTION FUNCTIONS ----------------------------------------------------

#' Batch Mean Centering (BMC)
#'
#' For known BEs, this method takes the batches, i.e., subgroup of samples
#' within a particular batch, and centers them to their mean.
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a matrix of corrected abundances.
#'
#' @keywords BECA Batch Mean Centering
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @return A matrix of batch-effect corrected counts
#'
#' @include mbecs_classes.R
mbecBMC <- function(input.obj, model.vars, type=c("clr","otu","tss")) {
    message("Applying Batch Mean-Centering (BMC) for batch-correction.")

    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="sxf", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # get unique batches
    input.batches <- unique(tmp.meta[[model.vars[1]]])

    corrected.cnts <- NULL
    for( batch.idx in input.batches ) {
        tmp <- scale(tmp.cnts[tmp.meta$sID[
            tmp.meta[[model.vars[1]]] %in% batch.idx], ],
            center=TRUE, scale=FALSE)

        corrected.cnts <- rbind.data.frame(corrected.cnts, tmp)
    }

    # reorder according to meta-table
    corrected.cnts <- corrected.cnts[tmp.meta$sID,]

    return(corrected.cnts)
}


#' Combat Batch Effects (ComBat)
#'
#' This method uses an non-/parametric empirical  Bayes framework to correct
#' for BEs. Described by Johnson et al. 2007 this method was initially conceived
#' to work with gene expression data and is part of the sva-package in R.
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a matrix of corrected abundances.
#'
#' @keywords BECA Batch Mean Centering
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @return A matrix of batch-effect corrected counts
#'
#' @include mbecs_classes.R
mbecBat <- function(input.obj, model.vars, type=c("clr","otu","tss")) {

    message("Applying ComBat (sva) for batch-correction.")
    ## ComBat requires 'fxs' orientation for inputs
    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="fxs", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # ToDo: tmp model
    if( length(model.vars) == 1 ) {
        corrected.cnts <- sva::ComBat(tmp.cnts, batch=tmp.meta[[model.vars[1]]],
                                      mod=1, par.prior=FALSE,
                                      prior.plots=FALSE)
    } else {
        tmp.mod <- stats::model.matrix( ~ tmp.meta[[model.vars[2]]]) #full model

        corrected.cnts <- sva::ComBat(tmp.cnts, batch=tmp.meta[[model.vars[1]]],
                                      mod=tmp.mod, par.prior=FALSE,
                                      prior.plots=FALSE)
    }
    return(corrected.cnts)
}


#' Remove Batch Effects (RBE)
#'
#' As part of the limma-package this method was designed to remove BEs from
#' Microarray Data. The algorithm fits the full-model to the data, i.e., all
#' relevant covariates whose effect should not be removed, and a model that only
#' contains the known BEs. The difference between these models produces a
#' residual matrix that (should) contain only the full-model-effect, e.g.,
#' treatment. As of now the mbecs-correction only uses the first input for
#' batch-effect grouping. ToDo: think about implementing a version for more
#' complex models.
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a matrix of corrected abundances.
#'
#' @keywords BECA Limma Remove Batch Effects
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @return A matrix of batch-effect corrected counts
#'
#' @include mbecs_classes.R
mbecRBE <- function(input.obj, model.vars, type=c("clr","otu","tss")) {

    message("Applying 'removeBatchEffect' (limma) for batch-correction.")
    ## removeBatchEffect requires 'fxs' orientation for inputs
    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="fxs", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    if( length(model.vars) == 1 ) {
        corrected.cnts <- limma::removeBatchEffect(
            tmp.cnts,
            batch=tmp.meta[[model.vars[1]]],
            design=matrix(1,ncol(tmp.cnts),1))
    } else {
        tmp.mod <- stats::model.matrix( ~ tmp.meta[[model.vars[2]]])

        corrected.cnts <- limma::removeBatchEffect(
            tmp.cnts, batch=tmp.meta[[model.vars[1]]], design=tmp.mod)
    }
    return(corrected.cnts)
}


#' Percentile Normalization (PN)
#'
#' This method was actually developed specifically to facilitate the integration
#' of microbiome data from different studies/experimental set-ups. This problem
#' is similar to the mitigation of BEs, i.e., when collectively analyzing two or
#' more data-sets, every study is effectively a batch on its own
#' (not withstanding the probable BEs within studies). The algorithm iterates
#' over the unique batches and converts the relative abundance of control
#' samples into their percentiles. The relative abundance of case-samples within
#' the respective batches is then transformed into percentiles of the associated
#' control-distribution. Basically, the procedure assumes that the control-group
#' is unaffected by any effect of interest, e.g., treatment or sickness, but
#' both groups within a batch are affected by that BE. The switch to percentiles
#' (kinda) flattens the effective difference in count values due to batch - as
#' compared to the other batches. This also introduces the two limiting aspects
#' in percentile normalization. It can only be applied to case/control designs
#' because it requires a reference group. In addition, the transformation into
#' percentiles removes information from the data.
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a matrix of corrected abundances.
#'
#' @keywords BECA Duvallet Percentile Normalisation
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'tss'.
#' @return A matrix of batch-effect corrected counts
#'
#' @include mbecs_classes.R
mbecPN <- function(input.obj, model.vars, type=c("clr","otu","tss")) {

    message("Applying Percentile Normalization (PN).")
    ## Percentile Normalization requires 'sxf' orientation for inputs
    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="sxf", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    if( !is(tmp.meta[,eval(model.vars[1])], "factor") ) {
        warning("Batch variable is not a factor! Converting to factor now.")
        tmp.meta[,eval(model.vars[1])] <- as.factor(
            tmp.meta[,eval(model.vars[1])])
    }


    if( !is(tmp.meta[,eval(model.vars[2])], "factor") ) {
        warning("Grouping variable is not a factor! Converting to factor now.")
        tmp.meta[,eval(model.vars[2])] <- as.factor(
            tmp.meta[,eval(model.vars[2])])
    }

    # check for case/control design
    if( nlevels(tmp.meta[,eval(model.vars[2])]) != 2 ) {
        warning("Grouping/Treatment contains ",
                nlevels(tmp.meta[,eval(model.vars[2])]),
                " different categories. Percentile normalization is designed to
                work with 2 classes only, i.e., case/control studies.")
        return(NULL)
    }
    # check/adjust for zero-values
    if( any(tmp.cnts == 0) ) {
        warning("Abundances contain zero values. Adding small uniform offset.")
        tmp.cnts <- apply(tmp.cnts, c(1,2),
                          function(x) ifelse(x == 0, stats::runif(1,0,10^-9),x))
    }
    # do normalisation
    corrected.cnts <- percentileNorm(tmp.cnts, tmp.meta[,eval(model.vars)])

    return(corrected.cnts)
}


#' Singular Value Decomposition (SVD)
#'
#' Basically perform matrix factorization and compute singular eigenvectors
#' (SEV). Assume that the first SEV captures the batch-effect and remove this
#' effect from the data. The interesting thing is that this works pretty well.
#' But since the SEVs are latent factors that are (most likely) confounded with
#' other effects it is not obvious to me that this is the optimal approach to
#' solve this issue.
#'
#' ToDo: IF I find the time to works on "my-own-approach" then this is the
#' point to start from!!!
#'
#' The input for this function is supposed to be an MbecData object that
#' contains total sum-scaled and cumulative log-ratio transformed abundance
#' matrices. Output will be a matrix of corrected abundances.
#'
#' @keywords Singular Value Decomposition
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @return A matrix of batch-effect corrected counts
#'
#' @include mbecs_classes.R
mbecSVD <- function(input.obj, model.vars, type=c("clr","otu","tss")) {

    message("Applying Singular Value Decomposition (SVD) for batch-correction.")
    ## SVD requires 'sxf' orientation for inputs
    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="sxf", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]
    # sd, mean, scale
    tmp.sd <- apply(tmp.cnts, 2, stats::sd)
    tmp.mean <- apply(tmp.cnts, 2, mean)
    # center and scale
    tmp.cnts.scale <- scale(tmp.cnts, center=TRUE, scale=TRUE)
    # produce square matrix
    tmp.cnts.square <- crossprod(tmp.cnts.scale)
    # apply singular value decomposition
    svd.res <- svd(tmp.cnts.square)

    # extract 1st singular vectors
    svd.u1 <- svd.res$u[ ,1]
    svd.v1 <- svd.res$v[ ,1]

    # deflate component 1 from the data
    tmp.t1 <- tmp.cnts.scale %*% svd.u1 / drop(sqrt(crossprod(svd.u1)))
    tmp.c1 <- crossprod(tmp.cnts.scale, tmp.t1) / drop(crossprod(tmp.t1))
    tmp.svd.defl  <- tmp.cnts.scale - tmp.t1 %*% t(tmp.c1)

    # add mean and standard deviation
    corrected.cnts <- tmp.cnts
    corrected.cnts[,] <- NA

    for( c.idx in seq_len(dim(corrected.cnts)[2]) ) {
        for( r.idx in seq_len(dim(corrected.cnts)[1]) ) {
            corrected.cnts[r.idx, c.idx] <- tmp.svd.defl[r.idx, c.idx] *
                tmp.sd[c.idx] + tmp.mean[c.idx]
        }
    }

    return(corrected.cnts)
}




#' Partial Least Squares Discriminant Analysis
#'
#' This function estimates latent dimensions from the explanatory matrix
#' \code{X}. The latent dimensions are maximally associated with the outcome
#' matrix \code{Y}. It is a built-in function of \code{PLSDA_batch} and has been
#' adjusted to work in the MBECS-package. To that end, the function
#' \code{mixOmics::explained_variance} was replaced with a computation based on
#' \code{vegan::cca} since this is already used in the MBECS package.
#' Additionally, the matrix deflation function was replaced with own code. The
#' near zero-variance correction function is taken from the caret -package. The
#' credit for algorithm and implementation goes to
#' 'https://github.com/EvaYiwenWang/PLSDAbatch' and the associated publication
#' that is referenced in the documentation and vignette.
#'
#' @keywords PLSDA batch correction
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars Vector of covariate names. First element relates to batch.
#' @param type Which abundance matrix to use, one of 'otu, tss, clr'. DEFAULT is
#' 'clr'.
#' @return A matrix of batch-effect corrected counts
#'
#' @include mbecs_classes.R
mbecPLSDA <- function(input.obj, model.vars, type=c("clr","otu","tss")) {

    . <- Freq <- weight <- NULL

    ## ToDo: Streamline this
    message("Applying PLSDA for batch-correction.")

    type <- match.arg(type)
    tmp <- mbecGetData(input.obj, orientation="sxf", type=eval(type))
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    ### first set-up the parameters that we do not want to deal with right now
    ncomp.trt = 2
    ncomp.bat = 2
    keepX.trt = rep(ncol(tmp.cnts), ncomp.trt)
    keepX.bat = rep(ncol(tmp.cnts), ncomp.bat)
    max.iter = 500
    tol = 1e-06
    near.zero.var = TRUE
    balance = TRUE
    ### first set-up the parameters that we do not want to deal with right now


    n = nrow(tmp.cnts)
    # make sure covariates are factors
    tmp.meta <- mbecHelpFactor(tmp.meta, model.vars)

    ## this will keep the row-names.. which is nice
    Y.bat.mat <- stats::model.matrix(~ 0 + tmp.meta[[model.vars[1]]], tmp.meta)
    ## we also want to set colnames
    colnames(Y.bat.mat) = levels(tmp.meta[[model.vars[1]]])
    #levels in batch-factor
    q.bat = ncol(Y.bat.mat)

    # basically remove all the features that exhibit variance close to zero
    if (near.zero.var == TRUE) {
        nzv = caret::nearZeroVar(tmp.cnts)
        if (length(nzv > 0)) {
            warning("Zero- or near-zero variance predictors.
                    \nReset predictors matrix\n
                    to not near-zero variance predictors.\nSee $nzv\n
                    for problematic predictors.")
            tmp.cnts = tmp.cnts[, -nzv$Position, drop = FALSE]
            if (ncol(tmp.cnts) == 0) {
                stop("No more predictors")
            }
        }
    }
    # basically remove all the features that exhibit variance close to zero

    # limit number of batch associated dimensions to number of features
    n.taxa = ncol(tmp.cnts) # number of features
    if (ncomp.bat > n.taxa) {
        warning("Reset maximum number of variates 'ncomp.bat' to ncol(X) = ",
                n.taxa, ".")
        ncomp.bat = n.taxa
    }
    if (length(keepX.bat) != ncomp.bat) {
        stop("length of 'keepX.bat' must be equal to ", ncomp.bat,
             ".")
    }
    if (any(keepX.bat > n.taxa)) {
        stop("each component of 'keepX' must be lower than or equal to ",
             n.taxa, ".")
    }

    ## ToDo: take care of this weight mess
    tmp.meta$weight <- 1

    # if treatment is given, this part will preserve associated variation
    if( length(model.vars) >= 2 ) {
        ### build result vector/matrices
        ## this will keep the row-names.. which is nice
        Y.trt.mat <- stats::model.matrix(~ 0 + tmp.meta[[model.vars[2]]], tmp.meta)
        ## we also want to set colnames
        colnames(Y.trt.mat) = levels(tmp.meta[[model.vars[2]]])
        #levels in batch-factor
        q.trt = ncol(Y.trt.mat)


        if (ncomp.trt > n.taxa) {
            warning("Reset max. number of variates 'ncomp.trt' to ncol(X) = ",
                    n.taxa, ".")
            ncomp.trt = n.taxa
        }
        if (length(keepX.trt) != ncomp.trt) {
            stop("length of 'keepX.trt' must be equal to ", ncomp.trt,
                 ".")
        }
        if (any(keepX.trt > n.taxa)) {
            stop("each component of 'keepX' must be lower than or equal to ",
                 n.taxa, ".")
        }

        ## UNBALANCED DESIGN
        # for unbalanced design, adjust the weight matrix
        if (!balance) {

            tmp.meta <- base::table(tmp.meta[, eval(model.vars[1])],
                                    tmp.meta[, eval(model.vars[2])],
                                    dnn = eval(model.vars)) %>%
                as.data.frame() %>%
                dplyr::right_join(., tmp.meta,
                                  by =eval(model.vars)) %>%
                dplyr::mutate(weight = 1/Freq) %>%
                dplyr::mutate(weight = sqrt(weight/min(weight)))

            # get a weighted matrix
            tmp.cnts.scale <- scale(tmp.meta$weight * tmp.cnts, center = TRUE,
                                    scale = TRUE)
            tmp.cnts.mean <- attributes(tmp.cnts.scale)$`scaled:center`
            tmp.cnts.sd <- attributes(tmp.cnts.scale)$`scaled:scale`

            Y.bat.scale = scale(tmp.meta$weight * Y.bat.mat, center = TRUE,
                                scale = TRUE)
            Y.trt.scale = scale(tmp.meta$weight * Y.trt.mat, center = TRUE,
                                scale = TRUE)
        }
        # in case of a balanced design - no weighting
        else {
            tmp.cnts.scale <- scale(tmp.cnts, center = TRUE, scale = TRUE)
            tmp.cnts.mean <- attributes(tmp.cnts.scale)$`scaled:center`
            tmp.cnts.sd <- attributes(tmp.cnts.scale)$`scaled:scale`
            Y.bat.scale = scale(Y.bat.mat, center = TRUE, scale = TRUE)
            Y.trt.scale = scale(Y.trt.mat, center = TRUE, scale = TRUE)
        }
        plsda_trt <-  externalPLSDA(X = tmp.cnts.scale, Y = Y.trt.scale,
                            ncomp = ncomp.trt, keepX = keepX.trt, tol = tol,
                            max.iter = max.iter)
        tmp.cnts.notrt <- plsda_trt$defl_data$X

        ### TRY OUT THE caret plsda function
        letest <- caret::plsda(x = tmp.cnts.scale, y = Y.bat.mat,
                               ncomp = ncomp.trt, scale = TRUE, tol = tol,
                               max.iter = max.iter,
                               near.zero.var = FALSE, logratio = "CLR", multilevel = NULL,
                               all.outputs = TRUE)


    }
    else {
        tmp.cnts.scale <- scale(tmp.cnts, center = TRUE, scale = TRUE)
        tmp.cnts.mean <- attributes(tmp.cnts.scale)$`scaled:center`
        tmp.cnts.sd <- attributes(tmp.cnts.scale)$`scaled:scale`
        Y.bat.scale = scale(Y.bat.mat, center = TRUE, scale = TRUE)
        plsda_trt = NULL
        tmp.cnts.notrt = tmp.cnts.scale
    }

    plsda_bat <- externalPLSDA(X = tmp.cnts.notrt, Y = Y.bat.scale, ncomp = ncomp.bat,
                       keepX = keepX.bat, tol = tol, max.iter = max.iter)
    bat_loadings <- plsda_bat$loadings$a
    tmp.cnts.temp <- tmp.cnts.scale
    for (h in seq_len(ncomp.bat)) {
        a.bat = bat_loadings[, h]
        t.bat = tmp.cnts.temp %*% a.bat
        tmp.cnts.temp <- mbecDeflate(tmp.cnts.temp, t.bat)
    }
    tmp.cnts.nobat <- tmp.cnts.temp
    tmp.cnts.nobat.final <- t(t(tmp.cnts.nobat) * tmp.cnts.sd + tmp.cnts.mean)
    tmp.cnts.nobat.final <- tmp.cnts.nobat.final/tmp.meta$weight
    if (length(model.vars) >= 2) {
        tmp.cnts.notrt.final <- t(t(tmp.cnts.notrt) *
                                      tmp.cnts.sd + tmp.cnts.mean)
        tmp.cnts.notrt.final <- tmp.cnts.notrt.final/tmp.meta$weight
    }
    else {
        tmp.cnts.notrt.final = NULL
    }

    return(tmp.cnts.nobat.final)
}




