# CORRECTION FUNCTIONS ----------------------------------------------------

#' Run Correction Pipeline
#'
#' Run all corerction algorithms selected by method and add corrected counts
#' as matrices to the data-set.
#'
#' @keywords Batch-Effect Correction Pipeline
#' @param input.obj phyloseq object or numeric matrix (correct orientation is handeled internally)
#' @param model.vars two covariates of interest to select by first variable selects panels and second one determines coloring
#' @param method algorithms to use
#' @param nc.features (OPTIONAL) A vector of features names to be used as negative controls in RUV-3. If not supplied, the algorithm will use an 'lm' to find pseudo-negative controls
#' @return an object of class MbecDataSet
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This call will use 'ComBat' for batch effect correction and store the new
#' # counts in a list-obj in the output.
#' study.obj <- mbecRunCorrections(input.obj=datadummy,
#' model.vars=c("group","batch"), method=c("bat","bmc"))
#'
#' # This call will use 'Percentile Normalization' for batch effect correction
#' # and replace the old count matrix.
#' study.obj <- mbecRunCorrections(datadummy, model.vars=c("group","batch"),
#' method=c("pn"))
mbecRunCorrections <- function(input.obj, model.vars=c("group","batch"),
                           method=c("ruv3","bmc","bat","rbe","fab","pn","svd"),
                           nc.features=NULL) {

  input.obj <- mbecProcessInput(input.obj, required.col=eval(model.vars))

  for( m.idx in method ) {
    input.obj <- mbecCorrection(input.obj = input.obj, method = eval(m.idx), update = FALSE)
  }

  return(input.obj)

}

#' Batch Effect Correction
#'
#' Either corrects or accounts for (known) batch effects with one of several algorithms.
#' The methods 'lm, lmm, sva, ruv-2 and ruv-4" assess batch effects and return STH. I GUESS...
#' The remaining methods return the batch-corrected counts.
#'
#' ASSESSMENT METHODS
#'
#' Linear (Mixed) Models: A simple linear mixed model with covariates 'treatment' an 'batch', or
#' respective variables in your particular data-set, will be fitted to each feature and the
#' significance for the treatment variable extracted.
#'
#' Surrogate variable Analysis (SVA): Surrogate Variable Analysis (SVA): Two step approach that
#' (1.) identify the number of latent factors to be estimated by fitting a full-model with effect
#' of interest and a null-model with no effects. The function num.sv then calculates the number of
#' latent factors. In the next (2.) step, the sva function will estimate the surrogate variables.
#' And adjust for them in full/null-model . Subsequent F-test gives significance values for each
#' feature - these P-values and Q-values are accounting for surrogate variables (estimated BEs).
#'
#' Remove unwanted Variation 2 (RUV-2): Estimates unknown BEs by using negative control variables
#' that, in principle, are unaffected by treatment/study/biological effect (aka the effect of
#' interest in an experiment). These variables are generally determined prior to the experiment.
#' An approach to RUV-2 without the presence of negative control variables is the estimation of
#' pseudo-negative controls. To that end an lm or lmm (depending on whether or not the study design
#' is balanced) with treatment is fitted to each feature and the significance calculated.
#' The features that are not significantly affected by treatment are considered as pseudo-negative
#' control variables. Subsequently, the actual RUV-2 function is applied to the data and returns
#' the p-values for treatment, considering unwanted BEs (whatever that means).
#'
#' Remove Unwanted Variation 4 (RUV-4): The updated version of RUV-2 also incorporates the residual
#' matrix (w/o treatment effect) to estimate the unknown BEs. To that end it follows the same
#' procedure in case there are no negative control variables and computes pseudo-controls from the
#' data via l(m)m. As RUV-2, this algorithm also uses the parameter 'k' for the number of latent
#' factors. RUV-4 brings the function 'getK()' that estimates this factor from the data itself.
#' The calculated values are however not always reliable. A value of k=0 fo example can occur and
#' should be set to 1 instead. The output is the same as with RUV-2.
#'
#' CORRECTION METHODS
#'
#' Remove Unwanted Variation 3 (RUV-3): This algorithm requires negative control-features, i.e.,
#' OTUs that are known to be unaffected by the batch effect, as well as technical replicates.
#' The algorithm will check for the existence of a replicate column in the covariate data. If the
#' column is not present, the execution stops and a warning message will be displayed.
#'
#' Batch Mean Centering (BMC): For known BEs, this method takes the batches, i.e., subgroup of
#' samples within a particular batch, and centers them to their mean.
#'
#' Combat Batch Effects (ComBat): This method uses an non-/parametric empirical  Bayes framework to
#' correct for BEs. Described by Johnson et al. 2007 this method was initially conceived to work
#' with gene expression data and is part of the sva-package in R.
#'
#' Remove Batch Effects (RBE): As part of the limma-package this method was designed to remove BEs
#' from Microarray Data. The algorithm fits the full-model to the data, i.e., all relevant
#' covariates whose effect should not be removed, and a model that only contains the known BEs.
#' The difference between these models produces a residual matrix that (should) contain only the
#' full-model-effect, e.g., treatment. As of now the mbecs-correction only uses the first input
#' for batch-effect grouping. ToDo: think about implementing a version for more complex models.
#'
#' FAbatch: Implemented in the bapred-package this method is based on the approach described in
#' Hornung et al. 2017. "It is a combination of two commonly used approaches: location-and-scale
#' adjustment and data cleaning by adjustment for distortions due to latent factors."
#' HOWEVER, since it can't handle the zero inflated counts of compositional microbiome-data -->
#' maybe just remove this - or keep in case of non-microbiome data?!
#'
#' Percentile Normalization (PN): this method was actually developed specifically to facilitate
#' the integration of microbiome data from different studies/experimental set-ups. This problem is
#' similar to the mitigation of BEs, i.e., when collectively analyzing two or more data-sets,
#' every study is effectively a batch on its own (not withstanding the probable BEs within studies).
#' The algorithm iterates over the unique batches and converts the relative abundance of control
#' samples into their percentiles. The relative abundance of case-samples within the respective
#' batches is then transformed into percentiles of the associated control-distribution. Basically,
#' the procedure assumes that the control-group is unaffected by any effect of interest, e.g.,
#' treatment or sickness, but both groups within a batch are affected by that BE. The switch to
#' percentiles (kinda) flattens the effective difference in count values due to batch - as compared
#' to the other batches. This also introduces the two limiting aspects in percentile normalization.
#' It can only be applied to case/control designs because it requires a reference group.
#' In addition, the transformation into percentiles removes information from the data.
#'
#' Singular Value Decomposition (SVD): Basically perform matrix factorization and compute singular
#' eigenvectors (SEV). Assume that the first SEV captures the batch-effect and remove this effect
#' from the data. The interesting thing is that this works pretty well (with the test-data anyway)
#' But since the SEVs are latent factors that are (most likely) confounded with other effects
#' it is not obvious to me that this is the optimal approach to solve this issue.
#' ToDo: IF I find the time to works on "my-own-approach" then this is the point to start from!!!
#'
#' The function returns an MbecData-object with either updated counts or a list that contains
#' matrices for all applied BE-correction methods - this facilitates the comparison of multiple
#' BE-correction approaches to determine the one that is best suited for the data-set. Input can
#' be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data. The
#' covariate table requires an 'sID' column that contains sample IDs equal to the sample naming in
#' the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is handeled internally)
#' @param model.vars two covariates of interest to select by first variable selects panels and second one determines coloring
#' @param method algorithm to use
#' @param update either update the input counts or create a new matrix of corrected counts within the input
#' @param nc.features (OPTIONAL) A vector of features names to be used as negative controls in RUV-3. If not supplied, the algorithm will use an 'lm' to find pseudo-negative controls
#' @return an object of class MbecDataSet
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This call will use 'ComBat' for batch effect correction and store the new
#' # counts in a list-obj in the output.
#' study.obj <- mbecCorrection(input.obj=datadummy,
#' model.vars=c("group","batch"), method="bat", update=FALSE)
#'
#' # This call will use 'Percentile Normalization' for batch effect correction
#' # and replace the old count matrix.
#' study.obj <- mbecCorrection(datadummy, model.vars=c("group","batch"),
#' method="pn", update=TRUE)
mbecCorrection <- function(input.obj, model.vars=c("group","batch"),
                           method=c("lm","lmm","sva","ruv2","ruv4","ruv3","bmc",
                                    "bat","rbe","fab","pn","svd"),
                           update=TRUE, nc.features=NULL) {

  input.obj <- mbecProcessInput(input.obj, required.col=eval(model.vars))

  ## 00. Check if 'method' was chosen correctly.
  method <- match.arg(method, choices = c("lm","lmm","sva","ruv2","ruv4","ruv3",
                                          "bmc","bat","rbe","fab","pn","svd"))
  # START WITH ASSESSMENT METHODS
  if( method == "lm" ) {
    message("Applying Linear Model (LM) to account for batch-effects.")
    res.assess <- mbecLM(input.obj = input.obj, method = "lm",
                          model.vars=model.vars)
  } else if( method == "lmm" ) {
    message("Applying Linear Mixed Model (LMM) to account for batch-effects.")
    res.assess <- mbecLM(input.obj = input.obj, method = "lmm",
                          model.vars=model.vars)
  } else if( method == "sva" ) {
    res.assess <- mbecSVA(input.obj, model.vars)
  } else if( method == "ruv2" ) {
    res.assess <- mbecRUV2(input.obj, model.vars, nc.features)
  } else if( method == "ruv4" ) {
    res.assess <- mbecRUV4(input.obj, model.vars, nc.features)
    ## CORRECTION METHODS FROM HERE
  } else if( method == "ruv3" ) {
    res.correct <- mbecRUV3(input.obj, model.vars, nc.features)
  } else if( method == "bmc" ) {
    res.correct <- mbecBMC(input.obj, model.vars)
  } else if( method == "bat" ) {
    res.correct <- mbecBat(input.obj, model.vars)
  } else if( method == "rbe" ) {
    res.correct <- mbecRBE(input.obj, model.vars)
  } else if( method == "fab" ) {
    res.correct <- mbecFAB(input.obj, model.vars)
  } else if( method == "pn" ) {
    res.correct <- mbecPN(input.obj, model.vars)
  } else if( method == "svd" ) {
    res.correct <- mbecSVD(input.obj, model.vars)
  }
  ## Prepare return-object
  message("Calling 'mbecSetData' - log is: ", attr(input.obj, "log"))
  # If 'update' is TRUE the counts will be updated - otherwise the corrected counts will be placed in a list.
  return.obj <- mbecSetData(input.obj = input.obj, new.cnts = res.correct,
                            log=method, type=method, update=update)
  return(return.obj)
  # Point of no return. *badum tsss*
}


#' Surrogate variable Analysis (SVA)
#'
#' Two step approach that (1.) identify the number of latent factors to be
#' estimated by fitting a full-model with effect of interest and a null-model
#' with no effects. The function num.sv then calculates the number of latent
#' factors. In the next (2.) step, the sva function will estimate the surrogate
#' variables. And adjust for them in full/null-model . Subsequent F-test gives
#' significance values for each feature - these P-values and Q-values are
#' accounting for surrogate variables (estimated BEs).
#'
#' The output will be SOMETHING
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @return SOMETHING
#' @include mbecs_classes.R
mbecSVA <- function(input.obj, model.vars) {

  message("Applying Surrogate Variable Analysis (SVA) to account for batch-effects.")
  tmp <- mbecGetData(input.obj, orientation="fxs")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  tmp.mod <- stats::model.matrix( ~ tmp.meta[[model.vars[1]]])
  tmp.mod0 <- stats::model.matrix( ~ 1, data=tmp.meta[[model.vars[1]]])
  tmp.sva.n <- sva::num.sv(dat = as.matrix(tmp.cnts), mod = tmp.mod)
  tmp.sva <- sva::sva(as.matrix(tmp.cnts), tmp.mod, tmp.mod0, n.sv = tmp.sva.n)

  message("Number of significant surrogate variables is: ", tmp.sva$n.sv)
  # do sth. that I have to understand yet
  tmp.mod.bat <- cbind(tmp.mod, tmp.sva$sv)
  tmp.mod0.bat <- cbind(tmp.mod0, tmp.sva$sv)
  tmp.sva.trt_p <- sva::f.pvalue(as.matrix(tmp.cnts), tmp.mod.bat, tmp.mod0.bat)
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
#' The output will be SOMETHING
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-3. If not supplied, the algorithm will use an 'lm'
#' to find pseudo-negative controls
#' @return SOMETHING
#' @include mbecs_classes.R
mbecRUV2 <- function(input.obj, model.vars, nc.features=NULL) {
  tmp <- mbecGetData(input.obj, orientation="sxf")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  ## check for neg. controls or generate pseudo negatives
  if( !is.null(nc.features) ) {
    tmp.nc <- colnames(tmp.cnts) %in% nc.features
    names(tmp.nc) <- colnames(tmp.cnts)
  } else {
    message("No negative control features provided. Using pseudo-negative controls.")
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lm", model.vars=model.vars)
    tmp.nc <- tmp.group.p > 0.05
  }

  message("Applying Remove Unwanted Variantion v2 (RUV-II) for batch-correction.")
  tmp.ruv2 <- ruv::RUV2(Y = tmp.cnts, X = tmp.meta[[model.vars[1]]],
                        ctl = tmp.nc, k = 3) # k is subjective

  tmp.ruv2.trt_p <- tmp.ruv2$p
  tmp.ruv2.trt_adjp <- stats::p.adjust(tmp.ruv2.trt_p, method="fdr")

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
#' The output will be SOMETHING
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-3. If not supplied, the algorithm will use an 'lm'
#' to find pseudo-negative controls
#' @return SOMETHING
#' @include mbecs_classes.R
mbecRUV4 <- function(input.obj, model.vars, nc.features=NULL) {
  tmp <- mbecGetData(input.obj, orientation="sxf")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  ## check for neg. controls or generate pseudo negatives
  if( !is.null(nc.features) ) {
    tmp.nc <- colnames(tmp.cnts) %in% nc.features
    names(tmp.nc) <- colnames(tmp.cnts)
  } else {
    message("No negative control features provided. Using pseudo-negative controls.")
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lm", model.vars=model.vars)
    tmp.nc <- tmp.group.p > 0.05
  }

  message("Applying Remove Unwanted Variantion v4 (RUV-IV) for batch-correction.")
  tmp.k <- ruv::getK(Y = tmp.cnts, X = tmp.meta$group, ctl = tmp.nc)
  tmp.k <- ifelse(tmp.k$k !=0, tmp.k$k, 1)

  tmp.ruv4 <- ruv::RUV4(Y = tmp.cnts, X = tmp.meta$group, ctl = tmp.nc, k = tmp.k)
  tmp.ruv4.trt_p <- tmp.ruv4$p
  tmp.ruv4.trt_adjp <- stats::p.adjust(tmp.ruv4.trt_p, method="fdr")

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
#' The function returns a matrix of corrected counts Correct orientation of
#' counts will be handled internally.
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @param nc.features (OPTIONAL) A vector of features names to be used as
#' negative controls in RUV-3. If not supplied, the algorithm will use an 'lm'
#' to find pseudo-negative controls
#' @return A matrix of batch-effect corrected counts
#' @include mbecs_classes.R
mbecRUV3 <- function(input.obj, model.vars, nc.features=NULL) {
  tmp <- mbecGetData(input.obj, orientation="sxf")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  ## check for neg. controls or generate pseudo negatives
  if( !is.null(nc.features) ) {
    tmp.nc <- colnames(tmp.cnts) %in% nc.features
    names(tmp.nc) <- colnames(tmp.cnts)
  } else {
    message("No negative control features provided. Using pseudo-negative controls.")
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lm", model.vars=model.vars)
    tmp.nc <- tmp.group.p > 0.05
  }

  message("Applying Remove Unwanted Variantion v3 (RUV-III) for batch-correction.")
  tmp.replicate <- tmp.meta$replicate
  if( length(tmp.replicate) == length(unique(tmp.replicate)) ) {
    stop("No technical replicates found. RUV-3 is not available for this data-set!")
  } else {
    tmp.replicate.matrix <- ruv::replicate.matrix(tmp.replicate)
  }

  corrected.cnts <- ruv::RUVIII(Y = tmp.cnts, M = tmp.replicate.matrix, ctl = tmp.nc)
  rownames(corrected.cnts) <- rownames(tmp.cnts)

  return(corrected.cnts)
}


#' Batch Mean Centering (BMC)
#'
#' For known BEs, this method takes the batches, i.e., subgroup of samples
#' within a particular batch, and centers them to their mean.
#'
#' The function returns a matrix of corrected counts Correct orientation of
#' counts will be handled internally.
#'
#' @keywords BECA Batch Mean Centering
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @return A matrix of batch-effect corrected counts
#' @include mbecs_classes.R
mbecBMC <- function(input.obj, model.vars) {

  message("Applying Batch Mean-Centering (BMC) for batch-correction.")
  tmp <- mbecGetData(input.obj, orientation="sxf")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  # get unique batches
  input.batches <- unique(tmp.meta[[model.vars[1]]])

  corrected.cnts <- NULL
  for( batch.idx in input.batches ) {
    tmp <- scale(tmp.cnts[tmp.meta$sample[tmp.meta[[model.vars[1]]] %in% batch.idx], ], center = TRUE, scale = FALSE)
    corrected.cnts <- rbind.data.frame(corrected.cnts, tmp)
  }

  return(corrected.cnts)
}


#' Combat Batch Effects (ComBat)
#'
#' This method uses an non-/parametric empirical  Bayes framework to correct
#' for BEs. Described by Johnson et al. 2007 this method was initially conceived
#' to work with gene expression data and is part of the sva-package in R.
#'
#' The function returns a matrix of corrected counts Correct orientation of
#' counts will be handled internally.
#'
#' @keywords BECA Batch Mean Centering
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars First variable determines Batch, second (optional) variable
#' describes study-group. ToDo: more vars and own models
#' @return A matrix of batch-effect corrected counts
#' @include mbecs_classes.R
mbecBat <- function(input.obj, model.vars) {

  message("Applying ComBat (sva) for batch-correction.")
  ## ComBat requires 'fxs' orientation for inputs
  tmp <- mbecGetData(input.obj, orientation="fxs")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  # ToDo: tmp model
  if( length(model.vars) == 1 ) {
    corrected.cnts <- sva::ComBat(tmp.cnts, batch = tmp.meta[[model.vars[1]]],
                                  mod = 1, par.prior = FALSE,
                                  prior.plots = FALSE)
  } else {
    tmp.mod <- stats::model.matrix( ~ tmp.meta[[model.vars[2]]]) # full model

    corrected.cnts <- sva::ComBat(tmp.cnts, batch = tmp.meta[[model.vars[1]]],
                                  mod = tmp.mod, par.prior = FALSE,
                                  prior.plots = FALSE)
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
#' The function returns a matrix of corrected counts Correct orientation of
#' counts will be handled internally.
#'
#' @keywords BECA Limma Remove Batch Effects
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars First variable determines Batch, second (optional) variable
#' describes study-group. ToDo: more vars and own models
#' @return A matrix of batch-effect corrected counts
#' @include mbecs_classes.R
mbecRBE <- function(input.obj, model.vars) {

  message("Applying 'removeBatchEffect' (limma) for batch-correction.")
  ## removeBatchEffect requires 'fxs' orientation for inputs
  tmp <- mbecGetData(input.obj, orientation="fxs")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  if( length(model.vars) == 1 ) {
    corrected.cnts <- limma::removeBatchEffect(
      tmp.cnts,
      batch = tmp.meta[[model.vars[1]]],design = matrix(1,ncol(tmp.cnts),1))
  } else {
    tmp.mod <- stats::model.matrix( ~ tmp.meta[[model.vars[2]]]) # full model

    corrected.cnts <- limma::removeBatchEffect(
      tmp.cnts, batch = tmp.meta[[model.vars[1]]],design = tmp.mod)
  }
  return(corrected.cnts)
}


#' FAbatch
#'
#' Implemented in the bapred-package this method is based on the approach
#' described in Hornung et al. 2017. "It is a combination of two commonly used
#' approaches: location-and-scale adjustment and data cleaning by adjustment
#' for distortions due to latent factors." HOWEVER, since it can't handle the
#' zero inflated counts of compositional microbiome-data --> maybe just remove
#' this - or keep in case of non-microbiome data?!
#'
#' The function returns a matrix of corrected counts Correct orientation of
#' counts will be handled internally.
#'
#' @keywords BECA Hornung FABatch
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @return A matrix of batch-effect corrected counts
#' @include mbecs_classes.R
mbecFAB <- function(input.obj, model.vars) {

  message("Applying FAbatch for batch-correction.")
  ## FAbatch requires 'sxf' orientation for inputs
  tmp <- mbecGetData(input.obj, orientation="sxf")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  # covariates have to be 'numeric' factors ...
  corrected.cnts <- bapred::fabatch(x=as.matrix(tmp.cnts),
                                    y=as.factor(as.numeric(tmp.meta[[model.vars[1]]])),
                                    batch = as.factor(as.numeric(tmp.meta[[model.vars[2]]])))
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
#' The function returns a matrix of corrected counts Correct orientation of
#' counts will be handled internally.
#'
#' @keywords BECA Duvallet Percentile Normalisation
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @return A matrix of batch-effect corrected counts
#' @include mbecs_classes.R
mbecPN <- function(input.obj, model.vars) {

  message("Applying Percentile Normalization (PN) for batch-correction.")
  ## Percentile Normalization requires 'sxf' orientation for inputs
  tmp <- mbecGetData(input.obj, orientation="sxf")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]
  # check for case/control design
  if( nlevels(tmp.meta[,eval(model.vars[2])]) != 2 ) {
    warning("Grouping/Treatment contains ",
            nlevels(tmp.meta[,eval(model.vars[2])]),
            " different categories. Percentile normalization is designed to work
            with 2 classes only, i.e., case/control studies.")
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


#'Singular Value Decomposition (SVD)
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
#' The function returns a matrix of corrected counts Correct orientation of
#' counts will be handled internally.
#'
#' @keywords BECA Duvallet Percentile Normalisation
#' @param input.obj phyloseq object or numeric matrix (correct orientation is
#' handeled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @return A matrix of batch-effect corrected counts
#' @include mbecs_classes.R
mbecSVD <- function(input.obj, model.vars) {

  message("Applying Singular Value Decomposition (SVD) for batch-correction.")
  ## SVD requires 'sxf' orientation for inputs
  tmp <- mbecGetData(input.obj, orientation="sxf")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]
  # sd, mean, scale
  tmp.sd <- apply(tmp.cnts, 2, stats::sd)
  tmp.mean <- apply(tmp.cnts, 2, mean)
  # center and scale
  tmp.cnts.scale <- scale(tmp.cnts, center = TRUE, scale = TRUE)
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








