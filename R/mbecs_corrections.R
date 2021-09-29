# CORRECTION FUNCTIONS ----------------------------------------------------

####################################################################################################
#' Batch Effect Correction
#'
#' Either corrects or accounts for (known) batch effects with one of several algorithms.
#' The methods 'lm, lmm, sva, ruv-2 and ruv-4" assess batch effects and return STH. I GUESS...
#' The remaining methods return the batch-corrected counts.
#'
#' Assessment Methods
#'
#' Linear (Mixed) Models:
#'
#' Surrogat3e variable Analysis (SVA):
#'
#' Remove unwanted Variation 2 (RUV-2):
#'
#' Remove Unwanted Variation 4 (RUV-4):
#'
#' Correction Methods
#' Remove Unwanted Variation 3 (RUV-3):
#'
#' Batch Mean Centering (BMC):
#'
#' Combat Batch Effects (ComBat):
#'
#' Remove Batch Effects (RBE):
#'
#' FAbatch: maybe just remove this - or keep in case of non-microbiome data
#'
#' Percentile Normalization (PN):
#'
#' Support Vector Decomposition (SVD):
#'
#' The function returns either a plot-frame or the finished ggplot object. Input for th data-set can
#' be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data. The
#' covariate table requires an 'sID' column that contains sample IDs equal to the sample naming in
#' the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Batch-Effect Correction Assessment
#' @param input.obj, mbecData object or numeric matrix (correct orientation is handled internally)
#' @param model.vars two covariates of interest to select by first variable selects panels and second one determines coloring
#' @param method, algorithm to use
#' @return an object of class MbecDataSet
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' This will return p-value for the linear model fit
#' \dontrun{val.score <- mbecLM(input.obj, model.vars=c("group","batch"), method="lm")}
mbecCorrection <- function(input.obj, model.vars=c("group","batch"),
                            method=c("lm","lmm","sva","ruv2","ruv4","ruv3","bmc","bat","rbe","fab","pn","svd"), update=TRUE) {

  ## ToDo:
  # - implement the whole 'sID' thingy
  # - adjust input-check
  model.vars <- c("group","batch")
  # - add logging functionality
  # - think about percentile normalization --> maybe the percentiles can be used to improve lm/lmm?!

  ## VALIDATE input and change to 'MbecData' if needed
  input.obj <- mbecProcessInput(input.obj, required.col=eval(model.vars))

  ## 00. Check if 'method' was chosen correctly.
  method <- match.arg(method, choices = c("lm","lmm","sva","ruv2","ruv4","ruv3","bmc","bat","rbe","fab","pn","svd"))

  if( method == "lm" ) {
    message("Applying Linear Model (LM) to account for batch-effects.")
    # fit a linear model for treatment and batch to every feature and return estimate of treatment significance for each
    # call helper function to handle this
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lm", model.vars=model.vars)

    ## RESULT is: tmp.group.p

  } else if( method == "lmm" ) {
    message("Applying Linear Mixed Model (LMM) to account for batch-effects.")
    # Use 'lmm' for unbalanced treatment x batch designs
    # fit a mixed-model with random batch-effect for each feature and extract adjusted p-values for treatment
    # call helper function to handle this
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lmm", model.vars=model.vars)

    ## RESULT is: tmp.group.p

  } else if( method == "sva" ) {
    message("Applying Surrogate Variable Analysis (SVA) to account for batch-effects.")

    ## needs fxs orientation
    tmp <- mbecGetData(input.obj, orientation="fxs")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # 'model.vars[1]' needs to be the experimental grouping; 'model.vars[2]' is the batch vector
    # first estimate the number of surrogate variables, i.e., number of latent factors to account for
    tmp.formula <- as.formula(paste(paste("x", model.vars[1], sep=" ~ "), paste(f.terms[], collapse=" + "), sep=" + "))

    tmp.mod <- stats::model.matrix( ~ tmp.meta[[model.vars[1]]])
    tmp.mod0 <- stats::model.matrix( ~ 1, data=tmp.meta[[model.vars[1]]])


    tmp.mod <- stats::model.matrix( ~ tmp.meta$group)
    tmp.mod0 <- stats::model.matrix( ~ 1, data = tmp.meta$group)
    tmp.sva.n <- sva::num.sv(dat = as.matrix(tmp.cnts), mod = tmp.mod)
    tmp.sva <- sva::sva(as.matrix(tmp.cnts), tmp.mod, tmp.mod0, n.sv = tmp.sva.n)

    message(paste("Number of significant surrogate variables is: ", tmp.sva$n.sv, sep=""))
    # do sth. that I have to understand yet
    tmp.mod.bat <- cbind(tmp.mod, tmp.sva$sv)
    tmp.mod0.bat <- cbind(tmp.mod0, tmp.sva$sv)
    tmp.sva.trt_p <- sva::f.pvalue(as.matrix(tmp.cnts), tmp.mod.bat, tmp.mod0.bat)
    tmp.sva.trt_adjp <- stats::p.adjust(tmp.sva.trt_p, method = 'fdr')

    # RESULT is: tmp.sva.trt_adjp

  } else if( method == "ruv2" ) {
    message("Applying Remove Unwanted Variantion v2 (RUV-II) for batch-correction.")
    ## needs sxf orientation
    tmp <- mbecGetData(input.obj, orientation="sxf")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # RUV-2 uses only negative controls
    ## ToDo: implement sth. to check if there are actually some negative controls available
    # without negative controls - fit lm to every feature and check whether treatment has a significant impact
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lm", model.vars=model.vars)

    ## everything that is not significntly impacted by the batch is considered a normal control
    tmp.nc <- tmp.group.p > 0.05

    tmp.ruv2 <- ruv::RUV2(Y = tmp.cnts, X = tmp.meta[[model.vars[1]]],
                          ctl = tmp.nc, k = 3) # k is subjective

    ## testing things
    # fit = variance_adjust(tmp.ruv2)
    tmp.ruv2.trt_p <- tmp.ruv2$p
    tmp.ruv2.trt_adjp <- p.adjust(tmp.ruv2.trt_p, method="fdr")

    # RESULT is: tmp.ruv2.trt_adjp

  } else if( method == "ruv4" ) {
    message("Applying Remove Unwanted Variantion v4 (RUV-IV) for batch-correction.")
    ## needs sxf orientation
    tmp <- mbecGetData(input.obj, orientation="sxf")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # RUV-4 requires negative controls
    ## ToDo: implement sth. to check if there are actually some negative controls available
    # without negative controls - fit lm to every feature and check whether treatment has a significant impact
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lm", model.vars=model.vars)

    ## everything that is not significntly impacted by the batch is considered a negative control
    tmp.nc <- tmp.group.p > 0.05

    # determine appropriate k-value
    tmp.k <- ruv::getK(Y = tmp.cnts, X = tmp.meta$group, ctl = tmp.nc)
    # set 'k' to 1 if algorithm chose 0 as value
    tmp.k <- ifelse(tmp.k$k !=0, tmp.k$k, 1)

    tmp.ruv4 <- ruv::RUV4(Y = tmp.cnts, X = tmp.meta$group, ctl = tmp.nc, k = tmp.k)
    tmp.ruv4.trt_p <- tmp.ruv4$p
    tmp.ruv4.trt_adjp <- p.adjust(tmp.ruv4.trt_p, method="fdr")

    ## RESULT is: tmp.ruv4.trt_adjp

  } else if( method == "ruv3" ) {
    message("Applying Remove Unwanted Variantion v3 (RUV-III) for batch-correction.")
    ## start with RUV-3 - requires negative controls and technical sample replicates
    # - get negative controls by fitting lm to every feature and selecting those that are non-significant for batch
    # - check if technical replicates are available, requires 'replicate' column

    ## check and prepare inputs
    ## needs sxf orientation
    tmp <- mbecGetData(input.obj, orientation="sxf")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # check if technical replicates are available - use 'replicate' column and look for replicate values
    tmp.replicate <- tmp.meta$replicate
    if( length(tmp.replicate) == length(unique(tmp.replicate)) ) {
      message("No technical replicates found. RUV-3 is not available for this data-set!")
      stop()
    } else {
      ## ToDo: should probably check if further requirements for replicates are satisfied
      # construct replicate matrix
      tmp.replicate.matrix <- ruv::replicate.matrix(tmp.replicate)
    }

    ## ToDo: implement sth. to check if there are actually some negative controls available
    # without negative controls - fit lm to every feature and check whether treatment has a significant impact
    tmp.group.p <- mbecLM(input.obj = input.obj, method = "lm")

    ## everything that is not significntly impacted by the batch is considered a normal control
    tmp.nc <- tmp.group.p > 0.05

    # all required inputs available --> apply RUV-3
    corrected.cnts <- ruv::RUVIII(Y = tmp.cnts, M = tmp.replicate.matrix, ctl = tmp.nc)
    rownames(corrected.cnts) <- rownames(tmp.cnts)

  } else if( method == "bmc" ) {
    message("Applying Batch Mean-Centering (BMC) for batch-correction.")
    ## check and prepare inputs
    ## needs sxf orientation
    tmp <- mbecGetData(input.obj, orientation="sxf")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # get unique batches
    input.batches <- unique(tmp.meta[[model.vars[2]]])

    corrected.cnts <- NULL
    for( batch.idx in input.batches ) {
      tmp <- scale(tmp.cnts[tmp.meta$sample[tmp.meta[[model.vars[2]]] %in% batch.idx], ], center = T, scale = F)
      corrected.cnts <- rbind.data.frame(corrected.cnts, tmp)
    }

  } else if( method == "bat" ) {
    message("Applying ComBat (sva) for batch-correction.")
    ## ComBat requires 'fxs' orientation for inputs
    tmp <- mbecGetData(input.obj, orientation="fxs")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # ToDo: tmp model
    tmp.mod <- model.matrix( ~ tmp.meta[[model.vars[1]]]) # full model

    corrected.cnts <- sva::ComBat(tmp.cnts, batch = tmp.meta[[model.vars[2]]],
                                  mod = tmp.mod, par.prior = F, prior.plots = F)

  } else if( method == "rbe" ) {
    message("Applying 'removeBatchEffect' (limma) for batch-correction.")
    ## removeBatchEffect requires 'fxs' orientation for inputs
    tmp <- mbecGetData(input.obj, orientation="fxs")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # ToDo: tmp model
    tmp.mod <- model.matrix( ~ tmp.meta[[model.vars[1]]]) # full model

    corrected.cnts <- limma::removeBatchEffect(tmp.cnts, batch = tmp.meta[[model.vars[2]]],
                                               design = tmp.mod)

  } else if( method == "fab" ) {
    message("Applying FAbatch for batch-correction.")
    ## FAbatch requires 'sxf' orientation for inputs
    tmp <- mbecGetData(input.obj, orientation="sxf")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

    # covariates have to be 'numeric' factors ...
    corrected.cnts <- bapred::fabatch(x=as.matrix(tmp.cnts),
                                      y=as.factor(as.numeric(tmp.meta[[model.vars[1]]])),
                                      batch = as.factor(as.numeric(tmp.meta[[model.vars[2]]])))

  } else if( method == "pn" ) {
    message("Applying Percentile Normalization (PN) for batch-correction.")
    ## Percentile Normalization requires 'sxf' orientation for inputs
    tmp <- mbecGetData(input.obj, orientation="sxf")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]
    # check for case/control design
    if( nlevels(tmp.meta[,eval(model.vars[1])]) != 2 ) {
      warning("Grouping/Treatment contains ", nlevels(tmp.meta[,eval(model.vars[1])]), " classes. Percentile normalization is designed to work with 2 classes only, i.e., case/control studies.")
    }
    # check/adjust for zero-values
    if( any(tmp.cnts == 0) ) {
      warning("Abundances contain zero values. Adding small uniform offset.")
      tmp.cnts <- apply(tmp.cnts, c(1,2), function(x) ifelse(x == 0, runif(1,0,10^-9),x))
    }

    # do normalisation
    corrected.cnts <- percentileNorm(tmp.cnts, tmp.meta[,eval(model.vars)])

  } else if( method == "svd" ) {
    message("Applying Singular Value Decomposition (SVD) for batch-correction.")
    ## SVD requires 'sxf' orientation for inputs
    tmp <- mbecGetData(input.obj, orientation="sxf")
    tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]
    # sd, mean, scale
    tmp.sd <- apply(tmp.cnts, 2, sd)
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

    for( c.idx in 1:dim(corrected.cnts)[2] ) {
      for( r.idx in 1:dim(corrected.cnts)[1] ) {
        corrected.cnts[r.idx, c.idx] <- tmp.svd.defl[r.idx, c.idx] * tmp.sd[c.idx] + tmp.mean[c.idx]
      }
    }
  }

  ## Prepare return-object
  message("Calling 'mbecSetData' - log is: ", attr(input.obj, "log"))
  # If 'update' is TRUE the counts will be updated - otherwise the corrected counts will be placed in a list.
  return.obj <- mbecSetData(input.obj = input.obj, new.cnts = corrected.cnts,
                            log=method, type=method, update=update)

  return(return.obj)
  # Point of no return. *badum tsss*
}



