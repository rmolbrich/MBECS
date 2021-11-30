# HELPER FUNCTIONS --------------------------------------------------------


#' Capitalize Word Beginning
#'
#' Change the first letter of the input to uppercase. Used in plotting functions
#' to make covariates, i.e., axis-lables look nicer.
#'
#' @keywords uppercase
#' @param input A string that needs to be Capitalized
#' @return Input with first letter capitalized
mbecUpperCase <- function(input=character()) {
  return(gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",input,perl = TRUE))
}


#' Linear (Mixed) Model Feature to Batch Fit
#'
#' Helper function that fits lm/lmm with covariates 'treatment' and 'batch' to every feature in the
#' data-set. Returns the fdr corrected significance value for the "treatment" variable. The method
#' 'lm' will fit the linear model\code{y ~ model.vars[1] + model.vars[2]} and the linear mixed model
#' will consider the second term as random effect, i.e., \code{y ~ model.vars[1] + (1|model.vars[2])}.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input for th data-set can
#' be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data. The
#' covariate table requires an 'sID' column that contains sample IDs equal to the sample naming in
#' the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Significance Linear Mixed Model Batch
#' @param input.obj, mbecData object or numeric matrix (correct orientation is handeled internally)
#' @param model.vars two covariates of interest to select by first variable selects panels and second one determines coloring
#' @param method, either 'lm' or 'lmm' for linear (mixed) models or 'auto' to detect correct method (not implemented yet)
#' @return vector of fdr corrected p-values that show significance of treatment for every feature
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return p-value for the linear model fit of every feature.
#' val.score <- mbecLM(input.obj=datadummy, model.vars=c("group","batch"),
#' method="lm")
mbecLM <- function(input.obj, method=c("lm","lmm"), model.vars=c("group","batch")) {
  # ToDo: standard model is '~group+batch' but maybe an alternative mode is nice
  #       alternative correction methods
  #       auto mode selection procedure --> detect unbalanced design?!

  ## check and prepare inputs
  tmp <- mbecGetData(input.obj, orientation="sxf", required.col=eval(model.vars))
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  if( method == "auto" ) {
    # ToDo: sth.
    message("This feature is supposed to detect unbalanced designs and select lmm instead of lm. Alas, it has not been implemented and doesn't do jack-shit so far.")

  } else if( method == "lm" ) {
    # fit lm to every feature with treatment and batch as as model parameters
    # then extract p-value for treatment
    tmp.group.p <- apply(tmp.cnts, 2, FUN = function(y){
      nc.lm <- stats::lm(y ~ get(model.vars[1]) + get(model.vars[2]), data = tmp.meta)
      nc.lm.summary <- summary(nc.lm)
      # extract p-value of group (treatment)
      p <- nc.lm.summary$coefficients[2,4]
    })

  } else if( method == "lmm" ) {
    # overkill, but just keep it for now
    f.terms <- paste("(1|",model.vars[2],")", sep="")

    tmp.group.p <- apply(tmp.cnts, 2, FUN = function(x) {

      tmp.formula <- stats::as.formula(paste(paste("x", model.vars[1], sep=" ~ "), paste(f.terms[], collapse=" + "), sep=" + "))
      nc.lmm <- eval(bquote(lmerTest::lmer(.(tmp.formula), data = tmp.meta)))
      nc.lmm.summary <- summary(nc.lmm)
      p <- nc.lmm.summary$coefficients[2,5]

    })
  }

  # correct for multiple testing
  tmp.group.p <- stats::p.adjust(tmp.group.p, method = 'fdr')

  return(tmp.group.p)
}


# TRANSFORMATION FUNCTIONS ------------------------------------------------


#' Log-Ratio Transformation
#'
#' Wrapper to help perform log-ratio transformations ,adapted from packages 'mixOmics' and
#' 'robCompositions' to work on matrices and Phyloseq objects alike.
#'
#' The function returns an MbecData object with tranformed counts and covariate information. Input for the data-set
#' can be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data.
#' The covariate table requires an 'sID' column that contains sample IDs equal to the sample naming
#' in the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Log Ratio Transformation
#' @param input.obj either pyhloseq-object (OTU orientation is handled) or numeric matrix (samples x features)
#' @param method one of 'CLR' or 'ILR'
#' @param offset optional offset in case of sparse matrix
#' @param required.col (OPTIONAL) A vector of column names in the meta-data that need to be present. Sanity check for subsequent steps.
#' @return MbecDataObject with transformed counts
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return the cumulative log-ratio transformed counts in an MbecData object
#' \dontrun{mbec.LRT <- LRTransform(input.obj=list(counts, covariates),
#' method="CLR", offset=0)}
#'
#' # This will return the inverse log-ratio transformed counts in an MbecData object
#' \dontrun{mbec.LRT <- LRTransform(input.obj=list(counts, covariates),
#' method="ILR", offset=0)}
LRTransform <- function(input.obj, method = c("none", "CLR", "ILR"), offset = 0, required.col=NULL) {
  ## 00. Check if 'method' was chosen correctly and get optional arguments
  method <- match.arg(method)

  ## VALIDATE input and change to 'MbecData' if needed
  if( !is.null(required.col) ) {
    input.obj <- mbecProcessInput(input.obj, required.col=eval(required.col))
  } else {
    input.obj <- mbecProcessInput(input.obj)
  }

  ## needs sxf orientation
  tmp <- mbecGetData(input.obj, orientation="sxf")

  if (method == "ILR") {
    if (!is(input.obj, "ilr")) {
      tmp.cnts = ilr.transfo(tmp[[1]], offset = offset)
    }

  } else if (method == "CLR") {

    tmp.cnts <- clr.transfo(tmp[[1]], offset = offset)
  }

  # rebuild sample AND feature names for reassembly
  colnames(tmp.cnts) <- colnames(tmp[[1]])
  rownames(tmp.cnts) <- rownames(tmp[[1]])

  phyloseq::otu_table(input.obj) = phyloseq::otu_table(tmp.cnts, taxa_are_rows = FALSE)

  return(input.obj)
}


#' Percentile Normalization
#'
#' Wrapper to help perform percentile normalization on a matrix of counts. Takes counts and a
#' data-frame of grouping variables and returns a matrix of transformed counts. This is designed
#' (by the Developers of the procedure) to work with case/control experiments by taking the
#' untreated group as reference and adjusting the other groupings of TRT x Batch to it.
#'
#' The function returns an MbecData object with tranformed counts and covariate information. Input for the data-set
#' can be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data.
#' The covariate table requires an 'sID' column that contains sample IDs equal to the sample naming
#' in the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Log Ratio Transformation
#' @param cnts numeric matrix (samples x features)
#' @param meta data-frame of covariate columns, first column contains study groups, second column contains batches
#' @return numeric matrix of corrected counts
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return a matrix of normalised counts, according to the covariate
#' # information in meta
#' mtx.pn_counts <- percentileNorm(cnts=datadummy$cnts,
#' meta=datadummy$meta[,c("group","batch")])
percentileNorm <- function(cnts, meta) {

  ref.group <- levels(meta[,1])[1]
  message("Group ",ref.group, " is considered control group, i.e., reference for normalization procedure. To change reference please 'relevel()' grouping factor accordingly.")

  norm.cnts <- cnts; norm.cnts[,] <- NA

  # for every batch
  for( b.idx in levels(meta[,2]) ) {
    # for every feature
    for( f.idx in seq_len(ncol(cnts)) ) {
      # which are the control-group values
      ctrl.group.vec <- cnts[which((meta[,1] %in% ref.group) & (meta[,2] %in% b.idx)), f.idx]
      # for every sample in the batch
      for( s.idx in which(meta[,2] %in% b.idx) ) {
        # call 'poscore' and get normalized value
        norm.cnts[s.idx, f.idx] <- poscore(ctrl.group.vec, cnts[s.idx, f.idx], "mean")

      }
    }
  }
  return(norm.cnts)
}


#' Percentile of Score
#'
#' Helper function that calculates percentiles of scores for batch-correction method 'pn'
#' (percentile normalization). R-implementation of Claire Duvallet's 'percentileofscore()' for
#' python.
#'
#' Calculates the number of values that bigger than reference (left) and the number of values that
#' are smaller than the reference (right). Percentiles of scores are given in the interval \eqn{I:[0,100]}.
#' Depending on type of calculation, the score will be computed as follows:
#'
#' \code{rank = (right + left + ifelse(right > left, 1, 0)) * 50/n}
#'
#' \code{weak = right / n*100}
#'
#' \code{strict = left / n*100}
#'
#' \code{mean = (right + left) * 50/n)}
#'
#' @keywords Log Ratio Transformation
#' @param cnt.vec, vector of cnts that acts as reference for score calculation
#' @param cnt, value to calculate score for
#' @param type, one of 'rank', 'weak', 'strict' or 'mean' to determine how score is calculated
#' @return a score for given counts in relation to reference counts
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


# EXTERNAL FUNCTIONS ------------------------------------------------------

### NOT MINE - reference or do sth. else
# 1 - ilr transform of the data, isoLMR function from robCompositions package, with changes
# https://github.com/matthias-da/robCompositions/blob/master/R/isomLR.R
# ---

#' KA changed the function to add a min value when many zeroes in data (prob with log and division by 0 otherwise)
#' @param fast if TRUE, it is approx. 10 times faster but numerical problems may occur for high dimensional data
#' @noRd
ilr.transfo = function(x, fast = TRUE, offset = 0) {
  if(any(x==0) & offset ==0)
    stop("make sure you use pseudo counts before normalisation to avoid 0 values with log ratio transformation")
  # ilr transformation
  x.ilr = matrix(NA, nrow = nrow(x), ncol = ncol(x)-1)
  D = ncol(x)
  # KA added: a little something to avoid 0 values
  if (fast)
  {
    for (i in seq_len(ncol(x.ilr)) )
    {
      #x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, (i+1) : D, drop = FALSE]), 1, prod) + offset)^(1 / (D-i))) / (x[,i]+ offset)) ToDo: remove once it works
      x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, seq.int(from=(i+1),to=D,by=1), drop = FALSE]), 1, prod) + offset)^(1 / (D-i))) / (x[,i]+ offset))
      #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop = FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
    }
  } else {
    for (i in seq_len(ncol(x.ilr)) )
    {
      x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(apply(as.matrix(x[, seq.int(from=(i+1),to=D,by=1)]), 1, function(x){exp(log(x))})/(x[, i]+ offset) + offset)
      #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(apply(as.matrix(x[,(i+1):D]), 1, function(x){exp(log(x))})/(x[,i]))
    }
  }
  ### ToDo: take care of this class-mess!
  class(x.ilr) = c(class(x.ilr), 'ilr')
  return(as.matrix(x.ilr))
}

# 2 - back transformation from ilr to clr space
clr.backtransfo = function(x) {
  # construct orthonormal basis
  V = matrix(0, nrow = ncol(x), ncol = ncol(x)-1)
  for( i in seq_len(ncol(V)) )
  {
    V[seq_len(i), i] = 1/i
    V[i+1, i] = (-1)
    V[, i] = V[, i] * sqrt(i/(i+1))
  }
  rownames(V) = colnames(x)
  return(V)

}

# CLR transformation
clr.transfo = function(x, offset = 0) {
  if (any(is.na(x) | x < 0)) {
    stop('\nFor CLR transformation, data must be non-negative with no missing values\n', call. = FALSE)
  }
  if(any(x==0) & offset ==0)
    stop("make sure you use pseudo counts before normalisation to avoid 0 values with log ratio transformation")

  # KA added
  #offset = min(x[which(x != 0)])*0.01


  #if (dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
  if (dim(x)[2] == 1)
  {
    res = list(x.clr = x, gm = rep(1, dim(x)[1]))
  } else{
    geometricmean = function (x) {
      #       if (any(na.omit(x == 0)))
      #         0
      #       else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
      #     }
      # KA changed to
      exp(mean(log(x + offset)))
    }
    gm = apply(x, 1, geometricmean)
    # KA changed
    x.clr = log((x + offset) / (gm))
    res = x.clr #list(x.clr = x.clr, gm = gm)
  }
  ### ToDo: take care of this class-mess!
  class(res) = c(class(res), "clr")
  return(res)
}



#' Centered Log-Ratio Transformation
#'
#' Internal function that performs CLR-transformation on input-matrix.
#' Formula is: clr(mtx) = ln( mtx / geometric_mean(mtx_samples))
#'
#' @keywords Log Ratio Transformation
#' @param input.mtx A matrix of counts (samples x features)
#' @param offset An (OPTIONAL) offset in case of sparse matrix. Function will
#' add an offset of 1/#features if matrix is sparse and offset not provided.
#' @return A matrix of transformed counts of same size and orientation as the
#' input.
mbecCLR <- function(input.mtx, offset = 0) {
  if( dim(input.mtx)[2] < 2 ) {
    message("No basis for transformation. Matrix contains less than 2 features, returning unchanged.")
    return(input.mtx)
  }
  # 1. stop for negative values and NAs
  if( any(input.mtx < 0 | is.na(input.mtx)) ) {
    stop("Examine your data for NAs and negative values, CLR transformation requires complete positive values.\n")
  }
  if( any(input.mtx==0) & offset==0 ) {
    message("Found zeros, function will add a small pseudo-count (1/#features) to enable log-ratio transformation.")
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







