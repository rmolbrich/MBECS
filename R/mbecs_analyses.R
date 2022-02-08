# ANALYSIS FUNCTIONS ------------------------------------------------------

#' Relative Log Expression Plot
#'
#' Takes two covariates, i.e., group and batch, and computes the RLE-plot over
#' the grouping of the first covariate, colored by the second covariate.
#' Effectively illustrating the relative expression between samples from
#' different batches within the respective study groups. Other covariates can be
#' chosen as input and the function will check for factors and convert if
#' necessary. Categorical factors, e.g., group membership, sex and batch,
#' produce the best result.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input
#' is an MbecData-object. If cumulative log-ratio (clr) and total sum-scaled
#' (tss) abundance matrices are part of the input, i.e., 'mbecTransform()' was
#' used, they can be selected as input by using the 'type' argument with either
#' "otu", "clr" or "tss". If batch effect corrected matrices are available, they
#' can be used by specifying the 'type' argument as "cor" and using the 'label'
#' argument to select the appropriate matrix by its denominator, e.g., for batch
#' correction method ComBat this would be "bat", for RemoveBatchEffects from the
#' limma package this is "rbe". Default correction method-labels are "ruv3",
#' "bmc","bat","rbe","pn","svd".
#'
#' The combination of 'type' and 'label' argument basically accesses the
#' attribute 'cor', a list that stores all matrices of corrected counts.
#' This list can also be accessed via getter and setter methods. Hence, the user
#' can supply their own matrices with own names.
#'
#' @keywords RLE relative log expression
#' @param input.obj MbecData-object
#' @param model.vars two covariates of interest to select by. First relates to
#' 'batch' and the second to relevant grouping.
#' @param type Which abundance matrix to use for the calculation.
#' @param label Which corrected abundance matrix to use for analysis.
#' @param return.data logical if TRUE returns the data.frame required for
#' plotting. Default (FALSE) will return plot object.
#' @return Either a ggplot2 object or a formatted data-frame to plot from.
#' @export
#'
#' @examples
#' # This will return the data.frame for plotting.
#' data.RLE <- mbecRLE(input.obj=dummy.mbec, type="clr",
#' model.vars=c('group','batch'), return.data=TRUE)
#'
#' # This will return the ggplot2 object for display, saving and modification.
#' plot.RLE <- mbecRLE(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' type="clr", return.data=FALSE)
mbecRLE <- function(input.obj, model.vars = c("batch","group"), type="clr",
                    label=character(), return.data = FALSE) {

  tmp <- mbecGetData(input.obj=input.obj, orientation="fxs",
                     required.col=eval(model.vars), type=eval(type),
                     label=label)
  tmp.cnts <- tmp[[1]]
  tmp.meta <- tmp[[2]] %>% tibble::rownames_to_column(., var="specimen")

  if( !(type == "cor") ) {
    label <- type
  }

  tmp.long <- NULL
  for(g.idx in unique(tmp.meta[, eval(model.vars[2])])) {
    message("Calculating RLE for group: ", g.idx)

    tmp.cnts.group <-
      dplyr::select(
        tmp.cnts,tmp.meta$specimen[tmp.meta[,eval(model.vars[2])] %in% g.idx])

    feature.med <- apply(tmp.cnts.group, 1, stats::median)

    tmp.group.long <- apply(tmp.cnts.group, 2,
                            function(sample.col) sample.col - feature.med) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(cols = everything(),
                          names_to = "specimen", values_to = "values")
    tmp.long <- rbind.data.frame(tmp.long, tmp.group.long)
  }

  tmp.long <- dplyr::left_join(tmp.long, tmp.meta,
                               by = "specimen") %>%
    dplyr::mutate(plot.order = paste(get(model.vars[2]),
                                     get(model.vars[1]), sep = "_")) %>%
    dplyr::arrange(plot.order) %>%
    dplyr::mutate(specimen = factor(specimen, levels = unique(specimen)))

  if (return.data) {
    return(tmp.long)
  }
  return(mbecRLEPlot(tmp.long, model.vars, label=label))
}


#' Principal Component Analysis Plot
#'
#' Takes two covariates, i.e., group and batch, and computes the ordination-plot
#' for user-selected principal components. Covariates determine sample-shape and
#' color and can be switched to shift the emphasis on either group. In addition
#' to the ordination-plot, the function will show the distribution of
#' eigenvalues (colored by the second covariate) on their respective principal
#' components.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input
#' is an MbecData-object. If cumulative log-ratio (clr) and total sum-scaled
#' (tss) abundance matrices are part of the input, i.e., 'mbecTransform()' was
#' used, they can be selected as input by using the 'type' argument with either
#' "otu", "clr" or "tss". If batch effect corrected matrices are available, they
#' can be used by specifying the 'type' argument as "cor" and using the 'label'
#' argument to select the appropriate matrix by its denominator, e.g., for batch
#' correction method ComBat this would be "bat", for RemoveBatchEffects from the
#' limma package this is "rbe". Default correction method-labels are "ruv3",
#' "bmc","bat","rbe","pn","svd".
#'
#' The combination of 'type' and 'label' argument basically accesses the
#' attribute 'cor', a list that stores all matrices of corrected counts.
#' This list can also be accessed via getter and setter methods. Hence, the user
#' can supply their own matrices with own names.
#'
#' @keywords PCA principal component analysis
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct
#' orientation is handled internally)
#' @param model.vars two covariates of interest to select by first variable
#' selects color (batch) and second one determines shape (group)
#' @param pca.axes numeric vector which axes to plot, first is X and second is Y
#' @param type Which abundance matrix to use for the calculation.
#' @param label Which corrected abundance matrix to use for analysis.
#' @param return.data logical if TRUE returns the data.frame required for
#' plotting. Default (FALSE) will return plot object.
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return the data.frame for plotting.
#' data.PCA <- mbecPCA(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), pca.axes=c(1,2), return.data=TRUE)
#'
#' # This will return the ggplot2 object for display, saving and modification.
#' # Selected PCs are PC3 on x-axis and PC2 on y-axis.
#' plot.PCA <- mbecPCA(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), pca.axes=c(3,2), return.data=FALSE)
setGeneric("mbecPCA", signature = "input.obj",
           function(input.obj, model.vars = c("batch", "group"),
                    pca.axes = c(1, 2), type="clr", label=character(),
                    return.data = FALSE) standardGeneric("mbecPCA"))

## In this form it works for 'phyloseq' and 'mbecData' objects
.mbecPCA <- function(input.obj, model.vars=c("batch", "group"),
                     pca.axes=c(1,2), type="clr", label=character(),
                     return.data=FALSE) {

  tmp <- mbecGetData(input.obj=input.obj, orientation = "sxf",
                     required.col=eval(model.vars), type=eval(type),
                     label=label)

  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  if( !(type == "cor") ) {
    label <- type
  }

  # calculate IQR and sort counts in decreasing order
  iqr <- apply(tmp.cnts, 2, stats::IQR)
  tmp.cnts <- tmp.cnts[, order(iqr, decreasing = TRUE)]

  PCA <- stats::prcomp(tmp.cnts, scale = FALSE)

  axes.number <- dim(PCA$x)[2]
  axes.names <- paste("PC", seq_len(axes.number), sep = "")

  plot.df <- PCA$x %>%
    data.frame(stringsAsFactors = FALSE) %>%
    dplyr::rename_at(seq_len(axes.number), ~axes.names) %>%
    tibble::rownames_to_column(var = "sID") %>%
    dplyr::left_join(tmp.meta, by = "sID")

  metric.df <- data.frame(var.explained=round((100 * PCA$sdev^2)/
                                      (sum(PCA$x^2/max(1,nrow(PCA$x) - 1))), 2),
                                      row.names = axes.names) %>%
    dplyr::mutate(axis.min = floor(apply(PCA$x, 2, function(col) min(col)))) %>%
    dplyr::mutate(axis.max = ceiling(apply(PCA$x, 2, function(col) max(col))))

  for (idx in seq_along(model.vars) ) {
    if (!is.factor(plot.df[, eval(model.vars[idx])])) {
      warning("Grouping variables need to be factors.
              Coercing to factor now, adjust beforehand to get best results.")
      plot.df[,eval(model.vars[idx])] <- factor(plot.df[,eval(model.vars[idx])])
    }
  }

  # No plotting, just return data.
  if (return.data) {
    return(list(plot.df, metric.df, pca.axes))
  }

  return(mbecPCAPlot(plot.df, metric.df, model.vars, pca.axes, label=label))
}


#' Principal Component Analysis Plot for MbecData
#'
#' Takes two covariates, i.e., group and batch, and computes the ordination-plot
#' for user-selected principal components. Covariates determine sample-shape and
#' color and can be switched to shift the emphasis on either group. In addition
#' to the ordination-plot, the function will show the distribution of
#' eigenvalues (colored by the second covariate) on their respective principal
#' components.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input
#' is an MbecData-object. If cumulative log-ratio (clr) and total sum-scaled
#' (tss) abundance matrices are part of the input, i.e., 'mbecTransform()' was
#' used, they can be selected as input by using the 'type' argument with either
#' "otu", "clr" or "tss". If batch effect corrected matrices are available, they
#' can be used by specifying the 'type' argument as "cor" and using the 'label'
#' argument to select the appropriate matrix by its denominator, e.g., for batch
#' correction method ComBat this would be "bat", for RemoveBatchEffects from the
#' limma package this is "rbe". Default correction method-labels are "ruv3",
#' "bmc","bat","rbe","pn","svd".
#'
#' The combination of 'type' and 'label' argument basically accesses the
#' attribute 'cor', a list that stores all matrices of corrected counts.
#' This list can also be accessed via getter and setter methods. Hence, the user
#' can supply their own matrices with own names.
#'
#' @keywords PCA principal component analysis
#' @param input.obj MbecData object
#' @param model.vars two covariates of interest to select by first variable
#' selects color (batch) and second one determines shape (group).
#' @param pca.axes numeric vector which axes to plot, first is X and second is Y
#' @param type Which abundance matrix to use for the calculation.
#' @param label Which corrected abundance matrix to use for analysis.
#' @param return.data logical if TRUE returns the data.frame required for
#' plotting. Default (FALSE) will return plot object.
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return the data.frame for plotting.
#' data.PCA <- mbecPCA(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), pca.axes=c(1,2), return.data=TRUE)
#'
#' # This will return the ggplot2 object for display, saving and modification.
#' # Selected PCs are PC3 on x-axis and PC2 on y-axis.
#' plot.PCA <- mbecPCA(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), pca.axes=c(3,2), return.data=FALSE)
setMethod("mbecPCA", "MbecData", function(input.obj,
                                          model.vars = c("batch", "group"),
                                          pca.axes = c(1, 2), type="clr",
                                          label=character(),
                                          return.data=FALSE) {
  .mbecPCA(input.obj, model.vars=model.vars, pca.axes=pca.axes, type=type,
           label=label, return.data=return.data)
})


#' Feature Differential Abundance Box-Plot
#'
#' Displays the abundance of a selected feature, grouped/colored by a covariate,
#' i.e., batch, in a box-plot. Includes the density-plot, i.e., the distribution
#' of counts for each sub-group. Selection methods for features are "TOP" and
#' "ALL" which select the top-n or all features respectively. The default value
#' for the argument 'n' is 10. If 'n' is supplied with a vector of feature
#' names, e.g., c("OTU1","OTU5", "OTU10"), of arbitrary length, the argument
#' 'method' will be ignored and only the given features selected for plotting.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input
#' is an MbecData-object. If cumulative log-ratio (clr) and total sum-scaled
#' (tss) abundance matrices are part of the input, i.e., 'mbecTransform()' was
#' used, they can be selected as input by using the 'type' argument with either
#' "otu", "clr" or "tss". If batch effect corrected matrices are available, they
#' can be used by specifying the 'type' argument as "cor" and using the 'label'
#' argument to select the appropriate matrix by its denominator, e.g., for batch
#' correction method ComBat this would be "bat", for RemoveBatchEffects from the
#' limma package this is "rbe". Default correction method-labels are "ruv3",
#' "bmc","bat","rbe","pn","svd".
#'
#' The combination of 'type' and 'label' argument basically accesses the
#' attribute 'cor', a list that stores all matrices of corrected counts.
#' This list can also be accessed via getter and setter methods. Hence, the user
#' can supply their own matrices with own names.
#'
#' @keywords Box abundance density
#' @param input.obj MbecData object
#' @param method One of 'ALL' or 'TOP' for 'n' most variable features, DEFAULT
#' is 'ALL'.
#' @param n Number of OTUs to display for 'TOP' method, or vector of specific
#' feature names to select.
#' @param model.var Covariate to group by, default is "batch".
#' @param type Which abundance matrix to use for the calculation.
#' @param label Which corrected abundance matrix to use for analysis.
#' @param return.data logical if TRUE returns the data.frame required for
#' plotting. Default (FALSE) will return plot object.
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return the plot-frame of all features in the data-set.
#' data.Box <- mbecBox(input.obj=dummy.mbec, method='ALL', model.var='batch',
#' type='clr', return.data=TRUE)
#'
#' # This will return the ggplot2 object of the top 5 most variable features.
#' plot.Box <- mbecBox(input.obj=dummy.mbec, method='TOP', n=5,
#' model.var='batch', type='otu', return.data=FALSE)
mbecBox <- function(input.obj, method = c("ALL", "TOP"), n = 10,
                    model.var = "batch", type="clr", label=character(),
                    return.data = FALSE) {

  input.obj <- mbecProcessInput(input.obj = input.obj, required.col = model.var)
  # needs sxf orientation
  tmp <- mbecGetData(input.obj = input.obj, orientation = "sxf",
                     required.col = eval(model.var), type=eval(type),
                     label=label)
  tmp[[2]] <- tibble::rownames_to_column(tmp[[2]], var = "specimen")

  if( !(type == "cor") ) {
    label <- type
  }

  otu.idx <- colnames(tmp[[1]])

  tmp <- tmp[[1]] %>%
    tibble::rownames_to_column(var = "specimen") %>%
    dplyr::left_join(tmp[[2]], by = c(specimen = "specimen"))

  if (method[1] == "TOP") {
    iqr <- apply(tmp[, otu.idx], 2, stats::IQR)
    iqr <- iqr[order(iqr, decreasing = TRUE)]
    otu.idx <- names(iqr)[seq_len(min(length(otu.idx), n))]

    tmp <- tmp %>%
      dplyr::select(c(dplyr::all_of(otu.idx), "specimen", eval(model.var)))

  } else if (length(method) >= 2) {
    message("'Method' parameter contains multiple elements -
            using to select features.")
    otu.idx <- method
    tmp <- tmp %>%
      dplyr::select(c(dplyr::all_of(otu.idx), "specimen", eval(model.var)))

  }  # else is 'select-all-mode'

  if (return.data) {
    return(list(tmp, otu.idx))
  }

  return(mbecBoxPlot(tmp, otu.idx, model.var, label=label))
}


#' Feature Differential Abundance Heatmap
#'
#' Shows the abundance value of selected features in a heatmap. By default, the
#' function expects two covariates group and batch to depict clustering in these
#' groups. More covariates can be included. Selection methods for features are
#' "TOP" and "ALL" which select the top-n or all features respectively. The
#' default value for the argument 'n' is 10. If 'n' is supplied with a vector
#' of feature names, e.g., c("OTU1","OTU5", "OTU10"), of arbitrary length, the
#' argument method' will be ignored and only the given features selected for
#' plotting.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input
#' is an MbecData-object. If cumulative log-ratio (clr) and total sum-scaled
#' (tss) abundance matrices are part of the input, i.e., 'mbecTransform()' was
#' used, they can be selected as input by using the 'type' argument with either
#' "otu", "clr" or "tss". If batch effect corrected matrices are available, they
#' can be used by specifying the 'type' argument as "cor" and using the 'label'
#' argument to select the appropriate matrix by its denominator, e.g., for batch
#' correction method ComBat this would be "bat", for RemoveBatchEffects from the
#' limma package this is "rbe". Default correction method-labels are "ruv3",
#' "bmc","bat","rbe","pn","svd".
#'
#' The combination of 'type' and 'label' argument basically accesses the
#' attribute 'cor', a list that stores all matrices of corrected counts.
#' This list can also be accessed via getter and setter methods. Hence, the user
#' can supply their own matrices with own names.
#'
#' @keywords Heat abundance clustering
#' @param input.obj MbecData object
#' @param model.vars Covariates of interest to show in heatmap.
#' @param center Flag to activate centering, DEFAULT is TRUE.
#' @param scale Flag to activate scaling, DEFAULT is TRUE.
#' @param method One of 'ALL' or 'TOP' or a vector of feature names.
#' @param n Number of features to select in method TOP.
#' @param type Which abundance matrix to use for the calculation.
#' @param label Which corrected abundance matrix to use for analysis.
#' @param return.data Logical if TRUE returns the data.frame required for
#' plotting. Default (FALSE) will return plot object.
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return the plot-frame of all features in the data-set.
#' data.Heat <- mbecHeat(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' center=TRUE, scale=TRUE, method='ALL', return.data=TRUE)
#'
#' # This will return the ggplot2 object of the top 5 most variable features.
#' plot.Heat <- mbecHeat(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' center=TRUE, scale=TRUE, method='TOP', n=5, return.data=FALSE)
mbecHeat <- function(input.obj, model.vars = c("batch", "group"), center = TRUE,
                     scale = TRUE, method = "TOP", n = 10, type="clr",
                     label=character(), return.data = FALSE) {

  cols <- pals::tableau20(20)[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19)]

  tmp <- mbecGetData(input.obj=input.obj, orientation="sxf",
                     required.col=eval(model.vars), type=eval(type),
                     label=label)
  tmp.cnts <- tmp[[1]]
  tmp.meta <- tmp[[2]]
  otu.idx <- colnames(tmp[[1]])

  if( !(type == "cor") ) {
    label <- type
  }

  for (g.idx in c(seq_along(model.vars))) {
    if (!is.factor(tmp.meta[, eval(model.vars[g.idx])])) {
      warning("Grouping variables need to be factors. Coercing variable: ",
              eval(model.vars[g.idx]),
              " to factor now, adjust beforehand to get best results.")
      tmp.meta[, eval(model.vars[g.idx])] <-
        factor(tmp.meta[, eval(model.vars[g.idx])])
    }
  }
  tmp.cnts <- base::scale(tmp.cnts, center=eval(center), scale=eval(scale))
  tmp.cnts <- base::scale(t(tmp.cnts), center=eval(center), scale=eval(scale))
  if (method[1] == "TOP") {
    iqr <- apply(tmp.cnts[otu.idx, ], 1, stats::IQR)
    iqr <- iqr[order(iqr, decreasing = TRUE)]
    otu.idx <- names(iqr)[seq_len(min(length(otu.idx), n))]
    tmp.cnts <- tmp.cnts[otu.idx, ]
  } else if (length(method) >= 2)
  {
    message("'Method' parameter contains multiple elements -
            using to select features.")
    tmp.cnts <- tmp.cnts[method, ]

  }  # else is 'select-all-mode'
  if (return.data) {
    return(list(tmp.cnts, tmp.meta))
  }
  return(mbecHeatPlot(tmp.cnts, tmp.meta, model.vars, label=label))
}


#' Mosaic Sample Group Allocation
#'
#' Depicts the dispersion of samples over two (preferentially categorical)
#' covariates of interest. Effectively showing, the un-/evenness within and
#' between covariates to inform the choice of methods for the subsequent steps
#' in an analysis.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input
#' for the data-set can be an MbecData-object.
#'
#' @keywords Mosaic sample allocation
#' @param input.obj MbecData object
#' @param model.vars Two covariates of interest to the sample allocation.
#' @param return.data Logical if TRUE returns the data.frame required for
#' plotting. Default (FALSE) will return plot object.
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' # This will return the plot-df of the samples grouped by group and batch.
#' data.Mosaic <- mbecMosaic(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), return.data=TRUE)
#'
#' # Return the ggplot2 object of the samples grouped by group and batch
#' plot.Mosaic <- mbecMosaic(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), return.data=FALSE)
mbecMosaic <- function(input.obj, model.vars = c("batch", "group"),
                       return.data = FALSE) {

  cols <- pals::tableau20(20)

  tmp <- mbecGetData(input.obj, orientation = "sxf",
                     required.col = eval(model.vars))
  tmp.meta <- tmp[[2]]

  if (length(model.vars) < 2) {
    message("Only one variable specified for Mosaic-plot, two are required!")
  } else if (length(model.vars) > 2) {
    message("More than 2 variables specified. Mosaic will take the first two.")
  }

  for (g.idx in eval(model.vars)) {
    if (!is.factor(tmp.meta[, eval(g.idx)])) {
      warning("Grouping variables need to be factors. Coercing variable: '",
              eval(g.idx),
              "' to factor now, adjust beforehand to get best results.")
      tmp.meta[, eval(g.idx)] <- as.factor(tmp.meta[, eval(g.idx)])
    }
  }
  n.observations <- base::dim(tmp.meta)[1]
  study.summary <- base::table(tmp.meta[, eval(model.vars[1])],
                               tmp.meta[, eval(model.vars[2])]) %>%
    as.data.frame() %>%
    dplyr::mutate(Freq.scaled = Freq/n.observations)

  if (return.data)
  {
    return(study.summary)
  }  # else return the plots

  return(mbecMosaicPlot(study.summary, model.vars))

}


# VARIANCE CALCULATION ----------------------------------------------------


#' Estimate Explained Variance
#'
#' The function offers a selection of methods/algorithms to estimate the
#' proportion of variance that can be attributed to covariates of interest.
#' This shows, how much variation is explained by the treatment effect, which
#' proportion is introduced by processing in batches and the leftover variance,
#' i.e., residuals that are not currently explained. Covariates of interest
#' (CoI) are selected by the user and the function will incorporate them into
#' the model building for the respective algorithm. The user can select from
#' five different approaches to adapt to the characteristics of the data-set,
#' e.g., LMMs are a better choice than LMs for a very unbalanced study design.
#' Available approaches are: Linear Model (lm), Linear Mixed Model (lmm),
#' Redundancy Analysis (rda), Principal Variance Component Analysis (pvca) or
#' Silhouette Coefficient (s.coef).
#'
#' Linear Model (lm): An additive model of all covariates is fitted to each
#' feature respectively and the proportion of variance is extracted for each
#' covariate (OTU_x ~ covariate_1 + covariate_2 + ...).
#'
#' Linear Mixed Model (lmm): All but the first covariate are considered mixed
#' effects. A model is fitted to each OTU respectively and the proportion of
#' variance extracted for each covariate.
#' (OTU_x ~ covariate_1 + (1|covariate_2) + (1|...)).
#'
#' partial Redundancy Analysis (rda): Iterates over given covariates, builds a
#' model of all covariates that includes one variable as condition/constraint
#' and then fits it to the feature abundance matrix. The difference in explained
#' variance between the full- and the constrained-model is then attributed to
#' the constraint. (cnts ~ group + Condition(batch) vs. cnts ~ group + batch)
#'
#' Principal Variance Component Analysis (pvca): Algorithm - calculate the
#' correlation of the fxs count-matrix - from there extract the eigenvectors and
#' eigenvalues and calculate the proportion of explained variance per
#' eigenvector (i.e. principal component) by dividing the eigenvalues by the
#' sum of eigenvalues. Now select as many PCs as required to fill a chosen
#' quota for the total proportion of explained variance. Iterate over all PCs
#' and fit a linear mixed model that contains all covariates as random effect
#' and all unique interactions between two covariates. Compute variance
#' covariance components form the resulting model --> From there we get the
#' Variance that each covariate(variable) contributes to this particular PC.
#' Then just standardize variance by dividing it through the sum of variance for
#' that model. Scale each PCs results by the proportion this PC accounted for in
#' the first place. And then do it again by dividing it through the total amount
#' of explained variance, i.e. the cutoff to select the number of PCs to take
#' (obviously not the cutoff but rather the actual values for the selected PCs).
#' Finally take the average over each random variable and interaction term and
#' display in a nice plot.
#'
#' Silhouette Coefficient (s.coef): Calculate principal components and get
#' sample-wise distances on the resulting (sxPC) matrix. Then iterate over all
#' the covariates and calculate the cluster silhouette (which is basically
#' either zero, if the cluster contains only a single element, or it is the
#' distance to the closest different cluster minus the distance of the sample
#' within its own cluster divided (scaled) by the maximum distance). Average
#' over each element in a cluster for all clusters and there is the
#' representation of how good the clustering is. This shows how good a
#' particular covariate characterizes the data, i.e., a treatment variable for
#' instance may differentiate the samples into treated and untreated groups
#' which implies two clusters. In an ideal scenario, the treatment variable,
#' i.e., indicator for some biological effect would produce a perfect
#' clustering. In reality, the confounding variables, e.g., batch, sex or age,
#' will also influence the ordination of samples. Hence, the clustering
#' coefficient is somewhat similar to the amount of explained variance metric
#' that the previous methods used. If used to compare an uncorrected data-set to
#' a batch-corrected set, the expected result would be an increase of clustering
#' coefficient for the biological effect (and all other covariates - because a
#' certain amount of uncertainty was removed from the data) and a decrease for
#' the batch effect.
#'
#' The function returns a data-frame for further analysis - the report functions
#' (mbecReport and mbecReportPrelim) will automatically produce plots. Input for
#' the data-set can be an MbecData-object, a phyloseq-object or a list that
#' contains counts and covariate data. The covariate table requires an 'sID'
#' column that contains sample IDs equal to the sample naming in the counts
#' table. Correct orientation of counts will be handled internally.
#'
#' @keywords Model Evaluation Variance
#' @param input.obj MbecData object
#' @param model.vars Vector of covariates to include in model-construction, in
#' case parameter 'model.form' is not supplied.
#' @param method Select method of modeling: Linear Model (lm), Linear Mixed
#' Model (lmm), Redundancy Analysis (rda), Principal Variance Component Analysis
#' (pvca) or Silhouette Coefficient (s.coef).
#' @param model.form string that describes a model formula, i.e.,
#' 'y ~ covariate1 + (1|covariate2)'.
#' @param type Which abundance matrix to use for the calculation.
#' @param label Which corrected abundance matrix to use for analysis.
#' @param no.warning (OPTIONAL) True/False-flag that should turn of singularity
#' warnings, but it doesn't quite work
#' @param na.action (OPTIONAL) set NA handling, will take global option if not
#' supplied
#' @return Data.frame that contains proportions of variance for given covariates
#' in every feature.
#' @include mbecs_classes.R
#' @export
#'
#' @examples
#' # This will return a data-frame that contains the variance attributable to
#' # group and batch according to linear additive model.
#' df.var.lm <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c("batch", "group"), method='lm', type='clr')
#' # This will return a data-frame that contains the variance attributable to
#' # group and batch according to principal variance component analysis.
#' df.var.pvca <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c("batch", "group"), method='pvca')
mbecModelVariance <- function(input.obj, model.vars=character(),
                              method=c("lm","lmm", "rda", "pvca", "s.coef"),
                              model.form=NULL, type="clr", label=character(),
                              no.warning = TRUE, na.action = NULL) {
  oldw <- getOption("warn")
  if (no.warning) {
    options(warn = -1)
  }
  on.exit(options(warn = oldw))

  # handle optional parameters
  if (is.null(na.action)) {
    na.action <- getOption("na.action")
  }

  method <- match.arg(method, choices = c("lm","lmm", "rda", "pvca", "s.coef"))
  type <- match.arg(type, choices = c("otu","clr","tss","ass","cor"))

  ## PVCA stuff
  pct_threshold <- 0.5876  # threshold for explained variances

  tmp <- mbecGetData(input.obj = input.obj, orientation = "sxf",
                     required.col=eval(model.vars), type=eval(type),
                     label=label)
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  ## Adjust 'type' variable if selected matrix comes from 'cor' list
  if( type == "cor" && length(label) != 0 ) type <- label

  ## CHECK if grouping variables are factors
  for (g.idx in eval(model.vars)) {
    if (!is.factor(tmp.meta[, eval(g.idx)])) {
      warning("Grouping variables need to be factors. Coercing variable: ",
              eval(g.idx),
              " to factor now, adjust beforehand to get best results.")
      tmp.meta[, eval(g.idx)] <- as.factor(tmp.meta[, eval(g.idx)])
    }
  }

  if (method == "lm") {
    res <- mbecModelVarianceLM(model.form,model.vars,tmp.cnts,tmp.meta,type)
    return(res)

  } else if (method == "lmm") {
    res <- mbecModelVarianceLMM(model.form,model.vars,tmp.cnts,tmp.meta,type)
    return(res)

  } else if (method == "rda") {
    res <- mbecModelVarianceRDA(model.vars,tmp.cnts,tmp.meta,type)
    return(res)

  } else if (method == "pvca") {
    res <- mbecModelVariancePVCA(model.vars, tmp.cnts, tmp.meta,
                                 type, pct_threshold, na.action)
    return(res)

  } else if (method == "s.coef") {
    res <- mbecModelVarianceSCOEF(model.vars, tmp.cnts, tmp.meta, type)
    return(res)

  }
  message("Input doesn't apply to any available method. Nothing was done here.")
  return(NULL)
}


#' Estimate Explained Variance with Linear Models
#'
#' The function offers a selection of methods/algorithms to estimate the
#' proportion of variance that can be attributed to covariates of interest.
#' This shows, how much variation is explained by the treatment effect, which
#' proportion is introduced by processing in batches and the leftover variance,
#' i.e., residuals that are not currently explained. Covariates of interest
#' (CoI) are selected by the user and the function will incorporate them into
#' the model.
#'
#' Linear Model (lm): An additive model of all covariates is fitted to each
#' feature respectively and the proportion of variance is extracted for each
#' covariate (OTU_x ~ covariate_1 + covariate_2 + ...).
#'
#' @keywords LM Proportion of Variance
#' @param model.form Formula for linear model, function will create simple
#' additive linear model if this argument is not supplied.
#' @param model.vars Covariates to use for model building if argument
#' 'model.form' is not given.
#' @param tmp.cnts Abundance matrix in 'sample x feature' orientation.
#' @param tmp.meta Covariate table that contains at least the used variables.
#' @param type String the denotes data source, i.e., one of "otu","clr" or "tss"
#' for the transformed counts or the label of the batch corrected count-matrix.
#' @return Data.frame that contains proportions of variance for given covariates
#' in a linear modelling approach.
#' @include mbecs_classes.R
mbecModelVarianceLM <- function(model.form, model.vars, tmp.cnts, tmp.meta, type) {
  message("Fitting linear model to every feature and extract proportion of
          variance explained by covariates.")
  if (!is.null(model.form)) {
    message("Use provided model formula.")
    tmp.formula <- stats::as.formula(model.form)
  } else {
    message("Construct formula from covariates.")
    tmp.formula <- stats::as.formula(paste("y", " ~ ", paste(model.vars,
                                                            collapse = " + ")))
  }

  # Add a progress bar.
  features <- colnames(tmp.cnts)
  lm.pb <- utils::txtProgressBar(min=0, max=length(features), style=3,
                                 width=80, char="=")

  model.variances <- NULL
  for (x in seq_along(features)) {
    y <- tmp.cnts[[eval(x)]]

    model.fit <- stats::lm(tmp.formula, data = tmp.meta)
    model.variances <- rbind.data.frame(model.variances,
                                        mbecVarianceStats(model.fit))
    utils::setTxtProgressBar(lm.pb, x)
  }
  # close progress bar
  close(lm.pb)

  res <- dplyr::mutate(model.variances, type = eval(type))
  attr(res, "modelType") <- "anova"

  return(res)
}


#' Estimate Explained Variance with Linear Mixed Models
#'
#' The function offers a selection of methods/algorithms to estimate the
#' proportion of variance that can be attributed to covariates of interest.
#' This shows, how much variation is explained by the treatment effect, which
#' proportion is introduced by processing in batches and the leftover variance,
#' i.e., residuals that are not currently explained. Covariates of interest
#' (CoI) are selected by the user and the function will incorporate them into
#' the model.
#'
#' Linear Mixed Model (lmm): Only the first covariate is considered a mixed
#' effect. A model is fitted to each OTU respectively and the proportion of
#' variance extracted for each covariate.
#' (OTU_x ~ covariate_2.. + covariate_n + (1|covariate_1)
#'
#' @keywords LMM Proportion of Variance
#' @param model.form Formula for linear mixed model, function will create simple
#' additive linear mixed model if this argument is not supplied.
#' @param model.vars Covariates to use for model building if argument
#' 'model.form' is not given.
#' @param tmp.cnts Abundance matrix in 'sample x feature' orientation.
#' @param tmp.meta Covariate table that contains at least the used variables.
#' @param type String the denotes data source, i.e., one of "otu","clr" or "tss"
#' for the transformed counts or the label of the batch corrected count-matrix.
#' @return Data.frame that contains proportions of variance for given covariates
#' in a linear mixed modelling approach.
#' @include mbecs_classes.R
mbecModelVarianceLMM <- function(model.form, model.vars, tmp.cnts, tmp.meta,
                                 type) {

  message("Fitting linear-mixed model to every feature and extract proportion
          of variance explained by covariates.")
  # FixMe: maybe include the option to adjust this
  control <- lme4::lmerControl(calc.derivs=TRUE, check.rankX="stop.deficient")

  if (!is.null(model.form)) {
    message("Use provided model formula.")
    tmp.formula <- stats::as.formula(model.form)
  } else {
    message("Construct formula from covariates.")
    f.terms <- paste("(1|", model.vars, ")", sep = "")
    tmp.formula <- stats::as.formula(paste(paste("y", paste(model.vars[-1],
                                                            collapse = " + "),
                                                 sep = " ~ "),
                                           paste(f.terms[1], collapse = " + "),
                                           sep = " + "))
  }
  # Add a progress bar.
  features <- colnames(tmp.cnts)
  lmm.pb <- utils::txtProgressBar(min = 0, max = length(features), style = 3,
                                  width = 80, char="=")

  model.variances <- NULL
  for (x in seq_along(features)) {
    y <- tmp.cnts[[eval(x)]]  # meh

    model.fit <- lme4::lmer(tmp.formula, data = tmp.meta, control = control)
    model.variances <- rbind.data.frame(model.variances,
                                        mbecVarianceStats(model.fit))

    utils::setTxtProgressBar(lmm.pb, x)
  }
  close(lmm.pb)
  res <- dplyr::mutate(model.variances, type = eval(type))
  attr(res, "modelType") <- "linear-mixed model"

  return(res)
}

#' Estimate Explained Variance with Redundancy Analysis
#'
#' The function offers a selection of methods/algorithms to estimate the
#' proportion of variance that can be attributed to covariates of interest.
#' This shows, how much variation is explained by the treatment effect, which
#' proportion is introduced by processing in batches and the leftover variance,
#' i.e., residuals that are not currently explained. Covariates of interest
#' (CoI) are selected by the user and the function will incorporate them into
#' the model.
#'
#' partial Redundancy Analysis (rda): Iterates over given covariates, builds a
#' model of all covariates that includes one variable as condition/constraint
#' and then fits it to the feature abundance matrix. The difference in explained
#' variance between the full- and the constrained-model is then attributed to
#' the constraint. (cnts ~ group + Condition(batch) vs. cnts ~ group + batch)
#'
#' @keywords Partial Redundancy Analysis of Variance
#' @param model.vars Covariates to use for model building.
#' @param tmp.cnts Abundance matrix in 'sample x feature' orientation.
#' @param tmp.meta Covariate table that contains at least the used variables.
#' @param type String the denotes data source, i.e., one of "otu","clr" or "tss"
#' for the transformed counts or the label of the batch corrected count-matrix.
#' @return Data.frame that contains proportions of variance for given covariates
#' in a partial redundancy analysis approach.
#' @include mbecs_classes.R
mbecModelVarianceRDA <- function(model.vars, tmp.cnts, tmp.meta, type) {

  model.variances <- data.frame(matrix(nrow = length(model.vars),
                                      ncol = 1, dimnames = list(model.vars,
                                                                eval(type))))
  # Add a progress bar.
  rda.pb <- utils::txtProgressBar(min = 0, max = length(model.vars), style = 3,
                                  width = 80, char="=")

  for (condition.idx in seq_along(model.vars)) {
    tmp.formula <- stats::as.formula(
      paste("tmp.cnts"," ~ ",
            paste(model.vars[-eval(condition.idx)],
                  "+", collapse = " "), " Condition(",
            model.vars[eval(condition.idx)],")", sep = ""))

    tmp.rda.covariate <- vegan::rda(tmp.formula,
                                    data = tmp.meta)

    model.variances[condition.idx,eval(type)] <-
      summary(tmp.rda.covariate)$partial.chi *
      100/summary(tmp.rda.covariate)$tot.chi

    utils::setTxtProgressBar(rda.pb, condition.idx)
  }
  close(rda.pb)
  res <- data.frame(t(model.variances)) %>%
    dplyr::mutate(type = eval(type))
  attr(res, "modelType") <- "rda"

  return(res)
}


#' Estimate Explained Variance with Principal Variance Component Analysis
#'
#' The function offers a selection of methods/algorithms to estimate the
#' proportion of variance that can be attributed to covariates of interest.
#' This shows, how much variation is explained by the treatment effect, which
#' proportion is introduced by processing in batches and the leftover variance,
#' i.e., residuals that are not currently explained. Covariates of interest
#' (CoI) are selected by the user and the function will incorporate them into
#' the model.
#'
#' Principal Variance Component Analysis (pvca): Algorithm - calculate the
#' correlation of the fxs count-matrix - from there extract the eigenvectors and
#' eigenvalues and calculate the proportion of explained variance per
#' eigenvector (i.e. principal component) by dividing the eigenvalues by the
#' sum of eigenvalues. Now select as many PCs as required to fill a chosen
#' quota for the total proportion of explained variance. Iterate over all PCs
#' and fit a linear mixed model that contains all covariates as random effect
#' and all unique interactions between two covariates. Compute variance
#' covariance components form the resulting model --> From there we get the
#' Variance that each covariate(variable) contributes to this particular PC.
#' Then just standardize variance by dividing it through the sum of variance for
#' that model. Scale each PCs results by the proportion this PC accounted for in
#' the first place. And then do it again by dividing it through the total amount
#' of explained variance, i.e. the cutoff to select the number of PCs to take
#' (obviously not the cutoff but rather the actual values for the selected PCs).
#' Finally take the average over each random variable and interaction term and
#' display in a nice plot.
#'
#' @keywords Principal Variance Component Analysis
#' @param model.vars Covariates to use for model building.
#' @param tmp.cnts Abundance matrix in 'sample x feature' orientation.
#' @param tmp.meta Covariate table that contains at least the used variables.
#' @param type String the denotes data source, i.e., one of "otu","clr" or "tss"
#' for the transformed counts or the label of the batch corrected count-matrix.
#' @param pct_threshold Cutoff value for accumulated variance in principal
#' components.
#' @param na.action Set NA handling, will take global option if not supplied.
#' @return Data.frame that contains proportions of variance for given covariates
#' in a principal variance component analysis approach.
#' @include mbecs_classes.R
mbecModelVariancePVCA <- function(model.vars, tmp.cnts, tmp.meta, type,
                                  pct_threshold, na.action) {
  n.vars <- length(model.vars)
  s.names <- rownames(tmp.cnts)

  tmp.cnts <- apply(tmp.cnts, 2, scale, center = TRUE, scale = FALSE) %>%
    t()
  colnames(tmp.cnts) <- s.names
  tmp.cor <- stats::cor(tmp.cnts)

  eCor <- eigen(tmp.cor)
  eVal <- eCor$values
  eVec <- eCor$vectors
  rownames(eVec) <- s.names
  n.eVal <- length(eVal)
  sum.eVal <- sum(eVal)
  prop.PCs <- eVal/sum.eVal

  n.PCs <- max(3, min(which(vapply(cumsum(prop.PCs),
                                   function(x) x >= pct_threshold,
                                   FUN.VALUE = logical(1)))))
  # suppress tibble verbose output
  suppressMessages(lmm.df <- eVec %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::select_at(seq_len(eval(n.PCs))) %>%
    cbind(., tmp.meta))

  f.terms <- paste("(1|", model.vars, ")", sep = "")

  for( var.idx in seq_len((n.vars - 1)) ) {
    for( interaction.idx in seq.int(from=(var.idx + 1), to=(n.vars), by=1) ) {
      f.terms <- c(f.terms, paste("(1|", model.vars[var.idx], ":",
                                  model.vars[interaction.idx], ")", sep = ""))
    }
  }
  n.effects <- length(f.terms) + 1
  randomEffectsMatrix <- matrix(data = 0, nrow = n.PCs, ncol = n.effects)

  model.formula <- stats::as.formula(paste("lmm.df[,vec.idx]", " ~ ",
                                           paste(f.terms,
                                                 collapse = " + "), sep = ""))

  for (vec.idx in seq_len(n.PCs)) {
    randomEffects <- data.frame(lme4::VarCorr(Rm1ML <-
                                              lme4::lmer(model.formula,lmm.df,
                                                         REML = TRUE,
                                                         verbose = 0,
                                                         na.action=na.action)))

    randomEffectsMatrix[vec.idx, ] <- as.numeric(randomEffects[, 4])
  }

  names.effects <- randomEffects[, 1]
  randomEffectsMatrix.std <- randomEffectsMatrix/rowSums(randomEffectsMatrix)

  scaled.eVal <- eVal/sum(eVal)
  randomEffectsMatrix.wgt <- randomEffectsMatrix.std*scaled.eVal[seq_len(n.PCs)]

  model.variances <- colSums(randomEffectsMatrix.wgt)/
    sum(colSums(randomEffectsMatrix.wgt))
  names(model.variances) <- names.effects

  res <- data.frame(t(model.variances)) %>%
    dplyr::mutate(type = eval(type))
  attr(res, "modelType") <- "PVCA - lmm"

  return(res)
}


#' Estimate Explained Variance with Silhouette Coefficient
#'
#' The function offers a selection of methods/algorithms to estimate the
#' proportion of variance that can be attributed to covariates of interest.
#' This shows, how much variation is explained by the treatment effect, which
#' proportion is introduced by processing in batches and the leftover variance,
#' i.e., residuals that are not currently explained. Covariates of interest
#' (CoI) are selected by the user and the function will incorporate them into
#' the model.
#'
#' Silhouette Coefficient (s.coef): Calculate principal components and get
#' sample-wise distances on the resulting (sxPC) matrix. Then iterate over all
#' the covariates and calculate the cluster silhouette (which is basically
#' either zero, if the cluster contains only a single element, or it is the
#' distance to the closest different cluster minus the distance of the sample
#' within its own cluster divided (scaled) by the maximum distance). Average
#' over each element in a cluster for all clusters and there is the
#' representation of how good the clustering is. This shows how good a
#' particular covariate characterizes the data, i.e., a treatment variable for
#' instance may differentiate the samples into treated and untreated groups
#' which implies two clusters. In an ideal scenario, the treatment variable,
#' i.e., indicator for some biological effect would produce a perfect
#' clustering. In reality, the confounding variables, e.g., batch, sex or age,
#' will also influence the ordination of samples. Hence, the clustering
#' coefficient is somewhat similar to the amount of explained variance metric
#' that the previous methods used. If used to compare an uncorrected data-set to
#' a batch-corrected set, the expected result would be an increase of clustering
#' coefficient for the biological effect (and all other covariates - because a
#' certain amount of uncertainty was removed from the data) and a decrease for
#' the batch effect.
#'
#' @keywords Silhouette Coefficient Assessment
#' @param model.vars Covariates to use for model building.
#' @param tmp.cnts Abundance matrix in 'sample x feature' orientation.
#' @param tmp.meta Covariate table that contains at least the used variables.
#' @param type String the denotes data source, i.e., one of "otu","clr" or "tss"
#' for the transformed counts or the label of the batch corrected count-matrix.
#' @return Data.frame that contains proportions of variance for given covariates
#' in a silhouette coefficient analysis approach.
#' @include mbecs_classes.R
mbecModelVarianceSCOEF <- function(model.vars, tmp.cnts, tmp.meta, type) {

  tmp.prcomp <- stats::prcomp(tmp.cnts,
                              center = TRUE, scale = FALSE)

  tmp.dist <- stats::dist(tmp.prcomp$x[,1:3],
                          method = "euclidian")

  avg.sil.df <- NULL
  for (var.elem in model.vars) {
    print(var.elem)

    tmp.sil <- cluster::silhouette(x=as.numeric(tmp.meta[,eval(var.elem)]),
                                  dist = tmp.dist)

    avg.sil.df <- rbind.data.frame(avg.sil.df,
                        data.frame(variable=var.elem,
                                   cluster=levels(tmp.meta[,eval(var.elem)]),
                                   sil.coefficient =
                                          c(summary(tmp.sil))$clus.avg.widths))
  }

  res <- avg.sil.df %>%
    dplyr::mutate(type = eval(type))
  attr(res, "modelType") <- "s.coef"

  return(res)
}



# VARIANCE STATISTICS -----------------------------------------------------



#' Wrapper for Model Variable Variance Extraction
#'
#' For a Linear (Mixed) Model, this function extracts the proportion of variance
#' that can be explained by terms and interactions and returns a named
#' row-vector.
#'
#' Linear Model: Perform an analysis of variance (ANOVA) on the model.fit and
#' return the Sum of squares for each term, scaled by the total sum of squares.
#'
#' Linear Mixed Model: employ helper function 'mbecMixedVariance' to extract
#' residuals, random effects and fixed effects components from the model. The
#' components are then transformed to reflect explained proportions of variance
#' for the model coefficients. The function implements transformation for
#' varying coefficients as well, but NO ADJUSTMENT for single or multiple
#' coefficients at this point.
#'
#' @keywords lm lmm proportion variance
#' @param model.fit A linear (mixed) model object of class 'lm' or 'lmerMod'.
#' @return A named row-vector, containing proportional variance for model terms.
#' @export
#'
#' @examples
#' # This will return the data.frame for plotting.
#' limo <- stats::lm(dummy.list$cnts[,1] ~ group + batch, data=dummy.list$meta)
#' vec.variance <- mbecVarianceStats(model.fit=limo)
mbecVarianceStats <- function(model.fit) {

  ### ToDo: implement glm versions of this

  # check validity of model fit
  mbecValidateModel(model.fit)

  # linear model
  if (is(model.fit, "lm")) {

    res <- mbecVarianceStatsLM(model.fit)

  } else if (is(model.fit, "lmerMod")) {

    res <- mbecVarianceStatsLMM(model.fit)

  } else if (is(model.fit, "glm")) {
    ### ToDon't
  }
  return(res)
}


#' Model Variable Variance Extraction from LM
#'
#' For a Linear Model, this function extracts the proportion of variance
#' that can be explained by terms and interactions and returns a named
#' row-vector.
#'
#' Linear Model: Perform an analysis of variance (ANOVA) on the model.fit and
#' return the Sum of squares for each term, scaled by the total sum of squares.
#'
#' @keywords lm proportion variance
#' @param model.fit A linear model object of class 'lm'.
#' @return A named row-vector, containing proportional variance for model terms.
mbecVarianceStatsLM <- function(model.fit) {

  vp <- stats::anova(model.fit) %>%
    data.frame() %>%
    dplyr::mutate(variance =
                    dplyr::select(., "Sum.Sq")/sum(dplyr::select(.,"Sum.Sq")),
                  .keep = "none") %>%
    t()

  return(vp)
}


#' Model Variable Variance Extraction from LMM
#'
#' For a Linear Mixed Model, this function extracts the proportion of variance
#' that can be explained by terms and interactions and returns a named
#' row-vector.
#'
#' Linear Mixed Model: employ helper function 'mbecMixedVariance' to extract
#' residuals, random effects and fixed effects components from the model. The
#' components are then transformed to reflect explained proportions of variance
#' for the model coefficients. The function implements transformation for
#' varying coefficients as well, but NO ADJUSTMENT for single or multiple
#' coefficients at this point.
#'
#' @keywords lmm proportion variance
#' @param model.fit A linear mixed model object of class 'lmerMod'.
#' @return A named row-vector, containing proportional variance for model terms.
mbecVarianceStatsLMM <- function(model.fit) {

  vc <- mbecMixedVariance(model.fit)

  lib.df <- data.frame(
    covariates = colnames(model.fit@frame),
    row.names = vapply(colnames(model.fit@frame),
                       function(covariate)
                         paste(paste(covariate,
                                     levels(model.fit@frame[, eval(covariate)]),
                                     sep = ""), collapse = ","),
                       FUN.VALUE = character(1)))

  total.var.LUT <- unlist(lapply(vc, function(var.comp) {
    if (length(var.comp) > 1) {
      weights <- (table(
        model.fit@frame[[lib.df[eval(paste(names(var.comp),
                                           collapse = ",")),
                                ]]])/nrow(model.fit@frame))
      var.comp %*% weights
    } else {
      names(var.comp) <- NULL
      var.comp
    }
  }), use.names = TRUE, recursive = FALSE)

  vp <- list()

  for (effect in names(vc)) {

    # value of this effect divided by the total variance
    tmp.t.var <- sum(total.var.LUT[-which(names(total.var.LUT) %in%
                                            eval(effect))])
    vp[[eval(effect)]] <- vc[[eval(effect)]]/(tmp.t.var +
                                                vc[[eval(effect)]])
  }

  vp <- unlist(vp)
  names(vp) <- gsub(".(Intercept)", replacement = "",
                    names(vp), fixed = TRUE)
  vp <- data.frame(t(vp))

  return(vp)
}


#' Mixed Model Variance-Component Extraction
#'
#' A helper function that extracts the variance components of linear mixed
#' models, i.e., residuals, random-effects, fixed-effects, scales them to
#' sample-size and returns a list of components.
#'
#' Uses 'lme4::VarCorr' to extract Residuals and random-effects components.
#' Standard Deviation of Residuals is stored as 'sc' attribute in the output of
#' 'VarCorr'.
#'
#' Uses 'lme4::fixef' to extract fixed-effects components, i.e., parameter
#' estimates. The attribute 'pp' of the model contains the dense model matrix
#' for fixed-effects parameters (X). The fixed effects variance, sigma2f, is the
#' variance of the matrix-multiplication beta times X (parameter vector by
#' model matrix)
#'
#' @keywords lmm proportion variance
#' @param model.fit A linear mixed model object of class 'lmerMod'.
#' @return A named list, containing proportional variance for model terms that
#' describe mixed effects.
#' @export
#'
#' @examples
#' # This will return the variance of random/mixed components.
#' limimo <- lme4::lmer(dummy.list$cnts[,1] ~ group + (1|batch),
#' data=dummy.list$meta)
#' list.variance <- mbecMixedVariance(model.fit=limimo)
mbecMixedVariance <- function(model.fit) {

  rVC <- lme4::VarCorr(model.fit)

  randomVar <- c(Residuals = attr(rVC, "sc")^2, lapply(rVC, diag))

  n.scaling <- (stats::nobs(model.fit) - 1)/stats::nobs(model.fit)

  fixedVar <- t(t(model.fit@pp$X) * lme4::fixef(model.fit)) %>%
    as.data.frame() %>%
    dplyr::select(!"(Intercept)") %>%
    dplyr::mutate(total.var = apply(., 1, sum)) %>%
    apply(2, function(effect) stats::var(effect) * n.scaling)

  fixedVar <- lapply(
    split(utils::head(fixedVar, -1)/sum(utils::head(fixedVar,-1)) *
            utils::tail(fixedVar, 1), names(utils::head(fixedVar, -1))),
                     unname)

  return(c(randomVar, fixedVar))
}


#' Validate Linear (Mixed) Models
#'
#' A helper function that calculates the collinearity between model variables
#' and stops execution, if the maximum value is bigger than the allowed
#' threshold.
#'
#' ToDo: maybe some additional validation steps and more informative output.
#'
#' @keywords collinearity model validation
#' @param model.fit lm() or lmm() output.
#' @param colinearityThreshold Cut-off for model rejection, I=[0,1].
#' @return No return values. Stops execution if validation fails.
#' @export
#'
#' @examples
#' # This will just go through if colinearity threshold is met.
#' limimo <- lme4::lmer(dummy.list$cnts[,1] ~ group + (1|batch),
#' data=dummy.list$meta)
#' mbecValidateModel(model.fit=limimo, colinearityThreshold=0.999)
mbecValidateModel <- function(model.fit, colinearityThreshold = 0.999) {
  ## ToDo: health & Safety
  if( colinScore(model.fit) > colinearityThreshold) {
    warning("Some covariates are strongly correlated.
            Please re-evaluate the variable selection and try again.")
  }
}


#' Variable Correlation Linear (Mixed) Models
#'
#' Takes a fitted model and computes maximum correlation between covariates as
#' return value. Return value contains actual correlation-matrix as 'vcor'
#' attribute.
#'
#' ToDo: maybe some additional validation steps and more informative output.
#'
#' @keywords collinearity model validation
#' @param model.fit lm() or lmm() output
#' @return Maximum amount of correlation for given model variables.
#' @export
#'
#' @examples
#' # This will return the maximum colinearity score in the given model
#' limimo <- lme4::lmer(dummy.list$cnts[,1] ~ group + (1|batch),
#' data=dummy.list$meta)
#' num.max_corr <- colinScore(model.fit=limimo)
colinScore <- function(model.fit) {
  # get variance-covariance matrix
  V <- stats::vcov(model.fit)
  # 'NA' values or non-square matrix
  if (any(is.na(V)) || nrow(V) == 0) {
    score <- ifelse(any(is.na(V)), 1, 0)
    attr(score, "vcor") <- NA
  } else {
    # scale to correlation
    C <- stats::cov2cor(as.matrix(V))
    # get largest correlation
    score <- max(abs(C[lower.tri(C)]))
    attr(score, "vcor") <- C
  }
  return(score)
}
