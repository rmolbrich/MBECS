# ANALYSIS FUNCTIONS ------------------------------------------------------


#' Relative Log Expression Plot
#'
#' Takes two covariates, i.e., group and batch, and computes the RLE-plot over the grouping of the
#' first covariate, colored by the second covariate. Effectively illustrating the relative
#' expression between samples from different batches within the respective study groups. Other
#' covariates can be chosen as input and the function will check for factors and convert if
#' necessary. Categorical factors, e.g., group membership, sex and batch, produce the best result.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input for th data-set can
#' be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data. The
#' covariate table requires an 'sID' column that contains sample IDs equal to the sample naming in
#' the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords RLE relative log expression
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param model.vars two covariates of interest to select by first variable selects panels and second one determines coloring
#' @param return.data logical if TRUE returns the data.frame required for plotting (NO plotting here bucko)
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#'
#' @examples
#' This will return the data.frame for plotting.
#' \dontrun{p.RLE <- mbecRLE(input.obj=list(counts, covariates),
#' model.vars=c("treatment","batches"), return.data=TRUE)}
#'
#' This will return the ggplot2 object for display, saving and modification.
#' \dontrun{p.RLE <- mbecRLE(input.obj=phyloseq, model.vars=c("treatment","sex"),
#' return.data=FALSE)}
mbecRLE <- function(input.obj, model.vars=c("group","batch"), return.data=FALSE) {

  cols <- pals::tableau20(20)

  ## check for correct inputs
  tmp <- mbecGetData(input.obj=input.obj, orientation="fxs", required.col=eval(model.vars))
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]] %>% rownames_to_column(., var = "specimen")

  ## SPLIT - COMPUTE - MERGE
  tmp.long <- NULL
  for( g.idx in unique(tmp.meta[,eval(model.vars[1])]) ) {
    message(paste("Calculating RLE for group: ",g.idx,sep=""))

    tmp.cnts.group <- dplyr::select(tmp.cnts, tmp.meta$specimen[tmp.meta[,eval(model.vars[1])] %in% g.idx])

    # Median per feature in this group
    feature.med = apply(tmp.cnts.group, 1, median)

    tmp.group.long <- apply(tmp.cnts.group, 2, function(sample.col) sample.col - feature.med) %>% # subtract feature-median from sample
      as.data.frame(.) %>%
      tidyr::pivot_longer(., cols = colnames(.), names_to = "specimen", values_to = "values") # re-arrange to long format

    tmp.long <- rbind.data.frame(tmp.long, tmp.group.long)
  }

  # the factor levels of samples need to be in the same order as batches per group
  # order by group_batch and then set factor levels for samples and batches in that order
  tmp.long <- dplyr::left_join(tmp.long, tmp.meta, by = "specimen") %>% # merge with sample data
    dplyr::mutate(plot.order = paste(get(model.vars[1]), get(model.vars[2]), sep="_")) %>%
    dplyr::arrange(., plot.order) %>%
    dplyr::mutate(specimen = factor(specimen, levels = unique(specimen)))

  if( return.data ) {
    return(tmp.long)
  } # else create the plot

  rle.plot <- ggplot2::ggplot(tmp.long, ggplot2::aes(x = specimen, y = values, fill = get(model.vars[2]))) +
    ggplot2::stat_boxplot(color="black",notch = TRUE,
                 outlier.colour = "#E42032", outlier.fill = "white",outlier.shape = 1, outlier.stroke = .5) +
    #facet_wrap(~Strain, ncol=2) +
    ggplot2::facet_grid(cols=vars(get(model.vars[1])), scales="free", space="free_x", drop=T) +
    ggplot2::scale_fill_manual(values = cols) +
    theme_rle() +
    ggplot2::guides(fill=guide_legend(title=element_blank()))

  return(rle.plot)
}



#' Principal Component Analysis Plot
#'
#' Takes two covariates, i.e., group and batch, and computes the ordination-plot for user-selected
#' principal components. Covariates determine sample-shape and color and can be switched to shift
#' the emphasis on either group. In addition to the ordination-plot, the function will show the
#' distribution of eigenvalues (colored by the second covariate) on their respective principal
#' components.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input for th data-set
#' can be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data.
#' The covariate table requires an 'sID' column that contains sample IDs equal to the sample naming
#' in the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords PCA principal component analysis
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param model.vars two covariates of interest to select by first variable selects shape and second one determines coloring
#' @param pca.axes numeric vector which axes to plot, first is X and second is Y
#' @param return.data logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' This will return the data.frame for plotting.
#' \dontrun{p.PCA <- mbecPCA(input.obj=list(counts, covariates),
#' model.vars=c("treatment","batches"), pca.axes=c(1,2), return.data=TRUE)}
#'
#' This will return the ggplot2 object for display, saving and modification. Selected PCs are PC3 on
#' x-axis and PC2 on y-axis.
#' \dontrun{p.PCA <- mbecPCA(input.obj=list(counts, covariates),
#' model.vars=c("treatment","batches"), pca.axes=c(3,2), return.data=FALSE)}
setGeneric("mbecPCA", signature="input.obj",
           function(input.obj, model.vars=c("group","batch"), pca.axes=c(1,2), return.data=FALSE)
             standardGeneric("mbecPCA")
)

## In this form it works for 'phyloseq' and 'mbecData' objects
.mbecPCA <- function(input.obj, model.vars=c("group","batch"), pca.axes=c(1,2), return.data=FALSE) {

  tmp <- mbecGetData(input.obj, orientation="sxf", required.col=eval(model.vars))
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  # calculate IQR and sort counts in decreasing order
  iqr <- apply(tmp.cnts,2,stats::IQR)
  tmp.cnts <- tmp.cnts[,order(iqr,decreasing=T)]

  PCA <- stats::prcomp(tmp.cnts, scale = F)

  axes.number <- dim(PCA$x)[2]
  axes.names <- paste("PC",1:axes.number, sep="")

  plot.df <- PCA$x %>%
    data.frame(., stringsAsFactors = FALSE) %>%
    dplyr::rename_at(1:axes.number, ~ axes.names) %>%
    tibble::rownames_to_column(., var = "sample") %>%
    dplyr::left_join(., tmp.meta, by = "sample")

  metric.df <- data.frame("var.explained"=round((100*PCA$sdev^2) / (sum(PCA$x^2 / max(1, nrow(PCA$x)-1))),2),row.names = axes.names) %>%
    dplyr::mutate("axis.min"=floor(apply(PCA$x, 2, function(col) min(col))))%>%
    dplyr::mutate("axis.max"=ceiling(apply(PCA$x, 2, function(col) max(col))))

  for( g.idx in c(1:length(model.vars)) ) {
    if( !is.factor(plot.df[,eval(model.vars[g.idx])]) ) {
      warning("Grouping variables need to be factors. Coercing to factor now, adjust beforehand to get best results.")
      plot.df[,eval(model.vars[g.idx])] <- factor(plot.df[,eval(model.vars[g.idx])])
    }
  }

  # No plotting, just return data.
  if( return.data ) {
    return(list(plot.df, metric.df, pca.axes))
  }

  # Release the Ploten
  ### PCA-PLOT
  # change first letter of group denominator to uppercase for pretty-plotting
  var.one <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.vars[1],perl = TRUE)
  var.two <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.vars[2],perl = TRUE)

  var.shape <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.vars[1],perl = TRUE)
  var.color <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.vars[2],perl = TRUE)

  title <- paste("PCA:",var.shape, "-", var.color)

  if( length(model.vars) >= 2 ){
    pMain <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[1]+1])), y = get(colnames(plot.df[pca.axes[2]+1])), colour = get(model.vars[2]), shape = get(model.vars[1]))) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::labs(colour = var.color, shape = var.shape) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      ggplot2::ylim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      ggplot2::xlab(paste0(colnames(plot.df[pca.axes[1]+1]), ': ', metric.df$var.explained[pca.axes[1]], '% expl.var')) +
      ggplot2::ylab(paste0(colnames(plot.df[pca.axes[2]+1]), ': ', metric.df$var.explained[pca.axes[2]], '% expl.var')) +
      theme_pca()

    pTop <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[1]+1])), fill = get(model.vars[2]), linetype = get(model.vars[2]))) +
      ggplot2::geom_density(size = 0.2, alpha = 0.5) + ylab('Density') +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      theme_pca() +
      ggplot2::labs(title = title) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = element_text(size = rel(0.8)),
                     plot.title = element_text(hjust = 0.5, size = rel(1.5)))

    pRight <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x=get(colnames(plot.df[pca.axes[2]+1])), fill = get(model.vars[2]), linetype = get(model.vars[2]))) +
      ggplot2::geom_density(size = 0.2,alpha = 0.5) +  coord_flip() + ylab('Density') +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      theme_pca() +
      ggplot2::theme(axis.title.x = element_text(size = rel(0.8)),
                     axis.title.y = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     plot.title = ggplot2::element_blank())

  }else{
    pMain <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[1]+1])), y = get(colnames(plot.df[pca.axes[2]+1])), colour = get(model.vars[1]))) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::labs(colour = var.color, shape = var.shape) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      ggplot2::ylim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      ggplot2::xlab(paste0(colnames(plot.df[pca.axes[1]+1]), ': ', metric.df$var.explained[pca.axes[1]], '% expl.var')) +
      ggplot2::ylab(paste0(colnames(plot.df[pca.axes[2]+1]), ': ', metric.df$var.explained[pca.axes[2]], '% expl.var')) +
      theme_pca()

    pTop <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[1]+1])), fill = get(model.vars[1]))) +
      ggplot2::geom_density(size = 0.2, alpha = 0.5) + ylab('Density') +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      theme_pca() +
      ggplot2::labs(title = title) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = element_text(size = rel(0.8)),
                     plot.title = element_text(hjust = 0.5, size = rel(1.5)))

    pRight <- ggplot(data = plot.df, ggplot2::aes(x=get(colnames(plot.df[pca.axes[2]+1])), fill = get(model.vars[1]))) +
      ggplot2::geom_density(size = 0.2,alpha = 0.5) +  coord_flip() + ylab('Density') +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      theme_pca() +
      ggplot2::theme(axis.title.x = element_text(size = rel(0.8)),
                     axis.title.y = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                     plot.title = ggplot2::element_blank())
  }

  g <- ggplot2::ggplotGrob(pMain)$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

  ret.plot <- gridExtra::grid.arrange(pTop + ggplot2::theme(legend.position = 'none'), legend, pMain +
                                        ggplot2::theme(legend.position = 'none'), pRight + ggplot2::theme(legend.position = 'none'),
                                      ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))

  return(ret.plot)

}

#' @include mbecs_classes.R
setMethod("mbecPCA", "MbecData",
          function(input.obj, model.vars=c("group","batch"), pca.axes=c(1,2), return.data=FALSE) {
            .mbecPCA(input.obj, model.vars=model.vars, pca.axes=pca.axes, return.data=return.data)
          }
)

setMethod("mbecPCA", "phyloseq",
          function(input.obj, model.vars=c("group","batch"), pca.axes=c(1,2), return.data=FALSE) {
            .mbecPCA(input.obj, model.vars=model.vars, pca.axes=pca.axes, return.data=return.data)
          }
)

setMethod("mbecPCA", "list",
          function(input.obj, model.vars=c("group","batch"), pca.axes=c(1,2), return.data=FALSE) {
            .mbecPCA(input.obj, model.vars=model.vars, pca.axes=pca.axes, return.data=return.data)
          }
)


#' Feature Differential Abundance Box-Plot
#'
#' Displays the abundance of a selected feature, grouped/colored by a covariate, i.e., batch, in a
#' box-plot. Includes the density-plot, i.e., the distribution of counts for each sub-group.
#' Selection methods for features are 'TOP' and 'ALL' which select the top-n or all features
#' respectively. The default value for n is 10 and can be changed with the accompanying parameter.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input for th data-set can
#' be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data. The
#' covariate table requires an 'sID' column that contains sample IDs equal to the sample naming in
#' the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Box abundance density
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param method one of 'ALL' or 'TOP' for 'n' most variable features, DEFAULT is 'ALL'
#' @param n number of OTUs to display for 'TOP' method
#' @param model.var covariate to group by, default is batch
#' @param return.data logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' This will return the plot-frame of all features i the data-set.
#' \dontrun{p.Box <- mbecBox(input.obj=list(counts, covariates), method="ALL", model.var="batch",
#' return.data=TRUE)}
#'
#' This will return the ggplot2 object of the top 15 most variable features.
#' \dontrun{p.Box <- mbecBox(input.obj=list(counts, covariates), method="TOP", n=15,
#' model.var="batch", return.data=FALSE)}
mbecBox <- function(input.obj, method=c("ALL","TOP"), n=10, model.var="batch", return.data=FALSE) {

  cols <- cols[c(1,3,5,7,9,11,13,15,17,19)]

  # ## 00. Init and test
  # method <- match.arg(method)
  # message(paste("Selection method is: ",method,n, sep=""))

  # needs sxf orientation
  tmp <- mbecGetData(input.obj, orientation="sxf", required.col=eval(model.var))
  tmp[[2]] <- rownames_to_column(tmp[[2]], var = "specimen")

  otu.idx <- colnames(tmp[[1]])

  tmp <- tmp[[1]] %>%
    tibble::rownames_to_column(var = "specimen") %>%
    dplyr::left_join(., tmp[[2]], by = c("specimen" = "specimen"))

  if( method[1] == "TOP" ) {
    # calculate IQR and order from largest to smallest
    iqr <- apply(tmp[,otu.idx],2,stats::IQR)
    iqr <- iqr[order(iqr, decreasing=TRUE)]
    otu.idx <- names(iqr)[1:min(length(otu.idx), n)]

    tmp <- tmp %>%
      dplyr::select(., c(dplyr::all_of(otu.idx), "specimen", eval(model.var)))

  } else if( length(method) >= 2 ) {
    message("'Method' parameter contains multiple elements - using to select features.")
    # calculate IQR and sort as well
    otu.idx <- method
    tmp <- tmp %>%
      dplyr::select(., c(dplyr::all_of(otu.idx), "specimen", eval(model.var)))

  } # else is 'select-all-mode'

  if( return.data ) {
    return(list(tmp, otu.idx))
  }

  ## Prepare the plots for selected features
  # upper case first letter for the legend box
  legend.title <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.var,perl = TRUE)
  ret.plot <- list()

  for( idx in otu.idx ) {
    p.box <- ggplot2::ggplot(data = tmp, ggplot2::aes(x = get(model.var), y = get(idx), fill = get(model.var))) + ggplot2::stat_boxplot(geom = "errorbar", width = 0.4) +
      ggplot2::geom_boxplot() + ggplot2::scale_fill_manual(values = cols) + theme_bw() +
      theme_box() +
      ggplot2::labs(fill = legend.title, y = 'value',title = idx)

    p.density <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(idx), fill = get(model.var))) +
      ggplot2::geom_density(alpha = 0.5) + ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(title = idx, x = 'Value', fill = legend.title) +
      theme_box()

    ## Put the plots in grid for plotting
    # modify legend
    g <- ggplot2::ggplotGrob(p.box)$grobs
    # extract legend
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

    # put plot into the list
    ret.plot[[eval(idx)]] <- gridExtra::arrangeGrob(p.box + ggplot2::theme(legend.position = 'none'),
                                                    p.density + ggplot2::theme(legend.position = 'none', plot.title = ggplot2::element_blank()),
                                                    legend,
                                                    ncol = 1, nrow = 3, heights = c(5,4.5,1))
  }

  return(ret.plot)
}


#' Feature Differential Abundance Heatmap
#'
#' Shows the abundance value of selected features in a heatmap. By default, the function expects two
#' covariates group and batch to depict clustering in these groups. More covariates can be included.
#' Selection methods for features are 'TOP' and 'ALL' which select the top-n or all features
#' respectively. The default value for n is 10 and can be changed with the accompanying parameter.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input for the data-set
#' can be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data.
#' The covariate table requires an 'sID' column that contains sample IDs equal to the sample naming
#' in the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Heat abundance clustering
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param model.vars covariates of interest to show in heatmap
#' @param center flag to activate centering
#' @param scale flag to activate scaling
#' @param method one of 'ALL' or 'TOP' or a vector of feature names
#' @param n number of features to select in method TOP
#' @param return.data logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' This will return the plot-frame of all features i the data-set.
#' \dontrun{p.Heat <- mbecHeat(input.obj=phyloseq.obj, model.vars=c("group","batch"), center=TRUE,
#' scale=TRUE, method="ALL", return.data=TRUE)}
#'
#' This will return the ggplot2 object of the top 15 most variable features.
#' \dontrun{p.Heat <- mbecHeat(input.obj=list(counts, covariates), model.vars=c("group","batch"),
#' center=TRUE, scale=TRUE, method="TOP", n=15, return.data=FALSE)}
mbecHeat <- function(input.obj, model.vars=c("group","batch"), center=TRUE, scale=TRUE, method="TOP", n=10, return.data=FALSE) {

  cols <- cols[c(1,3,5,7,9,11,13,15,17,19)]

  ## ToDo: adjust legend settings
  ## ToDo: for phyloseq select taxonomic level?!

  ## needs sxf orientation and the feature names
  tmp <- mbecGetData(input.obj, orientation="sxf", required.col=eval(model.vars))
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]
  otu.idx <- colnames(tmp[[1]])

  ## CHECK if model.vars are factors
  for( g.idx in c(1:length(model.vars)) ) {
    if( !is.factor(tmp.meta[,eval(model.vars[g.idx])]) ) {
      warning(paste0("Grouping variables need to be factors. Coercing variable: ", eval(model.vars[g.idx]), " to factor now, adjust beforehand to get best results."))
      tmp.meta[,eval(model.vars[g.idx])] <- factor(tmp.meta[,eval(model.vars[g.idx])])
    }
  }

  ## center & scale
  tmp.cnts <- base::scale(tmp.cnts, center = eval(center), scale = eval(scale))
  tmp.cnts <- base::scale(t(tmp.cnts), center = eval(center), scale = eval(scale))

  if( method[1] == "TOP" ) {
    # calculate IQR and order from largest to smallest
    iqr <- apply(tmp.cnts[otu.idx,],1,stats::IQR)
    iqr <- iqr[order(iqr, decreasing=TRUE)]
    otu.idx <- names(iqr)[1:min(length(otu.idx), n)]
    # select only wanted features
    tmp.cnts <- tmp.cnts[otu.idx,]
  } else if( length(method) >= 2 ) {
    message("'Method' parameter contains multiple elements - using to select features.")
    # calculate IQR and sort as well
    tmp.cnts <- tmp.cnts[method,]

  } # else is 'select-all-mode'

  if( return.data ) {
    return(list(tmp.cnts, tmp.meta))
  }

  p.title <- paste("Heatmap - Centered: ", center, " Scaled: ", scale, sep="")
  heat.plot <- pheatmap::pheatmap(tmp.cnts,
                                  scale = 'none',
                                  cluster_rows = F,
                                  cluster_cols = T,
                                  fontsize_row = 4, fontsize_col = 6,
                                  fontsize = 8,
                                  clustering_distance_rows = 'euclidean',
                                  clustering_method = 'ward.D',
                                  treeheight_row = 30,
                                  annotation_col = tmp.meta[,eval(model.vars)],
                                  #annotation_colors = ad.anno_metabo_colors,
                                  border_color = 'NA',
                                  main = p.title)

  return(heat.plot)
}


#' Mosaic Sample Group Allocation
#'
#' Depicts the dispersion of samples over two (preferentially categorical*) covariates of interest.
#' Effectively showing, the un-/evenness within and between covariates to inform the choice of
#' methods for the subsequent steps in an analysis.
#'
#' The function returns either a plot-frame or the finished ggplot object. Input for the data-set
#' can be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data.
#' The covariate table requires an 'sID' column that contains sample IDs equal to the sample naming
#' in the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Mosaic sample allocation
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param model.vars covariates of interest to the sample allocation
#' @param return.data logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' This will return the plot-frame of for the samples grouped by treatment and sex
#' \dontrun{p.Mosaic <- mbecMosaic(input.obj=phyloseq.obj, model.vars=c("treatment","sex"),
#' return.data=TRUE)}
#'
#' This will return the ggplot2 object of the samples grouped by group and batch
#' \dontrun{p.Mosaic <- mbecMosaic(input.obj=list(counts, covariates),
#' model.vars=c("group","batch"), return.data=FALSE)}
mbecMosaic <- function(input.obj, model.vars=c("group","batch"), return.data=FALSE) {

  cols <- pals::tableau20(20)
  ## ToDo: adjust legend settings
  ## ToDo: for phyloseq select taxonomic level?!

  ## needs sxf orientation
  tmp <- mbecGetData(input.obj, orientation="sxf", required.col=eval(model.vars))
  tmp.meta <- tmp[[2]]


  if( length(model.vars) < 2 ) {
    message("Only one variable specified for Mosaic-plot, two are required!")
  } else if( length(model.vars) > 2 ) {
    message("More than two variables specified. Mosaic will take the first two.")
  }

  ## CHECK if grouping variables are factors
  for( g.idx in eval(model.vars) ) {
    if( !is.factor(tmp.meta[,eval(g.idx)]) ) {
      warning(paste0("Grouping variables need to be factors. Coercing variable: ", eval(g.idx), " to factor now, adjust beforehand to get best results."))
      tmp.meta[,eval(g.idx)] <- as.factor(tmp.meta[,eval(g.idx)])
    }
  }

  # how many samples are there
  n.observations <- base::dim(tmp.meta)[1]

  study.summary <- base::table(tmp.meta[,eval(model.vars[1])], tmp.meta[,eval(model.vars[2])]) %>%
    as.data.frame(.) %>%
    mutate("Freq.scaled"=Freq / n.observations)

  # prepare plot-annotation
  vars.axes <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.vars,perl = TRUE)

  if( return.data ) {
    return(study.summary)
  } # else return the plots

  # split by batch
  plot.v2 <- ggplot2::ggplot(study.summary, ggplot2::aes(x = Var1, y= Freq.scaled, group = Var2, fill=Var1)) +
    ggplot2::facet_grid(cols=vars(Var2), scales="free", space="free_x", drop=T) +
    ggplot2::geom_bar(stat = "identity", width = 0.9) +
    ggplot2::guides(fill = guide_legend(title=eval(vars.axes[1]), reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ggplot2::ylab("Proportion of all observations") +
    theme_mosaic(legend_position = "bottom")

  # split by treatment
  plot.v1 <- ggplot2::ggplot(study.summary, ggplot2::aes(x = Var2, y= Freq.scaled, fill=Var2)) +
    ggplot2::facet_grid(cols=vars(Var1), scales="free", space="free_x", drop=T) +
    ggplot2::geom_bar(stat = "identity", width = 0.9) +
    ggplot2::guides(fill = ggplot2::guide_legend(title=eval(vars.axes[2]), reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ggplot2::ylab("Proportion of all observations") +
    theme_mosaic()

  mosaic.plot <- gridExtra::grid.arrange(plot.v2, plot.v1, ncol=1, nrow=2, heights=c(1,1))

  return(mosaic.plot)

}


# VARIANCE CALCULATION ----------------------------------------------------


#' Estimate Explained Variance
#'
#' The function offers a selection of methods/algorithms to estimate the proportion of variance that
#' can be attributed to covariates of interest. This shows, how much variation is explained by the
#' treatment effect, which proportion is introduced by processing in batches and the leftover
#' variance, i.e., residuals that are not currently explained. Covariates of interest (CoI) are
#' selected by the user and the function will incorporate them into the model building for the
#' respective algorithm. The user can select from five different approaches to adapt to the
#' characteristics of the data-set, e.g., LMMs are a better choice than LMs for a very unbalanced
#' study design. Available approaches are: Linear Model (lm), Linear Mixed Model (lmm),
#' Redundancy Analysis (rda), Principal Variance Component Analysis (pvca) or
#' Silhouette Coefficient (s.coef).
#'
#' Linear Model (lm): An additive model of all covariates is fitted to each feature respectively
#' and the proportion of variance is extracted for each covariate
#' (OTU_x ~ covariate_1 + covariate_2 + ...).
#'
#' Linear Mixed Model (lmm): All but the first covariate are considered mixed effects. A model is
#' fitted to each OTU respectively and the proportion of variance extracted for each covariate
#' (OTU_x ~ covariate_1 + (1|covariate_2) + (1|...)).
#'
#' partial Redundancy Analysis (rda): Iterates over given covariates, builds a model of all
#' covariates that includes one variable as condition/constraint and then fits it to the feature
#' abundance matrix. The difference in explained variance between the full- and the constrained-
#' model is then attributed to the constraint.
#' (cnts ~ group + Condition(batch) vs. cnts ~ group + batch)
#'
#' Principal Variance Component Analysis (pvca): Algorithm - calculate the correlation of the fxs
#' count-matrix - from there extract the eigenvectors and eigenvalues and calculate the proportion
#' of explained variance per eigenvector (i.e. principal component) by dividing the eigenvalues by
#' the sum of eigenvalues. Now select as many PCs as required to fill a chosen quota for the total
#' proportion of explained variance. Iterate over all PCs and fit a linear mixed model that contains
#' all covariates as random effect and all unique interactions between two covariates. Compute
#' variance covariance components form the resulting model --> From there we get the Variance that
#' each covariate(variable) contributes to this particular PC. Then just standardize variance by
#' dividing it through the sum of variance for that model. Scale each PCs results by the proportion
#' this PC accounted for in the first place. And then do it again by dividing it through the total
#' amount of explained variance, i.e. the cutoff to select the number of PCs to take (obviously
#' not the cutoff but rather the actual values for the selected PCs). Finally take the average
#' over each random variable and interaction term and display in a nice plot.
#'
#' Silhouette Coefficient (s.coef): Calculate principal components and get sample-wise distances on
#' the resulting (sxPC) matrix. Then iterate over all the covariates and calculate the cluster
#' silhouette (which is basically either zero, if the cluster contains only a single element, or it
#' is the distance to the closest different cluster minus the distance of the sample within its own
#' cluster divided (scaled) by the maximum distance). Average over each element in a cluster for all
#' clusters and there is the representation of how good the clustering is. This shows how good a
#' particular covariate characterizes the data, i.e., a treatment variable for instance may
#' differentiate the samples into treated and untreated groups which implies two clusters. In an
#' ideal scenario, the treatment variable, i.e., indicator for some biological effect would produce
#' a perfect clustering. In reality, the confounding variables, e.g., batch, sex or age, will also
#' influence the ordination of samples. Hence, the clustering coefficient is somewhat similar to the
#' amount of explained variance metric that the previous methods used. If used to compare an
#' uncorrected data-set to a batch-corrected set, the expected result would be an increase of
#' clustering coefficient for the biological effect (and all other covariates - because a certain
#' amount of uncertainty was removed from the data) and a decrease for the batch effect.
#'
#' The function returns a data-frame for further analysis - the report functions
#' (mbecReport and mbecReportPrelim) will automatically produce plots. Input for the data-set can
#' be an MbecData-object, a phyloseq-object or a list that contains counts and covariate data. The
#' covariate table requires an 'sID' column that contains sample IDs equal to the sample naming in
#' the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords Model Evaluation Variance
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handled internally)
#' @param model.vars vector of covariates to include in model-construction, in case parameter 'model.form' is not supplied
#' @param method select method of modeling: Linear Model (lm), Linear Mixed Model (lmm), Redundancy Analysis (rda), Principal Variance Component Analysis (pvca) or Silhouette Coefficient (s.coef)
#' @param model.form string that describes a model formula, i.e., "y ~ covariate1 + (1|covariate2)"
#' @param type creates a column with that string in the output df - to keep track of cnt-source
#' @return df that contains proportions of variance for given covariates in every feature
#' @include mbecs_classes.R
#'
#' @examples
#' This will return a data-frame that contains the variance attributable to group and batch
#' according to linear additive model.
#' \dontrun{df.var.lm <- mbecModelVariance(input.obj=phyloseq.obj, model.vars=c("group","batch"),
#' method="lm", type="RAW")}
#' This will return a data-frame that contains the variance attributable to group and batch
#' according to linear additive model.
#' \dontrun{df.var.pvca <- mbecModelVariance(input.obj=phyloseq.obj, model.vars=c("group","batch"),
#' method="pvca")}
mbecModelVariance <- function( input.obj, model.vars=character(), method=c("lm","lmm","rda","pvca"), model.form=NULL, type="NONE", no.warning=TRUE) {

  ### ToDo: selection cutoff for PCs in silhouette coefficient method?!
  ### ToDo: safety checks and logic to distinguish model types and also take care of this matrix-input issue
  ### ToDoAsWell: How2Make lm and lmm formulas.. or if  or whatever
  ## ToDo: check this out - itis part of lmm just put it her to remind me
  ## control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )
  oldw <- getOption("warn")
  if( no.warning ) {
    options(warn = -1)
  }
  on.exit(options(warn = oldw))

  ## PVCA stuff
  pct_threshold = .5876   # threshold for explained variances

  tmp <- mbecGetData(input.obj, orientation="sxf", required.col=eval(model.vars))
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  ## CHECK if grouping variables are factors
  for( g.idx in eval(model.vars) ) {
    if( !is.factor(tmp.meta[,eval(g.idx)]) ) {
      warning(paste0("Grouping variables need to be factors. Coercing variable: ", eval(g.idx), " to factor now, adjust beforehand to get best results."))
      tmp.meta[,eval(g.idx)] <- as.factor(tmp.meta[,eval(g.idx)])
    }
  }

  ## 1. implement for linear model
  if( method == "lm" ) {
    message("Fitting linear model to every feature and extract proportion of variance explained by covariates.")
    # if no formula is supplied a simple additive model will be constructed
    if( !is.null(model.form) ) {
      message("Use provided model formula.")
      # tmp.formula <- as.formula(deparse(model.form)) this would make a function work, but fails for string --> just string for now
      tmp.formula <- as.formula(model.form)
    } else {
      message("Construct formula from covariates.")
      tmp.formula = stats::as.formula(paste("y", " ~ ",
                                            paste(model.vars, collapse=" + ")))
    }
    ## so, for some reason this shit does not work with apply, because the iterator will not update in the formula
    ## maybe try again once the rest works
    ### DEBUG ###
    # model.variances <- sapply(colnames(tmp.cnts), FUN = function(x) {
    #   model.fit <- stats::lm(tmp.cnts[,eval(x)] ~ group + batch, data = tmp.meta)
    #   # data.frame(mbecVarianceStats(model.fit))
    #   data.frame(mbecVarianceStats(model.fit))
    # })
    ### DEBUG ###

    model.variances <- NULL
    for( x in colnames(tmp.cnts)) {
      y <- tmp.cnts[[eval(x)]]

      model.fit <- stats::lm(tmp.formula, data=tmp.meta)
      model.variances <- rbind.data.frame(model.variances, mbecVarianceStats(model.fit))
    }

    modelType = "anova"

  } else if( method == "lmm" ) {
    message("Fitting linear-mixed model to every feature and extract proportion of variance explained by covariates.")

    control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )

    # if no formula is supplied a simple additive model will be constructed
    if( !is.null(model.form) ) {
      message("Use provided model formula.")
      # tmp.formula <- as.formula(deparse(model.form)) this would make a function work, but fails for string --> just string for now
      tmp.formula <- as.formula(model.form)
    } else {
      message("Construct formula from covariates.")
      f.terms <- paste("(1|",model.vars,")", sep="")
      tmp.formula <- as.formula(paste(paste("y", model.vars[1], sep=" ~ "), paste(f.terms[-1], collapse=" + "), sep=" + "))

    }
    # now fit a model to every feature
    model.variances <- NULL
    for( x in colnames(tmp.cnts)) {
      # meh
      y <- tmp.cnts[[eval(x)]]

      model.fit <- lme4::lmer(tmp.formula, data=tmp.meta)
      model.variances <- rbind.data.frame(model.variances, mbecVarianceStats(model.fit))
    }

    modelType = "linear-mixed model"

  } else if( method == "rda" ) {
    # pre-analysis check-up
    # 1. response var: dimensional homogeneity (equal units of measurement) --> counts should be like that (else center and standardise)
    # 2. #explanatory variables << #samples --> else overdetermined system
    # 3. explanatory var: dimensional homogeneity for quantitative variables --> else center & standard (qualitative vars work anyways)
    # 4. Examine the distribution of each variable in you explanatory and response matrix as well as plots of each variable against other
    # variables in its own and any other matrix. If the relationships are markedly non-linear, apply transformations to linearise the
    # relationships and reduce the effect of outliers.
    # 5. If you wish to represent non-Euclidean relationships (e.g. Hellinger distances) between objects in an RDA ordination, you should
    # apply an ecologically-motivated transformation discussed on this page before analysis.

    ### perform as much RDAs as there are covariates, while conditioning on a single one each time.
    # set-up df for results
    model.variances = data.frame(matrix(nrow = length(model.vars), ncol = 1, dimnames = list(model.vars,eval(type))))

    # get significance and variance for whole model
    # tmp.formula = stats::as.formula(paste("tmp.cnts", " ~ ",
    #                                       paste(model.vars, collapse=" + "),sep=""))
    # tmp.rda.covariate <- vegan::rda(tmp.formula, data=tmp.meta)
    # tmp.sig <- anova(tmp.rda.covariate,permutations = how(nperm=999))
    # full.model.variance <- paste((tmp.sig$Variance[1] / sum(tmp.sig$Variance)), tmp.sig$`Pr(>F)`[1], sep="|")

    # iterate over the vector of covariate names (model.vars) to construct formulas for 'RDA' procedure
    for( condition.idx in 1:length(model.vars) ) {
      # the counts are always available in sxf format in variable 'tmp.cnts'
      # just iterate over all covariates and keep on to condition on
      tmp.formula = stats::as.formula(paste("tmp.cnts", " ~ ",
                                            paste(model.vars[-eval(condition.idx)], "+", collapse=" "),
                                            " Condition(", model.vars[eval(condition.idx)],")", sep=""))
      # compute RDA for this condition
      tmp.rda.covariate <- vegan::rda(tmp.formula, data=tmp.meta)

      # estimate significance of
      tmp.sig <- anova(tmp.rda.covariate,permutations = permute::how(nperm=999))$`Pr(>F)`[1]

      # calculate proportion of variance for this covariate as quotient of its eigenvalue and the sum of eigenvalues
      # partial.chi denotes the variance in the 'CONDITION'
      model.variances[condition.idx,eval(type)] <- summary(tmp.rda.covariate)$partial.chi*100/summary(tmp.rda.covariate)$tot.chi
      # model.variances[condition.idx,eval(type)] <- paste(summary(tmp.rda.covariate)$partial.chi*100/summary(tmp.rda.covariate)$tot.chi, tmp.sig, full.model.variance, sep="|")
    }

    modelType = "rda"
    # make-up some return format
    res <- data.frame(t(model.variances)) %>% mutate(type = eval(type))
    attr(res, "modelType") <- modelType

    return(res)


  } else if( method == "pvca" ) {

    #model.vars <- c("Time","Treatment","Batch")
    n.vars <- length(model.vars)
    s.names <- rownames(tmp.cnts)

    # center and/or scale the counts - if both arguments are set to false, it will just transpose the matrix for the subsequent steps
    tmp.cnts = apply(tmp.cnts, 2, scale, center = T, scale = FALSE) %>%
      t(.)

    #reset the sample-names
    colnames(tmp.cnts) <- s.names

    # calculate correlation - fxs
    tmp.cor <- cor(tmp.cnts)

    # get eigen-values/vectors of the correlation matrix
    eCor <- eigen(tmp.cor)
    eVal <- eCor$values
    eVec <- eCor$vectors
    rownames(eVec) <- s.names
    n.eVal <- length(eVal)
    sum.eVal <- sum(eVal)
    prop.PCs <- eVal/sum.eVal

    # figure out how many PCs to model by summing up the variance they account for until a threshold is met
    # takes at minimum 3 PCs regardless of threshold.
    # calculate cumulative sum for vector of variances - evaluate which position is larger than the cutoff -
    # select maximum out of 3 and the smallest PC that meets the cutoff
    n.PCs <- max(3, min(which(sapply(cumsum(prop.PCs), function(x) x >= pct_threshold))))

    # create a long-df that contains a column for all selected eigen-vectors and the required covariates (which obviously repeat for all vectors)
    # join selected eigenvectors and covariate data in a single df
    lmm.df <- eVec %>% tibble::as_tibble(., .name_repair = "unique") %>% select_at(1:eval(n.PCs)) %>% cbind(., tmp.meta)
    #lmm.df <- eVec %>% cbind(., tmp.meta)  %>% as_tibble(.) %>% select_at(1:eval(n.PCs))

    # # figure out how many effects there are - it is a random effect for every covariate - plus an interaction term for
    # # every combination of covariates - plus 1 for the residuals
    # n.effects <- length(model.vars) + (factorial(length(model.vars)) / (2 * (factorial(length(model.vars)-2)))) + 1

    # these are the basic random effects for all covariates
    f.terms <- paste("(1|",model.vars,")", sep="")

    # now add the random-interaction effects iterating through all but the last covariate and create a 'rid'
    # with all subsequent covariates - thus creating all combinations without repetition
    for(var.idx in 1:(n.vars-1)) {
      for(interaction.idx in (var.idx+1):(n.vars)) {
        # combine with the other formula terms and done
        f.terms <- c(f.terms, paste("(1|",model.vars[var.idx],":",model.vars[interaction.idx],")", sep=""))
      }
    }

    # OR just take the length of the f.terms vector +1 (for Residuals) to create a results.matrix
    n.effects <- length(f.terms) + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = n.PCs, ncol = n.effects)

    # now prepare the model formula, which only needs to be done once
    model.formula <- stats::as.formula(paste("lmm.df[,vec.idx]", " ~ ",paste(f.terms, collapse = " + "), sep=""))

    # for every selected eigen-vector, fit an lmm with random and interaction terms - calculate and extract the associated variances
    for( vec.idx in 1:n.PCs ) {
      randomEffects <- data.frame(lme4::VarCorr(Rm1ML <- lme4::lmer(model.formula, lmm.df, REML = TRUE, verbose = FALSE, na.action = na.omit)))
      # and put into result at the respective row
      randomEffectsMatrix[vec.idx,] <- as.numeric(randomEffects[,4])
    }

    # save effect-names and more importantly the actual order
    names.effects <- randomEffects[,1]

    # standardize variance by scaling every row by its sum
    randomEffectsMatrix.std <- randomEffectsMatrix / rowSums(randomEffectsMatrix)

    # weigh the values by the proportion of variance that the respective PCs accounted for
    # take the eigenvalues divided by sum of all eigenvalues and use them as weights for the
    # respective rows, i.e, effect-variances that belong to the PCs that the eigenvalues relate to
    scaled.eVal <- eVal / sum(eVal)
    randomEffectsMatrix.wgt <- randomEffectsMatrix.std * scaled.eVal[1:n.PCs]

    # at this point the total variance sums up to a value that should be close to the selected cutoff
    # for explained variance - so, now sum-up the variance for each effect, over the PCs and divide it by the
    # actual amount of variance that is associate with the selected PCs --> this puts it in the interval:[0,1]
    # for better comparability
    model.variances <- colSums(randomEffectsMatrix.wgt) / sum(colSums(randomEffectsMatrix.wgt))
    names(model.variances) <- names.effects

    modelType <- "PVCA - lmm"

  } else if( method == "s.coef" ) {

    ### prcomp and selection cutoff??!
    tmp.prcomp <- prcomp(tmp.cnts, center = TRUE, scale = FALSE)

    # get sample-wise distances
    tmp.dist <- dist(tmp.prcomp$x, method = "euclidian")

    # iterate over variables, calculate silhouette, compute avg. silhouette coefficient,
    # merge in neat little dataframe
    avg.sil.df <- NULL
    for( var.elem in model.vars ) {
      print(var.elem)

      tmp.sil = cluster::silhouette(x = as.numeric(tmp.meta[,eval(var.elem)]), dist = tmp.dist)

      avg.sil.df <- rbind.data.frame(avg.sil.df, data.frame("variable"=var.elem,
                                                            "cluster"=levels(tmp.meta[,eval(var.elem)]),
                                                            "sil.coefficient"= c(summary(tmp.sil))$clus.avg.widths))

    }

    res <- avg.sil.df %>% mutate(type = eval(type))
    attr(res, "modelType") <- "s.coef"

    return(res)
  }

  # transpose result and add column with correction/transformation type.. maybe more later though
  res <- mutate(model.variances, type = eval(type))
  attr(res, "modelType") <- modelType

  ### fancy return type from variancePart package - wait if we need that
  # res <- new("varPartResults", varPartMat, type=modelType, method="Variance explained (%)")

  return( res )
}


#' Model Variable Variance Extraction
#'
#' For a Linear (Mixed) Model, this function extracts the proportion of variance that can be
#' explained by terms and interactions and returns a named row-vector.
#'
#' Linear Model: Perform an analysis of variance (ANOVA) on the model.fit and return the
#' Sum of squares for each term, scaled by the total sum of squares.
#'
#' Linear Mixed Model: employ helper function 'mbecMixedVariance' to extract residuals,
#' random effects and fixed effects components from the model. The components are then
#' transformed to reflect explained proportions of variance for the model coefficients.
#' The function implements transformation for varying coefficients as well, but
#' NO ADJUSTMENT for single or multiple coefficients at this point.
#'
#' @keywords lm lmm proportion variance
#' @param model.fit linear (mixed) model object of class 'lm' or 'lmerMod'
#' @return a named row-vector, containing proportional variance for model terms
#' @export
#'
#' @examples
#' This will return the data.frame for plotting.
#' \dontrun{vec.variance <- mbecVarianceStats(model.fit=MyModel-obj)}
mbecVarianceStats <- function( model.fit ) {

  ### ToDo: implement glm versions of this

  # check validity of model fit
  mbecValidateModel( model.fit)

  # linear model
  if( class(model.fit) %in% "lm" ) {
    # get model coefficients
    # model.sum <- summary(model.fit)
    # w.cof2[i] <- model.sum$coefficients[3,1]

    vp = stats::anova(model.fit) %>%
      data.frame() %>%
      dplyr::mutate("variance" = dplyr::select(.,"Sum.Sq") / sum(dplyr::select(.,"Sum.Sq")), .keep="none") %>%
      t()

  } else if( class(model.fit) %in% "lmerMod" ) {  # linear-mixed model

    # this is the un-adjusted variance --> divide by total variance and done
    vc <- mbecMixedVariance(model.fit)

    # BUT maybe we want to adjust/condition for one or more variables..
    # and make this shit more complicated than it needs to be
    # adjust = NULL # some arbitrary adjustment for development
    #
    # # 1. make this adjustment validity check
    # # check for varying coefficient models - whatever that is.. copy and paste for now
    # if (max(sapply(varComp, length)) > 1 && !is.null(adjust)) {
    #   stop("The adjust and adjustAll arguments are not currently supported for varying coefficient models")
    # }
    # 2. for all covariates perform the variance calculation - for varying coefficients the total sum is somehow only this particular coefficient + all other variables - BUT the other coefficient level are ignored here
    # so, we need variances and scaled variance in subsequent steps - but only calculate once if possible

    # meh
    lib.df <- data.frame("covariates"=colnames(model.fit@frame),
                         row.names=sapply(colnames(model.fit@frame), function(covariate)
                           paste(paste(covariate, levels(model.fit@frame[,eval(covariate)]), sep=""), collapse = ",")))

    # and a look-up table to make variance calculations nice and easy
    total.var.LUT <- unlist(lapply(vc, function(var.comp) {
      if( length(var.comp) > 1 ) {
        weights = (table(model.fit@frame[[lib.df[eval(paste(names(var.comp), collapse = ",")),]]])/nrow(model.fit@frame))
        var.comp %*% weights
      } else {
        names(var.comp) <- NULL
        var.comp
      }
    } ), use.names = T, recursive = F)

    vp <- list()
    # now just sum total variance from LUT and leave out the current effect
    # the trick here is to  add the current effect later to the sum -->
    # in case it is a normal coefficient, it just works as a normal sum
    # but if it is a vector, i.e., varying coefficients, the separate addition produces a vector
    # instead of a scalar -> all the varying components can be updated in one fail-swoop *MUHAHAHAHA*
    for( effect in names(vc) ) {

      # value of this effect divided by the total variance
      tmp.t.var <- sum(total.var.LUT[-which(names(total.var.LUT) %in% eval(effect))])
      vp[[eval(effect)]] <- vc[[eval(effect)]] / (tmp.t.var + vc[[eval(effect)]])
    }

    vp <- unlist(vp)
    names(vp) <- gsub(".(Intercept)", replacement = "", names(vp), fixed=TRUE)
    vp <- data.frame(t(vp))

  } else if( class(model.fit) %in% "glm" ) {

    ### ToDon't
  }

  return(vp)
}


#' Mixed Model Variance-Component Extraction
#'
#' A helper function that extracts the variance components of linear mixed models, i.e., residuals,
#' random-effects, fixed-effects, scales them to sample-size and returns a list of components.
#'
#' Uses 'lme4::VarCorr' to extract Residuals and random-effects components. Standard Deviation of
#' Residuals is stored as 'sc' attribute in the output of 'VarCorr'.
#'
#' Uses 'lme4::fixef' to extract fixed-effects components, i.e., parameter estimates. The attribute
#' 'pp' of the model contains the dense model matrix for fixed-effects parameters (X). The fixed
#' effects variance, σ2f, is the variance of the matrix-multiplication β∗X (parameter vector by
#' model matrix)
#'
#' @keywords lmm proportion variance
#' @param model.fit linear mixed model object of class 'lmerMod'
#' @return a named list, containing proportional variance for model terms
#' @export
#'
#' @examples
#' This will return the data.frame for plotting.
#' \dontrun{list.variance <- mbecMixedVariance(model.fit=MyMixedModel-obj)}
mbecMixedVariance <- function(model.fit) {
  # remember: sd == sqrt(var)

  # VarCorr computes var, sd, cov for residuals and random effects
  rVC <- lme4::VarCorr(model.fit)
  # Standard Deviation of Residuals is stored as 'sc' attribute in rVC - and squared to get variances
  # iterate over the list(s) of random effects in rVC and squre sd to get variance - and don't forget to unlist()
  # so, technically this is the residualsRandomVar-vector
  # diag() is the same as attr(random.eff,"stddev")^2 and return the variances
  randomVar <- c("Residuals"=attr(rVC, "sc")^2,
                 lapply(rVC,diag))
  # scale to number of observations - aka (N-1)/N
  n.scaling <- (stats::nobs(model.fit) - 1) / stats::nobs(model.fit)
  # pp is class 'merPredD' and X is dense model matrix for the fixed-effects parameters
  # fixedEffect-value
  # with fixed effects parameters 'fixef()' extract the estimates for the fixed effects parameters
  # --> The fixed effects variance, σ2f, is the variance of the matrix-multiplication
  # β∗X (parameter vector by model matrix)
  #
  # extract both components and multiply each column by the corresponding value in the effects estimate
  # and then scale to sample-size, i.e., (N-1)/N
  fixedVar <- t(t(model.fit@pp$X) * lme4::fixef(model.fit)) %>%
    as.data.frame() %>%
    # 1. drop intercept
    dplyr::select(!"(Intercept)") %>%
    # 2. row-sums for all effects - can calculate total variance
    dplyr::mutate("total.var"=apply(., 1, sum)) %>%
    # 3. apply to get the scaled variance in each column
    apply(., 2, function(effect) var(effect) * n.scaling)

  # maybe fancy later, but for now this works
  fixedVar <- lapply(split(
    head(fixedVar, -1) / sum(head(fixedVar, -1)) * tail(fixedVar, 1),
    names(head(fixedVar, -1))),unname)

  # add attribute for total.var to every list element - to make the following steps less of a pita

  # concatenate and return the lists
  return(c(randomVar, fixedVar))
}


#' Validate Linear (Mixed) Models
#'
#' A helper function that calculates the collinearity between model variables and stops execution
#' if the maximum value is bigger than the allowed threshold.
#'
#' ToDo: maybe some additional validation steps and more informative output.
#'
#' @keywords collinearity model validation
#' @param model.fit lm() or lmm() output
#' @param colinearityThreshold cut-off for model rejection
#' @export
#'
#' @examples
#' This will return the data.frame for plotting.
#' \dontrun{mbecValidateModel(model.fit=MyMixedModel-obj, colinearityThreshold=0.999)}
mbecValidateModel <- function( model.fit, colinearityThreshold=0.999 ) {
  ## ToDo: health & Safety

  # implement for lmm as well
  if( class(model.fit) %in% "lm" ) {
    if( colinScore(model.fit) > colinearityThreshold ) {
      stop("Some covariates are strongly correlated. Try again.")
    }
  } else if( class(model.fit) %in% "lmerMod" ) {

    if( colinScore(model.fit) > colinearityThreshold ) {
      stop("Some covariates are strongly correlated. Try again.")
    }

    # ### DEVOPS
    #
    # # check that factors are random and continuous variables are fixed
    # ###################################################################
    #
    # # remove backticks with gsub manually
    # # solve issue that backticks are conserved is some but not all parts of lmer()
    #
    # # Simplified testing of random versus fixed effects
    # # allows (A|B) only where A is continuous
    #
    # # variables fit by regression
    # testVar = attr(attr(model.fit@frame, "terms"), "term.labels")
    # testVar = gsub("`", "", testVar)
    #
    # # get type for each variable
    # # keep only tested variables
    # varType = attr(attr(model.fit@frame, "terms"), "dataClasses")[-1]
    # varType = varType[testVar]
    #
    # # random effects
    # randVar = names(model.fit@flist)
    #
    # # fixed effects
    # # starting with all variables, remove random variables
    # fixedVar = setdiff(testVar, randVar)
    #
    # for( i in 1:length(varType) ){
    #
    #   # if factor is not random
    #   if( (showWarnings && ! dream) && varType[i] %in% c("factor", "character") && (! names(varType)[i] %in% randVar) ){
    #     stop(paste("Categorical variables modeled as fixed effect:", paste(names(varType)[i], collapse=', '), "\nThe results will not behave as expected and may be very wrong!!"))
    #   }
    #
    #   # If numeric/double is not fixed
    #   if( (showWarnings && ! dream) && varType[i] %in% c("numeric", "double") && (!names(varType)[i] %in% fixedVar) ){
    #     stop(paste("Continuous variable cannot be modeled as a random effect:", names(varType)[i]))
    #   }
    # }
    #
    # # show convergance message
    # if( showWarnings && !is.null(fit@optinfo$conv$lme4$messages) && (fit@optinfo$conv$lme4$messages != "boundary (singular) fit: see ?isSingular")){
    #   stop(fit@optinfo$conv$lme4$messages)
    # }

  }
}


#' Variable Correlation Linear (Mixed) Models
#'
#' Takes a fitted model and computes maximum correlation between covariates as return value.
#' Return value contains actual correlation-matrix as 'vcor' attribute.
#'
#' ToDo: maybe some additional validation steps and more informative output.
#'
#' @keywords collinearity model validation
#' @param model.fit lm() or lmm() output
#' @export
#'
#' @examples
#' This will return the data.frame for plotting.
#' \dontrun{num.max_corr <- colinScore(model.fit=MyMixedModel-obj)}
colinScore <- function(model.fit) {
  # get variance-covariance matrix
  V <- vcov(model.fit)
  # 'NA' values or non-square matrix
  if( any(is.na(V)) || nrow(V) == 0 ) {
    score <- ifelse(any(is.na(V)), 1,0)
    attr( score, "vcor") = NA
  } else {
    # scale to correlation
    C <- cov2cor(as.matrix(V))
    # get largest correlation
    score = max(abs(C[lower.tri(C)]))
    attr( score, "vcor") = C
  }
  return(score)
}


# VARIANCE PLOTTATION -----------------------------------------------------


#' Plot Proportion of Variance for L(M)M
#'
#' Covariate-Variances as modeled by linear (mixed) models will be displayed as box-plots.
#' It works with the output of 'mbecVarianceStats()' for methods 'lm' and 'lmm'. Format of this
#' output is a data.frame that contains a column for every model variable and as many rows as
#' there are features (OTUs, Genes, ..). Multiple frames may be used as input by putting them into
#' a list - IF the data.frames contain a column named 'type', this function will use 'facet_grid()'
#' to display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance linear mixed models
#' @param variance.obj, list or single output of 'mbecVarianceStats' with method lm
#' @export
#'
#' @examples
#' This will return a paneled plot that shows results for three variance assessments.
#' \dontrun{p.lmm <- mbecVarianceStatsPlot(variance.obj=list(df1, df2, df3))}
mbecVarianceStatsPlot <- function( variance.obj ) {

  plot.df <- variance.obj %>%
    bind_rows() %>% # this seems to work with single objects and lists
    gather(., "covariate", "variance", -type) %>%
    mutate(type = factor(type, levels = unique(type))) %>%
    mutate(variance = as.numeric(as.character(variance)))

  leplot <- ggplot(plot.df, aes(x = covariate, y = variance, fill = covariate)) +
    geom_boxplot() +
    facet_grid(cols = vars(type)) + ## this is the magic for comparative plotting
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          strip.text = element_text(size = 12), panel.grid = element_blank(),
          axis.text = element_text(size = 12), axis.title = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
    labs(x = 'Covariate', y = 'Proportion Variance', name = 'Covariate') + ylim(0,1)

  return(leplot)

}


#' Plot Proportion of Variance for pRDA
#'
#' Covariate-Variances as modeled by pRDA will be displayed as box-plots.
#' It works with the output of 'mbecVarianceStats()' for the method 'rda'. Format of this
#' output is a data.frame that contains a column for every model variable and as many rows as
#' there are features (OTUs, Genes, ..). Multiple frames may be used as input by putting them into
#' a list - IF the data.frames contain a column named 'type', this function will use 'facet_grid()'
#' to display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance partial Redundancy Analysis
#' @param rda.obj, list or single output of 'mbecVarianceStats' with method rda
#' @export
#'
#' @examples
#' This will return a paneled plot that shows results for three variance assessments.
#' \dontrun{p.rda <- mbecRDAStatsPlot(variance.obj=list(df1, df2, df3))}
mbecRDAStatsPlot <- function(rda.obj) {
  # first tidy-magic to create df for plotting
  leTest <- rda.obj %>%
    gather(., "covariate", "variance", -type) %>%
    mutate(type = factor(type, levels = unique(type))) %>%
    #separate(data=., col = "variance", into = c("variance","significance", "model.variance","model.significance"), sep="\\|") %>%
    mutate(variance = as.numeric(as.character(variance))) %>%
    mutate(variance.r = round(variance, 2))

  # now plot
  lePlot <- ggplot(data = leTest, aes(x = covariate, y = variance, fill = covariate)) +
    geom_bar(stat = "identity", position = 'dodge', colour = 'black') +
    # significance at the top
    # geom_text(data = leTest, aes(type, 100, label = significance),
    #           position = position_dodge(width = 0.9), size = 3) +
    # variance above the bars
    geom_text(data = leTest, aes(covariate, variance + 2.5, label = variance.r),
              position = position_dodge(width = 0.9), size = 3) +
    facet_grid(cols = vars(type)) + ## this is the magic for comparative plotting
    theme_bw() +
    labs(y = "Variance explained (%)") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          panel.grid = element_blank(), axis.text = element_text(size = 12),
          axis.title = element_text(size = 15), legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)) + ylim(0,100)

  return(lePlot)
  # FIN
}


#' Plot Proportion of Variance for PVCA
#'
#' Covariate-Variances as modeled by PVCA will be displayed as box-plots.
#' It works with the output of 'mbecVarianceStats()' for the method 'pvca'. Format of this
#' output is a data.frame that contains a column for every model variable and as many rows as
#' there are features (OTUs, Genes, ..). Multiple frames may be used as input by putting them into
#' a list - IF the data.frames contain a column named 'type', this function will use 'facet_grid()'
#' to display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance pvca
#' @param pvca.obj, list or single output of 'mbecVarianceStats' with method pvca
#' @export
#'
#' @examples
#' This will return a paneled plot that shows results for three variance assessments.
#' \dontrun{p.pvca <- mbecPVCAStatsPlot(variance.obj=list(df1, df2, df3))}
mbecPVCAStatsPlot <- function(pvca.obj) {
  # first tidy-magic to create df for plotting
  plot.df <- pvca.obj %>%
    gather(., "covariate", "variance", -type) %>%
    mutate(covariate=gsub("\\.",":",covariate)) %>%
    mutate(type = factor(type, levels = unique(type))) %>%
    mutate(variance = as.numeric(as.character(variance))) %>%
    mutate(variance.r = round(variance, 2)) %>%
    mutate(variance.p = round(variance*100, 2))

  # now plot
  lePlot <- ggplot(data = plot.df, aes(x = covariate, y = variance.p, fill = covariate)) +
    geom_bar(stat = "identity", position = 'dodge', colour = 'black') +
    geom_text(data = plot.df, aes(covariate, variance.p + 2.5, label = variance.p),
              position = position_dodge(width = 0.9), size = 3) + theme_bw() +
    facet_grid(cols=vars(type), scales="free", space="free_x", drop=T) +
    labs(x = "Random effects and Interactions", y = "Variance explained (%)") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          panel.grid = element_blank(), axis.text = element_text(size = 12),
          axis.title = element_text(size = 15), legend.title = element_text(size = 15),
          legend.text = element_text(size = 12)) + ylim(0,100)

  return(lePlot)
  # FIN
}


#' Plot Silhouette Coefficient
#'
#' The goodness of clustering assessed by the silhouette coefficient.
#' It works with the output of 'mbecVarianceStats()' for the method 's.coef'. Format of this
#' output is a data.frame that contains a column for every model variable and as many rows as
#' there are features (OTUs, Genes, ..). Multiple frames may be used as input by putting them into
#' a list - IF the data.frames contain a column named 'type', this function will use 'facet_grid()'
#' to display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance linear mixed models
#' @param scoef.obj, list or single output of 'mbecVarianceStats' with method s.coef
#' @export
#'
#' @examples
#' This will return a paneled plot that shows results for three variance assessments.
#' \dontrun{p.sc <- mbecSCOEFStatsPlot(variance.obj=list(df1, df2, df3))}
mbecSCOEFStatsPlot <- function(scoef.obj) {

  # first tidy-magic to create df for plotting
  plot.df <- scoef.obj %>%
    mutate(variable=gsub("\\.",":",variable)) %>%
    mutate(type = factor(type, levels = unique(type))) %>%
    mutate(sil.coefficient = as.numeric(as.character(sil.coefficient))) %>%
    mutate(sil.coefficient.r = round(sil.coefficient, 2))

  # now plot
  lePlot <- ggplot(plot.df, aes(x = variable, y = sil.coefficient, color = cluster, shape = variable)) +
    geom_point() + facet_grid(cols = vars(type)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          strip.text = element_text(size = 12), panel.grid = element_blank(),
          axis.text = element_text(size = 10), axis.title = element_text(size = 15),
          legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
    scale_color_manual(values = cols) +
    labs(x = 'Type', y = 'Silhouette Coefficient', name = 'Type')

  return(lePlot)
  # FIN
}


