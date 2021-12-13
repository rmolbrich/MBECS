# PRELIMINARY-ANALYSES PLOTS ----------------------------------------------

#' RLE plotting function
#'
#' Takes data.frame from mbecRLE and produces a ggplot2 object.
#'
#' @keywords RLE relative log expression
#' @param tmp.long 'mbecRLE' output
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @param cols vector of hex-colors to use
#' @return ggplot2 object
mbecRLEPlot <- function(tmp.long, model.vars, cols) {

  rle.plot <- ggplot2::ggplot(
    tmp.long,
    ggplot2::aes(x = specimen, y = values,fill = get(model.vars[2]))) +
    ggplot2::stat_boxplot(color = "black", notch = FALSE,
                          outlier.colour = "#E42032", outlier.fill = "white",
                          outlier.shape = 1, outlier.stroke = 0.5,
                          outlier.size = 0.5, outlier.alpha=0.5) +
    # facet_wrap(~Strain, ncol=2) +
    ggplot2::facet_grid(cols = ggplot2::vars(get(model.vars[1])), scales = "free",
                        space = "free_x", drop = TRUE) + ggplot2::scale_fill_manual(values = cols) +
    theme_rle() + ggplot2::guides(fill = ggplot2::guide_legend(title = ggplot2::element_blank()))

  return(rle.plot)
}



mbecBoxPlot <- function(tmp, otu.idx, model.var,
                        cols) {

  legend.title <- gsub("(^|[[:space:]])([[:alpha:]])",
                       "\\1\\U\\2", model.var, perl = TRUE)
  ret.plot <- list()

  for (idx in otu.idx) {
    p.box <- ggplot2::ggplot(data = tmp, ggplot2::aes(x = get(model.var),
                                                      y = get(idx), fill = get(model.var))) +
      ggplot2::stat_boxplot(geom = "errorbar",
                            width = 0.4) + ggplot2::geom_boxplot() +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::theme_bw() + theme_box() +
      ggplot2::labs(fill = legend.title, y = "value",
                    title = idx)

    p.density <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(idx),
                                                   fill = get(model.var))) + ggplot2::geom_density(alpha = 0.5) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(title = idx, x = "Value",
                    fill = legend.title) + theme_box()

    g <- ggplot2::ggplotGrob(p.box)$grobs

    legend <- g[[which(vapply(g, function(x) x$name,
                              FUN.VALUE = character(1)) == "guide-box")]]

    ret.plot[[eval(idx)]] <- gridExtra::arrangeGrob(p.box +
                                                      ggplot2::theme(legend.position = "none"),
                                                    p.density + ggplot2::theme(legend.position = "none",
                                                                               plot.title = ggplot2::element_blank()),
                                                    legend, ncol = 1, nrow = 3, heights = c(5,
                                                                                            4.5, 1))
  }
  return(ret.plot)
}



mbecHeatPlot <- function(center, scale, tmp.cnts, tmp.meta,
                         model.vars) {

  p.title <- paste("Heatmap - Centered: ", center, " Scaled: ",
                   scale, sep = "")
  heat.plot <- pheatmap::pheatmap(tmp.cnts, scale = "none",
                                  cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 4,
                                  fontsize_col = 6, fontsize = 8, clustering_distance_rows = "euclidean",
                                  clustering_method = "ward.D", treeheight_row = 30,
                                  annotation_col = tmp.meta[, eval(model.vars)], border_color = "NA",
                                  main = p.title,
                                  annotation_names_col = TRUE,
                                  show_colnames = FALSE)

  return(heat.plot)
}



mbecMosaicPlot <- function(study.summary,
                           model.vars) {

  vars.axes <- mbecUpperCase(model.vars)

  plot.v2 <- ggplot2::ggplot(study.summary,
                             ggplot2::aes(x = Var1,
                                          y = Freq.scaled,
                                          group = Var2, fill = Var1)) +
    ggplot2::facet_grid(cols = ggplot2::vars(Var2),
                        scales = "free",
                        space = "free_x",
                        drop = TRUE) +
    ggplot2::geom_bar(stat = "identity",
                      width = 0.9) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = eval(vars.axes[1]),
                                                 reverse = TRUE,
                                                 keywidth = 1, keyheight = 1)) +
    ggplot2::ylab("Proportion of all observations") +
    theme_mosaic(legend_position = "bottom")

  plot.v1 <- ggplot2::ggplot(study.summary,
                             ggplot2::aes(x = Var2,
                                          y = Freq.scaled,
                                          fill = Var2)) +
    ggplot2::facet_grid(cols = ggplot2::vars(Var1),
                        scales = "free",
                        space = "free_x",
                        drop = TRUE) +
    ggplot2::geom_bar(stat = "identity",
                      width = 0.9) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = eval(vars.axes[2]),
                                                 reverse = TRUE,
                                                 keywidth = 1, keyheight = 1)) +
    ggplot2::ylab("Proportion of all observations") +
    theme_mosaic()

  mosaic.plot <- gridExtra::grid.arrange(plot.v2,
                                         plot.v1, ncol = 1,
                                         nrow = 2, heights = c(1,
                                                               1))

  return(mosaic.plot)
}



# VARIANCE PLOTTATION -----------------------------------------------------


#' Plot Proportion of Variance for L(M)M
#'
#' Covariate-Variances as modeled by linear (mixed) models will be displayed as
#' box-plots. It works with the output of 'mbecVarianceStats()' for methods 'lm'
#' and 'lmm'. Format of this output is a data.frame that contains a column for
#' every model variable and as many rows as there are features
#' (OTUs, Genes, ..). Multiple frames may be used as input by putting them into
#' a list - IF the data.frames contain a column named 'type', this function will
#' use 'facet_grid()' to display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance linear mixed models
#' @param variance.obj, output of 'mbecVarianceStats' with method lm
#' @return A ggplot2 box-plot object.
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessments.
#' df.var.lm <- mbecModelVariance(input.obj=datadummy,
#' model.vars=c('group','batch'),
#' method='lm', type='RAW')
#' plot.lm <- mbecVarianceStatsPlot(variance.obj=df.var.lm)
mbecVarianceStatsPlot <- function(variance.obj) {

  plot.df <- variance.obj %>%
    dplyr::bind_rows() %>%
    tidyr::gather("covariate",
                  "variance", -type) %>%
    dplyr::mutate(type = factor(type,
                                levels = unique(type))) %>%
    dplyr::mutate(variance = as.numeric(as.character(variance)))

  leplot <- ggplot2::ggplot(plot.df,
                            ggplot2::aes(x = covariate,
                                         y = variance, fill = covariate)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(cols = ggplot2::vars(type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                       hjust = 1), strip.text = ggplot2::element_text(size = 12),
                   panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 15),
                   legend.title = ggplot2::element_text(size = 15),
                   legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::labs(x = "Linear (Mixed) Model",
                  y = "Proportion Variance",
                  name = "Covariate") +
    ggplot2::ylim(0, 1) +
    ggplot2::facet_grid(cols = ggplot2::vars(type), scales = "free",
                        space = "free_x", drop = TRUE)

  return(leplot)

}


#' Plot Proportion of Variance for pRDA
#'
#' Covariate-Variances as modeled by pRDA will be displayed as box-plots.
#' It works with the output of 'mbecVarianceStats()' for the method 'rda'.
#' Format of this output is a data.frame that contains a column for every model
#' variable and as many rows as there are features (OTUs, Genes, ..). Multiple
#' frames may be used as input by putting them into a list - IF the data.frames
#' contain a column named 'type', this function will use 'facet_grid()' to
#' display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance partial Redundancy Analysis
#' @param rda.obj, list or single output of 'mbecVarianceStats' with method rda
#' @return A ggplot2 box-plot object.
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for three variance
#' # assessments.
#' df.var.rda <- mbecModelVariance(input.obj=datadummy,
#' model.vars=c('group','batch'), method='rda', type='RAW')
#' plot.rda <- mbecRDAStatsPlot(rda.obj=df.var.rda)
mbecRDAStatsPlot <- function(rda.obj) {

  leTest <- rda.obj %>%
    tidyr::gather("covariate", "variance", -type) %>%
    dplyr::mutate(type = factor(type, levels = unique(type))) %>%
    dplyr::mutate(variance = as.numeric(as.character(variance))) %>%
    dplyr::mutate(variance.r = round(variance,
                                     2))

  lePlot <- ggplot2::ggplot(data = leTest, ggplot2::aes(x = covariate,
                                                        y = variance, fill = covariate)) + ggplot2::geom_bar(stat = "identity",
                                                                                                             position = "dodge", colour = "black") +
    ggplot2::geom_text(data = leTest, ggplot2::aes(covariate,
                                                   variance + 2.5, label = variance.r),
                       position = ggplot2::position_dodge(width = 0.9),
                       size = 3) + ggplot2::facet_grid(cols = ggplot2::vars(type)) +
    ggplot2::theme_bw() + ggplot2::labs(x = "RDA",
                                        y = "Variance explained (%)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                       hjust = 1), panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 15),
                   legend.title = ggplot2::element_text(size = 15),
                   legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::ylim(0, 100) +
    ggplot2::facet_grid(cols = ggplot2::vars(type), scales = "free",
                        space = "free_x", drop = TRUE)

  return(lePlot)
  # FIN
}


#' Plot Proportion of Variance for PVCA
#'
#' Covariate-Variances as modeled by PVCA will be displayed as box-plots.
#' It works with the output of 'mbecVarianceStats()' for the method 'pvca'.
#' Format of this output is a data.frame that contains a column for every model
#' variable and as many rows as there are features (OTUs, Genes, ..). Multiple
#' frames may be used as input by putting them into a list - IF the data.frames
#' contain a column named 'type', this function will use 'facet_grid()' to
#' display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance pvca
#' @param pvca.obj, output of 'mbecVarianceStats' with method pvca
#' @return A ggplot2 box-plot object.
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' df.var.pvca <- mbecModelVariance(input.obj=datadummy,
#' model.vars=c('group','batch'), method='pvca', type='RAW')
#' plot.pvca <- mbecPVCAStatsPlot(pvca.obj=df.var.pvca)
mbecPVCAStatsPlot <- function(pvca.obj) {

  plot.df <- pvca.obj %>%
    tidyr::gather("covariate", "variance",
                  -type) %>%
    dplyr::mutate(covariate = gsub("\\.",
                                   ":", covariate)) %>%
    dplyr::mutate(type = factor(type,
                                levels = unique(type))) %>%
    dplyr::mutate(variance = as.numeric(as.character(variance))) %>%
    dplyr::mutate(variance.r = round(variance,
                                     2)) %>%
    dplyr::mutate(variance.p = round(variance *
                                       100, 2))

  lePlot <- ggplot2::ggplot(data = plot.df,
                            ggplot2::aes(x = covariate, y = variance.p,
                                         fill = covariate)) + ggplot2::geom_bar(stat = "identity",
                                                                                position = "dodge", colour = "black") +
    ggplot2::geom_text(data = plot.df,
                       ggplot2::aes(covariate, variance.p +
                                      2.5, label = variance.p),
                       position = ggplot2::position_dodge(width = 0.9),
                       size = 3) + ggplot2::theme_bw() +
    ggplot2::facet_grid(cols = ggplot2::vars(type),
                        scales = "free", space = "free_x",
                        drop = TRUE) + ggplot2::labs(x = "PVCA - Random effects and Interactions",
                                                     y = "Variance explained (%)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                       hjust = 1), panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 15),
                   legend.title = ggplot2::element_text(size = 15),
                   legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::ylim(0, 100)+
    ggplot2::facet_grid(cols = ggplot2::vars(type), scales = "free",
                        space = "free_x", drop = TRUE)

  return(lePlot)
  # FIN
}


#' Plot Silhouette Coefficient
#'
#' The goodness of clustering assessed by the silhouette coefficient.
#' It works with the output of 'mbecVarianceStats()' for the method 's.coef'.
#' Format of this output is a data.frame that contains a column for every model
#' variable and as many rows as there are features (OTUs, Genes, ..). Multiple
#' frames may be used as input by putting them into a list - IF the data.frames
#' contain a column named 'type', this function will use 'facet_grid()' to
#' display side-by-side panels to enable easy comparison.
#'
#' @keywords plot proportion variance linear mixed models
#' @param scoef.obj, output of 'mbecVarianceStats' with method s.coef
#' @return A ggplot2 dot-plot object.
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' df.var.scoef <- mbecModelVariance(input.obj=datadummy,
#' model.vars=c('group','batch'), method='s.coef', type='RAW')
#' plot.scoef <- mbecSCOEFStatsPlot(scoef.obj=df.var.scoef)
mbecSCOEFStatsPlot <- function(scoef.obj) {
  ## ToDo: make my own colors - with black jack and hookers
  cols <- pals::tableau20(20)
  # first tidy-magic to create df for plotting
  plot.df <- scoef.obj %>%
    dplyr::mutate(variable = gsub("\\.", ":", variable)) %>%
    dplyr::mutate(type = factor(type, levels = unique(type))) %>%
    dplyr::mutate(sil.coefficient = as.numeric(as.character(sil.coefficient))) %>%
    dplyr::mutate(sil.coefficient.r = round(sil.coefficient, 2))

  # now plot
  lePlot <- ggplot2::ggplot(plot.df, ggplot2::aes(x = variable, y = sil.coefficient,
                                                  color = cluster, shape = variable)) + ggplot2::geom_point() + ggplot2::facet_grid(cols = ggplot2::vars(type)) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                                             hjust = 1), strip.text = ggplot2::element_text(size = 12), panel.grid = ggplot2::element_blank(),
                                         axis.text = ggplot2::element_text(size = 10), axis.title = ggplot2::element_text(size = 15),
                                         legend.title = ggplot2::element_text(size = 15), legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::scale_color_manual(values = cols) + ggplot2::labs(x = "Silhouette Coefficient", y = "Silhouette Coefficient",
                                                               name = "Type") +
    ggplot2::facet_grid(cols = ggplot2::vars(type), scales = "free",
                        space = "free_x", drop = TRUE)

  return(lePlot)
  # FIN
}

