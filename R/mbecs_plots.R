# PRELIMINARY-ANALYSES PLOTS ----------------------------------------------

#' RLE plotting function
#'
#' Takes data.frame from mbecRLE and produces a ggplot2 object.
#'
#' @keywords RLE Relative Log Expression Plot
#' @param rle.df 'mbecRLE'  data output
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring
#' @param label Name of the plot displayed as legend title.
#' @return ggplot2 object
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' rle.df <- mbecRLE(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' type="clr", return.data=TRUE)
#' plot.rle <- mbecRLEPlot(rle.df, c('group','batch'))
mbecRLEPlot <- function(rle.df, model.vars, label=NULL) {

  mbecCols <- pals::tableau20(20)
  n.tiles <- dim(unique(rle.df[,eval(model.vars[2])]))[1]

  if( n.tiles <= 3 ) {
    ncols = n.tiles
    nrows = 1
  } else {
    ncols = 3
    nrows = (n.tiles / 3)
  }

  rle.plot <- ggplot2::ggplot(
    rle.df,
    ggplot2::aes(x = specimen, y = values,fill = get(model.vars[1]))) +
    ggplot2::stat_boxplot(color = "black", notch = FALSE,
                          lwd=0.5, fatten=0.75,
                          outlier.colour = "#E42032", outlier.fill = "white",
                          outlier.shape = 1, outlier.stroke = 0.5,
                          outlier.size = 0.5, outlier.alpha=0.5) +
    ggplot2::facet_wrap(~get(model.vars[2]), ncol=ncols,nrow = nrows , scales = "free_x",
                        drop = TRUE) + ggplot2::theme(plot.margin=ggplot2::unit(c(0.2,0.2,0.05,0.2), "cm")) +
    # ggplot2::facet_grid(cols = ggplot2::vars(get(model.vars[2])), scales = "free",
    #                     space = "free_x", drop = TRUE) +
    ggplot2::scale_fill_manual(values = mbecCols) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = eval(label), fill = mbecUpperCase(eval(model.vars[1]))) +
    # adjustments for the legend
    ggplot2::theme(legend.position="bottom",
                   legend.text = ggplot2::element_text(color = "black", size=12),
                   legend.key = ggplot2::element_rect(size=12),
                   #axis.text.x=ggplot2::element_text(angle=35, hjust=1),
                   axis.title.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank(),
                   axis.ticks.x=ggplot2::element_blank())
    ggplot2::guides(fill = ggplot2::guide_legend(title = ggplot2::element_blank()))

  return(rle.plot)
}


#' Variability boxes plotting function
#'
#' Takes data.frame from mbecBox and produces a ggplot2 object.
#'
#' @keywords Abundance Variability Plot
#' @param tmp Count of selected features.
#' @param otu.idx Index of selected Otus in the data.
#' @param model.var Which covariate to group Otus by.
#' @param label Name of the plot displayed as legend title.
#' @return ggplot2 object
#' @export
#'
#' @examples
#' # This will return a list of the five most variable features grouped by the
#' # covariate 'batch'.
#' box.df <- mbecBox(input.obj=dummy.mbec, method='TOP', n=5,
#' model.var='batch', type="otu", return.data=TRUE)
#' plot.box <- mbecBoxPlot(box.df[[1]], box.df[[2]], 'batch')
mbecBoxPlot <- function(tmp, otu.idx, model.var, label=NULL) {

  mbecCols <- pals::tableau20(20)[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19)]
  x.angle = 0
  x.hjust = 0.5
  density.lwd = 0.2
  title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex =0.75

  legend.title <- paste(label, gsub("(^|[[:space:]])([[:alpha:]])",
                       "\\1\\U\\2", model.var, perl = TRUE), sep = " ")
  ret.plot <- list()

  for (idx in otu.idx) {
    p.box <-
      ggplot2::ggplot(data = tmp, ggplot2::aes(x = get(model.var),
                                               y = get(idx),
                                               fill = get(model.var))) +
      ggplot2::stat_boxplot(geom = "errorbar",
                            width = 0.4) + ggplot2::geom_boxplot() +
      ggplot2::scale_fill_manual(values = mbecCols) +
      ggplot2::theme_bw() +
        ggplot2::theme(
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = x.angle, hjust = x.hjust),
          panel.grid = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(size = 10),
          axis.title = ggplot2::element_text(size = 12),
          plot.title = ggplot2::element_text(hjust = 0.5),

          legend.position = 'bottom', legend.box = 'horizontal',
          legend.direction = 'horizontal',
          legend.key.height = ggplot2::unit(0.8, 'cm'),
          legend.key.width = ggplot2::unit(0.4, 'cm'),
          legend.title =
            ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
          legend.spacing.x = ggplot2::unit(0.4, 'cm'),
          legend.spacing.y = ggplot2::unit(0.4, 'cm'),
          legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))
        ) +
      ggplot2::labs(fill = legend.title, y = "value",
                    title = idx)

    p.density <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(idx),
                                                   fill = get(model.var))) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::scale_fill_manual(values = mbecCols) +
      ggplot2::labs(title = idx, x = "Value",
                    fill = legend.title) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = x.angle, hjust = x.hjust),
        panel.grid = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(hjust = 0.5),

        legend.position = 'bottom', legend.box = 'horizontal',
        legend.direction = 'horizontal',
        legend.key.height = ggplot2::unit(0.8, 'cm'),
        legend.key.width = ggplot2::unit(0.4, 'cm'),
        legend.title =
          ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
        legend.spacing.x = ggplot2::unit(0.4, 'cm'),
        legend.spacing.y = ggplot2::unit(0.4, 'cm'),
        legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex)))

    g <- ggplot2::ggplotGrob(p.box)$grobs

    legend <- g[[which(vapply(g, function(x) x$name,
                              FUN.VALUE = character(1)) == "guide-box")]]

    ret.plot[[eval(idx)]] <-
      gridExtra::arrangeGrob(p.box +
                               ggplot2::theme(legend.position = "none"),
                             p.density +
                               ggplot2::theme(legend.position = "none",
                                              plot.title =
                                                ggplot2::element_blank()),
                             legend, ncol = 1, nrow = 3, heights = c(5,4.5, 1))
  }
  return(ret.plot)
}


#' Heatmap plotting function
#'
#' Takes data.frame from 'mbecHeat()' and produces a ggplot2 object.
#'
#' @keywords RLE relative log expression
#' @param tmp.cnts Count values of selected features.
#' @param tmp.meta Covariate information for potting.
#' @param model.vars Two covariates of interest to select by first variable
#' selects panels and second one determines coloring.
#' @param label Name of the plot displayed as legend title.
#' @return ggplot2 object
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' heat.df <- mbecHeat(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' center=TRUE, scale=TRUE, method='TOP', n=5, return.data=TRUE)
#' plot.heat <- mbecHeatPlot(center=TRUE, scale=TRUE, tmp.cnts=heat.df[[1]],
#' tmp.meta=heat.df[[2]], model.vars=c('group','batch'))
mbecHeatPlot <- function(tmp.cnts, tmp.meta, model.vars, label=NULL) {

  p.title <- paste("Heatmap - ", label, sep = "")
  heat.plot <-
    pheatmap::pheatmap(tmp.cnts, scale = "none",
                       cluster_rows = FALSE, cluster_cols = TRUE,
                       fontsize_row = 4, fontsize_col = 6, fontsize = 8,
                       clustering_distance_rows = "euclidean",
                       clustering_method = "ward.D", treeheight_row = 30,
                       annotation_col = tmp.meta[, eval(model.vars)],
                       border_color = "NA", main = p.title,
                       annotation_names_col = TRUE, show_colnames = FALSE)
  return(heat.plot)
}


#' Mosaic plotting function
#'
#' Takes data.frame from mbecMosaic and produces a ggplot2 object.
#'
#' @keywords RLE relative log expression
#' @param study.summary 'mbecMosaic' output object.
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring.
#' @return ggplot2 object
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' mosaic.df <- mbecMosaic(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' return.data=TRUE)
#' plot.mosaic <- mbecMosaicPlot(study.summary=mosaic.df,
#' model.vars=c('group','batch'))
mbecMosaicPlot <- function(study.summary,
                           model.vars) {

  main_color = "#004B5A"
  x.angle = 0
  x.hjust = 0.5
  density.lwd = 0.2
  title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex =0.75

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
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_text(color = eval(main_color), size=12),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(color = "#7F7F7F"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), angle = 90),
                   legend.position = 'bottom', legend.box = 'horizontal',
                   legend.direction = 'horizontal',
                   legend.key.height = ggplot2::unit(0.2, 'cm'),
                   legend.key.width = ggplot2::unit(0.1, 'cm'),
                   legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
                   legend.spacing.x = ggplot2::unit(0.1, 'cm'),
                   legend.spacing.y = ggplot2::unit(0.1, 'cm'),
                   legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))) +
                ggplot2::theme(plot.margin=ggplot2::unit(c(0.2,0.2,0.05,0.2), "cm"))

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
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_text(color = eval(main_color), size=12),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "#7F7F7F"),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = ggplot2::rel(1), angle = 90),
      legend.position = 'bottom', legend.box = 'horizontal',
      legend.direction = 'horizontal',
      legend.key.height = ggplot2::unit(0.2, 'cm'),
      legend.key.width = ggplot2::unit(0.1, 'cm'),
      legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
      legend.spacing.x = ggplot2::unit(0.1, 'cm'),
      legend.spacing.y = ggplot2::unit(0.1, 'cm'),
      legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))) +
    ggplot2::theme(plot.margin=ggplot2::unit(c(0.05,0.2,0.2,0.2), "cm"))

  ## Function to extract legend
  g_legend <- function(a.gplot){
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }

  legend.v2 <- g_legend(plot.v2)
  legend.v1 <- g_legend(plot.v1)

  mosaic.plot <- gridExtra::grid.arrange(plot.v2 + ggplot2::theme(legend.position = "none"),
                                      plot.v1 + ggplot2::theme(legend.position = "none"),
                                      gridExtra::grid.arrange(legend.v1, legend.v2, ncol=2, nrow=1),
                                      ncol = 1, nrow = 3, widths = c(1), heights = c(1, 1,0.2),
                                      padding= -10)
  return(mosaic.plot)
}


#' PCA plotting function
#'
#' Takes data.frame from mbecPCA and produces a ggplot2 object.
#'
#' @keywords RLE relative log expression
#' @param plot.df Data.frame containing principal component data.
#' @param metric.df Data.frame containing covariate data.
#' @param model.vars two covariates of interest to select by first variable
#' selects panels and second one determines coloring.
#' @param pca.axes NMumerical two-piece vector that selects PCs to plot.
#' @param label Name of the plot displayed as legend title.
#' @return ggplot2 object
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' pca.df <- mbecPCA(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), pca.axes=c(1,2), return.data=TRUE)
#' plot.pca <- mbecPCAPlot(plot.df=pca.df[[1]], metric.df=pca.df[[2]],
#' model.vars=c('group','batch'), pca.axes=c(1,2))
mbecPCAPlot <- function(plot.df, metric.df, model.vars, pca.axes, label=NULL) {
  mbecCols <- pals::tableau20(20)[c(1,3,5,7,9,11,13,15,17)]
  x.angle = 0
  x.hjust = 0.5
  # density.lwd = 0.2
  # title.cex = 1.5
  legend.cex = 0.7
  legend.title.cex =0.75

  # meh
  # ks.table <- mbecPCTest(plot.df, pca.axes, model.vars)
  # plot.annotation.top <- paste(colnames(ks.table),
  #                              ks.table[1,], sep = ": ", collapse = " \n")
  # plot.annotation.right <- paste(colnames(ks.table),
  #                                ks.table[2,], sep = ": ", collapse = " \n")
  # meh

  var.color <- model.vars[1]; var.shape <- model.vars[2]
  label.col <- mbecUpperCase(model.vars[1])
  label.sha <- mbecUpperCase(model.vars[2])
  title <- paste("PCA:", label.sha, "-", label.col)
  if (length(model.vars) >= 2) {
    pMain <-
      ggplot2::ggplot(data = plot.df,
        ggplot2::aes(x = get(colnames(plot.df[pca.axes[1] + 1])),
                     y = get(colnames(plot.df[pca.axes[2] + 1])),
                     colour = get(var.color),shape = get(var.shape))) +
      ggplot2::scale_shape_manual(values=c(0,1,2,3,6,8,15,16,17,23,25,4,5,9)) +
      ggplot2::geom_point() + ggplot2::scale_color_manual(values = mbecCols) +
      ggplot2::labs(title = label, colour = label.col, shape = label.sha) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[1]],
                    metric.df$axis.max[pca.axes[1]]) +
      ggplot2::ylim(metric.df$axis.min[pca.axes[2]],
                    metric.df$axis.max[pca.axes[2]]) +
      ggplot2::xlab(paste0(colnames(plot.df[pca.axes[1] + 1]), ": ",
                           metric.df$var.explained[pca.axes[1]], "% expl.var")) +
      ggplot2::ylab(paste0(colnames(plot.df[pca.axes[2] + 1]), ": ",
                           metric.df$var.explained[pca.axes[2]], "% expl.var")) +
      ggplot2::theme_bw() +
      ggplot2::theme(

        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.position = 'right', legend.box = 'horizontal',
        legend.direction = 'vertical',
        legend.key.height = ggplot2::unit(0.2, 'cm'),
        legend.key.width = ggplot2::unit(0.1, 'cm'),
        legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
        legend.spacing.x = ggplot2::unit(0.1, 'cm'),
        legend.spacing.y = ggplot2::unit(0.1, 'cm'),
        legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex)))

    pTop <- ggplot2::ggplot(data = plot.df,
                            ggplot2::aes(x = get(colnames(plot.df[pca.axes[1] +
                              1])), fill = get(var.color),
                            linetype = get(var.shape))) + ggplot2::geom_density(size = 0.2,
                                                                                                                                                                     alpha = 0.5) + ggplot2::ylab("Density") + ggplot2::scale_fill_manual(values = mbecCols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      ggplot2::theme_bw() +
      ggplot2::theme(

        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.position = 'right', legend.box = 'horizontal',
        legend.direction = 'vertical',
        legend.key.height = ggplot2::unit(0.2, 'cm'),
        legend.key.width = ggplot2::unit(0.1, 'cm'),
        legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
        legend.spacing.x = ggplot2::unit(0.1, 'cm'),
        legend.spacing.y = ggplot2::unit(0.1, 'cm'),
        legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))) +
      ggplot2::labs(title = ggplot2::element_blank()) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.8)),
                     plot.title = ggplot2::element_text(hjust = 0.5,
                                                        size = ggplot2::rel(1.5))) #+
      # ggplot2::annotate("text", -Inf, Inf, label = plot.annotation.top, hjust = -0.05, vjust = 1.15,
      #                   colour = "#666666", size = ggplot2::rel(3))

    pRight <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[2] +
                                                                                      1])), fill = get(var.color), linetype = get(var.shape))) +
      ggplot2::geom_density(size = 0.2,alpha = 0.5) + ggplot2::coord_flip() +
      ggplot2::ylab("Density") + ggplot2::scale_fill_manual(values = mbecCols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      ggplot2::theme_bw() +
      ggplot2::theme(

        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        legend.position = 'right', legend.box = 'horizontal',
        legend.direction = 'vertical',
        legend.key.height = ggplot2::unit(0.2, 'cm'),
        legend.key.width = ggplot2::unit(0.1, 'cm'),
        legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.cex)),
        legend.spacing.x = ggplot2::unit(0.1, 'cm'),
        legend.spacing.y = ggplot2::unit(0.1, 'cm'),
        legend.text = ggplot2::element_text(size = ggplot2::rel(legend.cex))) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = ggplot2::rel(0.8)),
                                   axis.title.y = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                                   plot.title = ggplot2::element_blank()) #+
      # ggplot2::annotate("text", Inf, -Inf, label = plot.annotation.right, hjust = -0.05, vjust = 1.15,
      #                   colour = "#666666", size = ggplot2::rel(3))

  } else {
    pMain <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[1] +
                                                                                     1])), y = get(colnames(plot.df[pca.axes[2] + 1])), colour = get(var.color))) +
      ggplot2::geom_point() + ggplot2::scale_color_manual(values = mbecCols) +
      ggplot2::labs(title = label, colour = label.col, shape = label.sha) + ggplot2::xlim(metric.df$axis.min[pca.axes[1]],
                                                                           metric.df$axis.max[pca.axes[1]]) + ggplot2::ylim(metric.df$axis.min[pca.axes[2]],
                                                                                                                            metric.df$axis.max[pca.axes[2]]) + ggplot2::xlab(paste0(colnames(plot.df[pca.axes[1] +
                                                                                                                                                                                                       1]), ": ", metric.df$var.explained[pca.axes[1]], "% expl.var")) + ggplot2::ylab(paste0(colnames(plot.df[pca.axes[2] +
                                                                                                                                                                                                                                                                                                                 1]), ": ", metric.df$var.explained[pca.axes[2]], "% expl.var")) #+ theme_pca()

    pTop <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[1] +
                                                                                    1])), fill = get(var.color))) + ggplot2::geom_density(size = 0.2,
                                                                                                                                          alpha = 0.5) + ggplot2::ylab("Density") + ggplot2::scale_fill_manual(values = mbecCols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      #theme_pca() +
      ggplot2::labs(title = ggplot2::element_blank()) + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                                                  axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.8)), plot.title = ggplot2::element_text(hjust = 0.5,
                                                                                                                                                                     size = ggplot2::rel(1.5))) #+
      # ggplot2::annotate("text", -Inf, Inf, label = plot.annotation.top, hjust = -0.05, vjust = 1.15,
      #                   colour = "#666666", size = ggplot2::rel(3))

    pRight <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = get(colnames(plot.df[pca.axes[2] +
                                                                                      1])), fill = get(var.color))) + ggplot2::geom_density(size = 0.2,
                                                                                                                                            alpha = 0.5) + ggplot2::coord_flip() + ggplot2::ylab("Density") + ggplot2::scale_fill_manual(values = mbecCols) +
      ggplot2::xlim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      #theme_pca() +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = ggplot2::rel(0.8)),
                                   axis.title.y = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                                   plot.title = ggplot2::element_blank()) #+
      # ggplot2::annotate("text", Inf, -Inf, label = plot.annotation.right, hjust = -0.05, vjust = 1.15,
      #                   colour = "#666666", size = ggplot2::rel(3))
  }

  g <- ggplot2::ggplotGrob(pMain)$grobs
  legend <- g[[which(vapply(g, function(x) x$name, FUN.VALUE = character(1)) ==
                       "guide-box")]]
  ret.plot <- gridExtra::grid.arrange(pTop + ggplot2::theme(legend.position = "none"),
                                      legend, pMain + ggplot2::theme(legend.position = "none"), pRight + ggplot2::theme(legend.position = "none"),
                                      ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1.5, 3))

  return(ret.plot)





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
#' df.var.lm <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), method='lm', type='clr')
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
    ggplot2::geom_boxplot(lwd=0.5, fatten=0.75,
                          outlier.colour = "#E42032", outlier.fill = "white",
                          outlier.shape = 1, outlier.stroke = 0.5,
                          outlier.size = 0.5, outlier.alpha=0.5) +
    ggplot2::facet_grid(cols = ggplot2::vars(type)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
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
#' df.var.rda <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), method='rda', type='clr')
#' plot.rda <- mbecRDAStatsPlot(rda.obj=df.var.rda)
mbecRDAStatsPlot <- function(rda.obj) {

  leTest <- rda.obj %>%
    tidyr::gather("covariate", "variance", -type) %>%
    dplyr::mutate(type = factor(type, levels = unique(type))) %>%
    dplyr::mutate(variance = as.numeric(as.character(variance))) %>%
    dplyr::mutate(variance.r = round(variance,
                                     2))

  lePlot <-
    ggplot2::ggplot(data = leTest,
                    ggplot2::aes(x=covariate, y=variance, fill=covariate)) +
    ggplot2::geom_bar(stat = "identity",position = "dodge", colour = "black") +
    ggplot2::geom_text(data = leTest, ggplot2::aes(covariate, variance + 2.5,
                                                   label = variance.r),
                       position = ggplot2::position_dodge(width = 0.9),
                       size = 3) +
    ggplot2::facet_grid(cols = ggplot2::vars(type)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "RDA",y = "Variance explained (%)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,hjust = 1),
                   panel.grid = ggplot2::element_blank(),
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
#' df.var.pvca <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('batch','group'), method='pvca', type='clr')
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
#' df.var.scoef <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('batch','group'), method='s.coef', type='clr')
#' plot.scoef <- mbecSCOEFStatsPlot(scoef.obj=df.var.scoef)
mbecSCOEFStatsPlot <- function(scoef.obj) {
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
    ggplot2::scale_color_manual(values = cols) + ggplot2::labs(x = "Silhouette Coefficient", y = "Grouping",
                                                               name = "Type") +
    ggplot2::facet_grid(cols = ggplot2::vars(type), scales = "free",
                        space = "free_x", drop = TRUE)

  return(lePlot)
  # FIN
}

