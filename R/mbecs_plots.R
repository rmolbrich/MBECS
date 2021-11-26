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

  rle.plot <- ggplot2::ggplot(tmp.long, ggplot2::aes(x = specimen, y = values,
                                                     fill = get(model.vars[2]))) +
    ggplot2::stat_boxplot(color="black",notch = TRUE,
                          outlier.colour = "#E42032", outlier.fill = "white",
                          outlier.shape = 1, outlier.stroke = .5) +
    #facet_wrap(~Strain, ncol=2) +
    ggplot2::facet_grid(cols=ggplot2::vars(get(model.vars[1])), scales="free",
                        space="free_x", drop=TRUE) +
    ggplot2::scale_fill_manual(values = cols) +
    theme_rle() +
    ggplot2::guides(fill=ggplot2::guide_legend(title=ggplot2::element_blank()))

  return(rle.plot)
}



mbecBoxPlot <- function(tmp, otu.idx, model.var, cols) {

  legend.title <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.var,
                       perl = TRUE)
  ret.plot <- list()

  for( idx in otu.idx ) {
    p.box <- ggplot2::ggplot(data = tmp,
                             ggplot2::aes(x = get(model.var), y = get(idx),
                                          fill = get(model.var))) +
      ggplot2::stat_boxplot(geom = "errorbar", width = 0.4) +
      ggplot2::geom_boxplot() + ggplot2::scale_fill_manual(values = cols) +
      ggplot2::theme_bw() +
      theme_box() +
      ggplot2::labs(fill = legend.title, y = 'value',title = idx)

    p.density <- ggplot2::ggplot(tmp, ggplot2::aes(x = get(idx),
                                                   fill = get(model.var))) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(title = idx, x = 'Value', fill = legend.title) +
      theme_box()

    g <- ggplot2::ggplotGrob(p.box)$grobs

    legend <- g[[which(vapply(g, function(x) x$name,
                              FUN.VALUE=character(1)) == "guide-box")]]

    ret.plot[[eval(idx)]] <- gridExtra::arrangeGrob(p.box + ggplot2::theme(legend.position = 'none'),
                                                    p.density + ggplot2::theme(legend.position = 'none', plot.title = ggplot2::element_blank()),
                                                    legend,
                                                    ncol = 1, nrow = 3, heights = c(5,4.5,1))
  }
  return(ret.plot)
}



mbecHeatPlot <- function(center, scale, tmp.cnts, tmp.meta, model.vars) {

  p.title <- paste("Heatmap - Centered: ", center, " Scaled: ", scale, sep="")
  heat.plot <- pheatmap::pheatmap(tmp.cnts,
                                  scale = 'none',
                                  cluster_rows = FALSE,
                                  cluster_cols = TRUE,
                                  fontsize_row = 4, fontsize_col = 6,
                                  fontsize = 8,
                                  clustering_distance_rows = 'euclidean',
                                  clustering_method = 'ward.D',
                                  treeheight_row = 30,
                                  annotation_col = tmp.meta[,eval(model.vars)],
                                  border_color = 'NA',
                                  main = p.title)

  return(heat.plot)
}



mbecMosaicPlot <- function(study.summary, model.vars) {

  vars.axes <- mbecUpperCase(model.vars)

  plot.v2 <- ggplot2::ggplot(study.summary, ggplot2::aes(x = Var1, y= Freq.scaled, group = Var2, fill=Var1)) +
    ggplot2::facet_grid(cols=ggplot2::vars(Var2), scales="free", space="free_x", drop=TRUE) +
    ggplot2::geom_bar(stat = "identity", width = 0.9) +
    ggplot2::guides(fill = ggplot2::guide_legend(title=eval(vars.axes[1]), reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ggplot2::ylab("Proportion of all observations") +
    theme_mosaic(legend_position = "bottom")

  plot.v1 <- ggplot2::ggplot(study.summary, ggplot2::aes(x = Var2, y= Freq.scaled, fill=Var2)) +
    ggplot2::facet_grid(cols=ggplot2::vars(Var1), scales="free", space="free_x", drop=TRUE) +
    ggplot2::geom_bar(stat = "identity", width = 0.9) +
    ggplot2::guides(fill = ggplot2::guide_legend(title=eval(vars.axes[2]), reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ggplot2::ylab("Proportion of all observations") +
    theme_mosaic()

  mosaic.plot <- gridExtra::grid.arrange(plot.v2, plot.v1, ncol=1, nrow=2, heights=c(1,1))

  return(mosaic.plot)
}



