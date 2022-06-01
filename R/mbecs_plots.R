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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_wrap theme theme_bw
#' scale_fill_manual labs element_text element_blank element_rect guides
#' guide_legend
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' data(dummy.mbec)
#' rle.df <- mbecRLE(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' type="clr", return.data=TRUE)
#' plot.rle <- mbecRLEPlot(rle.df, c('group','batch'))
mbecRLEPlot <- function(rle.df, model.vars, label=NULL) {

    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")
    n.tiles <- dim(unique(rle.df[,eval(model.vars[2])]))[1]

    if( n.tiles <= 3 ) {
        ncols <- n.tiles
        nrows <- 1
    } else {
        ncols <- 3
        nrows <- ceiling((n.tiles / 3))
    }

    rle.plot <- ggplot(
        rle.df,
        aes_string(x="specimen", y="values",fill=model.vars[1])) +
        stat_boxplot(color="black", notch=FALSE,
                              lwd=0.5, fatten=0.75,
                              outlier.colour="#E42032", outlier.fill="white",
                              outlier.shape=1, outlier.stroke=0.5,
                              outlier.size=0.5, outlier.alpha=0.5) +
        facet_wrap(~get(model.vars[2]), ncol=ncols, nrow=nrows ,
                            scales="free_x", drop=TRUE) +
        theme(plot.margin=unit(c(0.2,0.2,0.05,0.2), "cm")) +
        scale_fill_manual(values = mbecCols) +
        theme_bw() +
        labs(title=eval(label),
                      fill=mbecUpperCase(eval(model.vars[1]))) +
        # adjustments for the legend
        theme(legend.position="bottom",
                       legend.text=element_text(color="black",
                                                         size=12),
                       legend.key=element_rect(size=12),
                       axis.title.x=element_blank(),
                       panel.background=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank()) +
    guides(fill=guide_legend(title=element_blank()))

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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_wrap theme theme_bw
#' scale_fill_manual scale_y_continuous labs element_text element_blank
#' element_rect guides guide_legend geom_boxplot unit
#'
#' @export
#'
#' @examples
#' # This will return a list of the five most variable features grouped by the
#' # covariate 'batch'.
#' data(dummy.mbec)
#' box.df <- mbecBox(input.obj=dummy.mbec, method='TOP', n=5,
#' model.var='batch', type="otu", return.data=TRUE)
#' plot.box <- mbecBoxPlot(box.df[[1]], box.df[[2]], 'batch')
mbecBoxPlot <- function(tmp, otu.idx, model.var, label=NULL) {

    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")
    x.angle <- 0
    x.hjust <- 0.5
    density.lwd <- 0.2
    title.cex <- 1.5
    legend.cex <- 0.7
    legend.title.cex <- 0.75


    ret.plot <- list()

    for (idx in otu.idx) {
        # set-up legend title
        legend.title <- paste(idx, "-",
                              gsub("(^|[[:space:]])([[:alpha:]])","\\1\\U\\2",
                                   model.var, perl = TRUE), sep =" ")

        if( !is.null(label) ) {
            legend.title <- paste(label, "-", legend.title, sep=" ")
        }

        p.box <-
            ggplot(data=tmp, aes_string(x=model.var,
                                                     y=idx,
                                                     fill=model.var)) +
            stat_boxplot(geom="errorbar", width=0.4) +
            geom_boxplot() +
            scale_fill_manual(values=mbecCols) +
            ggplot2::scale_y_continuous(labels = function(x)
                format(x, scientific = TRUE)) +
            theme_bw() +
            theme(
                panel.background=element_blank(),
                axis.line=element_blank(),
                axis.ticks=element_blank(),
                axis.text.x=element_text(angle=x.angle,
                                                    hjust=x.hjust),
                panel.grid=element_blank(),
                axis.title.x=element_blank(),
                axis.text=element_text(size=8, hjust=0.5),
                axis.title=element_text(size=12),
                plot.title=element_text(hjust=0.5),
                axis.text.y = element_text(angle = 90),
                legend.position='bottom', legend.box='horizontal',
                legend.direction='horizontal') +
            labs(fill=legend.title, y="Value", title=idx)

        p.density <- ggplot(tmp, aes_string(x=idx,
                                                       fill=model.var)) +
            geom_density(alpha=0.5) +
            scale_fill_manual(values=mbecCols) +
            ggplot2::scale_y_continuous(labels = function(x)
                format(x, scientific = TRUE)) +
            labs(title=element_blank(), y="Density", fill=legend.title) +
            theme_bw() +
            theme(
                panel.background=element_blank(),
                axis.line=element_blank(),
                axis.ticks=element_blank(),
                axis.text.x=element_text(angle=x.angle, hjust=x.hjust),
                panel.grid=element_blank(),
                axis.title.x=element_blank(),
                axis.text=element_text(size=8, hjust=0.5),
                axis.title=element_text(size=12),
                plot.title=element_text(hjust=0.5),
                axis.text.y = element_text(angle = 90),
                legend.position='bottom', legend.box='horizontal',
                legend.direction='horizontal')

        g <- ggplotGrob(p.box)$grobs

        # legend <- g[[which(vapply(g, function(x) x$name,
        #                           FUN.VALUE = character(1)) == "guide-box")]]

        ret.plot[[eval(idx)]] <-
            gridExtra::arrangeGrob(p.density + theme(legend.position="none"),
                                   p.box + theme(legend.position="bottom",
                                                 plot.title=element_blank()),
                                   ncol=1, nrow=2, heights=c(5, 5.5))
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
#'
#' @importFrom pheatmap pheatmap
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' data(dummy.mbec)
#' heat.df <- mbecHeat(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' center=TRUE, scale=TRUE, method='TOP', n=5, return.data=TRUE)
#' plot.heat <- mbecHeatPlot(tmp.cnts=heat.df[[1]],
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
                           #annotation_colors = mbecCols,
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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_grid theme theme_bw
#' scale_fill_manual labs element_text element_blank element_rect guides
#' guide_legend geom_bar unit ggplot_build element_line ggplot_gtable
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' data(dummy.mbec)
#' mosaic.df <- mbecMosaic(input.obj=dummy.mbec, model.vars=c('group','batch'),
#' return.data=TRUE)
#' plot.mosaic <- mbecMosaicPlot(study.summary=mosaic.df,
#' model.vars=c('group','batch'))
mbecMosaicPlot <- function(study.summary, model.vars) {

    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")
    # local variable references to shut up check()
    Var1 <- Var2 <- NULL

    main_color <- "#004B5A"
    x.angle <- 0
    x.hjust <- 0.5
    density.lwd <- 0.2
    title.cex <- 1.5
    legend.cex <- 0.7
    legend.title.cex <- 0.75

    vars.axes <- mbecUpperCase(model.vars)

    plot.v2 <- ggplot(
        study.summary, aes_string(x="Var1", y="Freq.scaled",
                                    group="Var2", fill="Var1")) +
        facet_grid(cols=vars(Var2), scales="free",
                            space="free_x", drop=TRUE) +
        geom_bar(stat="identity", width=0.9) +
        guides(
            fill=guide_legend(title=eval(vars.axes[1]), reverse=TRUE,
                                       keywidth=1, keyheight=1)) +
        scale_fill_manual(values=mbecCols) +
        ylab("Proportion of all observations") +
        theme(axis.text.x=element_blank(),
                       axis.text.y=element_text(color=eval(main_color),
                                                         size=12),
                       axis.ticks=element_blank(),
                       axis.line=element_line(color="#7F7F7F"),
                       axis.title.x=element_blank(),
                       axis.title.y=element_text(
                           size=rel(1), angle=90),
                       legend.position='bottom', legend.box='horizontal',
                       legend.direction='horizontal',
                       legend.key.height=unit(0.2, 'cm'),
                       legend.key.width=unit(0.1, 'cm'),
                       legend.title=element_text(
                           size=rel(legend.title.cex)),
                       legend.spacing.x=unit(0.1, 'cm'),
                       legend.spacing.y=unit(0.1, 'cm'),
                       legend.text=element_text(
                           size=rel(legend.cex))) +
        theme(plot.margin=unit(c(0.2,0.2,0.05,0.2), "cm"))

    plot.v1 <- ggplot(study.summary, aes_string(x="Var2",
                                                           y="Freq.scaled",
                                                           fill="Var2")) +
        facet_grid(cols=vars(Var1), scales="free",
                            space="free_x", drop=TRUE) +
        geom_bar(stat="identity", width=0.9) +
        guides(fill=guide_legend(title=eval(vars.axes[2]),
                                                   reverse=TRUE,
                                                   keywidth=1, keyheight=1)) +
        scale_fill_manual(values=mbecCols) +
        ylab("Proportion of all observations") +
        theme(axis.text.x=element_blank(),
                       axis.text.y=element_text(
                           color=eval(main_color), size=12),
                       axis.ticks=element_blank(),
                       axis.line=element_line(color="#7F7F7F"),
                       axis.title.x=element_blank(),
                       axis.title.y=element_text(
                           size=rel(1), angle=90),
                       legend.position='bottom', legend.box='horizontal',
                       legend.direction='horizontal',
                       legend.key.height=unit(0.2, 'cm'),
                       legend.key.width=unit(0.1, 'cm'),
                       legend.title=element_text(
                           size=rel(legend.title.cex)),
                       legend.spacing.x=unit(0.1, 'cm'),
                       legend.spacing.y=unit(0.1, 'cm'),
                       legend.text=element_text(
                           size=rel(legend.cex))) +
        theme(plot.margin=unit(c(0.05,0.2,0.2,0.2), "cm"))

    ## Function to extract legend
    g_legend <- function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        # leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        # legend <- tmp$grobs[[leg]]
        legend <- tmp$grobs[[which(
            vapply(tmp$grobs, function(x) x$name, FUN.VALUE = character(1)) ==
                                       "guide-box")]]
    }

    legend.v2 <- g_legend(plot.v2)
    legend.v1 <- g_legend(plot.v1)

    mosaic.plot <-
        gridExtra::grid.arrange(plot.v2 +
                                    theme(legend.position="none"),
                                plot.v1 +
                                    theme(legend.position="none"),
                                gridExtra::grid.arrange(legend.v1, legend.v2,
                                                        ncol=2, nrow=1),
                                ncol=1, nrow=3, widths=c(1),
                                heights=c(1, 1,0.2), padding= -10)
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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_grid theme theme_bw
#' scale_fill_manual labs element_text element_blank element_rect guides
#' guide_legend geom_bar unit ggplot_build scale_shape_manual xlim ylim
#' coord_flip ggplotGrob geom_point scale_color_manual xlab ylab rel
#' geom_density
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' data(dummy.mbec)
#' pca.df <- mbecPCA(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), pca.axes=c(1,2), return.data=TRUE)
#' plot.pca <- mbecPCAPlot(plot.df=pca.df[[1]], metric.df=pca.df[[2]],
#' model.vars=c('group','batch'), pca.axes=c(1,2))
mbecPCAPlot <- function(plot.df, metric.df, model.vars, pca.axes, label=NULL) {
    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")
    x.angle <- 0
    x.hjust <- 0.5
    legend.cex <- 0.7
    legend.title.cex <- 0.75

    var.color <- model.vars[1]
    var.shape <- model.vars[2]
    label.col <- mbecUpperCase(model.vars[1])
    label.sha <- mbecUpperCase(model.vars[2])
    title <- paste("PCA:", label.sha, "-", label.col)

    x.label <- paste0(colnames(plot.df[pca.axes[1] + 1]), ": ",
                      metric.df$var.explained[pca.axes[1]], "% expl.var")
    y.label <- paste0(colnames(plot.df[pca.axes[2] + 1]), ": ",
           metric.df$var.explained[pca.axes[2]], "% expl.var")

    if( !is.null(label) ) {
        x.label <- paste(label, "-", x.label, sep=" ")
    }

    if (length(model.vars) >= 2) {
        pMain <-
            ggplot(
                data=plot.df,
                aes_string(x=colnames(plot.df[pca.axes[1] + 1]),
                             y = colnames(plot.df[pca.axes[2] + 1]),
                             colour=var.color, shape=var.shape)) +
            scale_shape_manual(
                values=c(0,1,2,3,6,8,15,16,17,23,25,4,5,9)) +
            geom_point() +
            scale_color_manual(values=mbecCols) +
            labs(title=element_blank(), colour=label.col, shape=label.sha) +
            xlim(metric.df$axis.min[pca.axes[1]],
                          metric.df$axis.max[pca.axes[1]]) +
            ylim(metric.df$axis.min[pca.axes[2]],
                          metric.df$axis.max[pca.axes[2]]) +
            xlab(x.label) + ylab(y.label) + theme_bw() +
            theme(
                panel.background=element_blank(),
                axis.line=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                legend.position='right', legend.box='horizontal',
                legend.direction='vertical')

        pTop <-
            ggplot(data=plot.df,
                            aes_string(
                                x=colnames(plot.df[pca.axes[1]+1]),
                                fill=var.color, linetype=var.shape)) +
            geom_density(size=0.2, alpha=0.5) +
            ylab("Density") +
            scale_fill_manual(values=mbecCols) +
            xlim(metric.df$axis.min[pca.axes[1]],
                          metric.df$axis.max[pca.axes[1]]) +
            theme_bw() +
            theme(
                panel.background=element_blank(),
                axis.line=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                legend.position='right', legend.box='horizontal',
                legend.direction='vertical') +
            labs(title=element_blank()) +
            theme(axis.title.x=element_blank())

        pRight <-
            ggplot(data=plot.df,
                            aes_string(
                                x = colnames(plot.df[pca.axes[2]+1]),
                                fill=var.color, linetype=var.shape)) +
            geom_density(size=0.2,alpha=0.5) + coord_flip() +
            ylab("Density") +
            scale_fill_manual(values=mbecCols) +
            xlim(metric.df$axis.min[pca.axes[2]],
                          metric.df$axis.max[pca.axes[2]]) +
            theme_bw() +
            theme(
                panel.background=element_blank(),
                axis.line=element_blank(),
                axis.text=element_blank(),
                axis.ticks=element_blank(),
                legend.position='right', legend.box='horizontal',
                legend.direction='vertical') +
            labs(title=element_blank()) +
            theme(axis.title.y=element_blank())

    } else {
        pMain <-
            ggplot(data=plot.df,
                            aes_string(
                                x=colnames(plot.df[pca.axes[1] + 1]),
                                y=colnames(plot.df[pca.axes[2] + 1]),
                                colour=var.color)) +
            geom_point() +
            scale_color_manual(values=mbecCols) +
            theme_bw() +
            labs(title=element_blank(), colour=label.col, shape=label.sha) +
            xlim(metric.df$axis.min[pca.axes[1]],
                          metric.df$axis.max[pca.axes[1]]) +
            ylim(metric.df$axis.min[pca.axes[2]],
                          metric.df$axis.max[pca.axes[2]]) +
            xlab(x.label) + ylab(y.label)

        pTop <- ggplot(data=plot.df,
                                aes_string(
                                    x=colnames(plot.df[pca.axes[1]+1]),
                                    fill=var.color)) +
            geom_density(size = 0.2, alpha=0.5) +
            ylab("Density") +
            scale_fill_manual(values=mbecCols) +
            xlim(metric.df$axis.min[pca.axes[1]],
                          metric.df$axis.max[pca.axes[1]]) +
            theme_bw() +
            labs(title=element_blank()) +
            theme(axis.title.x=element_blank())

        pRight <-
            ggplot(data=plot.df,
                            aes_string(
                                x=colnames(plot.df[pca.axes[2]+1]),
                                fill=var.color)) +
            geom_density(size=0.2, alpha=0.5) +
            coord_flip() + ylab("Density") +
            scale_fill_manual(values=mbecCols) +
            xlim(metric.df$axis.min[pca.axes[2]],
                          metric.df$axis.max[pca.axes[2]]) +
            theme_bw() +
            labs(title=element_blank()) +
            theme(axis.title.y=element_blank())
    }

    g <- ggplotGrob(pMain)$grobs
    legend <-
        g[[which(vapply(g, function(x) x$name,
                        FUN.VALUE = character(1)) == "guide-box")]]
    ret.plot <-
        gridExtra::grid.arrange(pTop + theme(legend.position = "none"),
                                legend, pMain +
                                    theme(legend.position = "none"),
                                pRight +
                                    theme(legend.position = "none"),
                                ncol=2, nrow=2, widths=c(3, 1),
                                heights=c(1.5, 3))

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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_grid theme theme_bw
#' scale_fill_manual labs element_text element_blank element_rect guides
#' guide_legend geom_bar unit ggplot_build scale_shape_manual xlim ylim
#' coord_flip ggplotGrob vars
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessments.
#' data(dummy.mbec)
#' df.var.lm <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), method='lm', type='clr')
#' plot.lm <- mbecVarianceStatsPlot(variance.obj=df.var.lm)
mbecVarianceStatsPlot <- function(variance.obj) {

    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")

    # local variable references to shut up check()
    type <- variance <- NULL

    plot.df <- variance.obj %>%
        dplyr::bind_rows() %>%
        tidyr::gather("covariate",
                      "variance", -type) %>%
        dplyr::mutate(type=factor(type, levels=unique(type))) %>%
        dplyr::mutate(variance=as.numeric(as.character(variance)))

    leplot <-
        ggplot(plot.df, aes_string(x="covariate",
                                              y="variance", fill="covariate")) +
        geom_boxplot(lwd=0.5, fatten=0.75, outlier.colour="#E42032",
                              outlier.fill="white", outlier.shape=1,
                              outlier.stroke=0.5, outlier.size=0.5,
                              outlier.alpha=0.5) +
        facet_grid(cols=vars(type)) +
        scale_fill_manual(values = mbecCols) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1),
                       strip.text=element_text(size=12),
                       panel.grid=element_blank(),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=15),
                       legend.title=element_text(size=15),
                       legend.text=element_text(size=12)) +
        labs(x="Linear (Mixed) Model", y="Proportion of Variance",
                      name="Covariate") +
        ylim(0, 1) +
        facet_grid(cols=vars(type), scales="free",
                            space="free_x", drop=TRUE)

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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_grid theme theme_bw
#' scale_fill_manual labs element_text element_blank element_rect guides
#' guide_legend geom_bar unit ggplot_build scale_shape_manual xlim ylim
#' coord_flip ggplotGrob vars geom_text position_dodge
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for three variance
#' # assessments.
#' data(dummy.mbec)
#' df.var.rda <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('group','batch'), method='rda', type='clr')
#' plot.rda <- mbecRDAStatsPlot(rda.obj=df.var.rda)
mbecRDAStatsPlot <- function(rda.obj) {

    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")
    # local variable references to shut up check()
    type <- variance <- NULL

    leTest <- rda.obj %>%
        tidyr::gather("covariate", "variance", -type) %>%
        dplyr::mutate(type=factor(type, levels=unique(type))) %>%
        dplyr::mutate(variance=as.numeric(as.character(variance))) %>%
        dplyr::mutate(variance.r=round(variance, 2)) %>%
        dplyr::mutate(variance.offset=variance + 2.5)

    lePlot <-
        ggplot(data=leTest,
               aes_string(x="covariate", y="variance", fill="covariate")) +
        geom_bar(stat="identity",position="dodge", colour="black") +
        geom_text(data=leTest,
                  aes_string("covariate", "variance.offset",
                             label="variance.r"),
                  position = position_dodge(width=0.9), size=3) +
        facet_grid(cols=vars(type)) +
        scale_fill_manual(values = mbecCols) +
        theme_bw() +
        labs(x="RDA",y="Variance explained (%)") +
        theme(axis.text.x=element_text(angle=60,hjust=1),
                       panel.grid=element_blank(),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=15),
                       legend.title=element_text(size=15),
                       legend.text=element_text(size=12)) +
        ylim(0, 100) +
        facet_grid(cols=vars(type), scales="free",
                            space="free_x", drop=TRUE)

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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_grid theme theme_bw
#' scale_fill_manual labs element_text element_blank element_rect guides
#' guide_legend geom_bar unit ggplot_build scale_shape_manual xlim ylim
#' coord_flip ggplotGrob vars
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' data(dummy.mbec)
#' df.var.pvca <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('batch','group'), method='pvca', type='clr')
#' plot.pvca <- mbecPVCAStatsPlot(pvca.obj=df.var.pvca)
mbecPVCAStatsPlot <- function(pvca.obj) {

    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")
    # local variable references to shut up check()
    type <- covariate <- variance <- variance.p <- NULL

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
                                             100, 2)) %>%
        dplyr::mutate(variance.offset=variance.p + 2.5)

    lePlot <- ggplot(data = plot.df,
                              aes_string(x = "covariate", y = "variance.p",
                                           fill = "covariate")) +
        geom_bar(stat="identity", position="dodge", colour="black") +
        scale_fill_manual(values = mbecCols) +
        geom_text(data=plot.df,
                           aes_string(
                               "covariate", "variance.offset",
                               label="variance.p"),
                           position=position_dodge(width=0.9),
                           size = 3) + theme_bw() +
        facet_grid(cols=vars(type), scales="free",
                            space="free_x", drop=TRUE) +
        labs(x = "PVCA - Random effects and Interactions",
                      y = "Variance explained (%)") +
        theme(
            axis.text.x=element_text(angle=60,
                                              hjust=1),
            panel.grid=element_blank(),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=15),
                       legend.title=element_text(size=15),
                       legend.text=element_text(size=12)) +
        ylim(0, 100) +
        facet_grid(cols=vars(type), scales="free",
                            space="free_x", drop=TRUE)

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
#'
#' @importFrom ggplot2 ggplot aes_string stat_boxplot facet_grid theme theme_bw
#' scale_fill_manual labs element_text element_blank element_rect guides
#' guide_legend geom_bar unit ggplot_build scale_shape_manual xlim ylim
#' coord_flip ggplotGrob vars
#'
#' @export
#'
#' @examples
#' # This will return a paneled plot that shows results for the variance
#' # assessment.
#' data(dummy.mbec)
#' df.var.scoef <- mbecModelVariance(input.obj=dummy.mbec,
#' model.vars=c('batch','group'), method='s.coef', type='clr')
#' plot.scoef <- mbecSCOEFStatsPlot(scoef.obj=df.var.scoef)
mbecSCOEFStatsPlot <- function(scoef.obj) {

    # local variable references to shut up check()
    variable <- type <- sil.coefficient <- NULL

    mbecCols <- c("#9467bd","#BCBD22","#2CA02C","#E377C2","#1F77B4","#FF7F0E",
                  "#D62728","#8C564B","#E377C2","#7F7F7F","#17BECF")
    # first tidy-magic to create df for plotting
    plot.df <- scoef.obj %>%
        dplyr::mutate(variable=gsub("\\.", ":", variable)) %>%
        dplyr::mutate(type=factor(type, levels = unique(type))) %>%
        dplyr::mutate(sil.coefficient=
                          as.numeric(as.character(sil.coefficient))) %>%
        dplyr::mutate(sil.coefficient.r=round(sil.coefficient, 2))

    # now plot
    lePlot <-
        ggplot(plot.df, aes_string(
            x="variable", y="sil.coefficient", color="cluster",
            shape="variable")) +
        geom_point() + facet_grid(cols=vars(type)) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1),
                       strip.text=element_text(size=12),
                       panel.grid=element_blank(),
                       axis.text=element_text(size=10),
                       axis.title=element_text(size=15),
                       legend.title=element_text(size=15),
                       legend.text=element_text(size=12)) +
        scale_color_manual(values=mbecCols) +
        labs(x="Silhouette Coefficient", y="Grouping", name="Type") +
        facet_grid(cols=vars(type), scales="free",
                            space="free_x", drop=TRUE)

    return(lePlot)
    # FIN
}

