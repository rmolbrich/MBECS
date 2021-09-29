# HELPER FUNCTIONS --------------------------------------------------------

#' DEPRECATED
#' Helper function to put together, display and save plots that were create with MBEC.
#' @param data.obj, any output of MBEC's plotting functions.
#' @param p.type, which kind of plot to expect
#' @param f.title, string filename activates plotting if provided
#' @export
mbecPlots <- function(data.obj, p.type=NULL, grouping="batch", f.title=NULL, save.plot=FALSE) {

  ## Because 'f.title' is file- and plot-title, we need to create sth. if no user input is provided

  ### ToDo: make titleing great again

  if( is.null(f.title) ) {
    p.title <- paste(p.type, grouping[1], sep = "_")
  } else {
    p.title <- f.title
  }

  cols <- cols[c(1,3,5,7,9,11,13,15,17,19,2,4,6,8,10,12,14,16,18,20)]

  # ToDo: different groupings for RLE plot ?!
  if( p.type == "b_rle" ) {
    ### RLE-PLOT
    ret.plot <- ggplot(data.obj, aes(x = sample, y = values, fill = get(grouping))) +
      stat_boxplot(color="black",notch = TRUE,
                   outlier.colour = "#E42032", outlier.fill = "white",outlier.shape = 1, outlier.stroke = .5) +
      #facet_wrap(~Strain, ncol=2) +
      facet_grid(cols=vars(group), scales="free", space="free_x", drop=T) +
      scale_fill_manual(values = cols) +
      theme_new +
      guides(fill=guide_legend(title=element_blank())) +
      theme(axis.title.x = element_blank(),
            panel.background = element_blank(),
            #axis.text.x=element_blank(),
            axis.ticks.x=element_blank())


  } else if( p.type == "b_box" ) {
    ### BOX-PLOT
    ## ToDo: get rid of this
    x.angle = 0
    x.hjust = 0.5
    density.lwd = 0.2
    title.cex = 1.5
    legend.cex = 0.7
    legend.title.cex =0.75
    batch.legend.title <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",grouping,perl = TRUE)

    ret.plot <- list()

    for( idx in data.obj[[2]] ) {
      ## prepare the plots
      p.box <- ggplot(data = data.obj[[1]], aes(x = get(grouping), y = get(idx), fill = get(grouping))) + stat_boxplot(geom = "errorbar", width = 0.4) +
        geom_boxplot() + scale_fill_manual(values = cols) + theme_bw() +
        theme(axis.text.x = element_text(angle = x.angle, hjust = x.hjust), panel.grid = element_blank(),
              axis.title.x = element_blank(), axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              plot.title = element_text(hjust = 0.5)) +
        labs(fill = batch.legend.title, y = 'value',title = idx)

      p.density <- ggplot(data.obj[[1]], aes(x = get(idx), fill = get(grouping))) +
        geom_density(alpha = 0.5) + scale_fill_manual(values = cols) +
        labs(title = idx, x = 'Value', fill = batch.legend.title) +
        theme_bw() + theme(plot.title = element_blank(),
                           panel.grid = element_blank())

      ## Put the plots in grid for plotting
      # modify legend
      g <- ggplotGrob(p.box + theme(legend.position = 'bottom', legend.box = 'horizontal',
                                    legend.direction = 'horizontal',
                                    legend.key.height = unit(0.8, 'cm'),
                                    legend.key.width = unit(0.4, 'cm'),
                                    legend.title = element_text(size = rel(legend.title.cex)),
                                    legend.spacing.x = unit(0.4, 'cm'),
                                    legend.spacing.y = unit(0.4, 'cm'),
                                    legend.text = element_text(size = rel(legend.cex))))$grobs
      # extract legend
      legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

      # put plot into the list
      ret.plot[[eval(idx)]] <- arrangeGrob(p.box + theme(legend.position = 'none'),
                                           p.density + theme(legend.position = 'none'),
                                           legend,
                                           ncol = 1, nrow = 3, heights = c(5,5,1))
    }
  } else if( p.type == "b_pca") {
    ### PCA-PLOT
    ## ToDo: get rid of this
    x.angle = 0
    x.hjust = 0.5
    density.lwd = 0.2
    title.cex = 1.5
    legend.cex = 0.7
    legend.title.cex =0.75
    # change first letter of group denominator to uppercase for pretty-plotting
    batch.legend.title <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",grouping[1],perl = TRUE)
    trt.legend.title <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",grouping[2],perl = TRUE)
    title <- p.title

    plot.df <- data.obj[[1]]
    metric.df <- data.obj[[2]]
    plot.axes <- data.obj[[3]]

    if(length(grouping) >= 2){
      pMain <- ggplot(data = plot.df, aes(x = get(colnames(plot.df[plot.axes[1]+1])), y = get(colnames(plot.df[plot.axes[2]+1])), colour = get(grouping[1]), shape = get(grouping[2]))) +
        geom_point() +
        xlab(paste0(colnames(plot.df[plot.axes[1]+1]), ': ', metric.df$var.explained[plot.axes[1]], '% expl.var')) +
        ylab(paste0(colnames(plot.df[plot.axes[2]+1]), ': ', metric.df$var.explained[plot.axes[2]], '% expl.var')) +
        scale_color_manual(values = cols) + theme_bw() +
        xlim(metric.df$axis.min[plot.axes[1]], metric.df$axis.max[plot.axes[1]]) +
        ylim(metric.df$axis.min[plot.axes[2]], metric.df$axis.max[plot.axes[2]]) +
        labs(colour = batch.legend.title, shape = trt.legend.title)

      pTop <- ggplot(data = plot.df, aes(x = get(colnames(plot.df[plot.axes[1]+1])), fill = get(grouping[1]), linetype = get(grouping[2]))) +
        geom_density(size = density.lwd, alpha = 0.5) + ylab('Density') +
        theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)),
              plot.title = element_text(hjust = 0.5, size = rel(title.cex)),
              axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) +
        scale_fill_manual(values = cols) +
        xlim(metric.df$axis.min[plot.axes[1]], metric.df$axis.max[plot.axes[1]]) + labs(title = title)

      pRight <- ggplot(data = plot.df, aes(x=get(colnames(plot.df[plot.axes[2]+1])), fill = get(grouping[1]), linetype = get(grouping[2]))) +
        geom_density(size = density.lwd,alpha = 0.5) +  coord_flip() + ylab('Density') +
        theme(axis.title.x = element_text(size = rel(0.8)),
              axis.title.y = element_blank(), axis.line = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) +
        scale_fill_manual(values = cols) +
        xlim(metric.df$axis.min[plot.axes[2]], metric.df$axis.max[plot.axes[2]])

    }else{
      pMain <- ggplot(data = plot.df, aes(x = get(colnames(plot.df[plot.axes[1]+1])), y = get(colnames(plot.df[plot.axes[2]+1])), colour = get(grouping[1]))) +
        geom_point(shape = 16) +
        xlab(paste0(colnames(plot.df[plot.axes[1]+1]), ': ', metric.df$var.explained[plot.axes[1]], '% expl.var')) +
        ylab(paste0(colnames(plot.df[plot.axes[2]+1]), ': ', metric.df$var.explained[plot.axes[2]], '% expl.var')) +
        scale_color_manual(values = cols) + theme_bw() +
        xlim(metric.df$axis.min[plot.axes[1]], metric.df$axis.max[plot.axes[1]]) +
        ylim(metric.df$axis.min[plot.axes[2]], metric.df$axis.max[plot.axes[2]]) +
        labs(colour = batch.legend.title)

      pTop <- ggplot(data = plot.df, aes(x = PC1, fill = get(grouping[1]))) +
        geom_density(size = density.lwd, alpha=0.5) + ylab('Density') +
        theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)),
              plot.title = element_text(hjust = 0.5, size = rel(title.cex)),
              axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) + scale_fill_manual(values = cols) +
        xlim(metric.df$axis.min[plot.axes[1]], metric.df$axis.max[plot.axes[1]]) + labs(title = title)

      pRight <- ggplot(data = plot.df, aes(x=PC2, fill = get(grouping[1]))) +
        geom_density(size = density.lwd, alpha = 0.5) +  coord_flip() + ylab('Density') +
        theme(axis.title.x = element_text(size = rel(0.8)),
              axis.title.y = element_blank(), axis.line = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              panel.background = element_blank()) + scale_fill_manual(values = cols) +
        xlim(metric.df$axis.min[plot.axes[2]], metric.df$axis.max[plot.axes[2]])
    }


    g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'horizontal',
                                  legend.direction = 'vertical',
                                  legend.key.height = unit(0.2, 'cm'),
                                  legend.key.width = unit(0.1, 'cm'),
                                  legend.title = element_text(size = rel(legend.title.cex)),
                                  legend.spacing.x = unit(0.1, 'cm'),
                                  legend.spacing.y = unit(0.1, 'cm'),
                                  legend.text = element_text(size = rel(legend.cex))))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

    ret.plot <- gridExtra::grid.arrange(pTop + theme(legend.position = 'none'), legend, pMain +
                                          theme(legend.position = 'none'), pRight + theme(legend.position = 'none'),
                                        ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))

  } else if( p.type == "b_heat" ) {
    ### HEATMAP-PLOT

    ret.plot <- pheatmap::pheatmap(data.obj[[1]],
                                   scale = 'none',
                                   cluster_rows = F,
                                   cluster_cols = T,
                                   fontsize_row = 4, fontsize_col = 6,
                                   fontsize = 8,
                                   clustering_distance_rows = 'euclidean',
                                   clustering_method = 'ward.D',
                                   treeheight_row = 30,
                                   annotation_col = data.obj[[2]][,eval(grouping)],
                                   #annotation_colors = ad.anno_metabo_colors,
                                   border_color = 'NA',
                                   main = p.title)

  } else if( p.type == "b_hcl" ) {
    ### DENDROGRAM-PLOT

    mar.right <- (4 * dim(data.obj[[2]])[2]) + 1
    par(mar = c(4,1,1,mar.right))
    plot(data.obj[[1]], horiz = TRUE)
    dendextend::colored_bars(data.obj[[2]], data.obj[[1]], rowLabels = names(data.obj[[2]]), horiz = TRUE)
    legend("topleft", legend = useries, fill = tmp.cols, bg="white", cex=0.6)

    ret.plot <- recordPlot()

  } else if( p.type == "b_mosaic" ) {

    ret.plot <- ggplot2::ggplot(data = data.obj[[2]]) +
      ggmosaic::geom_mosaic(aes(x = ggmosaic::product(batch, group), fill = batch)) +
      scale_fill_manual(values = cols) + theme_bw() +
      labs(y="Batch", x="Experimental group", title = "Experimental groups to batches")
  }

  if( save.plot ) {

    if( p.type == "b_rle" ) {
      f.title <- paste("mbecRLE",f.title,sep="_")
    } else if( p.type == "b_box" ) {
      f.title <- paste("mbecBOX",f.title,sep="_")
    } else if( p.type == "b_pca" ) {
      f.title <- paste("mbecPCA",f.title,sep="_")
    } else if( p.type == "b_pca" ) {
      f.title <- paste("mbecHEAT",f.title,sep="_")
    } else if( p.type == "b_hcl" ) {
      f.title <- paste("mbecHCL",f.title,sep="_")
    } else if( p.type == "p_mosaic" ) {
      f.title <- paste("mbecMOSAIC",f.title,sep="_")
    }
    mbecSavePlot(ret.plot, f.title=f.title)

  }

  return(ret.plot)
}

#' DEPRECATED
#' Helper function to save plots as pdf (for now maybe..) in given path with given title.
#' @param plot.obj, output of 'mbecPlots()'
#' @param f.title, filename for output, defaults to 'mbecTest'
#' @param f.path, save directory, defaults to workdir
#' @export
mbecSavePlot <- function(plot.obj, f.title="mbecTest") {

  # ToDo:
  # a <- "old"
  # test <- function () {
  #   assign("a", "new", envir = .GlobalEnv)
  # }
  # test()
  # a  # display the new value

  if( any(class(plot.obj) %in% c("gtable","grob","ggplot","gg","pheatmap")) ) {

    pdf(file = file.path(PLOTDIR, paste0(f.title,".pdf")), width=9, height=9, fonts="Helvetica", onefile = T)
    grid.draw(plot.obj)
    dev.off()

  } else if( class(plot.obj) %in% c("list") ) {

    pdf(file = file.path(PLOTDIR, paste0(f.title,".pdf")), width=9, height=9, fonts="Helvetica", onefile = T)
    lapply(plot.obj, function(x) {grid.newpage(); grid.draw(x)})
    dev.off()

  } else if( class(plot.obj) %in% "recordedplot" ) {

    pdf(file = file.path(PLOTDIR, paste0(f.title,".pdf")), width=9, height=6, fonts="Helvetica", onefile = T)
    plot.obj
    dev.off()


  } else {
    stop(paste("Object of class: ", class(plot.obj)," is not supported by 'mbecSavePlots()'.", sep=""))
  }
}


#' Linear (Mixed) Model Feature to Batch Fit
#'
#' Helper function that fits lm/lmm with covariates 'treatment' and 'batch' to every feature in the
#' data-set. Returns the fdr corrected significance value for the "treatment" variable. The method
#' 'lm' will fit the linear model y ~ model.vars[1] + model.vars[2] and the linear mixed model
#' will consider the second term as random effect, i.e., y ~ model.vars[1] + (1|model.vars[2]).
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
#' This will return p-value for the linear model fit
#' \dontrun{val.score <- mbecLM(input.obj, model.vars=c("group","batch"), method="lm")}
mbecLM <- function(input.obj, model.vars=c("group","batch"), method=c("lm","lmm")) {
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

    tmp.group.p <- apply(tmp.cnts, 2, FUN = function(x){
      nc.lmm <- lmerTest::lmer(x ~ tmp.meta[,model.vars[1]] + (1|tmp.meta[,model.vars[2]]), data = tmp.meta)
      nc.lmm.summary <- summary(nc.lmm)
      p <- nc.lmm.summary$coefficients[2,5]
    })
  }

  # correct for multiple testing
  tmp.group.p <- p.adjust(tmp.group.p, method = 'fdr')

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
#' @param input.obj, either pyhloseq-object (OTU orientation is handled) or numeric matrix (samples x features)
#' @param method, one of 'CLR' or 'ILR'
#' @param offset, optional offset in case of sparse matrix
#' @return MbecDataObject with transformed counts
#' @export
#' @include mbecs_classes.R
#'
#' @examples
#' This will return the cumulative log-ratio tansformed counts in an MbecData object
#' \dontrun{mbec.LRT <- LRTransform(input.obj=list(counts, covariates),
#' method="CLR", offset=0)}
#'
#' This will return the inverse log-ratio transformed counts in an MbecData object
#' \dontrun{mbec.LRT <- LRTransform(input.obj=list(counts, covariates),
#' method="ILR", offset=0)}
LRTransform <- function(input.obj, meta.obj=NULL, method = c("none", "CLR", "ILR"), offset = 0) {
  ## 00. Check if 'method' was chosen correctly.
  method <- match.arg(method)

  ## VALIDATE input and change to 'MbecData' if needed
  input.obj <- mbecProcessInput(input.obj, meta.obj, required.col=eval(grouping))

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
#' This will return a matrix of normalised counts, according to the covariate information in meta
#' \dontrun{mtx.pn_counts <- percentileNorm(cnts=mtx_of_cnts, meta=grouping_info)}
percentileNorm <- function(cnts, meta) {

  ref.group <- levels(meta[,1])[1]
  message(paste0("Group ",ref.group, " is considered control group, i.e., reference for normalization procedure. To change reference please 'relevel()' grouping factor accordingly."))

  norm.cnts <- cnts; norm.cnts[,] <- NA

  # for every batch
  for( b.idx in levels(meta[,2]) ) {
    # for every feature
    for( f.idx in 1:ncol(cnts) ) {
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
#' are smaller than the reference (right). Percentiles of scores are given in the interval I:[0,100].
#' Depending on type of calculation, the score will be computed as follows:
#'
#' \dontrun{rank = (right + left + ifelse(right > left, 1, 0)) * 50/n}
#'
#' \dontrun{weak = right / n*100}
#'
#' \dontrun{strict = left / n*100}
#'
#' \dontrun{mean = (right + left) * 50/n)}
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
#' This will return a score for the supplied vector with default evaluation (strict).
#' \dontrun{val.score <- poscore(cnt.vec=ref.vec, cnt=adjust.vec, type="strict")}
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

# KA changed the function to add a min value when many zeroes in data (prob with log and division by 0 otherwise)
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
    for (i in 1 : ncol(x.ilr))
    {
      x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, (i+1) : D, drop = FALSE]), 1, prod) + offset)^(1 / (D-i))) / (x[,i]+ offset))
      #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop = FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
    }
  } else {
    for (i in 1 : ncol(x.ilr))
    {
      x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(apply(as.matrix(x[, (i+1):D]), 1, function(x){exp(log(x))})/(x[, i]+ offset) + offset)
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
  for( i in 1:ncol(V) )
  {
    V[1:i, i] = 1/i
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












