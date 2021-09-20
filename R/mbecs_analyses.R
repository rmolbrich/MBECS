# ANALYSIS FUNCTIONS ------------------------------------------------------

#' Relative Log Expression Plot
#'
#' Takes two covariates, i.e., group and batch, and computes the RLE-plot the grouping of the first covariate,
#' colored by the second covariate. Effectively illustrating the relative expression between samples from different
#' batches within the respective study groups. Other covariates can be chosen as input and the function will check for
#' factors and convert if necessary. Categorical factors, e.g., group membership, sex and batch, produce the best result.
#' The function returns either a plot-frame or the finished ggplot object. Input for th data-set can be an MbecData-object,
#' a phyloseq-object or a list that contains counts and covariate data. The covariate table requires an 'sID' column that
#' contains sample IDs equal to the sample naming in the counts table. Correct orientation of counts will be handled internally.
#'
#' @keywords RLE relative log expression
#' @param input.obj, list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param model.vars, two covariates of interest to select by first variable selects panels and second one determines coloring
#' @param return.data, logical if TRUE returns the data.frame required for plotting (NO plotting here bucko)
#' @export
#'
#' @examples
#' This will return the data.frame for plotting.
#' \dontrun{p.RLE <- mbecRLE(input.obj=list(counts, covariates), model.vars=c("treatment","batches"), return.data=TRUE)}
#'
#' This will return the ggplot2 object for display, saving and modification.
#' \dontrun{p.RLE <- mbecRLE(input.obj=phyloseq, model.vars=c("treatment","sex"), return.data=FALSE)}
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

  rle.plot <- ggplot(tmp.long, aes(x = specimen, y = values, fill = get(model.vars[2]))) +
    stat_boxplot(color="black",notch = TRUE,
                 outlier.colour = "#E42032", outlier.fill = "white",outlier.shape = 1, outlier.stroke = .5) +
    #facet_wrap(~Strain, ncol=2) +
    facet_grid(cols=vars(get(model.vars[1])), scales="free", space="free_x", drop=T) +
    scale_fill_manual(values = cols) +
    theme_rle() +
    guides(fill=guide_legend(title=element_blank()))

  return(rle.plot)
}



#' Creates nice PCA plot with axis-density graphs
#' @param input.obj, list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param pca.axes, numeric vector which axes to plot, first is X and second is Y
#' @param return.data, logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @include mbecs_classes.R
setGeneric("mbecPCA", signature="input.obj",
           function(input.obj, model.vars=c("group","batch"), pca.axes=c(1,2), return.data=FALSE)
             standardGeneric("mbecPCA")
)

## In this form it works for 'phyloseq' and 'mbecData' objects
.mbecPCA <- function(input.obj, model.vars=c("group","batch"), pca.axes=c(1,2), return.data=FALSE) {

  tmp <- mbecGetData(input.obj, orientation="sxf", required.col=eval(model.vars))
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  # calculate IQR and sort counts in decreasing order
  iqr <- apply(tmp.cnts,2,IQR)
  tmp.cnts <- tmp.cnts[,order(iqr,decreasing=T)]

  PCA <- prcomp(tmp.cnts, scale = F)

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
    pMain <- ggplot(data = plot.df, aes(x = get(colnames(plot.df[pca.axes[1]+1])), y = get(colnames(plot.df[pca.axes[2]+1])), colour = get(model.vars[2]), shape = get(model.vars[1]))) +
      geom_point() +
      scale_color_manual(values = cols) +
      labs(colour = var.color, shape = var.shape) +
      xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      ylim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      xlab(paste0(colnames(plot.df[pca.axes[1]+1]), ': ', metric.df$var.explained[pca.axes[1]], '% expl.var')) +
      ylab(paste0(colnames(plot.df[pca.axes[2]+1]), ': ', metric.df$var.explained[pca.axes[2]], '% expl.var')) +
      theme_pca()

    pTop <- ggplot(data = plot.df, aes(x = get(colnames(plot.df[pca.axes[1]+1])), fill = get(model.vars[2]), linetype = get(model.vars[2]))) +
      geom_density(size = 0.2, alpha = 0.5) + ylab('Density') +
      scale_fill_manual(values = cols) +
      xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      theme_pca() +
      labs(title = title) +
      theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)),
            plot.title = element_text(hjust = 0.5, size = rel(1.5)))

    pRight <- ggplot(data = plot.df, aes(x=get(colnames(plot.df[pca.axes[2]+1])), fill = get(model.vars[2]), linetype = get(model.vars[2]))) +
      geom_density(size = 0.2,alpha = 0.5) +  coord_flip() + ylab('Density') +
      scale_fill_manual(values = cols) +
      xlim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      theme_pca() +
      theme(axis.title.x = element_text(size = rel(0.8)),
            axis.title.y = element_blank(), axis.line = element_blank(),
            plot.title = element_blank())

  }else{
    pMain <- ggplot(data = plot.df, aes(x = get(colnames(plot.df[pca.axes[1]+1])), y = get(colnames(plot.df[pca.axes[2]+1])), colour = get(model.vars[1]))) +
      geom_point() +
      scale_color_manual(values = cols) +
      labs(colour = var.color, shape = var.shape) +
      xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      ylim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      xlab(paste0(colnames(plot.df[pca.axes[1]+1]), ': ', metric.df$var.explained[pca.axes[1]], '% expl.var')) +
      ylab(paste0(colnames(plot.df[pca.axes[2]+1]), ': ', metric.df$var.explained[pca.axes[2]], '% expl.var')) +
      theme_pca()

    pTop <- ggplot(data = plot.df, aes(x = get(colnames(plot.df[pca.axes[1]+1])), fill = get(model.vars[1]))) +
      geom_density(size = 0.2, alpha = 0.5) + ylab('Density') +
      scale_fill_manual(values = cols) +
      xlim(metric.df$axis.min[pca.axes[1]], metric.df$axis.max[pca.axes[1]]) +
      theme_pca() +
      labs(title = title) +
      theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)),
            plot.title = element_text(hjust = 0.5, size = rel(1.5)))

    pRight <- ggplot(data = plot.df, aes(x=get(colnames(plot.df[pca.axes[2]+1])), fill = get(model.vars[1]))) +
      geom_density(size = 0.2,alpha = 0.5) +  coord_flip() + ylab('Density') +
      scale_fill_manual(values = cols) +
      xlim(metric.df$axis.min[pca.axes[2]], metric.df$axis.max[pca.axes[2]]) +
      theme_pca() +
      theme(axis.title.x = element_text(size = rel(0.8)),
            axis.title.y = element_blank(), axis.line = element_blank(),
            plot.title = element_blank())
  }

  g <- ggplotGrob(pMain)$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

  ret.plot <- gridExtra::grid.arrange(pTop + theme(legend.position = 'none'), legend, pMain +
                                        theme(legend.position = 'none'), pRight + theme(legend.position = 'none'),
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

#' Produces box-plot + density for OTUs in phyloseq object.
#' @param input.obj, list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param method, one of 'ALL' or 'TOP' for 'n' most variable features, DEFAULT is 'ALL'
#' @param n, number of OTUs to display for 'TOP' method
#' @param model.var, covariate to group by, default is batch
#' @param return.data, logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @export
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
    iqr <- apply(tmp[,otu.idx],2,IQR)
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
    p.box <- ggplot(data = tmp, aes(x = get(model.var), y = get(idx), fill = get(model.var))) + stat_boxplot(geom = "errorbar", width = 0.4) +
      geom_boxplot() + scale_fill_manual(values = cols) + theme_bw() +
      theme_box() +
      labs(fill = legend.title, y = 'value',title = idx)

    p.density <- ggplot(tmp, aes(x = get(idx), fill = get(model.var))) +
      geom_density(alpha = 0.5) + scale_fill_manual(values = cols) +
      labs(title = idx, x = 'Value', fill = legend.title) +
      theme_box()

    ## Put the plots in grid for plotting
    # modify legend
    g <- ggplotGrob(p.box)$grobs
    # extract legend
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

    # put plot into the list
    ret.plot[[eval(idx)]] <- arrangeGrob(p.box + theme(legend.position = 'none'),
                                         p.density + theme(legend.position = 'none', plot.title = element_blank()),
                                         legend,
                                         ncol = 1, nrow = 3, heights = c(5,4.5,1))
  }

  return(ret.plot)
}


#' Produces box-plot + density for OTUs in phyloseq object.
#' @param input.obj, list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param grouping, covariate to group by, default is batch
#' @param center, flag to activate centering
#' @param scale, flag to activate scaling
#' @param method, one of 'ALL' or 'TOP' or a vector of feature names
#' @param return.data, logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @export
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
  tmp.cnts <- scale(tmp.cnts, center = eval(center), scale = eval(scale))
  tmp.cnts <- scale(t(tmp.cnts), center = eval(center), scale = eval(scale))

  if( method[1] == "TOP" ) {
    # calculate IQR and order from largest to smallest
    iqr <- apply(tmp.cnts[otu.idx,],1,IQR)
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


#' Produces mosaic plot to visualize sub-populations due to group-wise combinations, e.g., group x batch
#' @param input.obj, list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param grouping, two covariates of interest to display sample distribution for
#' @param return.data, logical if TRUE returns the data.frame required for plotting (NO plotting or saving here bucko)
#' @export
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
  n.observations <- dim(tmp.meta)[1]

  study.summary <- table(tmp.meta[,eval(model.vars[1])], tmp.meta[,eval(model.vars[2])]) %>%
    as.data.frame(.) %>%
    mutate("Freq.scaled"=Freq / n.observations)

  # prepare plot-annotation
  vars.axes <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",model.vars,perl = TRUE)

  if( return.data ) {
    return(study.summary)
  } # else return the plots

  # split by batch
  plot.v2 <- ggplot(study.summary, aes(x = Var1, y= Freq.scaled, group = Var2, fill=Var1)) +
    facet_grid(cols=vars(Var2), scales="free", space="free_x", drop=T) +
    geom_bar(stat = "identity", width = 0.9) +
    guides(fill = guide_legend(title=eval(vars.axes[1]), reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ylab("Proportion of all observations") +
    theme_mosaic(legend_position = "bottom")

  # split by treatment
  plot.v1 <- ggplot(study.summary, aes(x = Var2, y= Freq.scaled, fill=Var2)) +
    facet_grid(cols=vars(Var1), scales="free", space="free_x", drop=T) +
    geom_bar(stat = "identity", width = 0.9) +
    guides(fill = guide_legend(title=eval(vars.axes[2]), reverse = TRUE, keywidth = 1, keyheight = 1)) +
    ylab("Proportion of all observations") +
    theme_mosaic()

  mosaic.plot <- gridExtra::grid.arrange(plot.v2, plot.v1, ncol=1, nrow=2, heights=c(1,1))

  return(mosaic.plot)

}


# VARIANCE CALCULATION ----------------------------------------------------


#' Wrapper for calculation of variance statistics - Handles erros, validation, iteration and model selection
#' @param input.obj, list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param model.vars, vector of covariates to include in model-construction
#' @param m.formula, EVENTUALLY this will be used to set-up the model and chack covariate availability
#' @param method, WHEN 'm.formula' works, this will become obsolete
#' @param type, creates a column with that string in the output df - to keep track of cnt-source
#' @return df that contains proportions of variance for given covariates in every feature
#' @export
mbecModelVariance <- function( input.obj, model.vars=character(), method=c("lm","lmm","rda","pvca"), type="NONE") {

  ### ToDo: selection cutoff for PCs in silhouette coefficient method?!
  ### ToDo: safety checks and logic to distinguish model types and also take care of this matrix-input issue
  ### ToDoAsWell: How2Make lm and lmm formulas.. or if  or whatever
  ## ToDo: check this out - itis part of lmm just put it her to remind me
  ## control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )

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

    model.variances <- apply(tmp.cnts, 2, FUN = function(x){
      model.fit <- stats::lm(x ~ group + batch, data = tmp.meta)
      p<- mbecVarianceStats(model.fit)
    })
    modelType = "anova"

  } else if( method == "lmm" ) {
    message("Fitting linear-mixed model to every feature and extract proportion of variance explained by covariates.")

    control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )

    model.variances <- apply(tmp.cnts, 2, FUN = function(x) {
      model.fit <- lme4::lmer(x ~ group + (1|batch), data=tmp.meta, control=control)
      p <- mbecVarianceStats(model.fit)
    })

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
      tmp.sig <- anova(tmp.rda.covariate,permutations = how(nperm=999))$`Pr(>F)`[1]

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
    lmm.df <- eVec %>% as_tibble(.) %>% select_at(1:eval(n.PCs)) %>% cbind(., tmp.meta)
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
  res <- data.frame(t(model.variances)) %>% mutate(type = eval(type))
  attr(res, "modelType") <- modelType

  ### fancy return type from variancePart package - wait if we need that
  # res <- new("varPartResults", varPartMat, type=modelType, method="Variance explained (%)")

  return( res )
}


#' Helper function that takes a fitted model and extracts the
#' fraction of variances attributable to every model-variable
#' @param model.fit lm() or lmm() output
mbecVarianceStats <- function( model.fit ) {

  ### ToDo: implement lmm and glm versions of this


  # check validity of model fit
  mbecValidateModel( model.fit)


  if( class(model.fit) %in% "lm" ) {              # linear model
    # get model coefficients
    model.sum <- summary(model.fit)
    w.cof2[i] <- sum.res$coefficients[3,1]

    an = anova(model.fit)
    varFrac = an[['Sum Sq']] / sum( an[['Sum Sq']] )
    names(varFrac) = rownames(an)

  } else if( class(model.fit) %in% "lmerMod" ) {  # linear-mixed model

    # extract variance components
    vc = unlist(getVarianceComponents(model.fit))

    # create fractions
    varFrac = vc / sum(vc)

    # remove ".(Intercept)" string
    names(varFrac) = gsub("\\.\\(Intercept\\)", "", names(varFrac))

    ### ToDo: Variance calculation


  } else if( class(model.fit) %in% "glm" ) {

    ### ToDon't
  }



  return(varFrac)
}


#' Helper function that ensure validity of the model by testing colinearity
#' @param model.fit lm() or lmm() output
#' @param colinearityThreshold cut-off for model rejection
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
  }
}


#' Takes a fitted model and computes maximum correlation between covariates as return value.
#' Return value contains actual correlation-matrix as 'vcor' attribute.
#' Works for lm and lmm.
#' @param model.fit lm() or lmm() output
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



