# MBECS REPORT FUNCTIONS --------------------------------------------------


#' Wrapper function for the analyses required to construct preliminary or post-correction reports.
#' @param input.obj, list of phyloseq objects to compare, first element is considered uncorrected data
#' @param model.vars, required covariates to build models
#' @param return.plot, TRUE will return a list of all produced plots, FALSE will start rendering the report
#' @export
mbecReport <- function(input.obj, model.vars=c("group","batch"), return.plot=FALSE) {

  ## first determine whether this is a post- or a preliminary-report and then call the appropriate function

  # if list & !any are ps or mbecData objects --> prelim with cnts,meta input
  if( class(input.obj) == "list" ) {
    if( length(input.obj) == 2 & !any(sapply(input.obj, class) %in% c("MbecData", "phyloseq")) ) {
      # prelim-report from cnts + meta
      message("We have a preliminary report from cnts+meta-input!")
      tmp.report <- mbecReportPrelim(input.obj, model.vars, return.plot)

    } else {
      message("We have a comparative post-report!")
      # post-report from multiple phyloseq or MbecData objects
      # tmp.report <- mbecReportPost(input.obj, model.vars, return.plot)
    }
  } else {
    if( class(input.obj) %in% "phyloseq" ) {
      # no list and class phyloseq means preliminary report
      message("We have a preliminary report from phyloseq-input!")

      tmp.report <- mbecReportPrelim(input.obj, model.vars, return.plot)

    } else if( class(input.obj) %in% "MbecData" ) {
      # if only 'raw' MbecData --> call Prelim else Post
      if( length(attr(input.obj, "transformations")) >= 1 ) {
        message("We have a comparative post-report!")

        # tmp.report <- mbecReportPost(input.obj, model.vars, return.plot)

      } else {
        # if 'transformations' list is empty - only preliminary report
        # (technically, this can also be Post-correction.. I guess ToDo)
        tmp.report <- mbecReportPrelim(input.obj, model.vars, return.plot)
      }
    }
  }
}


#' Constructs an initial report of a single data-set without comparative analyses.
#' @param input.obj, list of phyloseq objects to compare, first element is considered uncorrected data
#' @param model.vars, required covariates to build models
#' @param return.plot, TRUE will return a list of all produced plots, FALSE will start rendering the report
#' @export
mbecReportPrelim <- function(input.obj, model.vars=c("group","batch"), return.plot = FALSE) {
  # only three situations here: list with cnts and meta, phyloseq or MbecData

  # Prepare the SIX exploratory plots
  prelim.report.list <- list()
  prelim.report.list[["mosaic"]] <- mbecMosaic( input.obj, model.vars=eval(model.vars) )
  prelim.report.list[["pca"]] <- mbecPCA(input.obj, model.vars=eval(model.vars), pca.axes = c(1,2))
  prelim.report.list[["rle"]] <- mbecRLE(input.obj, model.vars=eval(model.vars))
  prelim.report.list[["heat"]] <- mbecHeat(input.obj, method="TOP", n=5, model.vars=eval(model.vars))
  prelim.report.list[["box"]] <- mbecBox(input.obj, method="TOP", n=5, model.var=eval(model.vars)[2] )
  # Dendrogramm needs better function if possible
  #prelim.report.list[["dendro"]] <- mbecDendro(input.obj, method="TOP", n=5, model.var=eval(model.vars)[2] )

  # calculate the variance statistics
  prelim.report.list[["linmod"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="lm",
                                                      type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  prelim.report.list[["linmixmod"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="lmm",
                                                         type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  prelim.report.list[["rda"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="rda",
                                                   type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  prelim.report.list[["pvca"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="pvca",
                                                    type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  prelim.report.list[["scoef"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="s.coef",
                                                     type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  # to plot or not to plot .. and how
  # for now just call all stat-plot functions and replace the respective values in the list
  prelim.report.list[["linmod"]] <- mbecVarianceStatsPlot(prelim.report.list[["linmod"]])
  prelim.report.list[["linmixmod"]] <- mbecVarianceStatsPlot(prelim.report.list[["linmixmod"]])
  prelim.report.list[["rda"]] <- mbecRDAStatsPlot(prelim.report.list[["rda"]])
  prelim.report.list[["pvca"]] <- mbecPVCAStatsPlot(prelim.report.list[["pvca"]])
  prelim.report.list[["scoef"]] <- mbecSCOEFStatsPlot(prelim.report.list[["scoef"]])

  if( return.plot ) {
    return(prelim.report.list)
  }

  # construct the report

  input.obj <- mbecGetData(input.obj, orientation="fxs", required.col=eval(model.vars))


  rmarkdown::render("mbecReport_prelim_TMPLT.Rmd",
                    params = list(report.data=input.obj,
                                  report.vars=model.vars,
                                  report.list=prelim.report.list))
}



#' Covariate-Variances as modelled by lm/lmm will be displayed as box-plots. Both preliminary and post statistics.
#' @param variance.obj, list or single output of 'mbecVarianceStats' with method lm
#' @export
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


#' Plots variance statistics produced by 'mbecModelVariance' with method 'rda'
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


#' Plots variance statistics produced by 'mbecModelVariance' with method 'pvca'
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


#' Plots variance statistics produced by 'mbecModelVariance' with method silhouette-coefficient
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

