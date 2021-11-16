# MBECS REPORT FUNCTIONS --------------------------------------------------


#' MBECS Report Pipeline
#'
#' A wrapper for the generation of preliminary and comparative reports. The function will call the
#' appropriate sub-routines based on the given input. The subsequent sections will assume that the
#' covariates 'treatment' and 'batch' are the effects of interest (CoI). Other variables can be
#' chosen via the functions 'model.vars' parameter. Accepted inputs are MbecData-objects,
#' phyloseq-objects or a list that contains counts and covariate data. Multiple MbecData and
#' phyloseq objects in a list will produce a comparative report and single objects generate a
#' preliminary report.
#'
#' The reports are partitioned into three sections for both types of reports.
#'
#' Study Summary: This section provides a general overview of study design with a summary of
#' covariate variables, distribution of samples with respect to grouping, e.g., case/control or
#' un-treated/different treatment levels, and known batches as well as an ordination-plot for
#' the first two principal components.
#'
#' Visualization: This section aims to provide a visual assessment of the batch-effect. It provides
#' a plot of relative log-expression, a heatmap of the top 5 most variable features, a dendrogram
#' (hopefully as soon as it works) and box-plots of the expression of the top 5 most variable
#' features with respect to the batches.
#'
#' Variance Assessment: Diverse methods are applied to calculate the proportion of variance that
#' is explainable by the CoI. Linear (Mixed) Models with user supplied covariates are fit to the
#' data-set, partial Redundancy Analysis (CCA or pRDA) and Principal Variance Component Analysis
#' (PVCA) use different approaches to compute proportions of variance as well. Finally, the
#' Silhouette Coefficient is a representation for the goodness of fit. Which is how good grouping
#' covariates such as treatment or batch are represented in the clustering of samples.
#'
#' @keywords Comparative Preliminary Report
#' @param input.obj list(cnts, meta), phyloseq, MbecData object (correct orientation is handeled internally)
#' @param model.vars two covariates of interest to select by first variable selects panels and second one determines coloring
#' @param return.data logical if TRUE returns the data.frame required for plotting (NO plotting here bucko)
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
#'
#' @examples
#' # This will return the data.frame for plotting.
#' \dontrun{p.RLE <- mbecRLE(input.obj=list(counts, covariates),
#' model.vars=c("treatment","batches"), return.data=TRUE)}
#'
#' # This will return the ggplot2 object for display, saving and modification.
#' \dontrun{p.RLE <- mbecRLE(input.obj=phyloseq, model.vars=c("treatment","sex"),
#' return.data=FALSE)}
mbecReport <- function(input.obj, model.vars=c("group","batch"), return.data=FALSE) {

  ## first determine whether this is a post- or a preliminary-report and then call the appropriate function

  # if list & !any are ps or mbecData objects --> prelim with cnts,meta input
  if( is(input.obj, "list") ) {
    if( length(input.obj) == 2 & !any(lapply(input.obj, is) %in% c("MbecData", "phyloseq")) ) {
      # prelim-report from cnts + meta
      message("We have a preliminary report from cnts+meta-input!")
      tmp.report <- mbecReportPrelim(input.obj, model.vars, return.data)

    } else {
      message("We have a comparative post-report!")
      # post-report from multiple phyloseq or MbecData objects
      # tmp.report <- mbecReportPost(input.obj, model.vars, return.data)
    }
  } else {
    if( class(input.obj) %in% "phyloseq" ) {
      # no list and class phyloseq means preliminary report
      message("We have a preliminary report from phyloseq-input!")

      tmp.report <- mbecReportPrelim(input.obj, model.vars, return.data)

    } else if( class(input.obj) %in% "MbecData" ) {
      # if only 'raw' MbecData --> call Prelim else Post
      if( length(attr(input.obj, "transformations")) >= 1 ) {
        message("We have a comparative post-report!")

        # tmp.report <- mbecReportPost(input.obj, model.vars, return.data)

      } else {
        # if 'transformations' list is empty - only preliminary report
        # (technically, this can also be Post-correction.. I guess ToDo)
        tmp.report <- mbecReportPrelim(input.obj, model.vars, return.data)
      }
    }
  }

  return(tmp.report)
}


#' Constructs an initial report of a single data-set without comparative analyses.
#' @param input.obj, list of phyloseq objects to compare, first element is considered uncorrected data
#' @param model.vars, required covariates to build models
#' @param return.data, TRUE will return a list of all produced plots, FALSE will start rendering the report
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
mbecReportPrelim <- function(input.obj, model.vars=c("group","batch"), return.data = FALSE) {
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

  if( return.data ) {
    return(prelim.report.list)
  }

  # construct the report

  input.obj <- mbecGetData(input.obj, orientation="fxs", required.col=eval(model.vars))


  rmarkdown::render("mbecReport_prelim_TMPLT.Rmd",
                    params = list(report.data=input.obj,
                                  report.vars=model.vars,
                                  report.list=prelim.report.list))
}


