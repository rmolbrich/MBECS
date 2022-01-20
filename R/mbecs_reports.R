# MBECS REPORT FUNCTIONS --------------------------------------------------


#' Constructs an initial report of a single data-set without comparative analyses.
#'
#' Input can be of class MbecData, phyloseq or list(counts, meta-data). The
#' function will check if required covariates are present and apply
#' normalization with default parameters according to chosen type, i.e., 'clr'
#' (cumulative log-ratio) or 'tss' (total sum scaled).
#'
#' @param input.obj list of phyloseq objects to compare, first element is considered uncorrected data
#' @param model.vars required covariates to build models
#' @param type One of 'otu', 'tss' or 'clr' to determine the abundance matrix to use for evaluation.
#' @param return.data TRUE will return a list of all produced plots, FALSE will start rendering the report
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
mbecReportPrelim <- function(input.obj, model.vars=c("batch","group"), type="clr", return.data = FALSE) {
  # just start with input processing and transform into MbecData if required
  input.obj <- mbecProcessInput(input.obj, required.col=eval(model.vars))

  # for choice of abundance matrix other than 'otu' this checks if it is present
  # and calculates it if not
  type <- match.arg(type, choices = c("otu","clr","tss"))
  if( !(type == "otu") )
    input.obj <- mbecTransform(input.obj, method=eval(type))

  # Prepare the SIX exploratory plots
  prelim.report.list <- list()
  prelim.report.list[["mosaic"]] <- mbecMosaic(input.obj, model.vars=eval(model.vars), return.data = eval(return.data))
  prelim.report.list[["pca"]] <- mbecPCA(input.obj, model.vars=eval(model.vars), pca.axes = c(1,2), type=eval(type), return.data = eval(return.data))
  prelim.report.list[["rle"]] <- mbecRLE(input.obj, model.vars=eval(model.vars), type=eval(type), return.data = eval(return.data))
  prelim.report.list[["heat"]] <- mbecHeat(input.obj, method="TOP", n=10, model.vars=eval(model.vars), type=eval(type), return.data = eval(return.data))
  prelim.report.list[["box"]] <- mbecBox(input.obj, method="TOP", n=5, model.var=eval(model.vars)[1], type=eval(type), return.data = eval(return.data))
  # Dendrogramm needs better function if possible
  #prelim.report.list[["dendro"]] <- mbecDendro(input.obj, method="TOP", n=5, model.var=eval(model.vars)[2] )

  # calculate the variance statistics
  prelim.report.list[["linmod"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="lm",
                                                      type=eval(type))

  prelim.report.list[["linmixmod"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="lmm",
                                                          type=eval(type))

  prelim.report.list[["rda"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="rda",
                                                   type=eval(type))

  prelim.report.list[["pvca"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="pvca",
                                                    type=eval(type))

  prelim.report.list[["scoef"]] <- mbecModelVariance(input.obj, model.vars=model.vars, method="s.coef",
                                                     type=eval(type))

  # if only data is required - stop here and return the list of dfs
  if( return.data ) {
    return(prelim.report.list)
  }

  # to plot or not to plot .. and how
  # for now just call all stat-plot functions and replace the respective values in the list
  prelim.report.list[["linmod"]] <- mbecVarianceStatsPlot(prelim.report.list[["linmod"]])
  prelim.report.list[["linmixmod"]] <- mbecVarianceStatsPlot(prelim.report.list[["linmixmod"]])
  prelim.report.list[["rda"]] <- mbecRDAStatsPlot(prelim.report.list[["rda"]])
  prelim.report.list[["pvca"]] <- mbecPVCAStatsPlot(prelim.report.list[["pvca"]])
  prelim.report.list[["scoef"]] <- mbecSCOEFStatsPlot(prelim.report.list[["scoef"]])

  # construct the report
  ## ToDo: for next version include some configuration that allows users to add
  ## some file naming
  rmarkdown::render(system.file("rmd","mbecReport_prelim.Rmd",package="MBECS"),
                    output_dir = getwd(),
                    output_file = NULL,
                    params = list(report.data=mbecGetData(input.obj, orientation="fxs", required.col=eval(model.vars), type=eval(type)),
                                  report.vars=model.vars,
                                  report.list=prelim.report.list))
}


#' Constructs an initial report of a single data-set without comparative analyses.
#' @param input.obj, list of phyloseq objects to compare, first element is considered uncorrected data
#' @param model.vars, required covariates to build models
#' @param type One of 'otu', 'tss' or 'clr' to determine the abundance matrix to use for evaluation.
#' @param return.data, TRUE will return a list of all produced plots, FALSE will start rendering the report
#' @return either a ggplot2 object or a formatted data-frame to plot from
#' @export
mbecReportPost <- function(input.obj, model.vars=c("batch","group"), type="clr", return.data = FALSE) {

  n.cor <- length(input.obj@corrections)
  n.ass <- length(input.obj@assessments)

  if( n.cor == 0 && n.ass == 0) {
    stop("No corrections available.")
  }

  otu.idx <- mbecBox(input.obj, method = "TOP", n = 4, model.var = model.vars[1], type=eval(type), return.data = TRUE)[[2]]

  pre.list <- list()
  pre.list$mosaic <- mbecMosaic(input.obj, model.vars)
  pre.list$pca$pre <- mbecPCA(input.obj, model.vars, type=eval(type))
  pre.list$rle$pre <- mbecRLE(input.obj, model.vars=eval(model.vars), type=eval(type))
  pre.list$heat$pre <- mbecHeat(input.obj, model.vars, method = otu.idx, type = eval(type))[[4]] # to get the grob
  pre.list$box$pre <- mbecBox(input.obj, method = otu.idx, model.var = model.vars[1], type=eval(type))

  pre.list$linmod <- mbecModelVariance(input.obj, model.vars=model.vars, method="lm", type=eval(type))
  pre.list$linmixmod <- mbecModelVariance(input.obj, model.vars=model.vars, method="lmm", type=eval(type))
  pre.list$rda <- mbecModelVariance(input.obj, model.vars=model.vars, method="rda", type=eval(type))
  pre.list$pvca <- mbecModelVariance(input.obj, model.vars=model.vars, method="pvca", type=eval(type))
  pre.list$scoef <- mbecModelVariance(input.obj, model.vars=model.vars, method="s.coef", type=eval(type))

  # now iterate through all available assessment matrices and calculate metrics
  if( n.ass != 0 ) {
    for( ass.idx in names(input.obj@assessments) ) {
      # FixMe: this needs to do sth. at some point
    }
  }

  # "mosaic"    "pca"       "rle"       "heat"      "box"       "linmod"    "linmixmod" "rda"       "pvca"      "scoef"
  p.list <- pre.list
  # now iterate through all available corrected matrices and calculate metrics
  if( n.cor != 0 ) {
    for( cor.idx in names(input.obj@corrections) ) {

      p.list$pca[[eval(cor.idx)]] <- mbecPCA(input.obj, model.vars=eval(model.vars), pca.axes = c(1,2), type="cor", label = eval(cor.idx))
      p.list$box[[eval(cor.idx)]] <- mbecBox(input.obj, method=otu.idx, model.var=eval(model.vars)[1], type="cor", label = eval(cor.idx))
      p.list$heat[[eval(cor.idx)]] <- mbecHeat(input.obj, method=otu.idx, model.vars=eval(model.vars), type="cor", label = eval(cor.idx))[[4]]
      p.list$rle[[eval(cor.idx)]] <- mbecRLE(input.obj, model.vars=eval(model.vars), type="cor", label = eval(cor.idx))

      p.list$linmod <- rbind.data.frame(p.list$linmod, mbecModelVariance(input.obj, model.vars=model.vars, method="lm", type="cor", label = eval(cor.idx)))
      p.list$linmixmod <- rbind.data.frame(p.list$linmixmod, mbecModelVariance(input.obj, model.vars=model.vars, method="lmm", type="cor", label = eval(cor.idx)))
      p.list$rda <- rbind.data.frame(p.list$rda, mbecModelVariance(input.obj, model.vars=model.vars, method="rda", type="cor", label = eval(cor.idx)))
      p.list$pvca <- rbind.data.frame(p.list$pvca, mbecModelVariance(input.obj, model.vars=model.vars, method="pvca", type="cor", label = eval(cor.idx)))
      p.list$scoef <- rbind.data.frame(p.list$scoef, mbecModelVariance(input.obj, model.vars=model.vars, method="s.coef", type="cor", label = eval(cor.idx)))
    }
  }
  # FixMe: either do separate branches or just drop it
  if( return.data ) {
    return(p.list)
  }

  ## produce plot-panels for variance assessment analyses
  p.list[["linmod"]] <- mbecVarianceStatsPlot(p.list[["linmod"]])
  p.list[["linmixmod"]] <- mbecVarianceStatsPlot(p.list[["linmixmod"]])
  p.list[["rda"]] <- mbecRDAStatsPlot(p.list[["rda"]])
  p.list[["pvca"]] <- mbecPVCAStatsPlot(p.list[["pvca"]])
  p.list[["scoef"]] <- mbecSCOEFStatsPlot(p.list[["scoef"]])


  rmarkdown::render(system.file("rmd","mbecReport_post.Rmd",package="MBECS"),
                    output_dir = getwd(),
                    output_file = NULL,
                    params = list(report.data=mbecGetData(input.obj, orientation="fxs", required.col=eval(model.vars)),
                                  report.vars=model.vars,
                                  report.list=p.list))
}




