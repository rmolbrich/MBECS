#' MBECS: Evaluation and correction of batch effects in microbiome data-sets.
#'
#' The Microbiome Batch-Effect Correction Suite aims to provide a toolkit for
#' stringent assessment and correction of batch-effects in microbiome data sets.
#' To that end, the package offers wrapper-functions to summarize study-design
#' and data, e.g., PCA, Heatmap and Mosaic-plots, and to estimate the proportion
#' of variance that can be attributed to the batch effect. The function
#' \code{\link{mbecCorrection}} acts as a wrapper for various batch effects
#' correction algorithms (BECA) and in conjunction with the aforementioned
#' tools, it can be used to compare the effectiveness of correction methods on
#' particular sets of data.
#' All functions of this package are accessible on their own or within the
#' preliminary and comparative report pipelines respectively.
#'
#' @section Pipeline:
#'
#' \itemize{
#' \item{\code{\link{mbecProcessInput}}}
#' \item{\code{\link{mbecTransform}}}
#' \item{\code{\link{mbecReportPrelim}}}
#' \item{\code{\link{mbecCorrection}}}
#' \item{\code{\link{mbecRunCorrections}}}
#' \item{\code{\link{mbecReportPost}}}
#' }
#'
#' @section Exploratory functions:
#'
#' \itemize{
#' \item{\code{\link{mbecRLE}}}
#' \item{\code{\link{mbecPCA}}}
#' \item{\code{\link{mbecBox}}}
#' \item{\code{\link{mbecHeat}}}
#' \item{\code{\link{mbecMosaic}}}
#' }
#'
#' @section Variance functions:
#'
#' \itemize{
#' \item{\code{\link{mbecModelVariance}}}
#' \item{\code{\link{mbecVarianceStatsPlot}}}
#' \item{\code{\link{mbecRDAStatsPlot}}}
#' \item{\code{\link{mbecPVCAStatsPlot}}}
#' \item{\code{\link{mbecSCOEFStatsPlot}}}
#' }
#'
#' @importFrom magrittr "%>%"
#'
#' @docType package
#' @name MBECS
