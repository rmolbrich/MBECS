#' Mock-up microbiome abundance table and meta-data.
#'
#' A dataset containing pre-processed abundance table of
#' microbial communities and a matrix of covariate information.
#'
#' @format An mbecData object including tss and clr transformed counts:
#' \describe{
#'   \item{otu}{Compositional Abundance Data}
#'   \item{tss}{Compositional Abundance Data Sum-Scaled}
#'   \item{clr}{Compositional Abundance Data Log-Ratio Transformed}
#'   \item{meta}{Covariate Information}
#'   ...
#' }
"dummy.mbec"


#' Mock-up microbiome abundance table and meta-data.
#'
#' A dataset containing pre-processed abundance table of
#' microbial communities and a matrix of covariate information.
#'
#' @format A phyloseq object containing counts and meta-data:
#' \describe{
#'   \item{otu_table}{Compositional Abundance Data}
#'   \item{sam_data}{Covariate Information}
#'   ...
#' }
"dummy.ps"


#' Mock-up microbiome abundance table and meta-data.
#'
#' A dataset containing pre-processed abundance table of
#' microbial communities and a matrix of covariate information.
#'
#' @format A list object containing counts and meta-data:
#' \describe{
#'   \item{cnts}{Compositional Abundance Data}
#'   \item{meta}{Covariate Information}
#'   ...
#' }
"dummy.list"


