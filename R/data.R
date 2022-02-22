# DATA GENERATION ---------------------------------------------------------

#' Creates a dummy data-set with abundance matrix and meta-data.
#'
#' For given number of otus and samples this will create mockup microbiome data
#' that contains systematic and non-systematic batch effects. Comes with meta
#' data that denotes study groups and batches. The replicate column is fake
#' and only used to test RUV-implementations.
#'
#' 'Group' and 'batch' variables are actually taken into account in data
#' creation, but only to the degree that the random draws for values are
#' performed with different parameters respectively.
#'
#' THIS HAS ONLY A CONCEPTUAL SIMILARITY TO MICROBIOME DATA AT BEST AND IS IN
#' NO WAY USEFUL OTHER THAN TESTING PACKAGE FUNCTIONS AND VISUALIZING WORKFLOWS!
#'
#' @param n.otus Integer to determine number of features to "simulate".
#' @param n.samples Even integer to set number of samples to "simulate".
#' @return A list object that contains the made up abundance and the
#' accompanying meta-data.
#'
#' @export
#'
#' @examples
#' dummy.list <- mbecDummy(n.otus=100, n.samples=30)
mbecDummy <- function(n.otus=500, n.samples=40) {
    # make it binary
    replace1 <- function(x) {
        dplyr::if_else(x != 0, 1, 0)
    }

    n <- 2
    # create meta-data
    meta <- data.frame(
        "sID"=paste("S", 1:n.samples, sep=""),
        "group"=factor(c(rep("A", times=n.samples/n),
                         rep("B", times=n.samples/n))),
        "batch"=factor(rep(c("B1","B2"), times=n.samples/n)),
        "replicate"=NA,
        row.names = "sID")

    for( rep.idx in seq(n.samples) ) {
        #
        meta$replicate[rep.idx] <- (( rep.idx + (n - 1)) %% (n-(n-2)) ) + (n-1)

        if( rep.idx %% 4 == 0 ) n <- n +2
    }

    base_mtx <- Matrix::rsparsematrix(
        n.samples, n.otus, 0.314, symmetric = FALSE) %>%
        as.matrix() %>% as.data.frame() %>%
        dplyr::rename_with(~ gsub("V","OTU", .x)) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), replace1))

    rownames(base_mtx) <- rownames(meta)

    # prep mtx for systematic and non-systematic BE
    cnts.sys <- base_mtx
    cnts.nonsys <- base_mtx

    ## DEFINE distributions for NON-SYS, OTU, TREATMENT, BE
    # random impact of non-systematic BE
    nsys <- sample(0:1,size=n.otus,replace=TRUE)
    # some arbitrary distribution for abundance of all OTUs
    otu.dist <- round(stats::rexp(n.otus, .01))
    # some arbitrary distribution for treatment effect
    treat.dist <- stats::rlnorm(n = n.otus, meanlog=3, sdlog=0.5)
    # some arbitrary distribution for batch-effect
    batch_1 <- c(1:10,21:30)
    batch.dist <- stats::rnorm(n.otus, mean = 5, sd = 1)

    for( c.idx in 1:n.otus ) {
        ## SYSTEMATIC
        cnts.sys[cnts.sys[,c.idx] == 1, c.idx] <-
            abs(stats::rnorm(n=length(cnts.sys[cnts.sys[,c.idx] == 1, c.idx]),
                             mean=otu.dist[c.idx],  sd=1))

        # now add 'treatment-effect' to first twenty rows/samples
        # use normal distribution for all otus for slight variation
        cnts.sys[which(cnts.sys[1:20,c.idx] != 0), c.idx] <-
            cnts.sys[which(cnts.sys[1:20,c.idx] != 0), c.idx] +
            stats::rnorm(n=length(which(cnts.sys[1:20,c.idx] != 0)),
                         mean=treat.dist[c.idx], sd=1)

        # B1: now add systematic batch-effect to one batch
        # use normal distribution for all otus for slight variation
        cnts.sys[which(cnts.sys[batch_1,c.idx] != 0), c.idx] <-
            cnts.sys[which(cnts.sys[batch_1,c.idx] != 0), c.idx] +
            abs(stats::rnorm(n=length(which(cnts.sys[batch_1,c.idx] != 0)),
                             mean=stats::rnorm(1, mean = 5, sd = 1),
                             sd=abs(stats::rnorm(1, mean = 0, sd = 2))))

        ## NON-SYSTEMATIC
        cnts.nonsys[cnts.nonsys[,c.idx] == 1, c.idx] <-
            stats::rnorm(n=length(cnts.nonsys[cnts.nonsys[,c.idx] == 1, c.idx]),
                         mean=otu.dist[c.idx], sd=1)

        # now add 'treatment-effect' to first twenty rows/samples
        # use normal distribution for all otus for slight variation
        cnts.nonsys[which(cnts.nonsys[1:20,c.idx] != 0), c.idx] <-
            cnts.nonsys[which(cnts.nonsys[1:20,c.idx] != 0), c.idx] +
            stats::rnorm(n=length(which(cnts.nonsys[1:20,c.idx] != 0)),
                         mean=treat.dist[c.idx], sd=1)

        # B1: now add non-systematic batch-effect to one batch
        # use normal distribution for all otus for slight variation
        if( nsys[c.idx] == 0 ) {
            # just add some zero mean noise here
            cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] <-
                cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] +
                abs(stats::rnorm(n=length(which(cnts.nonsys[batch_1,c.idx] !=
                                                    0)),
                                 mean=stats::rnorm(1, mean = 0, sd = 1),
                                 sd=abs(stats::rnorm(1, mean = 0, sd = 2))))

        } else {
            cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] <-
                cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] +
                abs(stats::rnorm(n=length(which(cnts.nonsys[batch_1,c.idx] !=
                                                    0)),
                                 mean=stats::rnorm(1, mean = 50, sd = 1),
                                 sd=abs(stats::rnorm(1, mean = 0, sd = 2))))
        }
    }

    otu.cnts <- abs(floor(cnts.nonsys))
    dummy.list <- list("cnts"=otu.cnts, "meta"=meta)
    #usethis::use_data(dummy.list, overwrite = T)
    #dummy.mbec <- mbecTransform(dummy.list, method="tss")
    #dummy.mbec <- mbecTransform(dummy.mbec, method="clr")
    #usethis::use_data(dummy.mbec, overwrite = T)

    return(dummy.list)
}


# DATA-SETS ---------------------------------------------------------------

#' Mock-up microbiome abundance table and meta-data.
#'
#' An artificial data-set containing pre-processed abundance table of
#' microbial communities and a matrix of covariate information. The data was
#' created using the \code{\link{mbecDummy}} function for the sole purpose of
#' running examples and showing the package workflow. This particular object
#' was also processed with \code{\link{mbecTransform}} function in order to
#' generate "clr" and "tss" transformed abundance matrices.
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
#' An artificial data-set containing pre-processed abundance table of
#' microbial communities and a matrix of covariate information. The data was
#' created using the \code{\link{mbecDummy}} function for the sole purpose of
#' running examples and showing the package workflow. This particular object
#' was then converted using \code{\link{phyloseq}}.
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
#' An artificial data-set containing pre-processed abundance table of
#' microbial communities and a matrix of covariate information. The data was
#' created using the \code{\link{mbecDummy}} function for the sole purpose of
#' running examples and showing the package workflow.
#'
#' @format A list object containing counts and meta-data:
#' \describe{
#'   \item{cnts}{Compositional Abundance Data}
#'   \item{meta}{Covariate Information}
#'   ...
#' }
"dummy.list"


