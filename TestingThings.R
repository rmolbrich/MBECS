# TESTING -----------------------------------------------------------------


# for now take the ad.data dummy file and use it as 'code-example-execution-data'
data.dummy <- readRDS(file=file.path("dummy.matrix.rds"), refhook = NULL)

save(data.dummy, file = "data.dummy.Rda")
load(file.path("data","data.dummy.Rdata"))

#' \dontrun{p.RLE <- mbecRLE(input.obj=phyloseq, model.vars=c("treatment","sex"),
#' return.data=FALSE)}

p.RLE <- mbecRLE(input.obj=datadummy, model.vars=c("group","batch"), return.data=FALSE)
datadummy = data.dummy


p.Box <- mbecBox(input.obj=datadummy, method="TOP", n=10, model.var="batch", return.data=FALSE)

p.Heat <- mbecHeat(input.obj=datadummy, model.vars=c("group","batch"), center=TRUE,
                     scale=TRUE, method="TOP", n=10, return.data=T)

p.Mosaic <- mbecMosaic(input.obj=datadummy, model.vars=c("group","batch"), return.data=FALSE)




#' # This will return a score for the supplied vector with default evaluation
#' # (strict).
#' val.score <- poscore(cnt.vec=runif(100, min=0, max=100), cnt=42,
#' type="strict")

MbecData(input.obj=datadummy$cnts, meta.obj = datadummy$meta)
val.score <- poscore(cnt.vec=runif(100, min=0, max=100), cnt=42, type="strict")

https://support.bioconductor.org/api/email/M.Olbrich%40protonmail.com/

#' # This will return a paneled plot that shows results for the variance assessment.
#' df.var.pvca <- mbecModelVariance(input.obj=datadummy, model.vars=c("group","batch"),
#' method="s.coef", type="RAW")
#' plot.scoef <- mbecSCOEFStatsPlot(scoef.obj=df.var.scoef)



  ## Backup
  Authors@R:
  person(given = "Michael",
         family = "Olbrich",
         role = c("aut", "cre"),
         email = "M.Olbrich@protonmail.com",
         comment = c(ORCID = "0000-0003-2789-3382"))





# MAKE EXAMPLES GREAT AGAIN -----------------------------------------------

#' # This will check for the presence of variables 'group' and 'batch' in the meta-data and return
#' # an object of class 'MbecData'.
#' MbecData.obj <- mbecProcessInput(input.obj=datadummy,
#'     required.col=c("group","batch"))


MbecData.obj <- mbecProcessInput(input.obj=datadummy, required.col=c("group","batch"))

required.col <- c("group", "batch")

model.vars <- required.col
g.idx <- 1

stop("You need to supply a meta-frame that contains the columns: ", paste(required.col, collapse=", "), call. = FALSE)

ref.group <- "letest"

message("Group ",ref.group, " is considered control group, i.e., reference for normalization procedure. To change reference please 'relevel()' grouping factor accordingly.")

warning("Grouping variables need to be factors. Coercing variable: ", eval(model.vars[g.idx]), " to factor now, adjust beforehand to get best results.")

otu.idx <- c(1,2,3,4,5,6)
n <- 3


1:min(length(otu.idx), n)

seq_along(model.vars)

seq_len(min(length(otu.idx), n))

model.vars

1:length(model.vars)

1:(n.vars-1)

seq_len((n.vars-1))

var.idx <- 1

(var.idx+1):(n.vars)

seq.int(from=(var.idx+1), to=(n.vars), by=1)


corrected.cnts <- datadummy$cnts

seq_len(dim(corrected.cnts)[2])



# VAPPLY instead of sapply ------------------------------------------------

# Bioc-requirement:
# "Avoid class membership checks with class() / is() and == / !=; Use is(x, 'class') for S4 classes"
## doing it like this is apparently a NoNo
class(ps) == "phyloseq"
is(ps) == "phyloseq"

# supposedly we are doing this now.
is(ps, "phyloseq")
is(ps, "MbecData")

person(given = "Michael",
       family = "Olbrich",
       role = c("aut", "cre"),
       email = M.Olbrich@protonmail.com,
       comment = c(ORCID = "0000-0003-2789-3382"))

# if( length(input.obj) == 2 & !any(sapply(input.obj, class) %in% c("MbecData", "phyloseq")) )


# this works just fine
some.vec <- sapply(g, function(x) x$name)

legend <- g[[which(vapply(g, function(x) x$name, FUN.VALUE = character(1)) == "guide-box")]]

# hmm
n.PCs <- max(3, min(which(sapply(cumsum(prop.PCs), function(x) x >= pct_threshold))))

n.PCs <- max(3, min(which(vapply(cumsum(prop.PCs), function(x) x >= pct_threshold, FUN.VALUE=logical(1)))))
vapply(cumsum(prop.PCs), function(x) x >= pct_threshold, FUN.VALUE=logical(1))

# model variance fun
lib.df <- data.frame("covariates"=colnames(model.fit@frame),
                     row.names=vapply(colnames(model.fit@frame), function(covariate)
                       paste(paste(covariate, levels(model.fit@frame[,eval(covariate)]), sep=""), collapse = ","), FUN.VALUE = character(1)))

vapply(colnames(model.fit@frame), function(covariate)
  paste(paste(covariate, levels(model.fit@frame[,eval(covariate)]), sep=""), collapse = ","), FUN.VALUE = character(1))

covariate <- "batch"

paste(paste(covariate, levels(model.fit@frame[,eval(covariate)]), sep=""), collapse = ",")

# reports stuff - also fix class stuff because apparently class(bla)=="sth" is not good
!any(sapply(input.obj, class) %in% c("MbecData", "phyloseq"))

# how class should work
is(input.obj, "list")
bla <- sapply(input.obj, class)
is(input.obj$cnts, c("muh","ver"))

is(input.obj[[1]], c("sdf","fgh"))
is(input.obj[[1]], c("MbecData","phyloseq"))
is(input.obj[[1]], "MbecData")

any(is(input.obj[[1]]) %in% c("MbecData","phyloseq"))

# doesn't need apply at all because it's clear that we only have 2 elements
length(input.obj) == 2 & !any(lapply(input.obj, is) %in% c("MbecData", "phyloseq"))

length(input.obj) == 2 & !is(input.obj[[1]], c("MbecData", "phyloseq")) & !is(input.obj[[2]], c("MbecData", "phyloseq"))



bla <- sapply(input.obj, is)
vapply(input.obj, is, FUN.VALUE=c(character()))
blubb <- lapply(input.obj, is)




# TEST ModelVariance ------------------------------------------------------

model.vars<- c("batch","group")

df.lm <- mbecModelVariance(input.obj=datadummy, model.vars=model.vars, method="lm",
                  type="raw")


df.lmm <- mbecModelVariance(input.obj=datadummy, model.vars=model.vars, method="lmm",
                  type="raw")


df.rda <- mbecModelVariance(input.obj=datadummy, model.vars=model.vars, method="rda",
                            type="raw")

df.pvca <- mbecModelVariance(input.obj=datadummy, model.vars=model.vars, method="pvca")

df.scoef <- mbecModelVariance(input.obj=datadummy, model.vars=model.vars, method="s.coef",
                             type="raw")




# PRELIMINARY REPORT ------------------------------------------------------

le.test <- mbecReportPrelim(input.obj=datadummy, model.vars=c("group","batch"), return.data = FALSE)


# BECA TEST ---------------------------------------------------------------

# c("lm","lmm","sva","ruv2","ruv4")
# c("ruv3","bmc","bat","rbe","fab","pn","svd")

ps.bmc <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "bmc")
ps.bat <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "bat")
ps.rbe <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "rbe")
ps.fab <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "fab")
ps.pn <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "pn")
ps.svd <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "svd")


head(datadummy$cnts)


testdata <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","ps_dummy_56taxa_untransformed.rds"), refhook = NULL)

le.test <- mbecReportPrelim(input.obj=testdata, model.vars=c("group","batch"), return.data = FALSE)

correction.test <- mbecRunCorrections(input.obj=testdata, model.vars=c("group","batch"), method=c("bmc","bat","rbe","pn","svd"))

bla <- mbecReport(input.obj=correction.test, model.vars=c("group","batch"), return.data=FALSE)

test.bmc <- mbecBMC(input.obj=testdata, model.vars=c("group","batch"))
test.bat <- mbecBat(input.obj=testdata, model.vars=c("group","batch"))
test.rbe <- mbecRBE(input.obj=testdata, model.vars=c("group","batch"))
test.fab <- mbecFAB(input.obj=testdata, model.vars=c("group","batch"))
test.pn <- mbecPN(input.obj=testdata, model.vars=c("group","batch"))
test.svd <- mbecSVD(input.obj=testdata, model.vars=c("group","batch"))



# CLR 2.0 -----------------------------------------------------------------

clr(x) <- ln(x) - mean(ln(x))



testdata <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","ps_dummy_56taxa_untransformed.rds"), refhook = NULL)
tmp <- mbecGetData(testdata, orientation="sxf", required.col=c("group","batch"))
x <- tmp[[1]]; tmp.meta <- tmp[[2]]





test.mtx <- matrix(1:9, 3,3)
test.transform <- clr.transfo(test.mtx)

exp(mean(log(c(1,2,3)))) * c(0.3293169, 1.317268, 2.305218)

check.transform <- mbecCLR(test.mtx)

test.transform.x <- clr.transfo(x)
check.transform.x <- mbecCLR(x)




# ALR  --------------------------------------------------------------------

library(Compositional)
library(MASS)
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
y <- alr(x)
x1 <- alrinv(y)

# my own test
test.mtx <- matrix(1:9, 3,3)
alr.test <- alr(test.mtx)
x1 <- alrinv(alr.test)



y <- alr.test
x <- cbind(1, exp(y))
x/rowSums(x)


x <- test.mtx
x[, -1]/x[, 1]


alr()

function(x) {
  Rfast::Log(x[, -1]/x[, 1])
}



# ILR ---------------------------------------------------------------------

mbecILR <- function(input.mtx, offset = 0) {
  if( dim(input.mtx)[2] < 2 ) {
    message("No basis for transformation. Matrix contains less than 2 features,
            returning unchanged.")
    return(input.mtx)
  }
  # 1. stop for negative values and NAs
  if( any(input.mtx < 0 | is.na(input.mtx)) ) {
    stop("Examine your data for NAs and negative values, ILR transformation
         requires complete positive values.\n")
  }
  if( any(input.mtx==0) & offset==0 ) {
    message("Found zeros, function will add a small pseudo-count (1/#features) to enable log-ratio transformation.")
    offset <- 1/ncol(input.mtx)
  }
  # add the offset
  input.mtx <- input.mtx + offset

  return(input.mtx)
}




#
# ILR (Isometric log-ratio) transformation.
# `x` is an `n` by `k` matrix of positive observations with k >= 2.
#
ilr <- function(x, p=0) {
  y <- log(x)
  if (p != 0) y <- (exp(p * y) - 1) / p       # Box-Cox transformation
  y <- y - rowMeans(y, na.rm=TRUE)            # Recentered values
  k <- dim(y)[2]
  H <- contr.helmert(k)                       # Dimensions k by k-1
  H <- t(H) / sqrt((2:k)*(2:k-1))             # Dimensions k-1 by k
  # if(!is.null(colnames(x)))                   # (Helps with interpreting output)
  #   colnames(z) <- paste0(colnames(x)[-1], ".ILR")
  return(y %*% t(H))                          # Rotated/reflected values
}
#
# Specify a Dirichlet(alpha) distribution for testing.
#
alpha <- c(1,2,3,4)
#
# Simulate and plot compositional data.
#
n <- 1000
k <- length(alpha)
x <- matrix(rgamma(n*k, alpha), nrow=n, byrow=TRUE)
x <- x / rowSums(x)
colnames(x) <- paste0("X.", 1:k)
pairs(x, pch=19, col="#00000040", cex=0.6)
#
# Obtain the ILR.
#
y <- ilr(x, 0.5)
colnames(y) <- paste0("Y.", 1:(k-1))
#
# Plot the ILR.
#
pairs(y, pch=19, col="#00000040", cex=0.6)


bla <- ilr.transfo(x)
pairs(bla, pch=19, col="#00000040", cex=0.6)

hmm <- clr.backtransfo(bla)


clr.backtransfo = function(x) {
  # construct orthonormal basis
  V = matrix(0, nrow = ncol(x), ncol = ncol(x)-1)
  for( i in seq_len(ncol(V)) )
  {
    V[seq_len(i), i] = 1/i
    V[i+1, i] = (-1)
    V[, i] = V[, i] * sqrt(i/(i+1))
  }
  rownames(V) = colnames(x)
  return(V)

}

ilr.transfo = function(x, fast = TRUE, offset = 0) {
  # create output matrix
  x.ilr = matrix(NA, nrow = nrow(x), ncol = ncol(x)-1)
  D = ncol(x)
  if (fast) {
    for (i in seq_len(ncol(x.ilr)) ) {
      #x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, (i+1) : D, drop = FALSE]), 1, prod) + offset)^(1 / (D-i))) / (x[,i]+ offset)) ToDo: remove once it works
      x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, seq.int(from=(i+1),to=D,by=1), drop = FALSE]), 1, prod) + offset)^(1 / (D-i))) / (x[,i]+ offset))
      #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop = FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
    }
  } else {
    for (i in seq_len(ncol(x.ilr)) ) {
      x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(apply(as.matrix(x[, seq.int(from=(i+1),to=D,by=1)]), 1, function(x){exp(log(x))})/(x[, i]+ offset) + offset)
      #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(apply(as.matrix(x[,(i+1):D]), 1, function(x){exp(log(x))})/(x[,i]))
    }
  }
  ### ToDo: take care of this class-mess!
  class(x.ilr) = c(class(x.ilr), 'ilr')
  return(as.matrix(x.ilr))
}



# REAL LIFE DATA-SET ------------------------------------------------------

# load package
devtools::load_all()
# load data
rld <- readRDS(file=file.path("/Users/molbrich/Projects/Axel/Dog.Microbiome","data","RDS","phyloseq.OTU.rds"), refhook = NULL)
# remove the 'stool' samples
rld.nostool <- phyloseq::subset_samples(rld, Batch!="stool")
#rld.nostool <- phyloseq::subset_samples(rld, location=="palm")
# agglomerate to 'genus' level
rld.nostool <- tax_glom(rld.nostool, "Genus")
#test <- mbecProcessInput(rld, required.col = c("location","Batch"))
rld.clr <- LRTransform(rld.nostool, method="CLR", offset=0.01, required.col = c("location","Batch"))

## DEBUG ##
input.obj <- rld.clr
model.vars <- c("Batch","location") # for palm use group instead of batch

mbecTestModel(rld.clr, model.vars=c("Batch","location","diet"), mod)


### Preliminary check-up
mosaic <- mbecMosaic(rld.clr, model.vars=c("Batch","location") )
pca <- mbecPCA(rld.clr, model.vars=c("Batch"), pca.axes = c(1,2))
pca.loc <- mbecPCA(rld.clr, model.vars=c("Batch","location"), pca.axes = c(1,2))
# take a look
plot(pca)
rle <- mbecRLE(rld.clr, model.vars=c("Batch","location")) ## this is too large
heat <- mbecHeat(rld.clr, method="TOP", n=10, model.vars=c("Batch","location")) # include flag to switch if sIDs are shown
box <- mbecBox(rld.clr, method="TOP", n=5, model.var=c("Batch") )
# take a look
plot(box$Otu2)
plot(box$Otu8122)
plot(box$Otu6)
plot(box$Otu3)
plot(box$Otu36)

model.vars <- c("group","Batch") # for palm use group instead of batch
# calculate the variance statistics
var.linmod <- mbecModelVariance(rld.clr, model.vars=c("Batch","location"), method="lm",
                                                    type="prelim")

var.linmixmod <- mbecModelVariance(input.obj, model.vars=c("Batch"), method="lmm",
                                                        type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

var.rda <- mbecModelVariance(rld.clr, model.vars=c("Batch","location"), method="rda",
                                                 type="prelim")

var.pvca <- mbecModelVariance(rld.clr, model.vars=c("Batch","location"), method="pvca",
                              type="prelim")

var.scoef <- mbecModelVariance(rld.clr, model.vars=c("Batch","location"), method="s.coef",
                               type="prelim")

plot.linmod <- mbecVarianceStatsPlot(var.linmod)
plot.linmixmod <- mbecVarianceStatsPlot(var.linmixmod)
plot.rda <- mbecRDAStatsPlot(var.rda)
plot.pvca <- mbecPVCAStatsPlot(var.pvca)
plot.scoef <- mbecSCOEFStatsPlot(var.scoef)


mbecReportPrelim(rld.clr, model.vars = c("Batch","location"))






# REAL-LIFE DATA CORRECTIONS ----------------------------------------------
# load package
devtools::load_all()
# load data
rld <- readRDS(file=file.path("/Users/molbrich/Projects/Axel/Dog.Microbiome","data","RDS","phyloseq.OTU.rds"), refhook = NULL)
# remove the 'stool' samples
rld.nostool <- phyloseq::subset_samples(rld, Batch!="stool")
#rld.nostool <- phyloseq::subset_samples(rld, location=="palm")
# agglomerate to 'genus' level
rld.nostool <- tax_glom(rld.nostool, "Genus")
#test <- mbecProcessInput(rld, required.col = c("location","Batch"))
rld.clr <- LRTransform(rld.nostool, method="CLR", offset=0.01, required.col = c("location","Batch"))


## BMC
rld.clr.bmc <- mbecCorrection(rld.clr, model.vars=model.vars,method="bmc",update=TRUE, nc.features=NULL)
## BMC - var
var.linmod.bmc <- mbecModelVariance(rld.clr.bmc, model.vars=model.vars, method="lm",type="prelim")
var.rda.bmc <- mbecModelVariance(rld.clr.bmc, model.vars=model.vars, method="rda",type="prelim")
var.pvca.bmc <- mbecModelVariance(rld.clr.bmc, model.vars=model.vars, method="pvca",type="prelim")
var.scoef.bmc <- mbecModelVariance(rld.clr.bmc, model.vars=model.vars, method="s.coef",type="prelim")
## BMC - var - plot
plot.linmod.bmc <- mbecVarianceStatsPlot(var.linmod.bmc)
plot.rda.bmc <- mbecRDAStatsPlot(var.rda.bmc)
plot.pvca.bmc <- mbecPVCAStatsPlot(var.pvca.bmc)
plot.scoef.bmc <- mbecSCOEFStatsPlot(var.scoef.bmc)

## ComBat
rld.clr.bat <- mbecCorrection(rld.clr, model.vars=c("Batch","group"),method="bat",update=TRUE, nc.features=NULL)
## ComBat - var
var.linmod.bat <- mbecModelVariance(rld.clr.bat, model.vars=model.vars, method="lm",type="prelim")
var.rda.bat <- mbecModelVariance(rld.clr.bat, model.vars=model.vars, method="rda",type="prelim")
var.pvca.bat <- mbecModelVariance(rld.clr.bat, model.vars=model.vars, method="pvca",type="prelim")
var.scoef.bat <- mbecModelVariance(rld.clr.bat, model.vars=model.vars, method="s.coef",type="prelim")
## ComBat - var - plot
plot.linmod.bat <- mbecVarianceStatsPlot(var.linmod.bat)
plot.rda.bat <- mbecRDAStatsPlot(var.rda.bat)
plot.pvca.bat <- mbecPVCAStatsPlot(var.pvca.bat)
plot.scoef.bat <- mbecSCOEFStatsPlot(var.scoef.bat)

## RBE - HMMMMM
rld.clr.rbe <- mbecCorrection(rld.clr, model.vars=c("Batch","group"),method="rbe",update=TRUE, nc.features=NULL)
## RBE - var
var.linmod.rbe <- mbecModelVariance(rld.clr.rbe, model.vars=model.vars, method="lm",type="prelim")
var.rda.rbe <- mbecModelVariance(rld.clr.rbe, model.vars=model.vars, method="rda",type="prelim")
var.pvca.rbe <- mbecModelVariance(rld.clr.rbe, model.vars=model.vars, method="pvca",type="prelim")
var.scoef.rbe <- mbecModelVariance(rld.clr.rbe, model.vars=model.vars, method="s.coef",type="prelim")
## RBE - var - plot
plot.linmod.rbe <- mbecVarianceStatsPlot(var.linmod.rbe)
plot.rda.rbe <- mbecRDAStatsPlot(var.rda.rbe)
plot.pvca.rbe <- mbecPVCAStatsPlot(var.pvca.rbe)
plot.scoef.rbe <- mbecSCOEFStatsPlot(var.scoef.rbe)

# ToDo: sanity check as in ComBat - https://support.bioconductor.org/p/9139601/
rld.clr.svd <- mbecCorrection(rld.clr, model.vars=c("group","Batch"),method="svd",update=TRUE, nc.features=NULL)
## SVD - var
var.linmod.svd <- mbecModelVariance(rld.clr.svd, model.vars=model.vars, method="lm",type="prelim")
var.rda.svd <- mbecModelVariance(rld.clr.svd, model.vars=model.vars, method="rda",type="prelim")
var.pvca.svd <- mbecModelVariance(rld.clr.svd, model.vars=model.vars, method="pvca",type="prelim")
var.scoef.svd <- mbecModelVariance(rld.clr.svd, model.vars=model.vars, method="s.coef",type="prelim")
## SVD - var - plot
plot.linmod.svd <- mbecVarianceStatsPlot(var.linmod.svd)
plot.rda.svd <- mbecRDAStatsPlot(var.rda.svd)
plot.pvca.svd <- mbecPVCAStatsPlot(var.pvca.svd)
plot.scoef.svd <- mbecSCOEFStatsPlot(var.scoef.svd)



# put all in one object
leTest <- mbecSetData(rld.clr, new.cnts=mbecGetData(rld.clr.bmc, orientation="sxf")[[1]], log="bmc", type="bmc", update=FALSE)
leTest <- mbecSetData(leTest, new.cnts=mbecGetData(rld.clr.bat, orientation="sxf")[[1]], log="bat", type="bat", update=FALSE)
leTest <- mbecSetData(leTest, new.cnts=mbecGetData(rld.clr.rbe, orientation="sxf")[[1]], log="rbe", type="rbe", update=FALSE)
leTest <- mbecSetData(leTest, new.cnts=mbecGetData(rld.clr.svd, orientation="sxf")[[1]], log="svd", type="svd", update=FALSE)

# save the shit out of it

saveRDS(leTest, file=file.path("/Users/molbrich/Projects/Axel/Dog.Microbiome","data","RDS","phyloseq.OTU.BC.rds"), compress = TRUE)





# BETTER RLE HARDCODED --> ToDo: automate ---------------------------------


leTest <- mbecRLE(rld.clr, model.vars=c("location","Batch"), return.data = TRUE) ## this is too large

unique(tmp.long$location)

tmp.rle1 <- tmp.long[which(tmp.long$location %in% c("abdomen","axilla","ear","elbow")),]
tmp.rle2 <- tmp.long[which(tmp.long$location %in% c("flank","frontPaw","groin","hindPaw")),]
tmp.rle3 <- tmp.long[which(tmp.long$location %in% c("lips","palm","perineum","tail")),]




rle.plot1 <- ggplot2::ggplot(
  tmp.rle1,
  ggplot2::aes(x = specimen, y = values,fill = get(model.vars[2]))) +
  ggplot2::stat_boxplot(color = "black", notch = FALSE,
                        outlier.colour = "#E42032", outlier.fill = "white",
                        outlier.shape = 1, outlier.stroke = 0.5,
                        outlier.size = 0.5, outlier.alpha=0.5) +
  #ggplot2::facet_wrap(facets = ggplot2::vars(get(model.vars[1])), ncol=4, nrow=3) +
  ggplot2::facet_grid(cols = ggplot2::vars(get(model.vars[1])), scales = "free",
                      space = "free_x", drop = TRUE) +
  ggplot2::scale_fill_manual(values = cols) +
  theme_rle() + ggplot2::guides(fill = ggplot2::guide_legend(title = ggplot2::element_blank()))

rle.plot2 <- ggplot2::ggplot(
  tmp.rle2,
  ggplot2::aes(x = specimen, y = values,fill = get(model.vars[2]))) +
  ggplot2::stat_boxplot(color = "black", notch = FALSE,
                        outlier.colour = "#E42032", outlier.fill = "white",
                        outlier.shape = 1, outlier.stroke = 0.5,
                        outlier.size = 0.5, outlier.alpha=0.5) +
  #ggplot2::facet_wrap(facets = ggplot2::vars(get(model.vars[1])), ncol=4, nrow=3) +
  ggplot2::facet_grid(cols = ggplot2::vars(get(model.vars[1])), scales = "free",
                      space = "free_x", drop = TRUE) +
  ggplot2::scale_fill_manual(values = cols) +
  theme_rle() + ggplot2::guides(fill = ggplot2::guide_legend(title = ggplot2::element_blank()))

rle.plot3 <- ggplot2::ggplot(
  tmp.rle3,
  ggplot2::aes(x = specimen, y = values,fill = get(model.vars[2]))) +
  ggplot2::stat_boxplot(color = "black", notch = FALSE,
                        outlier.colour = "#E42032", outlier.fill = "white",
                        outlier.shape = 1, outlier.stroke = 0.5,
                        outlier.size = 0.5, outlier.alpha=0.5) +
  #ggplot2::facet_wrap(facets = ggplot2::vars(get(model.vars[1])), ncol=4, nrow=3) +
  ggplot2::facet_grid(cols = ggplot2::vars(get(model.vars[1])), scales = "free",
                      space = "free_x", drop = TRUE) +
  ggplot2::scale_fill_manual(values = cols) +
  theme_rle() + ggplot2::guides(fill = ggplot2::guide_legend(title = ggplot2::element_blank()))


rle.final <- gridExtra::arrangeGrob(rle.plot1 + ggplot2::theme(legend.position = "none"),
                                    rle.plot2 + ggplot2::theme(legend.position = "none"),
                                    rle.plot1, ncol = 1, nrow = 3, heights = c(1,1,1.5))

pdf(file = file.path("split_rle.pdf"), width=12, height=12, fonts="Helvetica", onefile = T)
plot(rle.final)
dev.off()


# RLD COMPARE CORRECTION METHODS ------------------------------------------

rmarkdown::render("mbecReport_prelim_TMPLT.Rmd",
                  params = list(report.data=mbecGetData(rld.clr, orientation="fxs", required.col=eval(model.vars)),
                                report.vars=model.vars,
                                report.list=prelim.data))


rld.corrected <- readRDS(file=file.path("/Users/molbrich/Projects/Axel/Dog.Microbiome","data","RDS","phyloseq.OTU.BC.rds"), refhook = NULL)

model.vars <- c("group","Batch")
input.obj <- rld.corrected

## Wrapper function for the analyses required to construct preliminary or post-correction reports.
#'@param input.obj, list of phyloseq objects to compare, first element is considered uncorrected data
#'@param formula, model formula or sth. similar to select groups to compare
mbecCompare <- function(input.obj, model.vars=c("group","batch")) {

  if( length(input.obj@transformations) != 0 ) {

    prelim.data <- mbecReportPrelim(rld.clr, model.vars, return.data = TRUE)


  }





  ## 3?! cases
  if( length(input.obj) > 1 ) {
    # handle list of 'MbecData' objects
    # list-length determines panel setup
    n.elem <- length(input.obj)
    n.otu <- 10

    ### GRAPHICAL REPRESENTATION ###
    # select TOP features from reference set to compare with corrected data sets
    box.topn <- mbecBox(input.obj[[1]], method="TOP",n=n.otu, model.var = model.vars[2], return.data = T)[[2]]

    ## 1. prepare plots to build panels from
    plot.list <- list()
    for( i.idx in 1:length(input.obj) ) {

      # if input is class 'phyloseq'
      if( class(input.obj[[i.idx]]) %in% "phyloseq" ) {
        f.title = "none"
      } else {
        f.title = attr(input.obj[[i.idx]], "type")
      }

      plot.list[["PCA"]][[i.idx]] <- mbecPCA(input.obj = input.obj[[i.idx]], model.vars=model.vars, pca.axes=c(1,2))
      plot.list[["RLE"]][[i.idx]] <- mbecRLE(input.obj[[i.idx]], return.data = F)
      plot.list[["HEAT"]][[i.idx]] <- mbecHeat(input.obj[[i.idx]], model.vars=model.vars)

      # temporary subset to features selected in 'box.topn' for box-plots
      tmp_subset <- prune_taxa(box.topn, input.obj[[i.idx]])
      plot.list[["BOX"]][[i.idx]] <- mbecBox(tmp_subset, method="ALL", model.var = model.vars[2], return.data = F)
    }

    ## 2. compute some numerical values
    stats.list <- NULL
    rda.array <- NULL
    pvca.array <- NULL
    silhouette.array <- NULL

    for( i.idx in seq_along(attr(input.obj, "transformations")) ) {
      # tmp workaround to  handle transformations
      tmp.mbec <- mbecSetData(input.obj, input.obj@transformations[[i.idx]],type=names(input.obj@transformations)[i.idx], update=TRUE)
      # variance per feature
      stats.list <- rbind.data.frame(stats.list, mbecModelVariance(tmp.mbec, model.vars=model.vars, method="lm",
                                                             type=ifelse(is.null(attr(tmp.mbec, "type")), "none", attr(tmp.mbec, "type"))))

      rda.array <- rbind.data.frame(rda.array, mbecModelVariance(tmp.mbec, model.vars=model.vars, method="rda",
                                                                 type=ifelse(is.null(attr(tmp.mbec, "type")), "none", attr(tmp.mbec, "type"))))

      pvca.array <- rbind.data.frame(pvca.array, mbecModelVariance(tmp.mbec, model.vars=model.vars, method="pvca",
                                                                   type=ifelse(is.null(attr(tmp.mbec, "type")), "none", attr(tmp.mbec, "type"))))

      silhouette.array <- rbind.data.frame(silhouette.array, mbecModelVariance(tmp.mbec, model.vars=model.vars, method="s.coef",
                                                                               type=ifelse(is.null(attr(tmp.mbec, "type")), "none", attr(tmp.mbec, "type"))))

    }

  } else if( length(names(attr(input.obj,"transformations"))) ) {
    # handle 'MbecData' object with a list of transformed counts
    ## get list of transformations + original 'raw' data

  } else {
    # either list of MbecData object with transformations
    # or single object without transformations
    message("Your particular input can not be processed at the moment.")
  }
}


# cut my life into pieces
plot.linmod.all <- mbecVarianceStatsPlot(stats.list)
plot.rda.all <- mbecRDAStatsPlot(rda.array)
plot.pvca.all <- mbecPVCAStatsPlot(pvca.array)
plot.scoef.all <- mbecSCOEFStatsPlot(silhouette.array)



pdf(file = file.path("Correction_variance_estimates_Genus_palm.pdf"), width=12, height=6, fonts="Helvetica", onefile = T)
plot(plot.linmod.all)
plot(plot.rda.all)
plot(plot.pvca.all)
plot(plot.scoef.all)
dev.off()









# INCLUDE A PROGRESS-BAR --------------------------------------------------

# https://ryouready.wordpress.com/2009/03/16/r-monitor-function-progress-with-a-progress-bar/

progressBarDemo <- function() {
  numberOfExecutions <- 20
  # Create the progress bar
  pb <- txtProgressBar(min = 0, max = numberOfExecutions)

  for(i in 1:numberOfExecutions) {
    Sys.sleep(0.2)
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }

  close(pb)
}

progressBarDemo()

#








# MAKE A CORRELATION PLOT -------------------------------------------------
library(reshape2)
library(ggplot2)

mydata <- mtcars[, c(1,3,4,5,6,7)]
cormat <- round(cor(mydata),2)
head(cormat)
melted_cormat <- melt(cormat)


## Use 'geom_raster()' because it is high performance version of 'geom_tile()'
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_raster()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri


# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  ggplot2::geom_raster()+
  ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  ggplot2::theme_minimal()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  ggplot2::coord_fixed()


## order by correlation
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_raster(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)


## aaaand make it pretty
ggheatmap +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


## prep correlation plot
# IF model.vars is NULL function takes all sample variables
mbecCORR <- function(input.obj, model.vars = NULL, return.data = FALSE) {

  tmp <- mbecGetData(input.obj = input.obj, orientation = "fxs",
                     required.col = eval(model.vars))

  tmp.meta <- tmp[[2]]

  # 1. split into numerical and categorical variables
  # 2. compute correlation for numerical vars
  # 3. compute canonical correlation for categorical variables




  ## rank correlation doesn't really work here
  loc.rank <- rank(tmp.meta$location)
  Bat.rank <- rank(tmp.meta$Batch)
  cor(loc.rank, Bat.rank)

  mA <- model.matrix(~ location - 1, tmp.meta)
  mB <- model.matrix(~ Batch - 1, tmp.meta)
  cor(mA, mB)

  cancor(mA, mB)$cor

  cancor(mB, mA)$cor



  cormat <- round(cor(test),2)








  vcov(tmp.meta)

  str(tmp.meta)

  test <- apply(tmp.meta, 2, factor)


  for (g.idx in unique(tmp.meta[, eval(model.vars[1])])) {
    message("Calculating RLE for group: ", g.idx)

    tmp.cnts.group <- dplyr::select(tmp.cnts, tmp.meta$specimen[tmp.meta[,
                                                                         eval(model.vars[1])] %in% g.idx])

    feature.med = apply(tmp.cnts.group, 1, stats::median)

    tmp.group.long <- apply(tmp.cnts.group, 2,
                            function(sample.col) sample.col - feature.med) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(cols = everything(),
                          names_to = "specimen", values_to = "values")
    tmp.long <- rbind.data.frame(tmp.long, tmp.group.long)
  }

  tmp.long <- dplyr::left_join(tmp.long, tmp.meta,
                               by = "specimen") %>%
    dplyr::mutate(plot.order = paste(get(model.vars[1]),
                                     get(model.vars[2]), sep = "_")) %>%
    dplyr::arrange(plot.order) %>%
    dplyr::mutate(specimen = factor(specimen, levels = unique(specimen)))

  if (return.data) {
    return(tmp.long)
  }
  return(mbecRLEPlot(tmp.long, model.vars, cols))
}




















# TOMFOOLERY --------------------------------------------------------------

tmp.meta <- tmp[[2]]

bla <- tibble::rownames_to_column(tmp.meta, var = "sID")

# calculate IQR and sort counts in decreasing order
iqr <- apply(tmp.cnts, 2, stats::IQR)
tmp.cnts <- tmp.cnts[, order(iqr, decreasing = TRUE)]

PCA <- stats::prcomp(tmp.cnts, scale = FALSE)

axes.number <- dim(PCA$x)[2]
axes.names <- paste("PC", seq_len(axes.number), sep = "")

plot.df <- PCA$x %>%
  data.frame(stringsAsFactors = FALSE) %>%
  dplyr::rename_at(seq_len(axes.number), ~axes.names) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::left_join(tmp.meta, by = "sample")

all(a %in% b)

a <- rownames(datadummy$meta)

b <- rownames(datadummy$cnts)
c <- colnames(datadummy$cnts)

identical(a,b)


test.phylo <- phyloseq(otu_table(datadummy$cnts, taxa_are_rows=F), sample_data(datadummy$meta))

test.meta <- datadummy$meta[c(1:20, 41:75, 21:40),]

identical(a,rownames(test.meta))

all(a %in% rownames(test.meta))

# test unordered names
test.phylo.2 <- phyloseq(otu_table(datadummy$cnts, taxa_are_rows=F), sample_data(test.meta))




# QUICK CAST TO DATAFRAME -------------------------------------------------
library("microbenchmark")

# Unit: seconds
# min       lq     mean   median       uq      max neval
# 11.43184 12.15854 12.52861 12.32509 12.76993 14.82869   100
microbenchmark(
  tmp.cnts <- data.frame(phyloseq::otu_table(input.obj), check.names = FALSE)
)

# Unit: microseconds
# min       lq     mean   median      uq     max neval
# 348.647 368.8625 446.1846 382.0095 467.176 981.608   100
microbenchmark(
  le.matrix <- as(phyloseq::otu_table(input.obj),"matrix")
)

# Unit: milliseconds
# min      lq     mean  median       uq      max neval
# 93.42146 100.912 148.5012 113.734 160.5625 998.5108   100
microbenchmark(
  le.frame <- data.frame(le.matrix)
)

# Unit: milliseconds
# min       lq     mean   median      uq      max neval
# 81.24113 85.49739 130.1809 91.09531 145.696 1079.763   100
microbenchmark(
  le.frame <- otuCast(input.obj)
)

otuCast <- function(leInput) {
  le.matrix <- as.data.frame(as(phyloseq::otu_table(leInput),"matrix"))

}



# SET-UP COLORMAP ---------------------------------------------------------

scales::show_col(palette())

###....



# TEST ESTIMATIBILITY OF MODEL-MATRIX -------------------------------------

## 1. a way to build model matrix out of multiple covariates

dd <- data.frame(a = gl(3,4), b = gl(4,1,12)) # balanced 2-way
options("contrasts") # typically 'treatment' (for unordered factors)
model.matrix(~ a + b, dd)

test.vars <- c("a","b")

paste("~ ", paste(test.vars, collapse=" + "))

model.matrix(as.formula(paste("~ ", paste(test.vars, collapse=" + "))), dd)

eTest <- stats::model.matrix(as.formula(paste("~ ", paste(model.vars[[-1]], collapse=" + "))), tmp.meta)

eTest <- stats::model.matrix(model.form, tmp.meta)


model.frame(model.formula,tmp.meta)[, -1]

labels(terms(model.formula))

terms.formula(model.formula)

limma::nonEstimable(tmp.mod)

test.formula <- model.formula[-2]

test2 <- test.formula[-2]

leTest <- stats::model.matrix( ~ tmp.meta[,model.vars])



mbecTestModel <- function(input.obj, model.vars, model.form=NULL) {
  # 1. extract covariate information
  tmp <- mbecGetData(input.obj, orientation="fxs")
  tmp.cnts <- tmp[[1]]; tmp.meta <- tmp[[2]]

  # 2. check if model-formula was supplied
  if( is.null(model.form) ) {
    # construct linear model from covariates
    message("Construct lm-formula from covariates.")

    model.form = stats::as.formula(paste("y", " ~ ", paste(model.vars, collapse = " + ")))

  }
  # if model.form is complete --> remove LHS
  if( length(model.form) == 3 ) model.form <- model.form[-2]

  # create model matrix from RHS
  model.mtx <- stats::model.matrix(model.form, tmp.meta)

  res.est <- limma::nonEstimable(model.mtx)

  if( !is.null(res.est) ) {
    message("There is a problem with the estimatibility of your model. Check out covariate: ", paste("'",res.est, "'", sep="", collapse=", "))
  }
   return(res.est)
}


model.form = stats::as.formula(paste("y", " ~ ", paste(model.vars[-4], collapse = " + ")))

model.form = stats::as.formula("~Batch + samplingSeason:housing + diet")


model.formula

gsub("Batch","",model.formula)

# LM Stuff
if (!is.null(model.form)) {
  message("Use provided model formula.")
  tmp.formula <- stats::as.formula(model.form)
} else {
  message("Construct formula from covariates.")
  tmp.formula = stats::as.formula(paste("y", " ~ ", paste(model.vars, collapse = " + ")))
}



tmp.group.p <- apply(tmp.cnts, 2, FUN = function(y){
  nc.lm <- stats::lm(y ~ get(model.vars[1]) + get(model.vars[2]), data = tmp.meta)
  nc.lm.summary <- summary(nc.lm)
  # extract p-value of group (treatment)
  p <- nc.lm.summary$coefficients[2,4]
})



# LMM Stuff
if (!is.null(model.form)) {
  message("Use provided model formula.")
  tmp.formula <- stats::as.formula(model.form)
} else {
  message("Construct formula from covariates.")
  f.terms <- paste("(1|", model.vars, ")", sep = "")

  tmp.formula <- stats::as.formula(paste(paste("y", paste(model.vars[-1], collapse = " + "), sep = " ~ "),
                                         paste(f.terms[1], collapse = " + "), sep = " + "))

}



####
f.terms <- paste("(1|",model.vars[2],")", sep="")

tmp.group.p <- apply(tmp.cnts, 2, FUN = function(x) {

  tmp.formula <- stats::as.formula(paste(paste("x", model.vars[1], sep=" ~ "), paste(f.terms[], collapse=" + "), sep=" + "))
  nc.lmm <- eval(bquote(lmerTest::lmer(.(tmp.formula), data = tmp.meta)))
  nc.lmm.summary <- summary(nc.lmm)
  p <- nc.lmm.summary$coefficients[2,5]

}
)






# CHECKS ------------------------------------------------------------------

dummies <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","data_dummies.rds"), refhook = NULL)

sponge <- LRTransform(dummies$sponge, method="CLR")
ad <- LRTransform(dummies$ad, method="CLR")
hd <- LRTransform(dummies$hd, method="CLR")

test.ds <- sponge

test.ds.corr <- mbecCorrection(test.ds, model.vars=c("batch","group"),method="bmc",update=FALSE, nc.features=NULL)
test.ds.corr <- mbecCorrection(test.ds.corr, model.vars=c("batch","group"),method="bat",update=FALSE, nc.features=NULL)
test.ds.corr <- mbecCorrection(test.ds.corr, model.vars=c("batch","group"),method="rbe",update=FALSE, nc.features=NULL)
test.ds.corr <- mbecCorrection(test.ds.corr, model.vars=c("batch","group"),method="pn",update=FALSE, nc.features=NULL)
test.ds.corr <- mbecCorrection(test.ds.corr, model.vars=c("batch","group"),method="svd",update=FALSE, nc.features=NULL)
test.ds.corr <- mbecCorrection(test.ds.corr, model.vars=c("batch","group"),method="ruv3",update=FALSE, nc.features=NULL)



test.ds.corr.sponge <- test.ds.corr


input.obj <- test.ds.corr.sponge


## For an mbecObject with full transformations list
mbecPipelineEval <- function(input.obj, model.vars=c("batch","group")) {

  # just start with input processing and transform into MbecData if required
  input.obj <- mbecProcessInput(input.obj, required.col=eval(model.vars))

  n.elem <- length(input.obj@transformations) # number batch corrected matrices
  n.otu <- 10 # number of features to select in heatmap and box plots
  # get the top most variable features for further comparisons
  box.topn <- mbecBox(input.obj, method="TOP",n=n.otu, model.var = model.vars[1], return.data = T)[[2]]

  # Prepare the exploratory plots
  post.report.list <- list()
  post.report.list[["mosaic"]] <- mbecMosaic(input.obj, model.vars=eval(model.vars) )
  post.report.list[["pca"]][["Base"]] <- mbecPCA(input.obj, model.vars=eval(model.vars), pca.axes = c(1,2))
  post.report.list[["rle"]][["Base"]] <- mbecRLE(input.obj, model.vars=eval(model.vars))
  post.report.list[["heat"]][["Base"]] <- mbecHeat(input.obj, method=box.topn, n=10, model.vars=eval(model.vars))
  post.report.list[["box"]][["Base"]] <- mbecBox(input.obj, method=box.topn, n=5, model.var=eval(model.vars)[1] )
  # variance calculations
  post.report.list[["stats"]] <- mbecModelVariance(input.obj, model.vars=eval(model.vars), method="lm",
                                                      type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  post.report.list[["rda"]] <- mbecModelVariance(input.obj, model.vars=eval(model.vars), method="rda",
                                                   type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  post.report.list[["pvca"]] <- mbecModelVariance(input.obj, model.vars=eval(model.vars), method="pvca",
                                                    type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  post.report.list[["scoef"]] <- mbecModelVariance(input.obj, model.vars=eval(model.vars), method="s.coef",
                                                     type=ifelse(is.null(attr(input.obj, "type")), "none", attr(input.obj, "type")))

  for( t.idx in names(input.obj@transformations) ) {
    post.report.list[["pca"]][[t.idx]] <- mbecPCA(mbecGetData(input.obj, transformation=eval(t.idx)), model.vars=eval(model.vars), pca.axes = c(1,2))
    post.report.list[["rle"]][[t.idx]] <- mbecRLE(mbecGetData(input.obj, transformation=eval(t.idx)), model.vars=eval(model.vars))
    post.report.list[["heat"]][[t.idx]] <- mbecHeat(mbecGetData(input.obj, transformation=eval(t.idx)), method=box.topn, n=10, model.vars=eval(model.vars))
    post.report.list[["box"]][[t.idx]] <- mbecBox(mbecGetData(input.obj, transformation=eval(t.idx)), method=box.topn, n=5, model.var=eval(model.vars)[1] )

    # variance per feature
    post.report.list[["stats"]] <- rbind.data.frame(post.report.list[["stats"]], mbecModelVariance(mbecGetData(input.obj, transformation=eval(t.idx)), model.vars=model.vars, method="lm",
                                                                 type=eval(t.idx)))

    post.report.list[["rda"]] <- rbind.data.frame(post.report.list[["rda"]], mbecModelVariance(mbecGetData(input.obj, transformation=eval(t.idx)), model.vars=model.vars, method="rda",
                                                               type=eval(t.idx)))

    post.report.list[["pvca"]] <- rbind.data.frame(post.report.list[["pvca"]], mbecModelVariance(mbecGetData(input.obj, transformation=eval(t.idx)), model.vars=model.vars, method="pvca",
                                                                 type=eval(t.idx)))

    post.report.list[["scoef"]] <- rbind.data.frame(post.report.list[["scoef"]], mbecModelVariance(mbecGetData(input.obj, transformation=eval(t.idx)), model.vars=model.vars, method="s.coef",
                                                                             type=eval(t.idx)))
  }

  # cut my life into pieces
  post.report.list[["stats"]] <- mbecVarianceStatsPlot(post.report.list[["stats"]])
  post.report.list[["rda"]] <- mbecRDAStatsPlot(post.report.list[["rda"]])
  post.report.list[["pvca"]] <- mbecPVCAStatsPlot(post.report.list[["pvca"]])
  post.report.list[["scoef"]] <- mbecSCOEFStatsPlot(post.report.list[["scoef"]])

  if( return.data ) {
    return(post.report.list)
  }

  # construct the report
  rep.data <- mbecGetData(input.obj, orientation="fxs", required.col=eval(model.vars))

  rmarkdown::render("mbecReport_post.Rmd",
                    params = list(report.data=rep.data,
                                  report.vars=model.vars,
                                  report.list=post.report.list))



  print(post.report.list$heat$Base)




}



gridExtra::marrangeGrob(grobs=post.report.list$pca, ncol = 2, nrow=2)









post.report.list[["pca"]][[t.idx]] <- mbecPCA(mbecGetData(input.obj, transformation=eval("bmc")), model.vars=eval(model.vars), pca.axes = c(1,2))






# METHOD VALIDATION -------------------------------------------------------

dummies <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","data_dummies.rds"), refhook = NULL)

sponge <- LRTransform(dummies$sponge, method="CLR")
ad <- LRTransform(dummies$ad, method="CLR")
hd <- LRTransform(dummies$hd, method="CLR")

test.ds <- sponge




# SUMMARY STATISTICS ------------------------------------------------------

# install.packages("devtools")
devtools::install_github("ropensci/skimr")











