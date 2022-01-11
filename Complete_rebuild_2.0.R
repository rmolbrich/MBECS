# TESTING NEW CLASS CONFIGURATION -----------------------------------------
# load package
devtools::load_all()

# A test-set that contains three studies (sponge, ad and hd)
dummies <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","data_dummies.rds"), refhook = NULL)

# 1. convert to MbecData
sponge <- mbecProcessInput(dummies$sponge, c("sample", "group", "batch"))

# 2. test getter-method
getTest <- mbecGetData(sponge, orientation="fxs", required.col=c("sample", "group", "batch"),
                         type="otu")

# 3. test setter-method
setTest <- mbecSetData(sponge, new.cnts=getTest[[1]], type="clr", label=NULL)
all(getTest[[1]] == setTest@clr)


# 4. test clr-transformation
sponge <- mbecTransform(sponge, method="clr")

getTest <- mbecGetData(sponge, orientation="fxs", required.col=c("sample", "group", "batch"),
                       type="clr")



# TEST ANALYSIS METHODS ---------------------------------------------------

# load package
devtools::load_all()

# A test-set that contains three studies (sponge, ad and hd)
dummies <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","data_dummies.rds"), refhook = NULL)

# convert to MbecData
test.obj <- mbecProcessInput(dummies$ad, c("sample", "group", "batch"))
# apply transformation
test.obj <- mbecTransform(dummies$ad, method="clr")

# RLE
mbecRLE(test.obj, model.vars = c("batch","group"), type="otu", return.data = FALSE)
mbecRLE(test.obj, model.vars = c("batch","group"), type="clr", return.data = FALSE)

# PCA
mbecPCA(test.obj, model.vars = c("batch","group"), type="otu", return.data = FALSE)
mbecPCA(test.obj, model.vars = c("batch","group"), type="clr", return.data = FALSE)

# BOX
p.Box <- mbecBox(input.obj=test.obj, method="TOP", n=4, model.var="batch", type="otu", return.data=FALSE)
gridExtra::marrangeGrob(p.Box, nrow=2, ncol=2)
p.Box <- mbecBox(input.obj=test.obj, method="TOP", n=4, model.var="batch", type="clr", return.data=FALSE)
gridExtra::marrangeGrob(p.Box, nrow=2, ncol=2)

# HEAT
p.Heat <- mbecHeat(input.obj=test.obj, model.vars=c("group","batch"), center=TRUE,
                   scale=TRUE, method="TOP", n=10, type="otu", return.data=F)
p.Heat <- mbecHeat(input.obj=test.obj, model.vars=c("group","batch"), center=TRUE,
                   scale=TRUE, method="TOP", n=10, type="clr", return.data=F)

# MOSAIC
p.Mosaic <- mbecMosaic(input.obj=test.obj, model.vars=c("group","batch"), return.data=FALSE)



# TEST VARIANCE CALCULATION METHODS ---------------------------------------


# load package
devtools::load_all()

# A test-set that contains three studies (sponge, ad and hd)
dummies <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","data_dummies.rds"), refhook = NULL)

# convert to MbecData
test.obj <- mbecProcessInput(dummies$ad, c("sample", "group", "batch"))
# apply transformation
test.obj <- mbecTransform(test.obj, method="clr")

### DEBUG ###
input.obj <- test.obj
model.vars <- c("batch","group") # for palm use group instead of batch
### DEBUG ###

# calculate the variance statistics
var.linmod <- mbecModelVariance(test.obj, model.vars=model.vars, method="lm",
                                type="clr")

var.linmixmod <- mbecModelVariance(test.obj, model.vars=model.vars, method="lmm",
                                   type="clr")

var.rda <- mbecModelVariance(test.obj, model.vars=model.vars, method="rda",
                             type="clr")

# FixMe - sth. is fishy here
var.pvca <- mbecModelVariance(test.obj, model.vars=model.vars, method="pvca",
                              type="clr")

var.scoef <- mbecModelVariance(test.obj, model.vars=model.vars, method="s.coef",
                               type="clr")

plot.linmod <- mbecVarianceStatsPlot(var.linmod)
plot.linmixmod <- mbecVarianceStatsPlot(var.linmixmod)
plot.rda <- mbecRDAStatsPlot(var.rda)
plot.pvca <- mbecPVCAStatsPlot(var.pvca)
plot.scoef <- mbecSCOEFStatsPlot(var.scoef)


mbecReportPrelim(test.obj, model.vars = c("batch","group"))





# TEST HELPER FUNCTIONS ---------------------------------------------------
# load package
devtools::load_all()

# A test-set that contains three studies (sponge, ad and hd)
dummies <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","data_dummies.rds"), refhook = NULL)

# convert to MbecData
test.obj <- mbecProcessInput(dummies$ad, c("sample", "group", "batch"))
# apply transformation
test.obj <- mbecTransform(test.obj, method="clr")
test.obj <- mbecTransform(test.obj, method="tss")

colSums(test.obj@tss)

### DEBUG ###
input.obj <- test.obj
model.vars <- c("batch","group") # for palm use group instead of batch
### DEBUG ###

mbecLM(test.obj, method="lmm", model.vars=c("batch","group"), type="otu", label=character())


# TEST CORRECTIONS --------------------------------------------------------
devtools::load_all()
dummies <- readRDS(file=file.path("/Users/molbrich/WorkBuschlab/Thesis/workspace","data","RDS","data_dummies.rds"), refhook = NULL)
test.obj <- mbecProcessInput(dummies$ad, c("sample", "group", "batch"))
test.obj <- mbecTransform(test.obj, method="clr")
test.obj <- mbecTransform(test.obj, method="tss")

test.bmc <- mbecBMC(input.obj=test.obj, model.vars=c("batch","group"), type = "otu")
test.bat <- mbecBat(input.obj=test.obj, model.vars=c("batch","group"), type = "tss")
test.rbe <- mbecRBE(input.obj=test.obj, model.vars=c("batch","group"), type = "clr")
test.pn <- mbecPN(input.obj=test.obj, model.vars=c("batch","group"), type = "tss")

test.svd <- mbecSVD(input.obj=test.obj, model.vars=c("batch","group"), type = "tss")

ps.bmc <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "bmc")
ps.bat <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "bat")
ps.rbe <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "rbe")
ps.fab <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "fab")
ps.pn <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "pn")
ps.svd <- mbecCorrection(input.obj = datadummy, model.vars=c("group","batch"), method = "svd")

# TESTING WITH REAL-LIFE DATA ---------------------------------------------
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



# FIX DEVTOOLS::CHECK ERRORS ----------------------------------------------

mbecBox(input.obj=datadummy, method='TOP', n=15,model.var='batch', type="otu", return.data=FALSE)



# SIMULATE DUMMY DATA -----------------------------------------------------


# Take dog-data and calculate probabillities of all taxa in every sample and every subgroup i.e. location.
# Thats basically the probability for any taxa to be sampled.. how to draw from this? --> for library size o, draw o times from a probability for n otus
# HowTo: aggregate the probability vectors to generate sth. to draw from?
# HowTo: simulate library size difference
# HowTo: account for sampling fraction

# jus use to create some artificial distribution with dummy meta data for now
mbecDummy <- function(n.otus=500, n.samples=40) {

  # build dataset
  # cnts <- matrix(1:16, nrow=4, ncol=4,
  #                dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4")))
  # cnts.norm <- matrix(rep.int(c(50,100), times=8), nrow=4, ncol=4,
  #                     dimnames=list(c("A","B","C","D"), c("F1","F2","F3","F4")))
  # meta <- data.frame("sID"=c("A","B","C","D"),
  #                    "group"=factor(c("X","X","Y","Y")),
  #                    "batch"=factor(c(1,2,1,2)), row.names = "sID")

  replace1 <- function(x) {
    dplyr::if_else(x != 0, 1, 0)
  }


  paste("R", c(1,2,1,2,3,4,3,4,5,6,5,6,7,8,7,8,9,10,9,10,11,12,11,12,13,14,13,14,15,16,15,16,17,18,17,18,19,20,19,20,21,22,21,22,23,24,23,24,25,26,25,26,27,28,27,28,29,30,29,30,31,32,31,32,33,34,33,34,35,36,35,36,37,38,37,38,39,40,39,40),sep="")



  # create meta-data
  meta <- data.frame("sID"=paste("S", 1:n.samples, sep=""),
                     "group"=factor(c(rep("A", times=n.samples/2),rep("B", times=n.samples/2))),
                     "batch"=factor(rep(c("B1","B2"), times=n.samples/2)),
                     "replicate"=factor(paste("R", c(1,2,1,2,3,4,3,4,5,6,5,6,7,8,7,8,9,10,9,10,11,12,11,12,13,14,13,14,
                                                     15,16,15,16,17,18,17,18,19,20,19,20),sep="")),
                     row.names = "sID")

  base_mtx <- Matrix::rsparsematrix(n.samples, n.otus, 0.314, symmetric = FALSE) %>%  # NxM sparse matrix
    as.matrix() %>% as.data.frame() %>%                              # transform to usable format
    dplyr::rename_with(., ~ gsub("V","OTU", .x)) %>%                        # rename columns to OTUxx
    dplyr::mutate(across(everything(), replace1))                           # replace non-zero values with '1'

  rownames(base_mtx) <- rownames(meta)

  # add groups and batches to the matrix
  #base_mtx$group <- factor(rep(c("case","control"),each=20), levels = c("control","case"))
  #base_mtx$batch <- factor(rep(c("b1","b2","b1","b2"),each=10), levels = c("b1","b2"))

  # prep mtx for systematic and non-systematic BE
  cnts.sys <- base_mtx
  cnts.nonsys <- base_mtx

  ## DEFINE distributions for NON-SYS, OTU, TREATMENT, BE
  # random impact of non-systematic BE
  nsys <- sample(0:1,size=n.otus,replace=TRUE)
  # some arbitrary distribution for abundance of all OTUs
  otu.dist <- round(rexp(n.otus, .01))
  # some arbitrary distribution for treatment effect
  treat.dist <- rlnorm(n = n.otus, meanlog=3, sdlog=0.5)
  # some arbitrary distribution for batch-effect (10% of max otu abundance as mean)
  batch_1 <- c(1:10,21:30)
  batch.dist <- rnorm(n.otus, mean = 5, sd = 1)


  for( c.idx in 1:n.otus ) {

    ## SYSTEMATIC
    cnts.sys[cnts.sys[,c.idx] == 1, c.idx] <- abs(rnorm(n=length(cnts.sys[cnts.sys[,c.idx] == 1, c.idx]),
                                                        mean=otu.dist[c.idx],
                                                        sd=1))

    # now add 'treatment-effect' to first twenty rows/samples - use normal distribution for all otus for slight variation
    cnts.sys[which(cnts.sys[1:20,c.idx] != 0), c.idx] <- cnts.sys[which(cnts.sys[1:20,c.idx] != 0), c.idx] +
      rnorm(n=length(which(cnts.sys[1:20,c.idx] != 0)),
            mean=treat.dist[c.idx],
            sd=1)

    # B1: now add systematic batch-effect to one batch - use normal distribution for all otus for slight variation
    cnts.sys[which(cnts.sys[batch_1,c.idx] != 0), c.idx] <- cnts.sys[which(cnts.sys[batch_1,c.idx] != 0), c.idx] +
      abs(rnorm(n=length(which(cnts.sys[batch_1,c.idx] != 0)),
                mean=rnorm(1, mean = 5, sd = 1),
                sd=abs(rnorm(1, mean = 0, sd = 2))))

    ## NON-SYSTEMATIC
    cnts.nonsys[cnts.nonsys[,c.idx] == 1, c.idx] <- rnorm(n=length(cnts.nonsys[cnts.nonsys[,c.idx] == 1, c.idx]),
                                                          mean=otu.dist[c.idx],
                                                          sd=1)

    # now add 'treatment-effect' to first twenty rows/samples - use normal distribution for all otus for slight variation
    cnts.nonsys[which(cnts.nonsys[1:20,c.idx] != 0), c.idx] <- cnts.nonsys[which(cnts.nonsys[1:20,c.idx] != 0), c.idx] +
      rnorm(n=length(which(cnts.nonsys[1:20,c.idx] != 0)),
            mean=treat.dist[c.idx],
            sd=1)

    # B1: now add non-systematic batch-effect to one batch - use normal distribution for all otus for slight variation
    if( nsys[c.idx] == 0 ) {
      # just add some zero mean noise here
      cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] <- cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] +
        abs(rnorm(n=length(which(cnts.nonsys[batch_1,c.idx] != 0)),
                  mean=rnorm(1, mean = 0, sd = 1),
                  sd=abs(rnorm(1, mean = 0, sd = 2))))

    } else {
      cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] <- cnts.nonsys[which(cnts.nonsys[batch_1,c.idx] != 0), c.idx] +
        abs(rnorm(n=length(which(cnts.nonsys[batch_1,c.idx] != 0)),
                  mean=rnorm(1, mean = 50, sd = 1),
                  sd=abs(rnorm(1, mean = 0, sd = 2))))
    }
  }


  otu.cnts <- abs(floor(cnts.nonsys))

  dummy.list <- list("cnts"=otu.cnts, "meta"=meta)

  usethis::use_data(dummy.list, overwrite = T)

  dummy.mbec <- mbecTransform(dummy.list, method="tss")
  dummy.mbec <- mbecTransform(dummy.mbec, method="clr")

  usethis::use_data(dummy.mbec, overwrite = T)

}


# Create the simulated data
m <- 50
n <- 20
nc <- 10 # negative controls without treatment effects
p <- 1
k <- 1
ctl <- rep(FALSE, n)
ctl[1:nc] <- TRUE
# treatment effect
X <- matrix(c(rep(0, floor(m/2)), rep(1, ceiling(m/2))), m, p)
beta <- matrix(rnorm(p*n, 5, 1), p, n) #treatment coefficients
beta[ ,ctl] <- 0
# batch effect
W <- as.matrix(rep(0, m), m, k)
W[c(1:12,38:50), 1] <-  1
alpha <- matrix(rnorm(k*n, 5, 1), k, n)
Y_alpha <- sapply(alpha, function(alpha){rnorm(m, mean =  alpha,
                                               abs(rnorm(1, mean = 0, sd = 2)))})
YY_alpha <- apply(Y_alpha, 2, function(x){x*W})

epsilon <- matrix(rnorm(m*n, 0, 1), m, n)
Y <- X%*%beta + YY_alpha + epsilon


# estimate batch coefficient for each OTU
w.cof <- c()
for(i in 1:ncol(Y)){
  res <- lm(Y[ ,i] ~ X + W)
  sum.res <- summary(res)
  w.cof[i] <- sum.res$coefficients[3,1]
}

par(mfrow = c(2,2))
hist(w.cof,col = 'gray')
plot(density(w.cof))
qqnorm(w.cof)
qqline(w.cof, col = 'red')
par(mfrow = c(1,1))











































