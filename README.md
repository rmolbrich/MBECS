# MBECS

The Microbiome Batch-Effect Correction Suite aims to provide a toolkit for 
stringent assessment and correction of batch-effects in microbiome data sets.
To that end, the package offers wrapper-functions to summarize study-design and 
data, e.g., PCA, Heatmap and Mosaic-plots, and to estimate the proportion of 
variance that can be attributed to the batch effect.
The `mbecsCorrection` function acts as a wrapper for various batch effect 
correcting algorithms (BECA) and in conjunction with the aforementioned tools, 
it can be used to compare the effectiveness of correction methods on particular 
sets of data.


## Installation

The `MBECS` package can be installed from Bioconductor. Note that Bioconductor 
follows a "release" and "development" schedule, where the release version is 
considered to be stable and updated every 6 months, and the development version 
contains latest updates.


### Release version

To install the stable release version, install `BiocManager` and the `MBECS` 
package as follows.

```
install.packages("BiocManager")
BiocManager::install("MBECS")
```

### Development version

To install the development version, there are two options.

(i) Install from the Bioconductor as `version = "devel"`. Information on how to 
use the development branch can be found 
[here](http://bioconductor.org/developers/how-to/useDevel/). 


```
install.packages("BiocManager")
BiocManager::install("MBECS", version = "devel")
```

(ii) To install the most current (but not necessarily stable) version, use the 
repository on GitHub:

```
# Use the devtools package to install from a GitHub repository.
install.packages("devtools")

# This will install the MBECS package from GitHub.
devtools::install_github("rmolbrich/MBECS")
```


## Workflow

This is an abridged version that shows the core functionality. For more 
detailed information about the packages functionality and the employed 
algorithms please refer to the package vignette.



### Get started

Load the package via the `library()` function.

```
library(MBECS)
```

The main application of this package is microbiome data. It is common practice 
to use the 
[`phyloseq`](https://bioconductor.org/packages/release/bioc/html/phyloseq.html) 
package for analyses of this type of data. The `MBECS` package extends the 
`phyloseq` class in order to 
provide its functionality. The user can utilize objects of class `phyloseq` or 
a `list` object that contains an abundance table as well as meta data. The 
package contains a dummy data-set of artificially generated data to illustrate 
this process.

Use the `data()`function to load the provided mockup data-sets at this point. 
The only purpose of this data is to illustrate package use. If your use your 
own data in the subsequent steps you can skip this one.

```
# List object 
data(dummy.list)
# Phyloseq object
data(dummy.ps)
# MbecData object
data(dummy.mbec)
```


#### Start from abundance table

For an input that consists of an abundance table and meta-data, both tables 
require sample names as either row or column names. They need to be passed in a 
`list` object with the abundance matrix as first element. The 
`mbecProcessInput()` function will handle the correct orientation and return an 
object of class `MbecData`.

```
# The dummy-list input object comprises two matrices:
names(dummy.list)
```

The optional argument `required.col` may be used to ensure that all covariate 
columns that should be there are available. For the dummy-data these are 
<span style="color: #CD3048;">"group"</span>, <span style="color: #CD3048;">
"batch"</span> and <span style="color: #CD3048;">"replicate"</span>.

```
mbec.obj <- mbecProcessInput(dummy.list, 
                             required.col = c("group", "batch", "replicate"))
```


#### Start from phyloseq object

The start is the same if the data is already of class `phyloseq`. The `dummy.ps`
object contains the same data as `dummy.list`, but it is of class `phyloseq`. 
Create an `MbecData` object from `phyloseq` input. 

The optional argument `required.col` may be used to ensure that all covariate 
columns that should be there are available. For the dummy-data these are 
<span style="color: #CD3048;">"group"</span>, <span style="color: #CD3048;">
"batch"</span> and <span style="color: #CD3048;">"replicate"</span>.

```
mbec.obj <- mbecProcessInput(dummy.ps, 
                             required.col = c("group", "batch", "replicate"))
```


### Apply transformations

The most common normalizing transformations in microbiome analysis are total 
sum scaling (TSS) and centered log-ratio transformation (CLR). Hence, the 
`MBECS` package offers these two methods. The resulting matrices will be stored 
in their respective `slots (tss, clr)` in the `MbecData` object, while the 
original abundance table will remain unchanged.

Use `mbecTransform()` to apply total sum scaling to the data. 

```
mbec.obj <- mbecTransform(mbec.obj, method = "tss")
```

Apply centered log-ratio transformation to the data. Due to the sparse nature 
of compositional microbiome data, the parameter `offset` may be used to add a 
small offset to the abundance matrix in order to facilitate the CLR 
transformation.

```
mbec.obj <- mbecTransform(mbec.obj, method = "clr", offset = 0.0001)
```


### Preliminary report

The function `mbecReportPrelim()` will provide the user with an overview of 
experimental setup and the significance of the batch effect. To that end it is 
required to declare the covariates that are related to batch effect and group 
effect respectively. In addition it provides the option to select the abundance 
table to use here. The CLR transformed abundances are the default and the 
function will calculate them if they are not present in the input. Technically, 
the user can start the analysis at this point because the function incorporates 
the functionality of the aforementioned processing functions.

The parameter `model.vars` is a character vector with two elements. The first 
denotes the covariate column that describes the batch effect and the second one 
should be used for the presumed biological effect of interest, e.g., the group 
effect in case/control studies. The `type` parameter selects which abundance 
table is to be used <span style="color: #CD3048;">"otu"</span>, 
<span style="color: #CD3048;">"clr"</span>, <span style="color: #CD3048;">"tss"
</span>.

```
mbecReportPrelim(input.obj=mbec.obj, model.vars=c("batch","group"), 
                 type="clr")
```

### Run corrections

The package acts as a wrapper for six different batch effect correcting 
algorithms (BECA).

-   Remove Unwanted Variation 3 (`ruv3`)
-   Batch Mean Centering (`bmc`)
-   ComBat (`bat`)
-   Remove Batch Effect (`rbe`)
-   Percentile Normalization (`pn`)
-   Support Vector Decomposition (`svd`)

The function `mbecCorrection()` will apply a single correction algorithm 
selected by the parameter `method` and return an object that contains the 
resulting corrected abundance matrix in its `cor slot` with the respective name.

```
mbec.obj <- mbecCorrection(mbec.obj, model.vars=c("batch","group"), 
                           method = "bat", type = "clr")
```

The function `mbecRunCorrections()` will apply all correction algorithms 
selected by the parameter `method` and return an object that contains all 
respective corrected abundance matrices in the `cor` slot. In this example 
there will be three in total, named like the methods that created them.

```
mbec.obj <- mbecRunCorrections(mbec.obj, model.vars=c("batch","group"),
                               method=c("ruv3","rbe","bmc","pn","svd"), 
                               type = "clr")
```



### Post report

The `mbecReportPost()` function will provide the user with a comparative report 
that shows how the chosen batch effect correction algorithms changed the 
data-set compared to the initial values.

The parameter `model.vars` is a character vector with two elements. The first 
denotes the covariate column that describes the batch effect and the second one 
should be used for the presumed biological effect of interest, e.g., the group 
effect in case/control studies. The `type` parameter selects which abundance 
table is to be used <span style="color: #CD3048;">"otu"</span>, 
<span style="color: #CD3048;">"clr"</span>, <span style="color: #CD3048;">"tss"
</span>.

```
mbecReportPost(input.obj=mbec.obj, model.vars=c("batch","group"), 
               type="clr")
```

### Retrieve corrrected data

Because the `MbecData` class extends the `phyloseq` class, all functions from 
`phyloseq` can be used as well. They do however only apply to the `otu_table` 
slot and will return an object of class `phyloseq`, i.e., any transformations 
or corrections will be lost. To retrieve an object of class `phyloseq` that 
contains the `otu_table` of corrected counts, for downstream analyses, the user 
can employ the `mbecGetPhyloseq()` function. As before, the arguments `type` and
`label` are used to specify which abundance table should be used in the 
returned object.

To retrieve the CLR transformed counts, set `type` accordingly.

```
ps.clr <- mbecGetPhyloseq(mbec.obj, type="clr")
```

If the batch-mean-centering corrected counts show the best results, select 
<span style="color: #CD3048;">"cor"</span> as `type` and set the `label` to 
<span style="color: #CD3048;">"bmc"</span>.

```
ps.bmc <- mbecGetPhyloseq(mbec.obj, type="cor", label="bmc")
```









