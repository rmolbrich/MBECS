MBECS - Microbiome Batch-Effect Correction Suite
================
Michael Olbrich
9/20/2021

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

# Study Summary

This section will contain all the information about this data-set.
Starting off with a literal summary, e.g., \#samples, \#samples per
group/treatment/other factor, \#covariates and \#factors, etc…

## Covariates

some sort of table to have an overview of the meta-data.

## Sample Distribution

Select two covariates of interest, e.g., ‘treatment’ and ‘batch’, and
display bars for the number of samples from one covariate grouping that
fall into the categories of the grouping.

## Sample Separation

Technically, the function ‘prcomp’ performs a singular value
decomposition to retrieve the eigenvectors for projection. Computes
covariance-matrix on counts, extract the orthogonal eigenvectors and
rank them by amount of variance they correspond to. Plot the samples on
a grid, using selected eigenvectors as axes. This shows
sample-relatedness (clustering) and can help identify the presence of
confounding factors.

# Visualization

Aka some more plots here.

## Relative Log Expression (RLE)

Separate samples by covariate of interest (CoI), e.g., treatment or
study group, calculate the median value for every feature count and
subtract it from all samples respectively. Express the remaining
variability in counts per sample in a box-plot (including outliers,
median values, whiskers, the whole shabang) and color them by batch
membership. Repeat for all factors in the CoI and use ‘facet\_grid’ to
display side-by-side.

## Heatmap

Center/Scale both features and samples. Select ‘ALL’, ‘TOPxx’ (by IQR
value) or list of features to display in a heatmap with covariates of
interest.

## Dendrogram

    ## Comes with the next update.

## BOX-plot

Select ‘ALL’ or ‘TOP’ xx features based on IQR, i.e., variability over
all samples. Produce box-plot showing expression with respect to a
covariate of interest (CoI)

# Variance Assessment

Several different approaches are used to estimate the amount of
variability attributable to covaraites of interest.

## Linear Model (LM)

This method fits a linear model to every feature respectively and
estimates the proportion of variance that the modelled covariates of
interest (coi) account for. The results are visualized in a box-plot
that shows the coi and the residual values. In case of the comparative
analysis between two correction methods or transformations, the plot
function will create panels that show the resulting boxes for the
respective transformations (count-matrices since meta is the same -
therefore same covariates and so on)

## Linear (Mixed) Model (LMM)

## Redundancy Analysis (pRDA)

pRDA / pCCA: (Legendre & Legendre (2012), Table 11.5 (p. 650)) The
Redundancy Analysis (RDA) A linear regression model is fitted to the
feature-matrix (i.e. counts) while conditioning on one COI at a time to
extract the proportion of explained variance for the variables. In this
case the result is a single number(value) for every covariate in every
transformation and thus the plot is one/or more panels of bar-plots.

Basically we take counts \~ group + Condition(batch) and subtract counts
\~ group and see how much variance batch accounts for - then repeat with
group as Condition

!! NEGATIVE COMPONENTS ARE NOT TO BE TRUSTED

Interpretation: Without a condition we can just see how much variance
(i.e. squared standard deviation in the distribution of counts) our
model (which is the constrained part) accounts for. Redundant terms
(i.e. parts of a variable that are already explained by sth. else as for
example Age and Birthdate). This is an indicator for the usefulness of
the model.

By performing pRDA with a condition the algorithm computes two models,
one with all variables and the other with the condition removed –&gt;
the difference in explained variance is attributed to the condition.
This is also why the sum of all conditions will be close to but not
exactly the total amount of explained variance!

Information concerning a number of constrained axes (RDA axes) and
unconstrained axes (PCA axes) are often presented in the results of an
RDA. - Inertia = sum(eigenvalues of all axes) –&gt; proportion for one
axis\_1 equals (eigenvalue\_1 / intertia) - The PCA axes represent the
unconstrained (i.e. residual uncharacterised factors)

## PrincipalVariance Component Analysis (PVCA)

Algorithm - calculate the correlation of the fxs count-matrix - from
there extract the eigenvectors and eigenvalues and calculate the
proportion of explained variance per eigenvector (i.e. principal
component) by dividing the eigenvalues by the sum of eigenvalues. Now
select as many PCs as required to fill a chosen quota for the total
proportion of explained variance. Iterate over all PCs and fit a linear
mixed model that contains all covariates as random effect and all unique
interactions between two covariates. Compute variance covariance
components form the resulting model –&gt; From there we get the Variance
that each covariate(variable) contributes to this particular PC. Then
just standardize variance by dividing it through the sum of variance for
that model. Scale each PCs results by the proportion this PC accounted
for in the first place. And then do it again by dividing it through the
total amount of explained variance, i.e. the cutoff to select the number
of PCs to take (but obviously not the cutoff but rather the actual
values for the selected PCs). Finally take the average over each random
variable and interaction term and display in a nice plot

## Silhouette Coefficient

Calculate principal components and get samplewise distances on the
resulting sxPC matrix. Then iterate over all the covariates and
calculate the cluster silhouette (which is basically either zero, if the
cluster contains only a single element, or it is the distance to the
closest different cluster minus the distance of the sample within its
own cluster divided (scaled) by the maximum distance). Average over each
element in a cluster for all clusters and there is the representation of
how good the clustering is.
