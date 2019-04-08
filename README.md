# confoundr

INTRODUCTION

This software implements three [diagnostics for confounding/selection-bias](https://www.ncbi.nlm.nih.gov/pubmed/27479649) that can be used in sequence. Built upon the framework of sequential exchangeability, these apply to any study of multivariate exposures e.g. time-varying exposures, direct effects, interaction, and censoring. The first two diagnostics pertain to the nature of confounding/selection-bias in the data, while the third is meant to examine residual confounding/selection-bias after applying certain adjustment methods. These tools are meant to help describe confounding/selection-bias in complex data _after_ investigators have selected covariates to adjust for (e.g., through subject-matter knowledge).

+ *Diagnostic 1* is a generalization of a "Table 1" for multivariate exposures (i.e. multiple exposures that are distinct in timing or character). It examines whether the prior covariate means are the same across exposure groups, among persons who have followed a particular exposure trajectory up to that point in time. Like a "Table 1" it is meant to describe whether exposure groups have different distributions of prior covariates. 

+ *Diagnostic 2* assesses whether there is feedback between exposures and confounders over time. If present, it indicates the use of [g-methods](https://www.ncbi.nlm.nih.gov/pubmed/28039382) to control for time-varying confounding for any definition of time-varying exposure that is explicitly or implicitly multivariate (e.g., ever use at time t). The diagnostic examines whether the covariate mean differs across prior exposure groups, after adjustment for covariates (that precede the exposure) through inverse probability weighting or propensity score stratification.

+ *Diagnostic 3* is meant to be applied after investigators have applied the earlier diagnostics and have chosen to use g-methods to adjust for confounding/selection-bias. The form of Diagnostic 3 is similar to that of Diagnostic 1 in that it is a generalized "Table 1" for weighted or stratified data. It can be applied to examine residual confounding/selection-bias when [inverse probability weights](https://www.ncbi.nlm.nih.gov/pubmed/10955408) are used to fit marginal structural models. It can also be applied to examine residual confounding/selection-bias when [propensity-score stratification](https://www.ncbi.nlm.nih.gov/pubmed/19817741) is used to implement the [parametric g-formula](https://www.ncbi.nlm.nih.gov/pubmed/23533091) or [marginal mean weighting through stratification](https://www.ncbi.nlm.nih.gov/pubmed/21843003).

CAPABILITIES

The tools can accommodate:
* Multivariate exposures that are binary or categorical (and continuous, when used in concert with modeling). 
* Varying temporal depth of covariate history.
* Unbalanced, sparse data with irregular measurement of exposures/covariates or missing data
* Artificial censoring rules.
* Requests for tables/plots at all times, specific times, or averages over selected dimensions of person-time.
* Data that are not time-indexed.
* Data that are supplied in "wide" or "long" format (e.g., from the `twang` and `CBPS` packages).

To install the package, use the following code:

```
install.packages("devtools",dependencies=TRUE)
library(devtools)
install_github("jwjackson/confoundr",
               dependencies=c("Depends","Imports"), 
               build = TRUE, 
               build_opts = c("--no-resave-data","--no-manual"))
```
To load the package, use the following code:

```
library(confoundr)
```

The package now contains the documentation and two vignettes.

For a cursory example with toy data, see:

```
vignette("quickdemo")
```


For a more involved example with simulated data based on a clinical trial, see:

```
vignette("selectionbias")
```


If you wish to use the functions directly, download the file Rfunctions_1_0_2.r' in the R directory. 

A PDF manual can be found in the INST directory.

For questions please [contact me](https://www.jhsph.edu/faculty/directory/profile/3410/john-w-jackson).

John W. Jackson, ScD
Assistant Professor
Department of Epidemiology
Johns Hopkins Bloomberg School of Public Health
