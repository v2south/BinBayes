# BinBayes.R

BinBayes.R is the software implementaion of Bayesian approach for the mixed effects analysis of accuracy studies using mixed binomial regression models. We will illustrated how to use BinBayes.R in the following two examples, including result interpretation and plotting.

## Introduction 

We assume that the user has both the [R](https://www.cran.r-project.org/) and [JAGS](http://mcmc-jags.sourceforge.net/) software packages installed and is familiar with the basic structure and syntax of the R language. In addition, we also require following three R packages: [<strong>coda</strong>](https://cran.r-project.org/web/packages/coda/index.html),[ <strong>lme4</strong>](https://cran.r-project.org/web/packages/lme4/index.html) and 
[<strong>rjags</strong>] (https://cran.r-project.org/web/packages/rjags/index.html). 
For installation of these pacakages, please see this [manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages) from CRAN. 



The R function BinBayes.R requires three input variables as follows:

*<strong> m_data </strong>  is a matrix or dataframe containing the data from your study. This dataframe should have four columns. The first column contains the subject identifier; the second column contains the item identifier; the third column contains the identifier for experimental condition; the fourth column holds a binary valued accuracy response. These columns should be labelled as subj, itemID, cond, Acc, respectively.
