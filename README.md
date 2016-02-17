# BinBayes.R

BinBayes.R is the software implementaion of Bayesian approach for the mixed effects analysis of accuracy studies using mixed binomial regression models. We will illustrated how to use BinBayes.R in the following two examples, including result interpretation and plotting.

## Introduction 

We assume that the user has both the [R](https://www.cran.r-project.org/) and [JAGS](http://mcmc-jags.sourceforge.net/) software packages installed and is familiar with the basic structure and syntax of the R language. In addition, we also require following three R packages: [<strong>coda</strong>](https://cran.r-project.org/web/packages/coda/index.html),[ <strong>lme4</strong>](https://cran.r-project.org/web/packages/lme4/index.html) and 
[<strong>rjags</strong>] (https://cran.r-project.org/web/packages/rjags/index.html). 
For installation of these pacakages, please see this [manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages) from CRAN. 



The R function BinBayes.R requires three input variables as follows:

* <strong> m_data </strong> is a matrix or dataframe containing the data from your study. This dataframe should have four columns. The first column contains the subject identifier; the second column contains the item identifier; the third column contains the identifier for experimental condition; the fourth column holds a binary valued accuracy response. These columns should be labelled as *subj*, *itemID*, *cond*, *Acc* respectively.

* <strong>link</strong> is a string that specifies the link function as ”Logit” or ”Probit”.
* <strong> model </strong>  is a string taking five possible values as follows:
  * "M1", baseline model with random subject and item effects with no effect of experimental condition.
  * "M2", model with random subject and item effects and with a fixed effect for the experimental condition.
  * "M3", model with random subject and item effects and where the effect of experimental condition varies across subjects.
  * "M4", model with random subject and item effects and where the effect of experimental condition varies across items.
  * "M5", model with random subject and item effects and where the effect of experimental condition varies across both subjects and items.


As illustrated below, the output of the BinBayes.R function will consist of an object having three components with the following names:
* <strong>bic</strong> is the value of the BIC for the fitted model.
* <strong>waic</strong>  is the value of the WAIC for the fitted model.
* <strong>post_summary</strong> is an mcmc list containing samples from the posterior distribution for all components of the fitted model.

## How to use BinBayes.R
```
> # Path is the file directory where you save the BinBayes.R
> path <- "/Users/Yin/Dropbox/Bayes Factor/BinBayes.R"
> # Load BinBayes.R
> source(path)
> # Read the data into R
> accuracy <- read.table("/Users/Yin/Dropbox/Bayes Factor/Prime3
Bayesian.txt",header=TRUE, na.strings=’.’,
colClasses=c(’factor’,’factor’,’factor’,’numeric’))
> # Remove cases with missing values
> accuracy<-na.omit(accuracy)
> accuracy
```
