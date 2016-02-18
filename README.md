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
The function BinBayes.R can be downloaed from [here](https://github.com/v2south/BinBayes/blob/master/BinBayes.R).
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
> accuracy
     subj itemID cond Acc
1     S01   i001   UD   1
2     S01   i002   RD   1
3     S01   i003   UC   1
4     S01   i004   RD   1
5     S01   i005   RD   1
6     S01   i006   UC   1
7     S01   i007   RD   1
8     S01   i008   UD   1
9     S01   i009   RC   1
10    S01   i010   UD   1
11    S01   i011   UC   1
12    S01   i012   RD   1
13    S01   i013   UD   1
.
.
.
> M4_result <- BinBayes(accuracy, "M4", "Logit")
> # BIC value
> M4_result$bic
[1] 2996.642
> # WAIC value
30
> M4_result$waic
[1] 2799.036
> # The following is part of the Post_Summary
> M4_result$post_summary[1:5,1:5]
[[1]]
          a[1]      a[2]       a[3]      a[4]        a[5]
[1,] -1.872848 0.7381711 -1.7706177 0.3784899 -0.11543930
[2,] -1.673983 0.8850012 -1.6350146 0.5198886  0.47122561
[3,] -1.555514 0.4032484 -1.2782973 0.4250902  0.68646389
[4,] -1.767328 1.3305617 -1.2789133 1.0805838  0.75447270
[5,] -1.704587 2.8370053  0.4086459 0.3240679 -0.01821947
```
Also, we can use the <strong>summary()</strong> function to get the posterior mean, posterior quantiles, and so on.
```
> summary(M4_result$post_summary)
                       Mean         SD     Naive SE Time-series SE
a[1]           -1.651674140 0.38186812 0.0027002154    0.004879634
a[2]            1.029455547 0.76992876 0.0054442184    0.007200982
a[3]           -0.447599745 0.51547315 0.0036449456    0.005530853
a[4]            0.546810962 0.67984946 0.0048072617    0.006687512
a[5]            0.136480949 0.60239446 0.0042595721    0.005966830
a[6]            0.535923970 0.66998439 0.0047375051    0.006276845
a[7]            0.519061256 0.68674737 0.0048560372    0.006519361
a[8]           -0.467438591 0.51890334 0.0036692007    0.005693154
a[9]            0.523479478 0.67846659 0.0047974833    0.006230641
a[10]           1.035360069 0.77143129 0.0054548430    0.007454912
a[11]           1.049869155 0.76432484 0.0054045928    0.007214576 
.
.
.
.
beta0           3.217671589 0.15527295 0.0010979456    0.004889324
sigma_a         1.042560073 0.10279647 0.0007268808    0.001817836
sigma_alpha_a   0.317394541 0.14507561 0.0010258394    0.020849972
sigma_b         0.443902891 0.08877508 0.0006277346    0.002598334

2. Quantiles for each variable:

                   2.5%       25%       50%       75%    97.5%
a[1]           -2.38914 -1.911754 -1.659339 -1.400410 -0.88743
a[2]           -0.34661  0.497792  0.975816  1.528980  2.64438
a[3]           -1.39199 -0.806694 -0.469248 -0.115705  0.60561
a[4]           -0.67710  0.069030  0.518771  0.983691  1.99242
a[5]           -0.95792 -0.286522  0.108670  0.525298  1.39232
a[6]           -0.66690  0.059904  0.502328  0.968129  1.95428
a[7]           -0.71599  0.040254  0.478461  0.961174  1.96379
a[8]           -1.42810 -0.827322 -0.484480 -0.131673  0.59941
a[9]           -0.69629  0.052312  0.487627  0.950006  1.95645
a[10]          -0.35858  0.494469  0.993932  1.532204  2.68696
a[11]          -0.31349  0.512905  1.005407  1.533181  2.68946
.
.
.
b[71]          -1.23225 -0.825627 -0.611229 -0.402596 -0.01648
b[72]          -0.47808 -0.069099  0.146281  0.371175  0.84990
beta0           2.91566  3.111778  3.214427  3.322229  3.52518
sigma_a         0.85684  0.970334  1.036779  1.107859  1.25747
sigma_alpha_a   0.09272  0.199065  0.304277  0.418378  0.62032
sigma_b         0.27230  0.384950  0.442368  0.501785  0.62098```
To get the 95% HPD interval, we can use the <strong>HPDinterval()</strong> functionas follows
```
> HPDinterval(M4_result$post_summary)
[[1]]
                   lower       upper
a[1]           -2.37513723 -0.876461878
a[2]           -0.42330027  2.542678120
a[3]           -1.45065339  0.535163576
a[4]           -0.71871493  1.922685809
```
To create the boxplots and density plots summarizing the posterior dis- tribution, we can first use the varnames() function in R to see all the variable names in the post summary component of the fitted model. We then ex- tract the corresponding variables to create posterior density plots and item effect boxplots for the parameters that we are interested in. For example:

```
> varnames(M4_result$post_summary)
[1] "a[1]"         "a[2]"         "a[3]"         "a[4]"         "a[5]"          
[6] "a[6]"         "a[7]"         "a[8]"         "a[9]"         "a[10]"         
.
.
.
[671] "b[67]"      "b[68]"        "b[69]"        "b[70]"        "b[71]"         
[676] "b[72]"      "beta0"        "sigma_a"      "sigma_alpha_a" "sigma_b"   
```

## Example 1
For this example, We were investigating the development of memory for visual scenes that occurs when one searches a scene for a particular object. We were specifically interested in what subjects might learn about other, non-target objects present in the scene while searching for the target object. In the first phase of the experiment, subject searched 80 scenes for a particular target object. In the test phase, they again searched the 80 scenes from the study phase as well as a set of 40 new scenes (new condition), looking for a specific target object in each case. For 40 of the scenes that had appeared in the first phase of the experiment, the target object was the same as in the first phase (studied condition), and for the other 40 scenes a new target was designated (alternate condition). In all 120 of these critical scenes, the target was present in the scene. Accuracy reflects whether or not the target was detected in the scene. Our primary interest was in whether there would be a benefit in the alternate condition, relative to the new condition, for having previously searched the scene (albeit for a different target). So the independent variable was scene type (studied, alternate, new) and the dependent variable was successful or failed detection of the target. There was also an additional set of scenes that did not contain the target that subjects were asked to search for, just to ensure that the task would be meaningful. Performance with these items was not analyzed. Dataset for this example could be downloaed from [here](https://github.com/v2south/BinBayes/blob/master/dataset/Scenes3_Bayesian.txt). 

## Example 2
For this example, We were investigating the influence of a semantic context on the identification of printed words shown either under clear (high contrast) or degraded (low contrast) conditions. The semantic context consisted of a prime word presented in advance of the target item. On critical trials, the target item was a word and on other trials the target was a nonword. The task was to classify the target on each trial as a word or a nonword (this is called a "lexical decision" task). Our interest is confined to trials with word targets. The prime word was either semantically related or unrelated to the target word (e.g., granite-STONE vs. attack-FLOWER), and the target word was presented either in clear or degraded form. Combining these two factors produced four conditions (related-clear, unrelated-clear, related-degraded, unrelated-degraded). For the current analysis, accuracy of response was the dependent measure.  Dataset for this example could be downloaed from [here](https://github.com/v2south/BinBayes/blob/master/dataset/Prime3_Bayesian.txt). 


## Future Work 
