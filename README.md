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
```
To get the 95% HPD interval, we can use the <strong>HPDinterval()</strong> functionas follows
```
> HPDinterval(M4_result$post_summary)
[[1]]
lower upper
a[1]           -2.37513723 -0.876461878
a[2]           -0.42330027  2.542678120
a[3]           -1.45065339  0.535163576
a[4]           -0.71871493  1.922685809
```
To create the boxplots and density plots summarizing the posterior dis- tribution, we can first use the varnames() function in R to see all the variable names in the post summary component of the fitted model. We then ex- tract the corresponding variables to create posterior density plots and item effect boxplots for the parameters that we are interested in. For example:

```
> varnames(M4_result$post_summary)
[1] "a[1]"
[6] "a[6]"
```

## Example 1
For this example, We were investigating the development of memory for visual scenes that occurs when one searches a scene for a particular object. We were specifically interested in what subjects might learn about other, non-target objects present in the scene while searching for the target object. In the first phase of the experiment, subject searched 80 scenes for a particular target object. In the test phase, they again searched the 80 scenes from the study phase as well as a set of 40 new scenes (new condition), looking for a specific target object in each case. For 40 of the scenes that had appeared in the first phase of the experiment, the target object was the same as in the first phase (studied condition), and for the other 40 scenes a new target was designated (alternate condition). In all 120 of these critical scenes, the target was present in the scene. Accuracy reflects whether or not the target was detected in the scene. Our primary interest was in whether there would be a benefit in the alternate condition, relative to the new condition, for having previously searched the scene (albeit for a different target). So the independent variable was scene type (studied, alternate, new) and the dependent variable was successful or failed detection of the target. There was also an additional set of scenes that did not contain the target that subjects were asked to search for, just to ensure that the task would be meaningful. Performance with these items was not analyzed. Dataset for this example could be downloaed from [here](https://github.com/v2south/BinBayes/blob/master/dataset/Scenes3_Bayesian.txt). 

## Example 2
For this example, We were investigating the influence of a semantic context on the identification of printed words shown either under clear (high contrast) or degraded (low contrast) conditions. The semantic context consisted of a prime word presented in advance of the target item. On critical trials, the target item was a word and on other trials the target was a nonword. The task was to classify the target on each trial as a word or a nonword (this is called a "lexical decision" task). Our interest is confined to trials with word targets. The prime word was either semantically related or unrelated to the target word (e.g., granite-STONE vs. attack-FLOWER), and the target word was presented either in clear or degraded form. Combining these two factors produced four conditions (related-clear, unrelated-clear, related-degraded, unrelated-degraded). For the current analysis, accuracy of response was the dependent measure.  Dataset for this example could be downloaed from [here](https://github.com/v2south/BinBayes/blob/master/dataset/Prime3_Bayesian.txt). 


## Future Work 
