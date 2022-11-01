# Additive Gene by Environment (GxE) Interaction Tests Under the Trend Effect of genotypes <br> <sub>Nilotpal Sanyal, Matthieu Rochemonteix, Summer Han</sub>

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation-of-cgen"
    id="toc-installation-of-cgen">Installation of CGEN</a>
-   <a href="#import-data" id="toc-import-data">Import Data</a>
-   <a href="#additive-gxe-tests-under-the-trend-of-effect-of-genotypes"
    id="toc-additive-gxe-tests-under-the-trend-of-effect-of-genotypes">Additive
    GxE tests under the trend of effect of genotypes</a>
    -   <a
        href="#without-g-e-independence-assumption-prospective-likelihood-indepfalse"
        id="toc-without-g-e-independence-assumption-prospective-likelihood-indepfalse">1.
        Without G-E independence assumption (Prospective likelihood)
        (<code>indep=FALSE</code>)</a>
    -   <a
        href="#under-g-e-independence-assumption-retrospective-likelihood-indeptrue"
        id="toc-under-g-e-independence-assumption-retrospective-likelihood-indeptrue">2.
        Under G-E independence assumption (Retrospective likelihood)
        (<code>indep=TRUE</code>)</a>
-   <a href="#references" id="toc-references">References</a>

# Introduction

This page serves as a tutorial for the additive gene by environment
(GxE) interaction tests under the trend effect of genotypes that have
been proposed by [Rochemonteix et al. (2020)](#ref1) and [Sanyal et al. (2021)](#ref2) and
implemented in the `additive.test` function of the R package `CGEN`.
Currently, these tests are available only for binary environmental
variables.

Specifically, under the trend effect of genotypes, we illustrate the
following tests:

-   Likelihood ratio tests (LRTs) without/with gene-environment
    independence assumption ([Rochemonteix et al. (2020)](#ref1)).

-   Wald tests without/with (UML/CML) gene-environment
    independence assumption, and a robust Wald test (EB) based on an empirical
    Bayes-type shrinkage estimator that combines estimates from the
    former Wald tests (UML and CML) ([Sanyal et al. (2021)](#ref2)).

The tests that use the gene-environment independence assumption are
based on the ‘retrospective likelihood’ function and constrained maximum 
likelihood (CML) estimates where the constraint is imposed by the independence 
assumption. When that assumption is not made, the use of the retrospective 
likelihood function along with the corresponding unconstrained maximum 
likelihood estimates (UML) is equivalent to the use of traditionally used 
‘propsective likelihood’ function and the estimates obtained therefrom.

We start by loading the package.

# Installation of CGEN

The R package `CGEN` can be installed from its Bioconductor repository
using:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CGEN")

library(CGEN)
```

# Import Data

For the illustration of the tests, we use the `Xdata2` dataset which
comes inside the `CGEN` package. This dataset contains sample covariate
and outcome data that have been taken from a lung cancer study. Let us
load the data and look at its contents.

``` r
data(Xdata2, package="CGEN") 

str(Xdata2)
#'data.frame':    11449 obs. of  8 variables:
# $ case.control: int  0 0 0 0 0 0 0 1 0 0 ...
# $ SNP         : int  0 0 0 1 1 0 2 1 1 0 ...
# $ smoking     : int  1 1 1 1 1 0 0 0 1 1 ...
# $ cov1        : num  25 25 15 35 5 25 5 15 15 5 ...
# $ cov2        : int  1 1 1 1 0 0 0 0 1 1 ...
# $ cov3        : int  4 8 3 4 4 4 8 0 8 4 ...
# $ cov4        : num  3.22 3.22 2.71 3.56 1.61 ...
# $ study       : int  4 4 3 5 2 4 2 3 3 2 ...
```

So, the data contain, for 11449 subjects, case-control type 0-1
outcomes, genotype values of a SNP (0,1,2), values of a binary (0-1)
environmental variable smoking, values of four covariates and a study
variable indicating which study subjects are taken from. We test for
additive SNP x smoking interaction in this data.

# Additive GxE tests under the trend of effect of genotypes

We use the `additive.test` function to conduct additive GxE interaction
tests under the assumption of an additive genetic model, i.e., the
presence of the trend effect of genotypes. The additive genetic model is
specified by the option `genetic.model = 0`. Further, the
gene-environment independence assumption is specified by setting
`indep=TRUE`. By default, `indep=FALSE` which corresponds to the
prospective likelihood-based analysis—which we illustrate first.

## 1. Without G-E independence assumption (Prospective likelihood) (`indep=FALSE`)

The following code performs tests for additive SNP x smoking interaction
under the trend effect of genotypes without any independence assumption. 
In the output, we get

-   LRT p value and test statistic under prospective likelihood.
-   Wald test p value and RERI estimates under prospective likelihood.

``` r
test_noindep <- additive.test(data = Xdata2, response.var = "case.control",
                snp.var = "SNP", exposure.var = "smoking", main.vars = c("cov1", "cov2",
                "cov3", "cov4", "study"), strata.var = "study", op = list(genetic.model = 0))

names(test_noindep)
#[1] "additive"   "model.info"

names(test_noindep$additive)
#[1] "pval.add.LRT" "pval.add.UML" "LRT.add"      "RERI.UML"    

# => p value of the LRT
test_noindep$additive$pval.add.LRT
#[1] 0.03197455

# => test statistic of the LRT
test_noindep$additive$LRT.add
#[1] 4.599861

# => p value of the Wald (UML) test
test_noindep$additive$pval.add.UML
#[1] 0.03924117

# => RERI, standard error and gradient estimates (UML)
test_noindep$additive$RERI.UML
#$value
#[1] -0.1530518
#
#$variance
#[1] 0.005511217
#
#$gradient
#     x11.UML     x2_1.UML x11:x2_1.UML 
#   0.5597640   -0.2059929    1.5068229
```

Additionally, the estimates of the additive risk model parameters, their standard
errors and covariance matrix, estimates of the odds ratio table, the fitted maximum 
likelihood value and some other information related to the model specification and 
fit are included in the `model.info` component of the output as pointed below.

``` r
names(test_noindep$model.info)
# [1] "response.var"     "snp.var"          "exposure.var"     "main.vars"        "strata.var"      
# [6] "op"               "GxEtable"         "method"           "parms.lm.UML"     "cov.lm.UML"      
#[11] "full.loglike.UML" "DF"               "OR.table.UML"    

# => Estimates of the additive risk model parameters
test_noindep$model.info$parms.lm.UML
#$parms.lm.UML
#             Estimate  Std. Error    z value     Pr(>|z|)
#Intercept -1.63045022 0.118039297 -13.812775 2.134487e-43
#cov1       0.02247682 0.003151272   7.132618 9.847750e-13
#cov2       0.12581973 0.039654108   3.172931 1.509086e-03
#cov3      -0.05452759 0.008103140  -6.729193 1.706071e-11
#cov4       0.21096033 0.050196679   4.202675 2.637791e-05
#study      0.20136331 0.027637750   7.285807 3.197506e-13
#x2_1       0.53813867 0.060343158   8.917973 4.748992e-19
#x11       -0.05439399 0.037462923  -1.451942 1.465178e-01
#x11:x2_1  -0.07374131 0.051971696  -1.418874 1.559356e-01
#x12       -0.11191773 0.079362236  -1.410214 1.584765e-01
#x12:x2_1  -0.16321845 0.112151620  -1.455337 1.455759e-01

# => Estimates of the covariance matrix of the additive risk model parameters
test_noindep$model.info$cov.lm.UML
#              Intercept          cov1          cov2          cov3          cov4         study
#Intercept  0.0152983986  2.378421e-04 -9.331093e-04 -1.761301e-04 -5.349101e-03 -6.958840e-04
#cov1       0.0002378421  9.557136e-06  1.989700e-06 -1.005618e-07 -1.142051e-04 -3.248776e-05
#cov2      -0.0009331093  1.989700e-06  1.572155e-03  6.407017e-06  1.153846e-05 -6.320078e-06
#cov3      -0.0001761301 -1.005618e-07  6.407017e-06  6.467444e-05 -6.432421e-06 -6.995111e-06
#cov4      -0.0053491011 -1.142051e-04  1.153846e-05 -6.432421e-06  2.929254e-03 -2.564966e-04
#study     -0.0006958840 -3.248776e-05 -6.320078e-06 -6.995111e-06 -2.564966e-04  7.311637e-04
#x2_1      -0.0013610431  2.190878e-06  1.196209e-04 -6.171406e-05  5.203345e-05 -1.785638e-04
#x11       -0.0012307572 -6.155937e-07  2.037488e-05  3.100335e-06 -1.440166e-06  6.894449e-06
#x11:x2_1   0.0013213374  6.863978e-07 -4.557594e-05 -2.898525e-06 -2.603386e-05 -8.483745e-06
#                   x2_1           x11      x11:x2_1
#Intercept -1.361043e-03 -1.230757e-03  1.321337e-03
#cov1       2.190878e-06 -6.155937e-07  6.863978e-07
#cov2       1.196209e-04  2.037488e-05 -4.557594e-05
#cov3      -6.171406e-05  3.100335e-06 -2.898525e-06
#cov4       5.203345e-05 -1.440166e-06 -2.603386e-05
#study     -1.785638e-04  6.894449e-06 -8.483745e-06
#x2_1       3.668713e-03  1.202591e-03 -2.307680e-03
#x11        1.202591e-03  1.464546e-03 -1.464767e-03
#x11:x2_1  -2.307680e-03 -1.464767e-03  2.720610e-03

# => Fitted maximum likelihood value
test_noindep$model.info$full.loglike.UML
#[1] -7305.029

# => Estimates of the odds ratio table
test_noindep$model.info$OR.table.UML
#          0        1
#0 1.0000000 1.712816
#1 0.9470589 1.506823
#2 0.8941178 1.300830
```


## 2. Under G-E independence assumption (Retrospective likelihood) (`indep=TRUE`)

Under the assumption that SNP and smoking are independent (specified by
`indep=TRUE`), the following code performs tests for additive SNP x
smoking interaction under the trend effect of genotypes. From the
output, we get

-   LRT p value and test statistic under retrospective likelihood.
-   Wald test p values and RERI estimates under retrospective
    likelihood.

``` r
test_indep <- additive.test(data = Xdata2, response.var = "case.control",
              snp.var = "SNP", exposure.var = "smoking", main.vars = c("cov1", "cov2",
              "cov3", "cov4", "study"), strata.var = "study", op = list(genetic.model = 0,
              indep = TRUE))

names(test_indep)
#[1] "additive"   "model.info"

names(test_indep$additive)
#[1] "pval.add.LRT" "pval.add.UML" "pval.add.CML" "pval.add.EB"  "LRT.add"      "RERI.UML"    
#[7] "RERI.CML"     "RERI.EB"     

# => p value of the LRT
test_indep$additive$pval.add.LRT
#[1] 0.4595864

# => test statistic of the LRT
test_indep$additive$LRT.add
#[1] 0.5469018

# => p value of the Wald UML test
test_indep$additive$pval.add.UML
#[1] 0.03924117

# => p value of the Wald CML test
test_indep$additive$pval.add.CML
#[1] 0.4661798

# => p value of the Wald EB test
test_indep$additive$pval.add.EB
#[1] 0.131338

# => RERI, standard error and gradient estimates for UML
test_indep$additive$RERI.UML
#$value
#[1] -0.1530518
#
#$variance
#[1] 0.005511217
#
#$gradient
#     x11.UML     x2_1.UML x11:x2_1.UML 
#   0.5597640   -0.2059929    1.5068229 
   
# => RERI, standard error and gradient estimates for CML
test_indep$additive$RERI.CML
#$value
#[1] -0.036034
#
#$variance
#[1] 0.002445213
#
#$gradient
#     x11.CML     x2_1.CML x11:x2_1.CML 
#    0.557127    -0.119424     1.473737
    
# => RERI, standard error and gradient estimates for EB
test_indep$additive$RERI.EB
#$value
#[1] -0.1194704
#
#$gradient
#     x11.UML     x2_1.UML x11:x2_1.UML      x11.CML     x2_1.CML x11:x2_1.CML 
#  0.62820373  -0.23117872   1.69105517  -0.06811736   0.01460142  -0.18018706 
#
#$variance
#[1] 0.006269466
```

Further, the estimates of the additive risk model parameters, their standard
errors and covariance matrix under UML and CML, estimates of the odds ratio table, 
the fitted maximum likelihood value and some other information related to the model 
specification and fit are included in the `model.info` component of the output as 
pointed below.

``` r
names(test_indep$model.info)
# [1] "response.var"     "snp.var"          "exposure.var"     "main.vars"        "strata.var"      
# [6] "op"               "GxEtable"         "method"           "parms.lm.UML"     "parms.lm.CML"    
#[11] "cov.lm.UML"       "cov.lm.CML"       "full.loglike.UML" "full.loglike.CML" "DF"              
#[16] "OR.table.CML"

# => UML estimates of the additive risk model parameters
test_indep$model.info$parms.lm.UML
#             Estimate  Std. Error    z value     Pr(>|z|)
#Intercept -1.63045022 0.118039297 -13.812775 2.134487e-43
#cov1       0.02247682 0.003151272   7.132618 9.847750e-13
#cov2       0.12581973 0.039654108   3.172931 1.509086e-03
#cov3      -0.05452759 0.008103140  -6.729193 1.706071e-11
#cov4       0.21096033 0.050196679   4.202675 2.637791e-05
#study      0.20136331 0.027637750   7.285807 3.197506e-13
#x2_1       0.53813867 0.060343158   8.917973 4.748992e-19
#x11       -0.05439399 0.037462923  -1.451942 1.465178e-01
#x11:x2_1  -0.07374131 0.051971696  -1.418874 1.559356e-01
#x12       -0.11191773 0.079362236  -1.410214 1.584765e-01
#x12:x2_1  -0.16321845 0.112151620  -1.455337 1.455759e-01

# => UML estimates of the covariance matrix of the additive risk model parameters
test_indep$model.info$cov.lm.UML
#              Intercept          cov1          cov2          cov3          cov4         study
#Intercept  0.0152983986  2.378421e-04 -9.331093e-04 -1.761301e-04 -5.349101e-03 -6.958840e-04
#cov1       0.0002378421  9.557136e-06  1.989700e-06 -1.005618e-07 -1.142051e-04 -3.248776e-05
#cov2      -0.0009331093  1.989700e-06  1.572155e-03  6.407017e-06  1.153846e-05 -6.320078e-06
#cov3      -0.0001761301 -1.005618e-07  6.407017e-06  6.467444e-05 -6.432421e-06 -6.995111e-06
#cov4      -0.0053491011 -1.142051e-04  1.153846e-05 -6.432421e-06  2.929254e-03 -2.564966e-04
#study     -0.0006958840 -3.248776e-05 -6.320078e-06 -6.995111e-06 -2.564966e-04  7.311637e-04
#x2_1      -0.0013610431  2.190878e-06  1.196209e-04 -6.171406e-05  5.203345e-05 -1.785638e-04
#x11       -0.0012307572 -6.155937e-07  2.037488e-05  3.100335e-06 -1.440166e-06  6.894449e-06
#x11:x2_1   0.0013213374  6.863978e-07 -4.557594e-05 -2.898525e-06 -2.603386e-05 -8.483745e-06
#                   x2_1           x11      x11:x2_1
#Intercept -1.361043e-03 -1.230757e-03  1.321337e-03
#cov1       2.190878e-06 -6.155937e-07  6.863978e-07
#cov2       1.196209e-04  2.037488e-05 -4.557594e-05
#cov3      -6.171406e-05  3.100335e-06 -2.898525e-06
#cov4       5.203345e-05 -1.440166e-06 -2.603386e-05
#study     -1.785638e-04  6.894449e-06 -8.483745e-06
#x2_1       3.668713e-03  1.202591e-03 -2.307680e-03
#x11        1.202591e-03  1.464546e-03 -1.464767e-03
#x11:x2_1  -2.307680e-03 -1.464767e-03  2.720610e-03

# => CML estimates of the additive risk model parameters
test_indep$model.info$parms.lm.CML
#                  Estimate  Std. Error     z value     Pr(>|z|)
#Intercept     -1.600452023 0.116190346 -13.7743976 3.634018e-43
#cov1           0.022386415 0.003142912   7.1228256 1.057366e-12
#cov2           0.124318932 0.039614500   3.1382179 1.699785e-03
#cov3          -0.054459013 0.008099091  -6.7240895 1.766938e-11
#cov4           0.210126279 0.050155451   4.1895003 2.795694e-05
#study          0.202221660 0.027600193   7.3268206 2.356770e-13
#x2_1           0.465720112 0.051206243   9.0949869 9.459789e-20
#x11           -0.087073182 0.032456121  -2.6827970 7.300930e-03
#x11:x2_1       0.009154455 0.036774022   0.2489381 8.034087e-01
#Allele_freq.1 -0.295012198 0.035331048  -8.3499418 6.829661e-17
#Allele_freq.2 -0.396743240 0.033161918 -11.9638209 5.497296e-33
#Allele_freq.3 -0.333456711 0.025535852 -13.0583744 5.693046e-39
#Allele_freq.4 -0.322738880 0.035995826  -8.9660085 3.074542e-19
#Allele_freq.5 -0.278420402 0.051568950  -5.3989930 6.701599e-08
#x12           -0.182457534 0.071408762  -2.5551141 1.061530e-02
#x12:x2_1       0.020031786 0.076765420   0.2609480 7.941326e-01

# => CML estimates of the covariance matrix of the additive risk model parameters
test_indep$model.info$cov.lm.CML
#                  Intercept          cov1          cov2          cov3          cov4         study
#Intercept      0.0149183836  2.363420e-04 -9.125414e-04 -1.731255e-04 -5.334602e-03 -6.826112e-04
#cov1           0.0002363420  9.541184e-06  1.974250e-06 -8.948280e-08 -1.141744e-04 -3.236240e-05
#cov2          -0.0009125414  1.974250e-06  1.569414e-03  6.289523e-06  1.105356e-05 -6.401576e-06
#cov3          -0.0001731255 -8.948280e-08  6.289523e-06  6.457902e-05 -6.436362e-06 -7.070684e-06
#cov4          -0.0053346020 -1.141744e-04  1.105356e-05 -6.436362e-06  2.928927e-03 -2.570816e-04
#study         -0.0006826112 -3.236240e-05 -6.401576e-06 -7.070684e-06 -2.570816e-04  7.303818e-04
#x2_1          -0.0007571200  3.053459e-06  8.075863e-05 -6.400035e-05  2.419504e-05 -1.825082e-04
#x11           -0.0008418266  9.247539e-07  8.709930e-09 -3.636213e-08 -1.515403e-05 -6.750915e-06
#x11:x2_1       0.0005868225 -3.481024e-07 -6.797736e-08  1.968512e-08  4.930166e-06 -3.399926e-06
#Allele_freq.1  0.0002239115 -9.120136e-07  7.422754e-08  1.249674e-07  1.562366e-05 -7.342697e-06
#Allele_freq.2  0.0002918164  8.344021e-07  2.731239e-08  3.441251e-08 -1.378131e-05 -3.833865e-06
#Allele_freq.3  0.0002535507 -6.584542e-07 -1.184019e-07 -1.608495e-07  1.067051e-05  6.117895e-06
#Allele_freq.4  0.0002856738 -6.334723e-07  2.534818e-07  2.813931e-07  1.127990e-05  1.337965e-05
#Allele_freq.5  0.0004026789 -2.774734e-07 -4.649136e-08 -6.133013e-08  4.875884e-06  1.892906e-05
#                       x2_1           x11      x11:x2_1 Allele_freq.1 Allele_freq.2 Allele_freq.3
#Intercept     -7.571200e-04 -8.418266e-04  5.868225e-04  2.239115e-04  2.918164e-04  2.535507e-04
#cov1           3.053459e-06  9.247539e-07 -3.481024e-07 -9.120136e-07  8.344021e-07 -6.584542e-07
#cov2           8.075863e-05  8.709930e-09 -6.797736e-08  7.422754e-08  2.731239e-08 -1.184019e-07
#cov3          -6.400035e-05 -3.636213e-08  1.968512e-08  1.249674e-07  3.441251e-08 -1.608495e-07
#cov4           2.419504e-05 -1.515403e-05  4.930166e-06  1.562366e-05 -1.378131e-05  1.067051e-05
#study         -1.825082e-04 -6.750915e-06 -3.399926e-06 -7.342697e-06 -3.833865e-06  6.117895e-06
#x2_1           2.611665e-03  5.905199e-04 -1.099144e-03 -1.306548e-04  3.069330e-07  7.026481e-05
#x11            5.905199e-04  1.056882e-03 -7.109619e-04 -2.621026e-04 -2.957803e-04 -3.300160e-04
#x11:x2_1      -1.099144e-03 -7.109619e-04  1.335450e-03  1.491739e-04 -1.038123e-06 -8.665475e-05
#Allele_freq.1 -1.306548e-04 -2.621026e-04  1.491739e-04  1.240445e-03  7.988810e-05  9.174099e-05
#Allele_freq.2  3.069330e-07 -2.957803e-04 -1.038123e-06  7.988810e-05  1.114959e-03  1.644835e-04
#Allele_freq.3  7.026481e-05 -3.300160e-04 -8.665475e-05  9.174099e-05  1.644835e-04  6.825037e-04
#Allele_freq.4  9.547743e-05 -3.953258e-04 -1.222033e-04  1.102479e-04  2.011962e-04  2.633754e-04
#Allele_freq.5 -2.032411e-05 -5.361886e-04  9.522093e-06  1.438327e-04  2.319524e-04  2.940179e-04
#              Allele_freq.4 Allele_freq.5
#Intercept      2.856738e-04  4.026789e-04
#cov1          -6.334723e-07 -2.774734e-07
#cov2           2.534818e-07 -4.649136e-08
#cov3           2.813931e-07 -6.133013e-08
#cov4           1.127990e-05  4.875884e-06
#study          1.337965e-05  1.892906e-05
#x2_1           9.547743e-05 -2.032411e-05
#x11           -3.953258e-04 -5.361886e-04
#x11:x2_1      -1.222033e-04  9.522093e-06
#Allele_freq.1  1.102479e-04  1.438327e-04
#Allele_freq.2  2.011962e-04  2.319524e-04
#Allele_freq.3  2.633754e-04  2.940179e-04
#Allele_freq.4  1.251823e-03  3.599500e-04
#Allele_freq.5  3.599500e-04  2.673354e-03

# => The fitted maximum likelihood value under UML
test_indep$model.info$full.loglike.UML
#[1] -7305.029

# => The fitted maximum likelihood value under CML
test_indep$model.info$full.loglike.CML
#[1] -18923.57

# => Estimates of the odds ratio table under CML
test_indep$model.info
#$OR.table.CML
#        0        1
#0 1.00000 1.593161
#1 0.91661 1.473737
#2 0.83322 1.354313

```


# References

<div id="ref1">
De Rochemonteix, M., Napolioni, V., Sanyal, N., Belloy, M.E., Caporaso, N.E., Landi, M.T., Greicius, M.D., Chatterjee, N. and Han, S.S., 2021. A likelihood ratio test for gene-environment interaction based on the trend effect of genotype under an additive risk model using the gene-environment independence assumption. American journal of epidemiology, 190(1), pp.129-141. https://doi.org/10.1093/aje/kwaa132. <br><br>

<div id="ref2">
Sanyal, N., Napolioni, V., de Rochemonteix, M., Belloy, M.E., Caporaso, N.E., Landi, M.T., Greicius, M.D., Chatterjee, N. and Han, S.S., 2021. A Robust Test for Additive Gene-Environment Interaction Under the Trend Effect of Genotype Using an Empirical Bayes-Type Shrinkage Estimator. American journal of epidemiology, 190(9), pp.1948-1960. https://doi.org/10.1093/aje/kwab124.
