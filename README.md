---
title: 'mixComp: An R package for estimating complexity of a mixture'
tags:
  - R
  - mixture distribution
  - complexity
authors:
  - name: Anja Weigel^[Co-first author] # note this makes a footnote saying 'Co-first author'
    affiliation: 2 # (Multiple affiliations must be quoted)
  - name: Fadoua Balabdaoui^[Co-first author] 
    affiliation: 1
  - name: Yulia Kulagina^[Corresponding author] 
    affiliation: 1
  - name: Lilian Mueller^[Contributor]
    affiliation: 2
  - name: Martin Maechler^[Contributor] 
    affiliation: 1
affiliations:
 - name: ETH Zurich, Seminar for Statistics, Switzerland
   index: 1
 - name: ETH Zurich, Switzerland
   index: 2
date: 19 April 2022
bibliography: refs.bib

---

# Summary

The **mixComp** package provides a number of methods for obtaining a consistent estimate of the complexity of a finite mixture (the focus is made on the univariate case). The considered approaches can be loosely grouped into three categories:
<ul>
  <li> methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; </li>
  <li> methods based on penalized minimum distance between the unknown probability density and a consistent estimator thereof. The distances considered in this survey are the Hellinger and the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973841.jpg">-distances; </li>
  <li> likelihood ratio test (LRT) - based techniques. </li>
</ul>
While not the primary goal, most methods simultaneously estimate the component weights and parameters. In this document, we give a brief overview of the methodology, and demonstrate the package's functionality in both real world examples and synthetically generated data. Moreover, we show how the package can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution.

# Statement of need

Mixture models occur in numerous settings including random
and fixed effects models, clustering, deconvolution, empirical Bayes problems and many others. In particular, they are often used to model the data
which are believed to come from a heterogeneous population, consisting of
several homogeneous subpopulations. In such cases the problem of finding
a good estimator for the number of components in the mixture arises naturally, 
whether of primary interest in itself or as a first step in estimating the
population distribution. Estimation of the order of a finite mixture model
is known to be a hard statistical task, and multiple techniques have been
suggested for solving it.

The **mixComp** package provides a variety of approaches for estimation of the complexity of
a finite mixture distribution, simultaneously estimating the component weights and parameters. 
It is applicable to parametric mixtures well beyond those whose component distributions are included in the **R** package **stats**, making it more customizable than most packages for model-based clustering. 
The estimation results can be printed out and plotted for further analysis and study.
The documentation contains multiple examples based on simulated as well as real data, several real-world 
datasets have been built-in to the package for convenience.

The use of the **mixComp** package might be of interest to researchers and practitioners who 
are studying phenomena that can be effectively modelled using mixture distributions. 
Among other things it can be used to identify settings and conditions, under which a certain 
method provides more accurate estimates than the others.

# Installation

To install from `CRAN`, use:
```
install.packages("mixComp")
```
# Section 1. Introduction to finite mixture models and mixComp

Mixture models have been used extensively in statistical applications and therefore have attracted a lot of attention from both theoretical and computational perspectives. Although the list of works on mixture models is too long to make an exhaustive inventory, we can refer to the following important papers and books: [[34]](#34), [[20]](#20),  [[21]](#21), [[36]](#36) and [[24]](#24).

The popularity of such models stems from the fact that they allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. To account for such heterogeneity, the (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method. This includes the well-known nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were well studied in the seminal work of [[20]](#20),  [[21]](#21). One remarkable property of the NPMLE of the mixing distribution is that it is, under some simple conditions, a discrete distribution function with at most <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973472.jpg"> number of jumps, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973472.jpg"> is the number of distinct observations in the random sample drawn from the mixture. This interesting feature is one reason, among others, for considering the smaller class of finite mixture models, i.e., mixture models with a discrete mixing distribution with a finite number of jumps. The model has the following simple interpretation: the population under study is assumed to consist of several homogeneous subpopulations. These subpopulations, typically referred to as the mixture's components, usually have a specific meaning depending on the problem at hand. In some very simple situations, the number of components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric maximum likelihood estimator (MLE) is known to occur at the fast rate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973529.jpg"> (under some regularity conditions). Also, the well-known expectation-maximization (EM) algorithm can be used to find the MLE of all the unknown parameters; see for example [[11]](#11). However, in many statistical applications such knowledge is rarely available and the number of components has to be estimated from the data. Although the mixture is still finite and the distribution of each component is assumed to belong to some parametric family, the estimation framework in this case is much harder than in the fully parametric one, where the number of components is known. In this paper, the terms *order*, *complexity*  and *number of components*  will  be used interchangeably to refer to this unknown number. The main goal of the package **mixComp** is to estimate the unknown complexity using several methods known from the statistical literature. These methods, which are discussed below in more detail, all come with theoretical guarantees for consistency as the sample size gets larger. Of course, consistency in this case means that an estimator is able to exactly recover the unknown complexity for large sample sizes. As expected, the performance of the methods varies according to the underlying mixture distribution and the sample size. This will be illustrated below through several synthetic as well as real datasets.

To describe the estimation problem, we start with some formal notation.  A distribution <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> is called a *finite mixture* if its density (we write density throughout and keep in mind that it may be taken with respect to the Lebesgue or the counting measure) is of the form

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973590.jpg">

where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973785.jpg"> is the mixture complexity, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973708.jpg"> are the component weights and the density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973751.jpg"> is the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975606.jpg">-th component of the mixture. As the scope of **mixComp** is limited to mixtures where the family of the component distributions is known, we replace <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973751.jpg"> by a parametric density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976027.jpg"> indexed by the (possibly multivariate, say <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg">-dimensional) parameter <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976173.jpg"> in the parameter space <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922033.jpg">.

Given some complexity <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">, the two relevant parameter spaces can therefore be defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976536.jpg">

and

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976599.jpg">

Throughout this document, it is assumed that the family of the component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976784.jpg"> is known, but the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976825.jpg">, the component weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977021.jpg"> and the mixture complexity <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973785.jpg"> are unknown, with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> being the parameter of interest. Assume now that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> is a finite mixture distribution with density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973590.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977574.jpg"> is an i.i.d. sample of size $n$ from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg">. The **mixComp** package aims to estimate the smallest such <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> on the basis of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">, either on its own or by simultaneously estimating the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977761.jpg"> and the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976173.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978071.jpg">.

In this setup, it seems natural to test for the number of components by comparing two consecutive models. Traditionally, the problem of choosing between nested models may be approached by applying the generalized likelihood ratio test and referring to the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978511.jpg"> distribution to assess significance, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978628.jpg"> is given by the number of constraints imposed on the alternative hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978742.jpg"> to arrive at the null hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978727.jpg">. However, in the context of mixture models, there are several issues hindering application of this classical theory. One of them is that there is no unique way of obtaining <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978727.jpg"> from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978742.jpg">. As an example, the two null hypotheses <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978884.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978911.jpg"> both yield the smaller model, showcasing the difficulties of applying the classical asymptotic theory of the likelihood ratio. This problem has been studied extensively in the literature and numerous alternative approaches to mixture complexity estimation have been suggested, laying the theoretical foundation for the subsequently described algorithms.

This document discusses various categories of functions found in the **mixComp** package, ranging from methods based on Hankel matrices, to techniques based upon distances between densities and likelihood ratio tests. The examples provided in the first sections all contain mixtures of 'standard' distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be performed using the functions from the **stats** package. The last section illustrates how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and valuating the density thereof.

Two main features distinguish this package from other mixture-related **R** [[29]](#29) packages: Firstly, it is focused on the estimation of the complexity rather than the component weights and parameters. While these are often estimated as a by-product, all methods contained in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture of interest. Secondly, it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** (see [[4]](#4)) and **flexmix** (see [[16]](#16), [[17]](#17), [[19]](#19)),   should both be mentioned at this point: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm. Notably, it also contains routines for selecting the number of components based on information criteria and parametric bootstrapping of the likelihood ratio test statistic values. However, they are limited to multinomial and (a variety of) normal mixtures as well as mixtures-of-regressions. Second, while **flexmix** was developed to deal with mixtures-of-regression, it sets itself apart from other packages by its extensibility, a design principle that we also aimed for when creating the  **mixComp** package. Other widely used packages dealing with mixture models are **mclust** [[31]](#31), which fits mixtures of Gaussians using the EM algorithm, **MixSim** [[25]](#25), which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [[22]](#22), which is used for grouped conditional data. Interested readers can find a comprehensive list of mixture-related packages on the CRAN Task View: Cluster Analysis and Finite Mixture Models website.

Before moving to the description of the different methods implemented in **mixComp** we would like to briefly mention pther theoretical work on the estimation of mixture complexity not currently included in the package. [[6]](#6) propose a method that is reminiscent of the ones described in Section 4. The main difference is that the authors consider distribution functions instead densities, i.e. they consider minimizing a penalized distance between the distribution function of the mixture and the empirical distribution function. The approach of [[13]](#13) is based on a minimum message length-like criterion, however, their method struggles to deal with mixtures with very different weights. [[38]](#38) propose a procedure based on alternating between splitting and merging the components in an EM-algorithm. This algorithm requires selecting two thresholds, the choice of which is somewhat unclear when working with a specific dataset. [[26]](#26) follow a Bayesian approach, taking the usual finite mixture model with Dirichlet weights and putting a prior distribution on the unknown number of components. 

The methods that were included in the package can be roughly devided into three categories: methods based on Hankel matrices, following the theory as described in [@hankel] and selected because of the fact that computation of $\textbf{w}$ and $\textbf{\theta}$ is not required, a method based on the likelihood ratio test statistic (LRTS) following [@lrt] since a likelihood ratio test seems like a natural approach in this setting and methods employing minimum distance calculations based on several works and included as a computationally more efficient alternative to the LRTS method for certain distributions and distances; see [@hell; @hellcont; @l2; @adap]. For example, when the distance is taken to be the Hellinger distance, such an approach is especially fast for discrete distributions. For a more fluid reading, the relevant theory will be portrayed at the beginning of each of the respective sections. The examples depicted in these first chapters all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be done by functions available from the **stats** package. The last chapter showcases how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and evaluating the density thereof.


# Section 2. Objects and functions defined in mixComp

Table 1 depicts five object classes defined in **mixComp**. The first two respectively represent a finite mixture distribution and a random sample drawn from such a distribution. The `Mix` object is printed as a matrix of component weights and parameters and is plotted as the density of the specified mixture, showing the overall as well as the component densities. The `rMix` object prints as the vector of observations and plots as a histogram, showcasing the individual components as well as the full sample. Both objects contain a number of attributes giving additional information, details of which can be found in the corresponding **R** help files. 

#### Table 1: Objects and functions defined in mixComp
|  Object class  | Created via                                  | Description                                     |
|:--------------:|:--------------------------------------------:|:-----------------------------------------------:|
| `Mix`          | `Mix`                                        | Represents a finite mixture                     |
| `rMix`         | `rMix`                                       | Randomly generated data from a finite mixture   |
| `datMix`       | `datMix` or `RtoDat`                         | Observed data from (presumably) a finite mixture|
| `hankDet`      | `nonparamHankel`                             | Vector of estimated Hankel matrix determinants  |
| `paramEst`     | `paramHankel(.scaled)`, `L2(.boot).disc`, `hellinger(.boot).disc`, `hellinger(.boot).cont` or `mix.lrt`  | Complexity estimate $\hat{p}$, together with estimates of the weights $\hat{\mathbf{w}}$ and the component parameters $\hat{\mathbf{\theta}}$>|

The generation of an object of class `Mix` hinges on four central arguments: a string `dist` specifying the name of the family of component densities (or kernels) $\{g(x;\theta):\theta \in \Theta \}$, a boolean`discrete` stating whether the distribution is discrete, a vector `w` giving the weights $w_i$ and a list `theta.list` (the component parameters can also be supplied via the `...` argument) containing the parameters of the component densities $\theta_i, i \in 1, \dots, p$. While the creation of `Mix` objects is mostly straightforward, two things should be noted in this regard: First, **mixComp** procedures will search for functions called `rdist` and `ddist` in the accessible namespaces. For most "standard" distributions, these functions are contained in the **stats** package and do not need to be user-written (compare with the Section 6). To make use of these functions, it is essential that the string `dist` is named correctly (e.g. to create a gaussian mixture on the basis of the **stats** package, `dist` has to be specified as `norm` instead of `normal`, `gaussian` etc. for the package to find the functions `dnorm` and `rnorm`). Second, the names of the list elements of `theta.list`(for the names of the `...` arguments) have to match the names of the formal arguments of the functions `ddist` and `rdist` exactly (e.g. for a gaussian mixture, the list elements have to be named `mean` and `sd`, as these are the formal arguments used by `rnorm` and `dnorm` functions of the **stats** package).

The following example creates two `Mix` objects, a 3-component mixture of normal distributions and a 3-component mixture of Poisson distributions. 

```{r mixobj}
set.seed(0)
# construct a Nix object:
normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
# plot the mixtures:
plot(normLocMix, main = "3-component normal mixture", cex.main = 0.9)
plot(poisMix, main = "3-component poisson mixture", cex.main = 0.9)
```
<p float="left">
  <img src="https://github.com/yuliadm/mixComp/blob/main/images/normMix.png" />
  <img src="https://github.com/yuliadm/mixComp/blob/main/images/poisMix.png" />
</p>

If required, random samples can be generated from these mixtures.
```{r rmix}
# generate random samples:
normLocRMix <- rMix(1000, obj = normLocMix)
poisRMix <- rMix(1000, obj = poisMix)
# plot the histograms of the random samples:
plot(normLocRMix, main = "Three component normal mixture", cex.main = 0.9)
plot(poisRMix, main = "Three component poisson mixture", cex.main = 0.9)
```
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/normRMix.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/poisRMix.png" />
</p>
  
The third object class shown in Table 1, called `datMix`, represents the data vector $\mathbf{X}$ based on which the mixture complexity is supposed to be estimated. These objects are most central to the package, as every procedure estimating the order of a mixture takes a `datMix` object as input. Apart from $\mathbf{X}$, it contains other "static" information needed for the estimation procedure (in contrast to "tuning parameters", which can be changed with every function call. An example of such a tuning parameter is the number of bootstrap replicates for a function employing a bootstrap procedure). A brief overview of which "static" attributes need to be supplied for each complexity estimation routine is given in Table 2. 

#### Table 2: Inputs needed for different functions contained in mixComp
| `R function`            |`dist`|`discrete`|`theta.bound.list`|  `MLE.function` | `Hankel.method` | `Hankel.function` |
|:-----------------------:|:----:|:--------:|:----------------:|:---------------:|:---------------:|:-----------------:|
|`nonparamHankel`         |  x   |          |                  |                 |        x        |         x         |
|`nonparamHankel(.scaled)`|  x   |    x     |         x        | o + i (optional)|        x        |         x         |
|`L2(.boot).disc`         |  x   |    x     |         x        |   i (optional)  |                 |                   |
|`hellinger(.boot).disc`  |  x   |    x     |         x        |   i (optional)  |                 |                   |
|`hellinger(.boot).cont`  |  x   |    x     |         x        |   i (optional)  |                 |                   |
|`mix.lrt`                |  x   |    x     |         x        | o + i (optional)|                 |                   |

In the table, "o" and "i" respectively stand for input used during optimization and initialization.

As a simple example of a given dataset to which mixture models have been applied extensively, take the Old Faithful dataset [@R; @faithful1; @faithful2]. In the context of mixture model estimation, the variable `waiting`, which gives the time in minutes between eruptions of the Old Faithful geyser in the Yellowstone National Park, is often considered to be the variable of interest. To estimate the number of components of the mixture distribution that provides a suitable approximation to the `waiting` data via **mixComp**, the raw data vector of observations has to be converted to a `datMix` object first. For the sake of exposition we specify all arguments of  the `datMix` function, starting with the vector of observations $\mathbf{X}$ and the string `dist`, specifying $\{g(x;\theta):\theta \in \Theta \}$ and the boolean `discrete`. As has often been done in the relevant literature, we assume that the data comes from a normal mixture.

```{r faithopts}
faithful.obs <- faithful$waiting
norm.dist <- "norm"
norm.discrete <- FALSE
```
Second, a named list of length $d$ containing the bounds of $\theta \in \Theta \subseteq \mathbf{R}^d$ has to be created. In this example, $\theta = \{\mu, \sigma\} \in \Theta = \mathbf{R} \times (0, \infty) \subseteq \mathbf{R}^2$.

```{r normlist}
# define the range for parameter values:
norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
```
Next, the argument `MLE.function` contains a single function if $d=1$ or a list of functions of length $d$ otherwise, specifying the $d$ functions needed to estimate the MLEs of $\theta$ based on $\textbf{X}$ if $p$ were equal to 1 (i.e., the MLEs of the component distribution). If this argument is supplied and the `datMix` object is handed over to a complexity estimation procedure relying on optimizing over a likelihood function, the `MLE.function` attribute will be used for the single component case. In case the objective function is either not a likelihood or corresponds to a mixture with more than 1 component, numerical optimization will be used based on **Rsolnp**'s function `solnp` [@solnp, @Rsolnp]. The initial values (for the parameters of a $j$-component mixture, say) supplied to the solver are then calculated as follows: the data $\textbf{X}$ is clustered into $j$ groups by the function `clara` (of the **cluster** package by [@cluster] and the data corresponding to each group is given to `MLE.function`. The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates. Specifying `MLE.function` is optional and if it is not, for example because the MLE solution does not exists in closed form, numerical optimization is used to find the relevant MLE's.

Presuming a normal mixture, one specifies $2$ functions, namely the MLE of the mean $\hat{\mu}_{MLE} = \frac{1}{n}\sum_{i = 1}^n X_i$ and the MLE of the standard deviation $\hat{\sigma}_{MLE} = \sqrt{\frac{1}{n}\sum_{i = 1}^n (X_i - \hat{\mu}_{MLE})^2}$.

```{r normfun}
# define the MLE functions for the mean and sd: 
MLE.norm.mean <- function(dat) mean(dat)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)
```
The last two arguments, `Hankel.method` and `Hankel.function`, need to be supplied if the mixture complexity is to be estimated based on the Hankel matrix of the moments of the mixing distribution. The reader is referred to the Section 3 for further information on how these arguments are to be specified (in this case, the simplifying assumption of unit variance is made. This would be a poor choice for the `waiting` data, so $p$ should not be estimated with one of the methods using these arguments, namely `nonparamHankel`, `paramHankel` and `paramHankel.scaled`, see Table 2). 

```{r normmom}
method <- "translation"
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}
```
Finally, all previously generated objects are combined to a `datMix` object.

```{r faithdatmix}
# construct a datMix object that summarizes all the necessary information:
faithful.dM <- datMix(faithful.obs, dist = norm.dist, discrete = norm.discrete,
                      theta.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = method,
                      Hankel.function = mom.std.norm)
```
In the preceding example, the data vector $\mathbf{X}$ was taken from an already existing dataset. As seen before, the `rMix` function can be used to generate a $n$-sized sample from a specific mixture. If this synthesized data is to be used in simulations (i.e. passed to one of the functions estimating the mixture complexity) an `rMix` object can be converted to a `datMix` object via the `RtoDat` function. Apart from `dist` and `discrete`, all `datMix` arguments have to be supplied to `RtoDat` likewise. 

Unlike the above mentioned objects whose creation precedes any type of mixture complexity estimation, objects of the bottom two classes (see Table 1) contain the results of the estimation procedures. Generally, the functions estimating the number of components differ in the types of families of component densities $\{g(x;\theta):\theta \in \Theta \}$ for which they allow and in whether they provide estimates of the weights $\hat{w}_i$ and the component parameters $\hat{\theta}_i, i \in 1, \dots, \hat{p}$, the latter determining the object class of the estimation result. These differences are shown in Table 3. The function `nonparamHankel` returns an object of class `hankDet`, which is a vector of determinants (scaled and/or penalized), each entry belonging to a certain complexity estimate. The link between these determinants and $p$ will be discussed in the Section 3. `paramEst` objects arise when using any other function estimating the mixture complexity, all of which additionally return estimates of the component weights and parameters. For both object classes, print and plot methods are available to summarize and visualize the estimation results. 

#### Table 3: Distribution restrictions and output types of different functions contained in mixComp
|     R function                    |Restrictions on the family of the component density |Estimation of $\mathbf{w}$ and $\mathbf{\theta}$|
|:---------------------------------:|:--------------------------------------------------:|:----------------------------------------------:|
|`nonparamHankel`                   |Compatible with `explicit`, `translation` or `scale`|                                                |
|`nonparamHankel(.scaled)`          |Compatible with `explicit`, `translation` or `scale`|                       x                        |  
|`L2(.boot).disc`                   |Discrete distributions                              |                       x                        | 
|`hellinger(.boot).disc`            |Discrete distributions                              |                       x                        | 
|`hellinger(.boot).cont`            |Continuous distributions                            |                       x                        | 
|`mix.lrt`                          |                                                    |                       x                        | 


# Section 3. Functions using Hankel matrices

$p$ is characterized by the smallest integer $j$ such that the determinant of the $(j+1) \times (j+1)$ Hankel matrix of the first $2j$ moments of the mixing distribution equals zero. Moreover, it can be shown that this determinant is zero for all $j \geq p$. Formally, for any vector $\mathbf{c} = (c_0, \dots, c_{2k}) \in \mathbb{R}^{2k+1}$ with $c_0 = 1$, the Hankel matrix of $\mathbf{c}$ is defined as the $(k+1)\times(k+1)$ matrix given by
$$H(\mathbf{c})_{i,j} = c_{i+j-2}, \quad \quad 1 \leq i,j \leq k+1.$$
Now, let $\textbf{c}^{2j+1} \in \mathbb{R}^{2j+1}$ be the vector containing the first $2j+1$ (raw) moments of the mixing distribution. For finite mixture models, this amounts to

$$c^{2j+1}_k = \sum_{i=1}^p w_i \theta^k_i, \quad \text{ for } k \in \{0,\dots, 2j\} \text{ and } \theta_i \in \mathbb{R}, i \in \{1,\dots,p\}.$$

Then, for all $j \geq 1$, $H(\textbf{c}^{2j+1})$ is non-negative and

$$ p = \min\{j \geq 1 : \det H(\textbf{c}^{2j+1}) = 0\}. $$

Making use of this fact, the first approach to estimating the order of a mixture that is implemented in **mixComp** relies on initially finding a consistent estimator of $\textbf{c}^{2j+1}$ based on $\textbf{X}$, say $\hat{\textbf{c}}^{2j+1}$, to then iteratively calculate the applicable Hankel matrix while increasing the assumed order $j$ until a sufficiently small value of $\det H(\hat{\textbf{c}}^{2j+1})$ is attained. However, since $\det H(\hat{\textbf{c}}^{2j+1})$ should be close to 0 for all $j \geq p$, this would lead to choosing $\hat{p}$ rather larger than the true value and it seems natural to introduce a penalty term. Therefore [@hankel] define the empirical penalized objective function as

$$J_n(j) \coloneqq \lvert \det H(\hat{\textbf{c}}^{2j+1}) \rvert + A(j)l(n),$$

with $l(n)$ being a positive function converging to $0$ as $n\to\infty$ and $A(j)$ being positive and strictly increasing. 
$$\hat{p} \coloneqq \argmin_{j \in \mathbb{N}} J_n(j)$$
is then a consistent estimator of $p$.

As an extension to simply adding a penalty term to the determinant, a scaling approach was considered by \citet{lilian}. Let $\hat{d}_j = \det H(\hat{\textbf{c}}^{2j+1})$, $d_j = \det H(\textbf{c}^{2j+1})$ and $j_m \geq p, j_m \in \mathbb{N}$. Since the estimated moments $\hat{\textbf{c}}^{2j+1}$ are asymptotically normal, one can apply the delta method giving

$$
\sqrt{n} \cdot
  \big(
    \hat{d}_1-d_1,
    \dots,
    \hat{d}_{p-1}-d_{p-1},
    \hat{d}_p-0,
    \dots,
    \hat{d}_{j_m}-0
  \big)^\top \quad \overset{\mathcal{D}}{\longrightarrow} \quad \mathcal{N}(0_{j_m \times 1}, \Sigma_ {j_m \times j_m}).
$$

Instead of inspecting the vector $(\hat{d}_1, \dots, \hat{d}_{j_m})$, one could therefore also base the complexity analysis on a vector of scaled determinants employing a nonparametric bootstrap procedure on $\textbf{X}$. 

To this end, let $\tilde{\Sigma} \in \mathbb{R}^{j_m \times j_m}$ denote the covariance matrix of the determinants $\hat{d}^{*b}_{j}$ calculated on the $b^{\text{th}}$ bootstrap sample for $b=1, \dots, B$ and $j = 1, \dots j_m$. Note that 

$$\tilde{\Sigma} \approx \frac{\Sigma}{n} \quad \text{ as }B \to \infty, n \to \infty$$


and write $\tilde{\Sigma}^{-1/2} = \sqrt{n} \cdot \hat{\Sigma}^{-1/2}$. Define the rescaled vector 
\begin{equation} \label{eq:scaled}
\big( y_1, \dots, y_p, \dots, y_{j_m} \big)^\top := \sqrt{n} \cdot {\hat{\Sigma}}^{-1/2} \big(
   \hat{d}_1,
    \dots,
    \hat{d}_p,
    \dots,
    \hat{d}_{j_m}
  \big)^\top.
\end{equation}

Note that in the case of the scaled version, the criterion to be minimized becomes 
$$
{J_n(j)}_{scaled} := \vert y_j \vert + A(j)l(n) \cdot \sqrt{n}.
$$
That is, the chosen penalty function should be multiplied by $\sqrt{n}$.

This approach was proposed to address the issue of determinants already being very small from the beginning (even for $j = 1$), which, in the simulations by [@lilian], made it hard to discern the "best" complexity estimate, a problem that was not reduced much by solely adding a penalty term.  

With this general framework in place, the computation now merely hinges on calculating $\hat{\textbf{c}}^{2j+1}$. The **mixComp** package offers three methods to do so. The method to use depends on the family of component densities $\{g(x;\theta):\theta \in \Theta \}$ and is linked to some function $f_j(\textbf{X})$ needed to estimate $\hat{\textbf{c}}^{2j+1}$. The calculation method and the relevant function are specified when creating the `datMix` object as arguments `Hankel.method` and `Hankel.function`.


#### 1. `Hankel.method = "explicit"`

This method can be applied when a closed form expression for estimates of the moments of the mixing distribution exists. `Hankel.function` then contains the function explicitly estimating $\textbf{c}^{2j+1}$. 

As an example, consider a mixture of geometric distributions, where it can be shown that

$$c^{2j+1}_j = 1 - \sum_{l = 0}^{j-1} f(l) = 1 - F(j-1),$$

with $F$ the true cumulative distribution function. Hence one may take

$$\hat{c}^{2j+1}_j = 1 - \hat{F}(j-1)$$

as an estimator, with $\hat{F}$ being the empirical distribution function.

```{r geommom}
# define the function for computing the moments:
explicit.geom <- function(dat, j){
  1 - ecdf(dat)(j - 1)
}
```

As a second example, consider what [@hankel, p. 283, equation (3)] called the "natural" estimator, i.e. using

\begin{equation}\label{eq:1}
\hat{c}^{2j+1}_j = f_j\left(\frac{1}{n} \sum_{i=1}^n \psi_j(X_i)\right)
\end{equation}

when

\begin{equation}\label{eq:2}
c^{2j+1}_j = f_j(\E[\psi_j(X_i)]).
\end{equation}

Note that the estimators of this form may also be supplied as `Hankel.method = "explicit"` with `Hankel.function` equal to the right-hand side of Equation \autoref{eq:1}. For example, the "natural" estimator is applicable in the case of Poisson mixtures. If $Y \sim Pois(\lambda)$, it is a well known fact that

$$\lambda^j = \E[Y(Y-1)\dots(Y-j+1)],$$

which, in combination with \autoref{eq:1} and \autoref{eq:2} suggests using

$$\hat{c}^{2j+1}_j = \frac{1}{n} \sum_{i=1}^n X_i(X_i-1)\dots(X_i-j+1)$$

as an estimator.

```{r geompois}
# define the function for computing the moments:
explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
         matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}
```

#### 2. `Hankel.method = "translation"`
 
In Example 3.1., [@hankel, p.284] describe how $\textbf{c}^{2j+1}$ can be estimated if the family of component distributions $(G_\theta)$ is given by $\text{d}G_\theta(x) = \text{d}G(x-\theta)$, where $G$ is a known probability distribution whose moments can be given explicitly. In this case, a triangular linear system can be solved for the estimated moments of the mixing distribution $\hat{\textbf{c}}^{2j+1}$ using the empirical moments of the mixture distribution and the known moments of $G$. The former can be estimated from the data vector $\textbf{X}$ whereas the latter has to be supplied by the user. Thus, `Hankel.function` contains a function of $j$ returning the $j$-th (raw) moment of $G$.

As an example, consider a mixture of normal distributions with unknown mean and unit variance. Then $G$ is the standard normal distribution, and its $j$-th moment $m_j$ is defined as

$$
m_j=
\begin{cases}
0 & \text{if } j \text{ odd},\\
(j-1)!! & \text{if } j \text{ even}.
\end{cases}
$$

```{r}
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}
```

#### 3. `Hankel.method = "scale"`

Similarly, in Example 3.2. [@hankel, p.285] describe how $\textbf{c}^{2j+1}$ can be estimated if the family of component distributions $(G_\theta)$ is given by \linebreak $\text{d}G_\theta(x) = \text{d}G(\frac{x}{\theta})$, where $G$ is a known probability distribution whose moments can be given explicitly. Likewise, a triangular linear system can be solved for $\hat{\textbf{c}}^{2j+1}$, using the empirical moments of the mixture distribution and the known moments of $G$. `Hankel.function` contains a function of $j$ returning the $j$-th moment of $G$. Note that squares have to be taken everywhere if for some integer $j$, $m_j = 0$ (compare with [@hankel, p.285]).

As an example, consider a mixture of normal distributions with zero mean and unknown variance. Then $G$ is again the standard normal distribution, and its $j$-th moment is defined as above.

Coming back to the overall goal of complexity estimation, the function `nonparamHankel` returns all estimated determinant values corresponding to complexities up to `j.max`, so that the user can pick the lowest $j$ generating a sufficiently small determinant. The function allows the inclusion of a penalty term as a function of the sample size `n` and the currently assumed complexity `j` which will be added to the determinant value (by supplying `pen.function`), and/or scaling of the determinants (by setting `scaled  = TRUE`). For scaling, a nonparametric bootstrap is used to calculate the covariance of the estimated determinants, with `B` being the size of the bootstrap sample. The inverse of the square root (i.e. the matrix $S$ such that $A = SS$, where $A$ is the (square) covariance matrix. The procedure uses **expm**'s `sqrtm` [@expm]) of this covariance matrix is then multiplied with the estimated determinant vector to get the scaled determinant vector.

We will initially apply this method to the two already generated datasets of 3-component Poisson and normal mixtures using the penalty $A(j)l(n) = \frac{j\log(n)}{\sqrt{n}}$ and scaling the determinants according to Equation \autoref{eq:scaled}.

First, for converting the previously simulated samples from 3-component Poisson and normal mixtures yielding the objects of class `rMix` to objects of class `datMix` one should apply the `RtoDat` function as follows:
```{r rtodat}
MLE.pois <- function(dat) mean(dat)

# create datMix objects:
pois.dM <- RtoDat(poisRMix, theta.bound.list = list(lambda = c(0, Inf)), 
                  MLE.function = MLE.pois, Hankel.method = "explicit",
                  Hankel.function = explicit.pois)


normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list, Hankel.method = "translation",
                     Hankel.function = mom.std.norm)
```

In the case of the scaled version of the method, the penalty should be multiplied by $\sqrt{n}$ as mentioned earlier. 
```{r nonph}
# define the penalty function:
pen <- function(j, n){
  j * log(n)
}

# apply the nonparamHankel function to the datMix objects:
set.seed(1)
poisdets_sca_pen <- nonparamHankel(pois.dM, j.max = 5, scaled = TRUE, 
                                   B = 1000, pen.function = pen)
normdets_sca_pen <- nonparamHankel(normLoc.dM, j.max = 5, scaled = TRUE, 
                                   B = 1000, pen.function = pen)
```

We can print and plot the results as suggested below.

```{r plotnonph, figures-side, fig.show="hold", out.width="50%"}
# print the results (for the Poisson mixture)
print(poisdets_sca_pen)
# plot results for both mixtures:
par(mar = c(5, 5, 1, 1))
plot(poisdets_sca_pen, main = "3-component Poisson mixture", cex.main = 0.9)
plot(normdets_sca_pen, main = "3-component Normal mixture", cex.main = 0.9)
```
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_art_1.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_art_2.png" />
</p>

Having created the data ourselves, we know that it comes from a 3-component Poisson mixture and a 3-component Gaussian mixture respectively. The resulting plots indicate that while theoretically sound, the scaled version of the Hankel method can struggle to correctly identify the number of components in practice.

As the preceding example shows, it can be quite difficult to determine the order estimate from the vector of estimated determinants alone. Thus, the package includes another option of estimating $p$ based on Hankel matrices, however, using a more "parametric" approach which goes hand in hand with estimating$\mathbf{w}$ and $\mathbf{\theta}$. The `paramHankel` procedure initially assumes the mixture to only contain a single component, setting $j = 1$, and then sequentially tests $p = j$ versus $p = j+1$ for $j = 1,2, \dots$, until the algorithm terminates. To do so, it determines the MLE for a $j$-component mixture $(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) = (\hat{w}_1, \dots, \hat{w}_j, \hat{\theta}_1, \dots, \hat{\theta}_j) \in W_j \times \Theta_j$, generates `B` parametric bootstrap samples of size $n$ from the distribution corresponding to $(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j)$ and calculates `B` determinants of the corresponding $(j+1) \times (j+1)$ Hankel matrices. The null hypothesis $H_0: p = j$ is rejected and $j$ is increased by $1$ if the determinant value based on the original data vector $\textbf{X}$ lies outside of the interval $[ql, qu]$, a range specified by the `ql` and `qu` empirical quantiles of the bootstrapped determinants. Otherwise, $j$ is returned as the order estimate $\hat{p}$, that is $\hat p$ is the first order for which the null hypothesis is not rejected.  

`paramHankel.scaled` functions similarly to `paramHankel` with the exception that the bootstrapped determinants are scaled by the empirical standard deviation of the bootstrap sample. To scale the original determinant, `B` nonparametric bootstrap samples of size $n$ are generated from the data, the corresponding determinants are calculated and their empirical standard deviation is used.

Applying `paramHankel.scaled` to the same Poisson and Normal mixtures results in the correct identification of the mixture complexity in both cases as can be seen in the plot:

```{r plotph, figures-side, fig.show="hold", out.width="50%"}
# apply papamHankel.scaled to datMix objects:
set.seed(1)
pois_sca_pen <- paramHankel.scaled(pois.dM)
norm_sca_pen <- paramHankel.scaled(normLoc.dM)
# plot the results for both mixtures:
par(mar=c(5, 5, 1, 1))
plot(pois_sca_pen,)
plot(norm_sca_pen)
```
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_art_1.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_art_2.png" />
</p>

Consider now, as a real-world example, the Children dataset whose content was taken from the Annual Report of the pension fund S.P.P. of 1952. The dataset initially appeared in work of [@thisted] and was subsequently analysed by many authors. It entails data on 4075 widows who recieved pension from the fund, with their number of children being our variable of interest. For example, there are 3062 widows without children, 587 widows with one child, etc. Many authors have noted that this data is not consistent with being a random sample from a Poisson distribution since the number of zeros found in the data is too large. Thisted approached this by fitting a mixture of two populations, one which is always zero and one which follows a Poisson distribution. **mixComp** includes this data stored as a dataframe. Here, we want to investigate 
how the Hankel matrix methods compare when fitting the data to a mixture of Poissons.

The estimation process starts with the construction of the `datMix` object.

```{r childex}
# convert the data to vector:
children.obs <- unlist(children)
# define the MLE function:
MLE.pois <- function(dat) mean(dat)
# construct a datMix object:
children.dM <- datMix(children.obs, dist = "pois", discrete = TRUE, 
                      Hankel.method = "explicit", 
                      Hankel.function = explicit.pois,
                      theta.bound.list = list(lambda = c(0, Inf)), 
                      MLE.function = MLE.pois)
```


First, we check the nonparametric method. We define the penalty $A(j)l(n)$ as $\frac{j\log(n)}{\sqrt{n}}$ and scale the determinants according to Equation \autoref{eq:scaled} (by multiplying the penalty by $\sqrt{n}$). The result suggests that the data comes from a 2-component mixture.

```{r childplotnph, fig.width = 5, fig.height = 4}
# define the penalty:
pen <- function(j, n) j * log(n)
# estimate the number of components:
set.seed(0)
(det_sca_pen <- nonparamHankel(children.dM, j.max = 5, scaled = TRUE, 
                              B = 1000, pen.function = pen))
#plot the results:
plot(det_sca_pen, main = "Non-parametric Hankel method for Children dataset",
     cex.main = 0.9)
```

Next, we check the fit of the parametric version. The printed result of `paramHankel.scaled` shows that this method also suggests 2 to be the number of components, with the first component corresponding to a Poisson distribution with$\lambda = 0.0306$. Note that the limit case $\lambda = 0$ results in a point mass at 0, and that this fit therefore nicely lines up with the idea of a component accounting for only the zero observations. The plot shows that this method yields a sensible fit overall.

```{r childplotph, fig.width = 5, fig.height = 4}
set.seed(0)
param_sca <- paramHankel.scaled(children.dM, j.max = 5, B = 1000, ql = 0.025, 
                          qu = 0.975)
plot(param_sca, breaks = 8, ylim = c(0, 0.8))
```
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_real.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_real.png" />
</p>


# Section 4. Functions using distances

Unlike the theory on Hankel matrices introduced in Section 3, many theoretical considerations rely on estimates of the weights $\mathbf{w}$ and the component parameters \mbox{\boldmath$\theta$}. As mentioned in the introduction, it is assumed that the family of component densities $\{g(x;\theta): \theta \in \Theta \}$ is known. To embed the subsequent algorithms in a more theoretical framework, consider the parametric family of mixture densities

$$\mathcal{F}_j = \{ f_{j, \mathbf{w},\bm{\theta}} : (\mathbf{w}, \mbox{\boldmath$\theta$}) \in W_j \times \Theta_j \}.$$ 


With $\{g(x;\theta): \theta \in \Theta \}$ set in advance, elements of $\mathcal{F}_j$ can be written as

$$f_{j,\mathbf{w},\bm{\theta}}(x) = \sum_{i = 1}^j w_i g(x; \theta_i).$$

Note that the support of $f$ will depend on the support of $g$ and $\mathcal{F}_j \subseteq \mathcal{F}_{j+1}$\footnote{This is obvious by setting $w_{j+1} = 0$.} for all $j$. Now take a specific mixture $f_0 = f_{p_0, \mathbf{w}_0,\bm{\theta}_0}$, where $(\mathbf{w}_0,\bm{\theta}_0) \in W_{p_0} \times \Theta_{p_0}$. Clearly, the mixture complexity is defined as
$$p_0 = \min\{j:f_0 \in \mathcal{F}_j\}.$$

The above suggests an estimation procedure based on initially finding the 'best' possible estimate (in a sense to be determined) $(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) \in W_j \times \Theta_j$ for a given value of $j$, in order to compare the thereby specified probability density/mass function 

$$\hat{f}_j(x) = f_{j, \hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j}(x),$$

with a non-parametric density/probability mass estimate $\tilde{f}_n(x)$. As the classes $\mathcal{F}_j$ and $\mathcal{F}_{j+1}$ are nested, the distance $D$ (to be defined below) between $\hat{f}_j$ and $\tilde{f}_n$ will only decrease with $j$. Thus, it makes sense to add some penalty term (increasing in $j$) to $D(\hat{f}_j, \tilde{f}_n)$ and find the first value of $j$ where the penalized distance for $j$ is smaller than that for $j+1$. Rearranging the terms gives rise to an algorithm starting at $j = 1$, involving some threshold $t(j,n)$ depending on the penalty, where, if $j$ is the first integer satisfying

\begin{equation}\label{eq:distances}
D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n) \leq t(j,n),
\end{equation}

then $j$ is taken as the estimate $\hat{p}$. If the inequality is not fulfilled, $j$ is increased by $1$ and the procedure is repeated. Consistency of estimators defined this way has been shown in a number of cases, amongst them those used by the **mixComp** algorithms, and the reader is referred to [@l2; @hell; @hellcont] for the proofs relevant to the results implemented in the package.

The preceding notation was held as broad as possible, since different distance measures $D$ and non-parametric estimators $\tilde{f}_n$ can be used. Those relevant to the package are mostly well-known, still, definitions can be found in the Appendix. Three procedures are implemented in **mixComp** based on the foregoing methodology: `L2.disc`, `hellinger.disc` and `hellinger.cont`.

#### 1. `L2.disc`

`L2.disc` employs the squared $L_2$ distance as the distance measure $D$ and is only to be used for discrete mixtures since the nonparametric estimate $\tilde{f}_n$ is defined as the empirical probability mass function. In this setting, the 'best' estimate $(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) \in W_j \times \Theta_j$ for a given $j$ corresponds to


$$(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) = \argmin_{(\mathbf{w}, \bm{\theta})} L_2^2(f_{j, \mathbf{w}, \bm{\theta}}, \tilde{f_n}) 
= \argmin_{(\mathbf{w}, \bm{\theta})} \left\{ \sum_{x=0}^{\infty} f_{j, \mathbf{w}, \bm{\theta}}^2(x) - \frac{2}{n} \sum_{i=1}^{n}f_{j, \mathbf{w}, \bm{\theta}}(X_i)\right\}.$$


As the squared $L_2$ distance might involve an infinite sum (for distributions with infinite support), the user has the option to determine the cut-off value using the `n.inf` argument, which is set to 1000 by default. The parameters $(\hat{\mathbf{w}}^{j+1}, \hat{\bm{\theta}}^{j+1})$ are obtained analogously. Once both parameter sets have been determined, the difference in their respective squared $L_2$ distances to $\tilde{f}_n$ is compared to a `threshold` (equaling $t(j,n)$ defined above. The threshold function can be entered directly or one of the predefined thresholds, called `LIC` or `SBC` and given respectively by

$$\frac{0.6}{n} \ln\left(\frac{j+1}{j}\right) \quad\quad \text{ and} \quad\quad \frac{0.6 \ln(n)}{n} \ln\left(\frac{j+1}{j}\right)$$

can be used. Note that, if a customized function is to be used, its arguments have to be named `j` and `n`. If the difference in squared distances is smaller than the selected threshold, the algorithm terminates and the true order is estimated as $j$, otherwise $j$ is increased by 1 and the procedure starts over. The reader is invited to consult [@l2] for further details.

#### 2. `hellinger.disc`

This second function presents an alternative estimation procedure for *discrete* mixtures, working much the same as `L2.disc`, however, using a different measure of distance and different thresholds. As the name suggests, it is based on the square of the Hellinger distance, causing the 'best' estimate $(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) \in W_j \times \Theta_j$ for a given $j$ to equal

$$(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) = \argmin_{(\mathbf{w}, \bm{\theta})} H^2(f_{j, \mathbf{w}, \bm{\theta}}, \tilde{f_n}) 
= \argmax_{(\mathbf{w}, \bm{\theta})} \sum_{x=0}^{X_{(n)}} \sqrt{f_{j, \mathbf{w}, \bm{\theta}}(x) \tilde{f}_n(x)},$$

with $X_{(n)} = \max_{i = 1}^n (X_i)$. The relevant theory can be found in [@hell]. In accordance with their work, the two predefined thresholds are given by

$$\text{AIC} = \frac{d+1}{n} \quad \quad \text{and} \quad \quad \text{SBC} = \frac{(d+1)\ln(n)}{2n}$$

(recall that $d$ is the number of component parameters, i.e., $\Theta \subseteq \mathbb{R}^d$). If a customized function is to be used, its arguments have to named `j` and `n` once more, so if the user wants to include the number of component parameters $d$, it has to be entered explicitly. 

#### 3. `hellinger.cont`
 
Unlike the two preceding functions, this procedure is applicable to *continuous* mixture models and uses a kernel density estimator (KDE) as $\tilde{f}_n$. Its `bandwidth` can be chosen by the user, or the adaptive KDE found in [@adap, p. 1720, equation (2)] may be used by specifying `bandwidth = "adaptive"`. The calculations are based on the continuous version of the squared Hellinger distance, where the 'best' estimate $(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) \in W_j \times \Theta_j$ for a given $j$ corresponds to

\begin{equation}\label{eq:hellcont}
(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) = \argmin_{(\mathbf{w}, \bm{\theta})} H^2(f_{j, \mathbf{w}, \bm{\theta}}, \tilde{f_n}) 
= \argmax_{(\mathbf{w}, \bm{\theta})} \int \sqrt{f_{j, \mathbf{w}, \bm{\theta}}(x)\tilde{f}_n(x)}\ dx.
\end{equation}

Since the computational burden of optimizing over an integral to find the 'best' weights and component parameters is immense, the algorithm approximates the objective function \autoref{eq:hellcont} by sampling $n_s = $ `sample.n` observations $Y_i$ from $\tilde{f}_n(x)$ and setting

$$
(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j) = \argmax_{(\mathbf{w}, \bm{\theta})} \sum_{i = 1}^{n_s} \sqrt{\frac{f_{j, \mathbf{w}, \bm{\theta}}(Y_i)}{\tilde{f}_n(Y_i)}}.
$$

This procedure uses the same thresholds as `hellinger.disc`.

As before, we initially show the fit these methods yield on the two artificial datasets. As can be seen in the resulting plots, both `hellinger.disc` and `hellinger.cont` correctly estimate that the data comes from 3-component mixtures.

```{r plothel, figures-side, fig.show="hold", out.width="50%"}
set.seed(0)
h_disc_pois <- hellinger.disc(pois.dM, threshold = "AIC")
h_cont_norm <- hellinger.cont(normLoc.dM, bandwidth = 0.5, sample.n = 5000, 
                      threshold = "AIC")
par(mar = c(5, 5, 1, 1))
plot(h_disc_pois)
plot(h_cont_norm)
```
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/dist_art_1.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/dist_art_2.png" />
</p>

For a real-world example, refer back to the `faithful` dataset and the corresponding `datMix` object which was created in Section 1. Fitting the distance methods to a continuous density requires a choice of bandwidth. While using the adaptive bandwidth is an option, if the user does not want to do so, it is recommended to use the function `kdensity` from the package **kdensity** [@kdensity] which automatically selects an optimal bandwidth (can be accessed via `kdensity(data)$bw`). If the user wants to compare different bandwidth values, it is advisable to look at the plots of the respective kernel density estimates using `kdensity` and to choose one that captures the shape of the data well without fitting to noise.

The following figures illustrate the above point by showing the KDE of the Old Faithful data with bandwidths 1, 4 and 8. Here, 4 seems to be an appropriate choice.

<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/bandwidth1.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/bandwidth4.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/bandwidth8.png" />
</p>

`hellinger.cont` fits a 2-component mixture to the data, which fits the data well and comprises similar parameter estimates to those found in the literature.

```{r faithplothel, fig.width = 5, fig.height = 4}
# estimate the number of components:
library(kdensity)
res <- hellinger.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw,
                      sample.n = 5000, threshold = "AIC")
plot(res)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/hell-cont-norm.png">

At this point, it is worth having a closer look at the thresholds. They each satisfy $t(j,n) \rightarrow 0$ as $n \rightarrow \infty$, the sole condition the authors require. Now, the consistency proofs for estimators defined via Equation \autoref{eq:distances} all rely on the fact that, as $n \rightarrow \infty$,

$$D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n) \rightarrow d_j > 0, \text{ for } j < p$$

and

$$D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n) \rightarrow 0, \text{ for } j \geq p,$$

where $p$ is the true complexity (compare with [@l2, p. 4253, Proof of the Theorem], [@hell, p. 4383, Proof] and [@hellcont, p. 1485, Proof of Theorem 1]. If however $t(j,n)$ goes to 0 faster than $D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n)$ for $j \geq p$, asymptotically, the decision rule outlined above will always lead to $j$ being rejected. Therefore, a second condition should be placed on $t(j,n)$, namely choosing it in accordance with 

$$D(\hat{f}_p, \tilde{f}_n) - D(\hat{f}_{p+1}, \tilde{f}_n) = o_p(t(j,n)).$$

Neither the $L_2$ Information Criterion  (LIC) nor the Akakike Information Criterion (AIC), nor in the continuous case, the Schwarz Bayesian Criterion (SBC), satisfy this condition, yet they are still part of the package for two reasons. First, since they were used in the original papers, they are included for the sake of completeness and reproducibility of original results. Second, consistency is an asymptotic property, and while the aforementioned thresholds do not fulfill it, they still perform well (and not rarely better than consistent thresholds) for smaller sample sizes. In the example above, the number of components is correctly identified under the non-consistent AIC threshold. Nonetheless, the user will get a warning when using one of non-consistent predefined thresholds.

The preceding example shows that $\hat{p}$ directly depends on the chosen threshold $t(j, n)$, as is also obvious from Equation \autoref{eq:distances}. While some thresholds can be motivated better than others from a theoretical perspective, the choice will ultimately always remain somewhat arbitrary. It would thus be desirable to have versions of the preceding functions which do not suffer from this drawback. `L2.boot.disc`, `hellinger.boot.disc` and `hellinger.boot.cont` all work similarly to their counterparts, with the exception that the difference in distances is not compared to a predefined threshold but a value generated by a bootstrap procedure.  At every iteration (of $j$), the procedure sequentially tests $p = j$ versus $p = j+1$ for $j = 1,2, \dots$, using a parametric bootstrap to generate \code{B} samples of size $n$ from a $j$-component mixture given the previously calculated 'best' parameter values $(\hat{\mathbf{w}}^j, \hat{\bm{\theta}}^j)$. For each of the bootstrap samples, again the 'best' estimates corresponding to densities with $j$ and $j+1$ components are calculated, as well as their difference in distances from $\tilde{f}_n$. The null hypothesis $H_0: p = j$ is rejected and $j$ is increased by $1$ if the original difference $D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n)$ lies outside of the interval $[ql, qu]$, specified by the `ql` and `qu` empirical quantiles of the bootstrapped differences. Otherwise, $j$ is returned as the order estimate $\hat{p}$. 

Since the bootstrap version returns a very similar result to the threshold version on the Old Faithful dataset, we introduce a new example here. Consider the so-called Shakespeare dataset which comprises the number of occurrences of the words that Shakespeare used in his writings. For example, the number of times Shakespeare used a word only once is 14 376, while the number of times the same word occurred exactly 10 times in his writing is 363. The same data have been considered in other papers, see e.g., [@sp68], [@Efron1976], [@CheeWang2016] and [@balabdkulagina]. In the last three papers, the underlying statistical question that the authors wanted to answer is: how many words did Shakespeare actually know? This problem is known under the name of 'species richness' and can be solved using a variety of approaches. The goal is to use the observed frequencies of species, here words, to estimate the unobserved number of words that Shakespeare knew and did not use in his writings. While there is a whole spectrum of methods for estimating species richness, we limit ourselves here to motivate fitting a finite mixture of geometrics to the data. It is known from [@steutel69] that the class of completely monotone probability mass functions defined on the set of non-negative integers, that is the class of $p$ such that $(-1)^k \nabla^k p(i) \ge 0$ for all integers $k \ge 0, i \ge 0$ coincides with the class of all mixtures of geometrics probability mass functions (here $\nabla p = p(i+1) - p(i)$ and $\nabla^{r+1} \equiv \nabla \circ \nabla^r $ for any integer $r \ge 1$). In \citet{bdF2019}, the monotone complete least squares estimator (LSE) was defined for such class, which is the first non-parametric estimator that was considered for an element in such a family. Complete monotonicity can be defined on any subset of the set of integers of the form $\{a, a+1, \ldots \} $ for $ a \ge 1$ since the change of variable $x \mapsto x-a$ brings us back to complete monotonicity on $\{0, 1, \ldots \}$. It can be clearly seen ([@balabdkulagina]) that the complete monotone estimator fits very well the empirical estimator of the word occurrences. This result strongly suggests that complete monotonicity is a very appropriate model. In the scope of this paper, we want to explore how fitting a finite mixture of geometric distributions with unknown number of components works for this dataset. This alternative approach is actually inspired by the fact that the complete monotone LSE is itself a finite mixture of geometrics (with a random number of components). Such a result is rather universal and its exact statement can be found in Proposition 2.3 in [@bdF2019]. 

Since we inherently do not observe the number of words Shakespeare did not use, the data start at 1. However, using $Y = X-1$ and assuming $Y$ is a geometric mixture with parameters $\{p, w_1, \dots, w_p, \theta_1, \dots, \theta_p\}$ leads to the following model for $X$:

\begin{equation}\label{eq:shakespeare}
f(x)  =  w_1  (1-\theta_1)^{x-1}  \theta_1  +  \ldots  +  w_{p} (1-\theta_{p})^{x-1}  \theta_{p}, \ \ x \in \{1,2,\ldots \, 100\}
\end{equation}

In accordance with the **R**-function `dgeom`, the parametrization we use for the probability mass function of a geometric distribution means that $\theta_i \in [0,1)$ is the success probability for the $i$-th component, $i~\in~\{1, \ldots, p\}$. Building on Equation \autoref{eq:shakespeare}, the clear appropriateness of the complete monotone model for the word frequencies in the Shakespeare data can be complemented by a more applied interpretation, its underyling assumption being that words in any language belong to different categories depending on the context in which they are used. As there is a finite number of words, this justifies the appropriateness of fitting a finite mixture model, whose components would correspond to the aforementioned categories. With the $i$-th component distribution given by $g(x;\theta_i) = (1-\theta_i)^{x-1}  \theta_i,\ x = 1,2,\dots$, this expression can be seen as the probability of a word belonging to category $i$ not appearing (in some new work) after having previously been used $x$ times.

The `datMix` object corresponding to the Shakespeare dataset is generated as follows: 

```
shakespeare.obs <- unlist(shakespeare) - 1
# define the MLE function:
MLE.geom <- function(dat) 1 / (mean(dat) + 1)

Shakespeare.dM <- datMix(shakespeare.obs, dist = "geom", discrete = TRUE, 
MLE.function = MLE.geom, theta.bound.list = list(prob = c(0, 1)))

# estimate the number of components and plot the results:
set.seed(0)
res <- hellinger.boot.disc(Shakespeare.dM, B = 50, ql = 0.025, qu = 0.975)
plot(res)
```

<img src="https://github.com/yuliadm/mixComp/blob/main/images/hell-boot-geom.png">

`hellinger.boot.disc` estimates that the data comes from a 3-component geometric mixture (thus clustering the english words Shakespeare used into three categories).

# Section 5. Functions using LRTS

As a third option of estimating the mixture complexity, the **mixComp** package provides an algorithm based on the likelihood ratio test statistic (LRTS), sequentially testing $p = j$ versus $p = j+1$ for $j = 1,2, \dots$, until the algorithm terminates. As noted in Section 1, it is not possible to use the generalized likelihood ratio test in the context of mixture models directly as the standard way of obtaining the asymptotic distribution of the LRTS under the null hypothesis cannot be applied. However, one can get estimates of the test's critical values by employing a bootstrap procedure.

Making use of this approach, the function \code{mix.lrt} iteratively increases the assumed order $j$ and finds the MLE for both, the density of a mixture with $j$ and $j+1$ components, giving $(\hat{\mathbf{w}}^{j}, \hat{\bm{\theta}}^{j}) \in W_j \times \Theta_j$ and $(\hat{\mathbf{w}}^{j+1}, \hat{\bm{\theta}}^{j+1}) \in W_{j+1} \times \Theta_{j+1}$. It then calculates the corresponding LRTS, defined as

$$\text{LRTS}= -2\ln\left(\frac{L_0}{L_1}\right) \quad \text{, with}$$

$$L_0 = L_{\textbf{X}}(\hat{\mathbf{w}}^{j}, \hat{\bm{\theta}}^{j}) \quad\quad \text{and} \quad\quad L_1 = L_{\textbf{X}}(\hat{\mathbf{w}}^{j+1}, \hat{\bm{\theta}}^{j+1})\text{,}$$

$L_{\textbf{X}}$ being the likelihood function given the data ${\textbf{X}}$.

Next, a parametric bootstrap is used to generate `B` samples of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-component mixture given the previously calculated MLE $(\hat{\mathbf{w}}^{j}, \hat{\bm{\theta}}^{j})$. For each of the bootstrap samples, the MLEs corresponding to the mixture densities with $j$ and $j+1$ components are calculated, as well as the LRTS. The null hypothesis $H_0: p = j$ is rejected and $j$ increased by $1$ if the LRTS based on the original data vector $\textbf{X}$ is larger than the chosen `quantile` of its bootstrapped counterparts. Otherwise, $j$ is returned as the order estimate $\hat{p}$. For further details, the reader is referred to [@lrt].

For the two artificial datasets, this method estimates 3-component mixtures with very similar parameters to the distance methods, so we go straight to a real-world example. Consider the Acidity dataset which comprises measurements of the acid neutralizing capacity (ANC) taken from 155 lakes in North-Central Wisconsin. The ANC indicates a lakes' capability to absorb acid, with low values potentially leading to a loss of biological resources. This dataset has been analysed as a mixture of normal distributions on the log scale by [@acidity1], [@acidity2] and [@acidity3]. While the former papers suggest the number of components to equal 2 (with 3 also being considered), the latter estimates $p$ to lie between 3 and 5. The `mix.lrt` method agrees with [@acidity1] and [@acidity2], returning a 2-component mixture with reasonable estimates for the component weights and parameters.

```{r lrtacid, fig.width = 5, fig.height = 4, results='hide', message=FALSE, warning=FALSE}
acidity.obs <- unlist(acidity)

acidity.dM <- datMix(acidity.obs, dist = "norm", discrete = FALSE, 
                     MLE.function = MLE.norm.list, 
                     theta.bound.list = norm.bound.list)

set.seed(0)
res <- mix.lrt(acidity.dM, B = 50, quantile = 0.95)
plot(res)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/lrt-norm.png">

# Section 6. Non-standard mixtures

In all preceding examples, the families of component densities $\{g(x;\theta):\theta \in \Theta \}$ belonged to one of the "standard" probability distributions included in the **stats** package, which provides the density/mass function, cumulative distribution function, quantile function and random variate generation for selected distributions. The function names are of the form `dxxx`, `pxxx`, `qxxx` and `rxxx` respectively. With some additional effort, it is possible to use the **mixComp** package on "non-standard" distributions -- the user merely has to provide functions evaluating the density and generating random numbers for the component distribution. In accordance with **R** conventions, the user-generated function `dxxx` has to take `x` and the distribution parameters as input and returns the value of the density function specified by the parameters at the point `x`. Likewise, `rxxx` requires `n` and the distribution parameters as input and returns `n` random numbers based on the distribution specified by the aforementioned parameters.

As an example, consider an artificial dataset generated by sampling $\mathbf{X}$ from a 3-component mixture of normals with means equal to 10, 11 and 13. Assume that the standard deviation of all components is known to be $0.5$, yet the number of components and their means are unknown. Then each of the components follows a $\mathcal{N}(\mu, 0.5)$ distribution, which shall be called `norm0.5`. The first step is always that of creating the `dxxx` and `rxxx` functions, since they will be called by the **mixComp** functions.

The following example creates the `Mix` and `rMix` objects
based on the density of a normal mixture with $\mathbf{w} = (0.3, 0.4, 0.3)$, $\bm{\mu} = (10, 11, 13)$ and \linebreak$\bm{\sigma} = (0.5, 0.5, 0.5)$ and plots the obtained mixture density and the corresponding random sample. 

```{r, figures-side, fig.show="hold", out.width="50%", results='hide', message=FALSE, warning=FALSE}
dnorm0.5 <- function(x, mean){
  dnorm(x, mean = mean,  sd = 0.5)
}
rnorm0.5 <- function(n, mean){
  rnorm(n, mean = mean,  sd = 0.5)
}
## create objects `Mix` and `rMix`:
set.seed(1)
norm0.5Mix <- Mix("norm0.5", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 11, 13))
norm0.5RMix <- rMix(1000, obj = norm0.5Mix)
## plot the results:
plot(norm0.5Mix)
plot(norm0.5RMix)
```
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/norm0.5Mix.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/norm0.5RMix.png" />
</p>

Below we will estimate of the mixture density using `mix.lrt` given a sample from the considered above 3-component normal mixture. We start by creating all necessary inputs:
```{r}
norm0.5.list <- vector(mode = "list", length = 1)
names(norm0.5.list) <- c("mean")
norm0.5.list$mean <- c(-Inf, Inf)

MLE.norm0.5 <- function(dat) mean(dat)

norm0.5.dM <- RtoDat(norm0.5RMix, theta.bound.list = norm0.5.list,
                     MLE.function = MLE.norm0.5)
```
Now the **mixComp** procedures can be used on the `datMix` object as usual. The results can be printed and plotted using `print` and `plot` functions.
```{r, results='hide', message=FALSE, warning=FALSE}
set.seed(1)
res <- mix.lrt(norm0.5.dM, B = 50, quantile = 0.95)

print(res)
plot(res)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/est-norm0.5.png">

# Section 7. Computational nuance for mixComp functions using the solnp() solver

Several functions in the **mixComp** package (namely, `nonparamHankel`, `paramHankel`, `hellinger.cont`, `hellinger.disc`, `L2.disc`, `mix.lrt`) make use of `solnp()`  function (**Rsolnp** library), which is a solver for general nonlinear programming problems. The above mentioned **mixComp** functions attempt to generate good starting values by clustering the data via `clara` function prior to applying  `solnp()`, which lead to convergence of the algorithm in most of the cases. However when running multiple simulations, `solnp()` might not converge for particular initial values with default control values. This may happen when very few observations are assigned to some of the clusters have, in which case the solver can get "stuck", not even resulting in bad exit status (codes 1 and 2 in returned by `solnp()` convergence value). This issue can be overcome by specifying the control parameters in the functions using the `solnp()` solver (e.g. by defining `L2.boot.disc(geom.dM, j.max = 5, B = 500, ql = 0.025, qu = 0.975, control = list("rho" = 0.1, tol = 0.0000001, trace = 0))`, see description of the `solnp()` function for further details on setting the control parameters) or by setting a time limit for the function execution and going to the next iteration whenever the time limit is exceeded. 


# Section 8. Summary and discussion

In this paper we presented the **R** package **mixComp**, a collection of routines developed to estimate a mixture's complexity from a data sample $\mathbf{X}$. Moreover, it provides the possibility of generating artificial data from specified mixtures as well as proper visualization tools for the complexity estimation result, plotting either the successive determinant values or the final fitted mixture. If estimates of the component weights and parameters are sought, $\hat{p}$ can be passed to one of the many **R** packages specialized in their calculation, or one of the **mixComp** functions returning weight and parameter estimates can be used. However, it should be noted that the theory on which this package is based concentrates on showing consistency for $\hat{p}$ -- other estimates obtained in the process are merely by-products. A primary goal of the package was to make it extendible, meaning that it can be utilized on mixtures of less well-known distributions, as long as sufficient information on the density and random variate generation is provided. We hope that this package will be useful for practitioners in the many areas where mixture models are applicable.

# Computational details

All computations and graphics in this paper have been done using **R** version 4.0.0 with the packages **boot** 1.3-24, **cluster** 2.1.0, **expm** 0.999-4, **matrixcalc** 1.0-3, **Rsolnp** 1.16 and **kdensity** 1.0.1. **R** itself
and all packages used are available from the Comprehensive **R** Archive Network (CRAN) at https://CRAN.R-project.org/.


# Appendix 
### Distance and non-parametric estimator definitions

Let $\textbf{X} = \{X_1, \dots, X_n\}$ be an i.i.d. sample of size $n$ from some continuous distribution with unknown density $f$. Its *kernel densty estimator* is defined as

$$\tilde{f}_n(x) \coloneqq \frac{1}{nc_n}\sum_{i=1}^n K\left( \frac{x-X_i}{c_n} \right),$$	

with kernel $K$ and bandwidth $c_n$.

As an extension, [@adap] defined the *adaptive kernel density estimator* as

$$\tilde{f}_n(x) \coloneqq \frac{1}{n}\sum_{i=1}^n \sum_{j=1}^k \frac{a_j(X_i)}{c_{n, j}} K\left( \frac{x-X_i}{c_{n, j}} \right),$$	

with kernel $K$, bandwidth $c_{n, j}$ and weights (positive and summing to $1$) $a_j(X_i)$.
	
Let $\textbf{X} = \{X_1, \dots, X_n\}$ be an i.i.d. sample of size $n$ from some discrete distribution with unknown probability mass function $f$. Its *empirical probability mass function* is defined as

$$\tilde{f}_n(x) \coloneqq \frac{1}{n}\sum_{i=1}^n \mathbbm{1}_{\{X_i = x\}}.$$	

Let $g, f$ be two probability mass functions defined on $\mathcal{X}$.
The squared $L_2$ distance between $g$ and $f$ is given by

$$L_2^2(g,f) \coloneqq \sum_{\mathcal{X}} (g(x)-f(x))^2.$$	
	

Let $g, f$ be two density or probability mass functions defined on $\mathcal{X}$.
The squared Hellinger distance between $g$ and $f$ is given by

$$H^2(g,f) \coloneqq \sum_{\mathcal{X}} \left( \sqrt{g(x)}-\sqrt{f(x)} \right)^2$$	
	
in the discrete case and by

$$H^2(g,f) \coloneqq \int_{\mathcal{X}} \left( \sqrt{g(x)}-\sqrt{f(x)} \right)^2 dx$$	

in the continuous case.

### Further results on the real-world datasets

The Table below showcases the results of all complexity estimation procedures when applied to the four real-world datasets that were discussed in this paper. The table also includes the results for the Shakespeare dataset borrowed from [@shakespeare]. It was not mentioned in the vignette but is included in the **mixComp** package. The following settings were used to calculate these results (default setting were used unless indicated otherwise):

* `set.seed(1)` was used for all complexity calculations,
* for `hellinger.cont` and `hellinger.boot.cont` the bandwidths suggested by **kdensity** were used,
* for `nonparamHankel`, we used $A(j)l(n) = \frac{j\log(n)}{\sqrt{n}}$  and `j.max = 6`,
* for the bootstrapped distance and the LRTS methods, we used `B = 50`.

#### Table A: Results for all **mixComp** methods used on real-world data
|  Method                   | Children | Old Faithful | Shakespeare | Acidity |
|:-------------------------:|:--------:|:------------:|:-----------:|:-------:|
| `nonparamHankel`          |    2     |      x       |      2      |    x    |
| `nonparamHankel`(scaled)  |    2     |      x       |      2      |    x    |
| `paramHankel`             |    2     |      x       |      2      |    x    |
| `paramHankel.scaled`      |    2     |      x       |      2      |    x    |
| `L2.disc`                 |    2     |      x       |      3      |    x    |
| `L2.boot.disc`            |    2     |      x       |      4      |    x    |
| `hellinger.disc`          |    2     |      x       |      3      |    x    |
| `hellinger.boot.disc`     |    2     |      x       |      3      |    x    |
| `hellinger.cont`          |    x     |      2       |      x      |    2    |
| `hellinger.boot.cont`     |    x     |      2       |      x      |    2    |
| `mix.lrt`                 |    2     |      2       |      3      |    2    |

# References



The **mixComp** package provides a number of methods for estimating the complexity of a finite mixture.
The considered approaches can be loosely grouped into three categories:
<ul>
  <li> methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; </li>
  <li> methods based on penalized minimum distance between the unknown probability density and a consistent estimator thereof. The distances considered in this survey are the Hellinger and the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973841.jpg">-distances; </li>
  <li> likelihood ratio test (LRT) - based techniques. </li>
</ul>
While not the primary goal, most methods simultaneously estimate the component weights and parameters. In this document, we give a brief overview of the methodology, and demonstrate the package's functionality in both real world examples and synthetically generated data. Moreover, we show how the package can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution.


The relevant pepers include but are not confined to (the full list of references is available in the vignette):
<ol>
<li> Dacunha-Castelle, Didier, and Elisabeth Gassiat. The estimation of the order of a mixture model. Bernoulli 3, no. 3 (1997): 279-299. <https://projecteuclid.org/download/pdf_1/euclid.bj/1177334456>. </li>
<li> Woo, Mi-Ja, and T. N. Sriram. Robust estimation of mixture complexity. Journal of the American Statistical Association 101, no. 476 (2006): 1475-1486. <doi:10.1198/016214506000000555>. </li>
<li> Woo, Mi-Ja, and T. N. Sriram. Robust estimation of mixture complexity for count data. Computational statistics & data analysis 51, no. 9 (2007): 4379-4392. <doi:10.1016/j.csda.2006.06.006>. </li>
<li> Umashanger, T., and T. N. Sriram. L2E estimation of mixture complexity for count data. Computational statistics & data analysis 53, no. 12 (2009): 4243-4254. <doi:10.1016/j.csda.2009.05.013>. </li>
<li> Karlis, Dimitris, and Evdokia Xekalaki. On testing for the number of components in a mixed Poisson model. Annals of the Institute of Statistical Mathematics 51, no. 1 (1999): 149-162. <doi:10.1023/A:1003839420071>. </li>
<li> Cutler, Adele, and Olga I. Cordero-Brana. Minimum Hellinger Distance Estimation for Finite Mixture Models. Journal of the American Statistical Association 91, no. 436 (1996): 1716-1723. <doi:10.2307/2291601>. </li>
</ol>
	
For illustrative purposes several datasets have been included in the package:
<ol>
<li> **accidents**, from Karlis, Dimitris, and Evdokia Xekalaki. On testing for the number of components in a mixed Poisson model. Annals of the Institute of Statistical Mathematics 51, no. 1 (1999): 149-162. <doi:10.1023/A:1003839420071>. </li>
<li> **acidity**, from Sybil L. Crawford, Morris H. DeGroot, Joseph B. Kadane & Mitchell J. Small (1992) Modeling Lake-Chemistry Distributions: Approximate Bayesian Methods for Estimating a Finite-Mixture Model, Technometrics, 34:4, 441-453. <doi:10.1080/00401706.1992.10484955>. </li>
<li> **children**, from Thisted, R. A. (1988). Elements of statistical computing: Numerical computation (Vol. 1). CRC Press. </li>
<li> **faithful**, from R package "datasets"; Azzalini, A. and Bowman, A. W. (1990). A look at some data on the Old Faithful geyser. Applied Statistics, 39, 357--365. <https://www.jstor.org/stable/2347385>. </li>
<li> **shakespeare**, from Efron, Bradley, and Ronald Thisted. "Estimating the number of unseen species: How many words did Shakespeare know?." Biometrika 63.3 (1976): 435-447. <doi:10.1093/biomet/63.3.435>. </li>
</ol>


## Installation

To install from `CRAN`, use:
```
install.packages("mixComp")
```

You can install the development version from GitHub with:
```
devtools::install_github("yuliamd/mixComp")
```

## Section 1. Introduction to finite mixture models and mixComp


Mixture models have been used extensively in statistical applications and therefore have attracted a lot of attention from both theoretical and computational perspectives. Although the list of works on mixture models is too long to make an exhaustive inventory, we can refer to the following important papers and books: [@Teicher63], [@LindsayI], [@LindsayII], [@Titterington] and [@McLachlan].

The popularity of such models stems from the fact that they allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. To account for such heterogeneity, the (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method. This includes the well-known nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were well studied in the seminal work of [@LindsayI], [@LindsayII]. One remarkable property of the NPMLE of the mixing distribution is that it is, under some simple conditions, a discrete distribution function with at most <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973472.jpg"> number of jumps, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973472.jpg"> is the number of distinct observations in the random sample drawn from the mixture. This interesting feature is one reason, among others, for considering the smaller class of finite mixture models, i.e., mixture models with a discrete mixing distribution with a finite number of jumps. The model has the following simple interpretation: the population under study is assumed to consist of several homogeneous subpopulations. These subpopulations, typically referred to as the mixture's components, usually have a specific meaning depending on the problem at hand. In some very simple situations, the number of components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric maximum likelihood estimator (MLE) is known to occur at the fast rate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973529.jpg"> (under some regularity conditions). Also, the well-known expectation-maximization (EM) algorithm can be used to find the MLE of all the unknown parameters; see for example [@Dempster]. However, in many statistical applications such knowledge is rarely available and the number of components has to be estimated from the data. Although the mixture is still finite and the distribution of each component is assumed to belong to some parametric family, the estimation framework in this case is much harder than in the fully parametric one, where the number of components is known. In this paper, the terms *order*, *complexity*  and *number of components*  will  be used interchangeably to refer to this unknown number. The main goal of the package **mixComp** is to estimate the unknown complexity using several methods known from the statistical literature. These methods, which are discussed below in more detail, all come with theoretical guarantees for consistency as the sample size gets larger. Of course, consistency in this case means that an estimator is able to exactly recover the unknown complexity for large sample sizes. As expected, the performance of the methods varies according to the underlying mixture distribution and the sample size. This will be illustrated below through several synthetic as well as real datasets.

To describe the estimation problem, we start with some formal notation.  A distribution <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> is called a *finite mixture* if its density (we write density throughout and keep in mind that it may be taken with respect to the Lebesgue or the counting measure) is of the form

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973590.jpg">

where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973785.jpg"> is the mixture complexity,

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973708.jpg"> 

are the component weights and the density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973751.jpg"> is the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975606.jpg">-th component of the mixture. As the scope of **mixComp** is limited to mixtures where the family of the component distributions is known, we replace <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973751.jpg"> by a parametric density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976027.jpg"> indexed by the (possibly multivariate, say <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg">-dimensional) parameter <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976173.jpg"> in the parameter space <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922033.jpg">.

Given some complexity <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">, the two relevant parameter spaces can therefore be defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976536.jpg">

and

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976599.jpg">

Throughout this document, it is assumed that the family of the component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976784.jpg"> is known, but the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976825.jpg">, the component weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977021.jpg"> and the mixture complexity <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973785.jpg"> are unknown, with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> being the parameter of interest. Assume now that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> is a finite mixture distribution with density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973590.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977574.jpg"> is an i.i.d. sample of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg">. The **mixComp** package aims to estimate the smallest such <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> on the basis of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">, either on its own or by simultaneously estimating the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977761.jpg"> and the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976173.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978071.jpg">.

In this setup, it seems natural to test for the number of components by comparing two consecutive models. Traditionally, the problem of choosing between nested models may be approached by applying the generalized likelihood ratio test and referring to the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978511.jpg">  distribution to assess significance, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978628.jpg"> is given by the number of constraints imposed on the alternative hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978742.jpg"> to arrive at the null hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978727.jpg">. However, in the context of mixture models, there are several issues hindering application of this classical theory. One of them is that there is no unique way of obtaining <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978727.jpg"> from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978742.jpg">. As an example, the two null hypotheses <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978884.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978911.jpg"> both yield the smaller model, showcasing the difficulties of applying the classical asymptotic theory of the likelihood ratio. This problem has been studied extensively in the literature and numerous alternative approaches to mixture complexity estimation have been suggested, laying the theoretical foundation for the subsequently described algorithms.

This document discusses various categories of functions found in the **mixComp** package, ranging from methods based on Hankel matrices, to techniques based upon distances between densities and likelihood ratio tests. The examples provided in the first sections all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be performed using the functions from the **stats** package. The last section illustrates how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and valuating the density thereof.

Two main features distinguish this package from other mixture-related `R` [@R] packages: Firstly, it is focused on the estimation of the complexity rather than the component weights and parameters. While these are often estimated as a by-product, all methods contained in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture of interest. Secondly, it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** [see @mixtools] and **flexmix** [see @flexmix1; @flexmix2; @flexmix3]  should both be mentioned at this point: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm. Notably, it also contains routines for selecting the number of components based on information criteria and parametric bootstrapping of the likelihood ratio test statistic values. However, they are limited to multinomial and (a variety of) normal mixtures as well as mixtures-of-regressions. Second, while **flexmix** was developed to deal with mixtures-of-regression, it sets itself apart from other packages by its extensibility, a design principle that we also aimed for when creating the  **mixComp** package. Other widely used packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], which is used for grouped conditional data. Interested readers can find a comprehensive list of mixture-related packages on the CRAN Task View: Cluster Analysis and Finite Mixture Models website.

Before moving to the description of the different methods implemented in **mixComp** we would like to briefly mention pther theoretical work on the estimation of mixture complexity not currently included in the package. [@chen] propose a method that is reminiscent of the ones described in Section 4. The main difference is that the authors consider distribution functions instead densities, i.e. they consider minimizing a penalized distance between the distribution function of the mixture and the empirical distribution function. The approach of [@figueiredo] is based on a minimum message length-like criterion, however, their method struggles to deal with mixtures with very different weights. [@xian] propose a procedure based on alternating between splitting and merging the components in an EM-algorithm. This algorithm requires selecting two thresholds, the choice of which is somewhat unclear when working with a specific dataset. [@miller] follow a Bayesian approach, taking the usual finite mixture model with Dirichlet weights and putting a prior distribution on the unknown number of components. 

The methods that were included in the package can be roughly devided into three categories: methods based on Hankel matrices, following the theory as described in [@hankel] and selected because of the fact that computation of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979158.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979228.jpg"> is not required, a method based on the likelihood ratio test statistic (LRTS) following [@lrt] since a likelihood ratio test seems like a natural approach in this setting and methods employing minimum distance calculations based on several works and included as a computationally more efficient alternative to the LRTS method for certain distributions and distances; see [@hell; @hellcont; @l2; @adap]. For example, when the distance is taken to be the Hellinger distance, such an approach is especially fast for discrete distributions. For a more fluid reading, the relevant theory will be portrayed at the beginning of each of the respective sections. The examples depicted in these first chapters all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be done by functions available from the **stats** package. The last chapter showcases how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and evaluating the density thereof.


## Section 2. Objects and functions defined in mixComp

Table 1 depicts five object classes defined in **mixComp**. The first two respectively represent a finite mixture distribution and a random sample drawn from such a distribution. The `Mix` object is printed as a matrix of component weights and parameters and is plotted as the density of the specified mixture, showing the overall as well as the component densities. The `rMix` object prints as the vector of observations and plots as a histogram, showcasing the individual components as well as the full sample. Both objects contain a number of attributes giving additional information, details of which can be found in the corresponding **R** help files. 

#### Table 1: Objects and functions defined in mixComp
|  Object class  | Created via                                  | Description                                     |
|:--------------:|:--------------------------------------------:|:-----------------------------------------------:|
| `Mix`          | `Mix`                                        | Represents a finite mixture                     |
| `rMix`         | `rMix`                                       | Randomly generated data from a finite mixture   |
| `datMix`       | `datMix` or `RtoDat`                         | Observed data from (presumably) a finite mixture|
| `hankDet`      | `nonparamHankel`                             | Vector of estimated Hankel matrix determinants  |
| `paramEst`     | `paramHankel(.scaled)`, `L2(.boot).disc`, `hellinger(.boot).disc`, `hellinger(.boot).cont` or `mix.lrt`  | Complexity estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg">, together with estimates of the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014622.jpg"> and the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014686.jpg">|

The generation of an object of class `Mix` hinges on four central arguments: a string `dist` specifying the name of the family of component densities (or kernels) <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014854.jpg">, a boolean`discrete` stating whether the distribution is discrete, a vector `w` giving the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977761.jpg"> and a list `theta.list` (the component parameters can also be supplied via the `...` argument) containing the parameters of the component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015122.jpg">. While the creation of `Mix` objects is mostly straightforward, two things should be noted in this regard: First, **mixComp** procedures will search for functions called `rdist` and `ddist` in the accessible namespaces. For most "standard" distributions, these functions are contained in the **stats** package and do not need to be user-written (compare with the Section 6). To make use of these functions, it is essential that the string `dist` is named correctly (e.g. to create a gaussian mixture on the basis of the **stats** package, `dist` has to be specified as `norm` instead of `normal`, `gaussian` etc. for the package to find the functions `dnorm` and `rnorm`). Second, the names of the list elements of `theta.list`(for the names of the `...` arguments) have to match the names of the formal arguments of the functions `ddist` and `rdist` exactly (e.g. for a gaussian mixture, the list elements have to be named `mean` and `sd`, as these are the formal arguments used by `rnorm` and `dnorm` functions of the **stats** package).

The following example creates two `Mix` objects, a 3-component mixture of normal distributions and a 3-component mixture of Poisson distributions. 

```{r mixobj}
set.seed(0)
# construct a Nix object:
normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
# plot the mixtures:
plot(normLocMix, main = "3-component normal mixture", cex.main = 0.9)
plot(poisMix, main = "3-component poisson mixture", cex.main = 0.9)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/normMix.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/poisMix.png">

If required, random samples can be generated from these mixtures.
```{r rmix}
# generate random samples:
normLocRMix <- rMix(1000, obj = normLocMix)
poisRMix <- rMix(1000, obj = poisMix)
# plot the histograms of the random samples:
plot(normLocRMix, main = "Three component normal mixture", cex.main = 0.9)
plot(poisRMix, main = "Three component poisson mixture", cex.main = 0.9)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/normRMix.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/poisRMix.png">

The third object class shown in Table 1, called `datMix`, represents the data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> based on which the mixture complexity is supposed to be estimated. These objects are most central to the package, as every procedure estimating the order of a mixture takes a `datMix` object as input. Apart from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">, it contains other "static" information needed for the estimation procedure (in contrast to "tuning parameters", which can be changed with every function call. An example of such a tuning parameter is the number of bootstrap replicates for a function employing a bootstrap procedure). A brief overview of which "static" attributes need to be supplied for each complexity estimation routine is given in Table 2. 

#### Table 2: Inputs needed for different functions contained in mixComp
| `R function`            |`dist`|`discrete`|`theta.bound.list`|  `MLE.function` | `Hankel.method` | `Hankel.function` |
|:-----------------------:|:----:|:--------:|:----------------:|:---------------:|:---------------:|:-----------------:|
|`nonparamHankel`         |  x   |          |                  |                 |        x        |         x         |
|`nonparamHankel(.scaled)`|  x   |    x     |         x        | o + i (optional)|        x        |         x         |
|`L2(.boot).disc`         |  x   |    x     |         x        |   i (optional)  |                 |                   |
|`hellinger(.boot).disc`  |  x   |    x     |         x        |   i (optional)  |                 |                   |
|`hellinger(.boot).cont`  |  x   |    x     |         x        |   i (optional)  |                 |                   |
|`mix.lrt`                |  x   |    x     |         x        | o + i (optional)|                 |                   |

In the table, "o" and "i" respectively stand for input used during optimization and initialization.

As a simple example of a given dataset to which mixture models have been applied extensively, take the Old Faithful dataset [@R; @faithful1; @faithful2]. In the context of mixture model estimation, the variable `waiting`, which gives the time in minutes between eruptions of the Old Faithful geyser in the Yellowstone National Park, is often considered to be the variable of interest. To estimate the number of components of the mixture distribution that provides a suitable approximation to the `waiting` data via **mixComp**, the raw data vector of observations has to be converted to a `datMix` object first. For the sake of exposition we specify all arguments of  the `datMix` function, starting with the vector of observations <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> and the string `dist`, specifying <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015320.jpg"> and the boolean `discrete`. As has often been done in the relevant literature, we assume that the data comes from a normal mixture.

```{r faithopts}
faithful.obs <- faithful$waiting
norm.dist <- "norm"
norm.discrete <- FALSE
```

Second, a named list of length <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg"> containing the bounds of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015725.jpg"> has to be created. In this example, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015826.jpg">.

```{r normlist}
# define the range for parameter values:
norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
```

Next, the argument `MLE.function` contains a single function if <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015658.jpg"> or a list of functions of length <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg"> otherwise, specifying the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg"> functions needed to estimate the MLE's of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979228.jpg"> based on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> if <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> were equal to 1 (i.e. the MLE's of the component distribution). If this argument is supplied and the `datMix` object is handed over to a complexity estimation procedure relying on optimizing over a likelihood function, the `MLE.function` attribute will be used for the single component case. In case the objective function is either not a likelihood or corresponds to a mixture with more than 1 component, numerical optimization will be used based on **Rsolnp**'s function `solnp` [@solnp, @Rsolnp]. The initial values (for the parameters of a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-component mixture, say) supplied to the solver are then calculated as follows: the data <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> is clustered into <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> groups by the function `clara` (of the **cluster** package by [@cluster] and the data corresponding to each group is given to `MLE.function`. The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates. Specifying `MLE.function` is optional and if it is not, for example because the MLE solution does not exists in closed form, numerical optimization is used to find the relevant MLE's.

Presuming a normal mixture, one specifies 2 functions, namely the MLE of the mean <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015978.jpg"> and the MLE of the standard deviation <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650016036.jpg">.

```{r normfun}
# define the MLE functions for the mean and sd: 
MLE.norm.mean <- function(dat) mean(dat)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)
```

The last two arguments, `Hankel.method` and `Hankel.function`, need to be supplied if the mixture complexity is to be estimated based on the Hankel matrix of the moments of the mixing distribution. The reader is referred to the Section 3 for further information on how these arguments are to be specified (in this case, the simplifying assumption of unit variance is made. This would be a poor choice for the `waiting` data, so <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> should not be estimated with one of the methods using these arguments, namely `nonparamHankel`, `paramHankel` and `paramHankel.scaled`, see Table 2). 

```{r normmom}
method <- "translation"
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}
```

Finally, all previously generated objects are combined to a `datMix` object.

```{r faithdatmix}
# construct a datMix object that summarizes all the necessary information:
faithful.dM <- datMix(faithful.obs, dist = norm.dist, discrete = norm.discrete,
                      theta.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = method,
                      Hankel.function = mom.std.norm)
```

In the preceding example, the data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> was taken from an already existing dataset. As seen before, the `rMix` function can be used to generate a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg">-sized sample from a specific mixture. If this synthesized data is to be used in simulations (i.e. passed to one of the functions estimating the mixture complexity) an `rMix` object can be converted to a `datMix` object via the `RtoDat` function. Apart from `dist` and `discrete`, all `datMix` arguments have to be supplied to `RtoDat` likewise. 

Unlike the above mentioned objects whose creation precedes any type of mixture complexity estimation, objects of the bottom two classes (see Table 1) contain the results of the estimation procedures. Generally, the functions estimating the number of components differ in the types of families of component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976784.jpg"> for which they allow and in whether they provide estimates of the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650016317.jpg"> and the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650016699.jpg">, the latter determining the object class of the estimation result. These differences are shown in Table 3. The function `nonparamHankel` returns an object of class `hankDet`, which is a vector of determinants (scaled and/or penalized), each entry belonging to a certain complexity estimate. The link between these determinants and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> will be discussed in the Section 3. `paramEst` objects arise when using any other function estimating the mixture complexity, all of which additionally return estimates of the component weights and parameters. For both object classes, print and plot methods are available to summarize and visualize the estimation results. 

#### Table 3: Distribution restrictions and output types of different functions contained in mixComp
|     R function                    |Restrictions on the family of the component density |Estimation of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979158.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979228.jpg">|
|:---------------------------------:|:--------------------------------------------------:|:----------------------------------------------:|
|`nonparamHankel`                   |Compatible with `explicit`, `translation` or `scale`|                                                |
|`nonparamHankel(.scaled)`          |Compatible with `explicit`, `translation` or `scale`|                       x                        |  
|`L2(.boot).disc`                   |Discrete distributions                              |                       x                        | 
|`hellinger(.boot).disc`            |Discrete distributions                              |                       x                        | 
|`hellinger(.boot).cont`            |Continuous distributions                            |                       x                        | 
|`mix.lrt`                          |                                                    |                       x                        | 


## Section 3. Functions using Hankel matrices

In 1997, [@hankel] showed that a mixture's order <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> is characterized by the smallest integer <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> such that the determinant of the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017103.jpg"> Hankel matrix of the first <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017156.jpg"> moments of the mixing distribution equals zero. Moreover, it can be shown that this determinant is zero for all <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017206.jpg">. Formally, for any vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017273.jpg"> with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017336.jpg"> equal to 1, the *Hankel matrix* of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017386.jpg"> is defined as the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017441.jpg"> matrix given by

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017499.jpg">

Now, let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017555.jpg"> be the vector containing the first <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017605.jpg"> (raw) moments of the mixing distribution. For finite mixture models, this amounts to

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017661.jpg">

Then, for all <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017720.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017786.jpg"> is non-negative and

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017845.jpg">

Making use of this fact, the first approach to estimating the order of a mixture that is implemented in **mixComp** relies on initially finding a consistent estimator of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018010.jpg"> based on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">, say <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg">, to then iteratively calculate the applicable Hankel matrix while increasing the assumed order <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> until a sufficiently small value of 
<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650203349.jpg"> is attained. However, since <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650203349.jpg"> should be close to 0 for all <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017206.jpg">, this would lead to choosing <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg"> rather larger than the true value and it seems natural to introduce a penalty term. Therefore [@hankel] define the empirical penalized objective function as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018549.jpg">

with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018991.jpg"> being a positive function converging to 0 as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019043.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019115.jpg"> being positive and strictly increasing. 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019176.jpg">

is then a consistent estimator of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg">.

As an extension to simply adding a penalty term to the determinant, a scaling approach was considered by [@lilian]. Let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019301.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019380.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019564.jpg">. Since the estimated moments <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg"> are asymptotically normal, one can apply the delta method giving

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019736.jpg">

Instead of inspecting the vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019799.jpg">, one could therefore also base the complexity analysis on a vector of scaled determinants, employing a nonparametric bootstrap procedure on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">.
To this end, let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019902.jpg"> denote the covariance matrix of the determinants <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021468.jpg"> calculated on the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021576.jpg">-th bootstrap sample for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021631.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021686.jpg">. Note that 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021736.jpg">

and write <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021787.jpg">. Define the rescaled vector 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021849.jpg">

Note that in the case of the scaled version,  the criterion to be minimized becomes 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021903.jpg">

That is, the chosen penalty function should be multiplied by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973529.jpg">.

This approach was proposed to address the issue of determinants already being very small from the beginning (even for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022011.jpg">), which, in the simulations by [@lilian], made it hard to discern the "best" complexity estimate, a problem that was not reduced much by solely adding a penalty term.  

With this general framework in place, the computation now merely hinges on calculating <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg">. The **mixComp** package offers three methods to do so. The method to use depends on the family of component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015320.jpg"> and is linked to some function <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022217.jpg"> needed to estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg">. The calculation method and the relevant function are specified when creating the `datMix` object as arguments `Hankel.method` and `Hankel.function`.


#### 1. `Hankel.method = "explicit"`

This method can be applied when a closed form expression for estimates of the moments of the mixing distribution exists. `Hankel.function` then contains the function explicitly estimating <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018010.jpg">. 

As an example, consider a mixture of geometric distributions, where it can be shown that

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022531.jpg">

with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> the true cumulative distribution function. Hence one may take

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022746.jpg">

as an estimator, with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022798.jpg"> being the empirical distribution function.

```{r geommom}
# define the function for computing the moments:
explicit.geom <- function(dat, j){
  1 - ecdf(dat)(j - 1)
}
```

As a second example, consider what [@hankel, p. 283, equation (3)] called the "natural" estimator, i.e. using

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022851.jpg">

when

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022913.jpg">

Note that the estimators of this form may also be supplied as `Hankel.method = "explicit"` with `Hankel.function`. For example, the "natural" estimator is applicable in the case of Poisson mixtures. If <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022964.jpg">, it is a well known fact that

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023020.jpg">

which then suggests using

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023088.jpg">

as an estimator.

```{r geompois}
# define the function for computing the moments:
explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
         matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}
```


#### 2. `Hankel.method = "translation"`
 
In Example 3.1., [@hankel, p.284] describes how <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018010.jpg"> can be estimated if the family of component distributions <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023420.jpg"> is given by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023472.jpg">, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is a known probability distribution whose moments can be given explicitly. In this case, a triangular linear system can be solved for the estimated moments of the mixing distribution <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg"> using the empirical moments of the mixture distribution and the known moments of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">. The former can be estimated from the data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> whereas the latter has to be supplied by the user. Thus, `Hankel.function` contains a function of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> returning the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-th (raw) moment of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">.

As an example, consider a mixture of normal distributions with unknown mean and unit variance. Then <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is the standard normal distribution, and its <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">th moment <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024376.jpg"> is defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024444.jpg">

```{r}
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}
```


#### 3. `Hankel.method = "scale"`

Similarly, example 3.2. in [@hankel, p.285] describes how <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018010.jpg"> can be estimated if the family of component distributions <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023420.jpg"> is given by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024662.jpg">, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is a known probability distribution whose moments can be given explicitly. Likewise, a triangular linear system can be solved for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg">, using the empirical moments of the mixture distribution and the known moments of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">. `Hankel.function` contains a function of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> returning the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-th moment of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">. Note that squares have to be taken everywhere if for some integer <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024922.jpg"> (compare with [@hankel, p.285]).

As an example, consider a mixture of normal distributions with zero mean and unknown variance. Then <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is again the standard normal distribution, and its <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-th moment is defined as above.

Coming back to the overall goal of complexity estimation, the function `nonparamHankel` returns all estimated determinant values corresponding to complexities up to `j.max`, so that the user can pick the lowest <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> generating a sufficiently small determinant. The function allows the inclusion of a penalty term as a function of the sample size `n` and the currently assumed complexity `j` which will be added to the determinant value (by supplying `pen.function`), and/or scaling of the determinants (by setting `scaled  = TRUE`). For scaling, a nonparametric bootstrap is used to calculate the covariance of the estimated determinants, with `B` being the size of the bootstrap sample. The inverse of the square root (i.e. the matrix <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650031611.jpg"> such that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650031632.jpg">, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650031644.jpg"> is the (square) covariance matrix. The procedure uses **expm**'s `sqrtm` [@expm]) of this covariance matrix is then multiplied with the estimated determinant vector to get the scaled determinant vector.

We will initially apply this method to the two already generated datasets of 3-component Poisson and normal mixtures using the penalty <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025161.jpg"> and scaling the determinants as described above.

First, for converting the previously simulated samples from 3-component Poisson and normal mixtures yielding the objects of class `rMix` to objects of class `datMix` one should apply the `RtoDat` function as follows:
```{r rtodat}
MLE.pois <- function(dat) mean(dat)

# create datMix objects:
pois.dM <- RtoDat(poisRMix, theta.bound.list = list(lambda = c(0, Inf)), 
                  MLE.function = MLE.pois, Hankel.method = "explicit",
                  Hankel.function = explicit.pois)


normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list, Hankel.method = "translation",
                     Hankel.function = mom.std.norm)
```

In the case of the scaled version of the method, the penalty should be multiplied by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973529.jpg"> as mentioned earlier. 
```{r nonph}
# define the penalty function:
pen <- function(j, n){
  j * log(n)
}

# apply the nonparamHankel function to the datMix objects:
set.seed(1)
poisdets_sca_pen <- nonparamHankel(pois.dM, j.max = 5, scaled = TRUE, 
                                   B = 1000, pen.function = pen)
normdets_sca_pen <- nonparamHankel(normLoc.dM, j.max = 5, scaled = TRUE, 
                                   B = 1000, pen.function = pen)
```

We can print and plot the results as suggested below.

```{r plotnonph, figures-side, fig.show="hold", out.width="50%"}
# print the results (for the Poisson mixture)
print(poisdets_sca_pen)
# plot results for both mixtures:
par(mar = c(5, 5, 1, 1))
plot(poisdets_sca_pen, main = "3-component Poisson mixture", cex.main = 0.9)
plot(normdets_sca_pen, main = "3-component Normal mixture", cex.main = 0.9)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_art_1.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_art_2.png">

The resulting plots indicate that while theoretically sound, the scaled version of the Hankel method can struggle to correctly identify the number of components in practice.

As the preceding example shows, it can be quite difficult to determine the order estimate from the vector of estimated determinants alone. Thus, the package includes another option of estimating <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> based on Hankel matrices, however, using a more "parametric" approach which goes hand in hand with estimating <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979158.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979228.jpg">. The `paramHankel` procedure initially assumes the mixture to only contain a single component, setting <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022011.jpg">, and then sequentially tests <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025663.jpg"> versus <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025642.jpg"> for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025687.jpg">, until the algorithm terminates. To do so, it determines the MLE for a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-component mixture <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025877.jpg">, generates `B` parametric bootstrap samples of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from the distribution corresponding to <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030786.jpg"> and calculates `B` determinants of the corresponding <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017103.jpg"> Hankel matrices. The null hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026196.jpg"> is rejected and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> increased by 1 if the determinant value based on the original data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> lies outside of the interval <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026269.jpg">, a range specified by the `ql` and `qu` empirical quantiles of the bootstrapped determinants. Otherwise, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is returned as the order estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg">, that is <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg"> is the first order for which the null hypothesis is not rejected. 

`paramHankel.scaled` functions similarly to `paramHankel` with the exception that the bootstrapped determinants are scaled by the empirical standard deviation of the bootstrap sample. To scale the original determinant, `B` nonparametric bootstrap samples of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> are generated from the data, the corresponding determinants are calculated and their empirical standard deviation is used.

Applying `paramHankel.scaled` to the same Poisson and Normal mixtures results in the correct identification of the mixture complexity in both cases as can be seen in the plot:

```{r plotph, figures-side, fig.show="hold", out.width="50%"}
# apply papamHankel.scaled to datMix objects:
set.seed(1)
pois_sca_pen <- paramHankel.scaled(pois.dM)
norm_sca_pen <- paramHankel.scaled(normLoc.dM)
# plot the results for both mixtures:
par(mar=c(5, 5, 1, 1))
plot(pois_sca_pen,)
plot(norm_sca_pen)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_art_1.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_art_2.png">

As another example, consider data generated from a three component geometric mixture, with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026503.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026570.jpg">. 

```{r geomex}
set.seed(1)
# define the geometric mixture:
geomMix <- Mix("geom", discrete=TRUE, w = c(0.1, 0.6, 0.3), prob = c(0.8, 0.2, 0.4))
# generate a random sample from the mixture:
geomRMix <- rMix(1000, obj = geomMix)
# construct the corresponding datMis object:
MLE.geom <- function(dat){
  1/(mean(dat)+1)
}
geom.dM <- RtoDat(geomRMix, Hankel.method = "explicit", 
                  Hankel.function = explicit.geom, 
                  theta.bound.list = list(prob = c(0, 1)), 
                  MLE.function = MLE.geom)
# etimate the number of components using paramHankel function:
(res <- paramHankel(geom.dM, j.max = 5, B = 1000, ql = 0.025, qu = 0.975))
```

Again, `paramHankel` correctly identifies the data as having been generated by a 3-component mixture.

Consider now, as a real-world example, the Children dataset whose content was taken from the Annual Report of the pension fund S.P.P. of 1952. The dataset initially appeared in work of [@thisted] and was subsequently analysed by many authors. It entails data on 4075 widows who recieved pension from the fund, with their number of children being our variable of interest. For example, there are 3062 widows without children, 587 widows with one child, etc. Many authors have noted that this data is not consistent with being a random sample from a Poisson distribution since the number of zeros found in the data is too large. Thisted approached this by fitting a mixture of two populations, one which is always zero and one which follows a Poisson distribution. **mixComp** includes this data stored as a dataframe. Here, we want to investigate 
how the Hankel matrix methods compare when fitting the data to a mixture of Poissons.

The estimation process starts with the construction of the `datMix` object.

```{r childex}
# convert the data to vector:
children.obs <- unlist(children)
# define the MLE function:
MLE.pois <- function(dat) mean(dat)
# construct a datMix object:
children.dM <- datMix(children.obs, dist = "pois", discrete = TRUE, 
                      Hankel.method = "explicit", 
                      Hankel.function = explicit.pois,
                      theta.bound.list = list(lambda = c(0, Inf)), 
                      MLE.function = MLE.pois)
```


First, we check the nonparametric method. We define the penalty <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026704.jpg"> as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026758.jpg"> (by multiplying <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026849.jpg"> by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973529.jpg">). The result suggests that the data comes from a 2-component mixture.

```{r childplotnph, fig.width = 5, fig.height = 4}
# define the penalty:
pen <- function(j, n) j * log(n)
# estimate the number of components:
set.seed(0)
(det_sca_pen <- nonparamHankel(children.dM, j.max = 5, scaled = TRUE, 
                              B = 1000, pen.function = pen))
#plot the results:
plot(det_sca_pen, main = "Non-parametric Hankel method for Children dataset",
     cex.main = 0.9)
```

Next, we check the fit of the parametric version. The printed result of `paramHankel.scaled` shows that this method also suggests 2 to be the number of components, with the first component corresponding to a Poisson distribution with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026963.jpg">. Note that the limit case <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027021.jpg"> results in a point mass at 0, and that this fit therefore nicely lines up with the idea of a component accounting for only the zero observations. The plot shows that this method yields a sensible fit overall.

```{r childplotph, fig.width = 5, fig.height = 4}
set.seed(0)
param_sca <- paramHankel.scaled(children.dM, j.max = 5, B = 1000, ql = 0.025, 
                          qu = 0.975)
plot(param_sca, breaks = 8, ylim = c(0, 0.8))
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_real.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_real.png">

## Section 4. Functions using distances

Unlike the theory on Hankel matrices introduced in Section 3, many theoretical considerations rely on estimates of the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979158.jpg"> and the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979228.jpg">. As mentioned in the introduction, it is assumed that the family of component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015320.jpg"> is known. To embed the subsequent algorithms in a more theoretical framework, consider the parametric family of mixture densities

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027370.jpg">

With <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015320.jpg"> set in advance, elements of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030256.jpg"> can be written as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027564.jpg">

Note that the support of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027645.jpg"> will depend on the support of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027658.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027735.jpg"> (This can be seen by setting <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027784.jpg">) for all <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">. Now take a specific mixture <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028187.jpg">, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028256.jpg">. Clearly, the mixture's complexity is defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028308.jpg">

The above suggests an estimation procedure based on initially finding the "best" possible estimate (in a sense to be determined) <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028366.jpg"> for a given value of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">, in order to compare the thereby specified pdf/pmf 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028449.jpg">

with a non-parametric density/probability mass estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030161.jpg">. As the classes <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030256.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030221.jpg"> are nested, the distance <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028687.jpg"> (to be defined below) between <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030378.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg"> will not increase with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">. Thus, it makes sense to add some penalty term (increasing in <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">) to <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030511.jpg"> and find the first value of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> where the penalized distance for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is smaller than that for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978303.jpg">. Rearranging the terms gives rise to an algorithm starting at <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022011.jpg">, involving some threshold <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030894.jpg"> depending on the penalty, where, if <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is the first integer satisfying

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030958.jpg">

then <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is taken as the estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg">. If the inequality is not fulfilled, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is increased by 1 and the procedure is repeated. Consistency of estimators defined this way has been shown in a number of cases, amongst them those used by the **mixComp** algorithms, and the reader is referred to [@l2; @hell; @hellcont] for the proofs relevant to the results implemented in the package.

The preceding notation was held as broad as possible, since different distance measures <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028687.jpg"> and non-parametric estimators <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg"> can be used. Those relevant to the package are mostly well-known, still, definitions can be found in the Appendix. Three procedures are implemented in **mixComp** based on the foregoing methodology: `L2.disc`, `hellinger.disc` and `hellinger.cont`.

#### 1. `L2.disc`

`L2.disc` employs the squared <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973841.jpg"> distance as the distance measure <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028687.jpg"> and is only to be used for *discrete* mixtures since the nonparametric estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg"> is defined as the empirical probability mass function. In this setting, the "best" estimate 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028366.jpg">

for a given <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> corresponds to

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650032377.jpg">

As the squared <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973841.jpg"> distance might involve an infinite sum (for distributions with infinite support), the user has the option to determine the cut-off value using the `n.inf` argument, which is set to 1000 by default. The parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650032650.jpg"> are obtained analogously. Once both parameter sets have been determined, the difference in their respective squared <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973841.jpg"> distances to <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg"> is compared to a `threshold` (equaling <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030894.jpg"> defined above. The threshold function can be entered directly or one of the predefined thresholds, called `LIC` or `SBC` and given respectively by

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650032536.jpg">

can be used. Note that, if a customized function is to be used, its arguments have to be named `j` and `n`. If the difference in squared distances is smaller than the selected threshold, the algorithm terminates and the true order is estimated as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">, otherwise <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is increased by 1 and the procedure starts over. The reader is invited to consult [@l2] for further details.

#### 2. `hellinger.disc`

This second function presents an alternative estimation procedure for *discrete* mixtures, working much the same as `L2.disc`, however, using a different measure of distance and different thresholds. As the name suggests, it is based on the square of the Hellinger distance, causing the "best" estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650032650.jpg"> for a given <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> to equal

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650047360.jpg">

with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650047443.jpg">. The relevant theory can be found in [@hell]. In accordance with their work, the two predefined thresholds are given by

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650047548.jpg">

(recall that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg"> is the number of dimensions, i.e. <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922033.jpg">). If a customized function is to be used, its arguments have to named `j` and `n` once more, so if the user wants to include the number of component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg">, it has to be entered explicitly. 

#### 3. `hellinger.cont`
 
Unlike the two preceding functions, this procedure is applicable to *continuous* mixture models and uses a kernel density estimator (KDE) as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg">. Its `bandwidth` can be chosen by the user, or the adaptive KDE found in [@adap, p. 1720, equation (2)] may be used by specifying `bandwidth = "adaptive"`. The calculations are based on the continuous version of the squared Hellinger distance, where the "best" estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028366.jpg"> for a given <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> corresponds to

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650047870.jpg">

Since the computational burden of optimizing over an integral to find the "best" weights and component parameters is immense, the algorithm approximates the objective function defined in the previous equation by sampling <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650047975.jpg"> `sample.n` observations <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650048345.jpg"> from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg"> and setting

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650048404.jpg">

Consider again the artificially created samples from the 3-component normal mixture with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650048978.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049033.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049112.jpg"> and the Poisson mixture with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049280.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049280.jpg">.
This procedure `hellinger.cont` uses the same thresholds as `hellinger.disc`.

```{r plothel, figures-side, fig.show="hold", out.width="50%"}
set.seed(0)
h_disc_pois <- hellinger.disc(pois.dM, threshold = "AIC")
h_cont_norm <- hellinger.cont(normLoc.dM, bandwidth = 0.5, sample.n = 5000, 
                      threshold = "AIC")
par(mar = c(5, 5, 1, 1))
plot(h_disc_pois)
plot(h_cont_norm)
```

For a real-world example, refer back to the `faithful` dataset and the corresponding `datMix` object which was created in Section 1. Fitting the distance methods to a continuous density requires a choice of bandwidth. While using the adaptive bandwidth is an option, if the user does not want to do so, it is recommended to use the function `kdensity` from the package **kdensity** [@kdensity] which automatically selects an optimal bandwidth. If the user wants to compare different bandwidth values, it is advisable to look at the plots of the respective kernel density estimates using `kdensity` and to choose one that captures the shape of the data well without fitting to noise.

<img src="https://github.com/yuliadm/mixComp/blob/main/images/bandwidth1.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/bandwidth4.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/bandwidth8.png">

`hellinger.cont` fits a 2-component mixture to the data, which fits the data well and comprises similar parameter estimates to those found in the literature.

```{r faithplothel, fig.width = 5, fig.height = 4}
# estimate the number of components:
library(kdensity)
res <- hellinger.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw,
                      sample.n = 5000, threshold = "AIC")
plot(res)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/hell-cont-norm.png">

At this point, it is worth having a closer look at the thresholds. They each satisfy <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049424.jpg"> as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019043.jpg">, the sole condition the authors require. Now, the consistency proofs for the estimators defined in this Section all rely on the fact that, as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019043.jpg">,

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049583.jpg">

and

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049649.jpg">

where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> is the true complexity (compare with [@l2, p. 4253, Proof of the Theorem], [@hell, p. 4383, Proof] and [@hellcont, p. 1485, Proof of Theorem 1]. If however <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049424.jpg"> faster than <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049788.jpg"> for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017206.jpg">, asymptotically, the decision rule outlined above will always lead to <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> being rejected. Therefore, a second condition should be placed on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030894.jpg">, namely choosing it in accordance with 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049991.jpg">

Neither the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973841.jpg"> Information Criterion  (LIC) nor the Akakike Information Criterion (AIC), nor in the continuous case, the Schwarz Bayesian Criterion (SBC), satisfy this condition, yet they are still part of the package for two reasons. First, since they were used in the original papers, they are included for the sake of completeness and reproducibility of original results. Second, consistency is an asymptotic property, and while the aforementioned thresholds do not fulfill it, they still perform well (and not rarely better than consistent thresholds) for smaller sample sizes. In the example above, the number of components is correctly identified under the non-consistent AIC threshold. Nonetheless, the user will get a warning when using one of non-consistent predefined thresholds.

The preceding example shows that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg"> directly depends on the chosen threshold <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030894.jpg">. While some thresholds can be motivated better than others from a theoretical perspective, the choice will ultimately always remain somewhat arbitrary. It would thus be desirable to have versions of the preceding functions which do not suffer from this drawback. `L2.boot.disc`, `hellinger.boot.disc` and `hellinger.boot.cont` all work similarly to their counterparts, with the exception that the difference in distances is not compared to a predefined threshold but a value generated by a bootstrap procedure.  At every iteration (of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">), the procedure sequentially tests <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025663.jpg"> versus <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025642.jpg"> for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025687.jpg">, using a parametric bootstrap to generate `B` samples of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-component mixture given the previously calculated "best" parameter values <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030786.jpg">. For each of the bootstrap samples, again the "best" estimates corresponding to densities with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978303.jpg"> components are calculated, as well as their difference in distances from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg">. The null hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026196.jpg"> is rejected and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> increased by 1 if the original difference <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049788.jpg"> lies outside of the interval <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026269.jpg">, specified by the `ql` and `qu` empirical quantiles of the bootstrapped differences. Otherwise, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is returned as the order estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg">. 

Since the bootstrap version returns a very similar result to the threshold version on the Old Faithful dataset, we introduce a new example here. Consider the so-called Shakespeare dataset which comprises the number of occurrences of the words that Shakespeare used in his writings. For example, the number of times Shakespeare used a word only once is 14 376, while the number of times the same word occurred exactly 10 times in his writing is 363. The same data have been considered in other papers, see e.g., [@sp68], [@Efron1976], [@CheeWang2016] and [@balabdkulagina]. In the last three papers, the underlying statistical question that the authors wanted to answer is: how many words did Shakespeare actually know? This problem is known under the name of 'species richness' and can be solved using a variety of approaches. The goal is to use the observed frequencies of species, here words, to estimate the unobserved number of words that Shakespeare knew and did not use in his writings. While there is a whole spectrum of methods for estimating species richness, we limit ourselves here to motivate fitting a finite mixture of geometrics to the data. It is known from [@steutel69] that the class of completely monotone probability mass functions defined on the set of non-negative integers, that is the class of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> such that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224248.jpg"> for all integers <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224316.jpg"> coincides with the class of all mixtures of geometrics probability mass functions (here <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224375.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224434.jpg"> for any integer <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224521.jpg">). In [@bdF2019], the monotone complete least squares estimator (LSE) was defined for such class, which is the first non-parametric estimator that was considered for an element in such a family. Complete monotonicity can be defined on any subset of the set of integers of the form <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224708.jpg"> since the change of variable <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224786.jpg"> brings us back to complete monotonicity on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650224850.jpg">. It can be clearly seen that the complete monotone estimator fits very well the empirical estimator of the word occurrences. This result strongly suggests that complete monotonicity is a very appropriate model. In the scope of this paper, we want to explore how fitting a finite mixture of geometric distributions with unknown number of components works for this dataset. This alternative approach is actually inspired by the fact that the complete monotone LSE is itself a finite mixture of geometrics (with a random number of components). Such a result is rather universal and its exact statement can be found in Proposition 2.3 in [@bdF2019]. 


Since we inherently do not observe the number of words Shakespeare did not use, the data start at 1. However, using <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650228576.jpg"> and assuming <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650228421.jpg"> is a geometric mixture with parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650228651.jpg"> leads to the following model for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650228501.jpg">:

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650228680.jpg">

In accordance with the R-function `dgeom`, the parametrization we use for the probability mass function of a geometric distribution means that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650228713.jpg"> is the success probability for the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975606.jpg">-th component, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650228824.jpg">. Building on the above equation, the clear appropriateness of the complete monotone model for the word frequencies in the Shakespeare data can be complemented by a more applied interpretation, its underyling assumption being that words in any language belong to different categories depending on the context in which they are used. As there is a finite number of words, this justifies the appropriateness of fitting a finite mixture model, whose components would correspond to the aforementioned categories. With the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975606.jpg">-th component distribution given by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650229293.jpg">, this expression can be seen as the probability of a word belonging to category <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975606.jpg"> not appearing (in some new work) after having previously been used <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650229825.jpg"> times.

The `datMix` object corresponding to the Shakespeare dataset is generated as follows: 


```
shakespeare.obs <- unlist(shakespeare) - 1
# define the MLE function:
MLE.geom <- function(dat) 1 / (mean(dat) + 1)

Shakespeare.dM <- datMix(shakespeare.obs, dist = "geom", discrete = TRUE, 
MLE.function = MLE.geom, theta.bound.list = list(prob = c(0, 1)))

# estimate the number of components and plot the results:
set.seed(0)
res <- hellinger.boot.disc(Shakespeare.dM, B = 50, ql = 0.025, qu = 0.975)
plot(res)
```

<img src="https://github.com/yuliadm/mixComp/blob/main/images/hell-boot-geom.png">

`hellinger.boot.disc` estimates that the data comes from a 3-component geometric mixture (thus clustering the english words Shakespeare used into three categories).



## Section 5. Functions using LRTS

As a third option of estimating the mixture complexity, the **mixComp** package provides an algorithm based on the likelihood ratio test statistic (LRTS), sequentially testing <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025663.jpg"> versus <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025642.jpg"> for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025687.jpg">, until the algorithm terminates. As noted in Section 1, it is not possible to use the generalized likelihood ratio test in the context of mixture models directly as the standard way of obtaining the asymptotic distribution of the LRTS under the null hypothesis cannot be applied. However, one can get estimates of the test's critical values by employing a bootstrap procedure.

Making use of this approach, the function `mix.lrt` iteratively increases the assumed order <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> and finds the MLE for both, the density of a mixture with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978303.jpg"> components, giving <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028366.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650051363.jpg">. It then calculates the corresponding LRTS, defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650051485.jpg">

with

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650051541.jpg">

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650051597.jpg"> being the likelihood function given the data <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">.

Next, a parametric bootstrap is used to generate `B` samples of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-component mixture given the previously calculated MLE <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030786.jpg">. For each of the bootstrap samples, the MLEs corresponding to densities of mixtures with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978303.jpg"> components are calculated, as well as the LRTS. The null hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026196.jpg"> is rejected and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> increased by 1 if the LRTS based on the original data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> is larger than the chosen `quantile` of its bootstrapped counterparts. Otherwise, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is returned as the order estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg">. For further details, the reader is referred to [@lrt].

For the two artificial datasets, this method estimates 3-component mixtures with very similar parameters to the distance methods, so we go straight to a real-world example. Consider the Acidity dataset which comprises measurements of the acid neutralizing capacity (ANC) taken from 155 lakes in North-Central Wisconsin. The ANC indicates a lakes' capability to absorb acid, with low values potentially leading to a loss of biological resources. This dataset has been analysed as a mixture of normal distributions on the log scale by [@acidity1], [@acidity2] and [@acidity3]. While the former papers suggest the number of components to equal 2 (with 3 also being considered), the latter estimates <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> to lie between 3 and 5. The `mix.lrt` method agrees with [@acidity1] and [@acidity2], returning a 2-component mixture with reasonable estimates for the component weights and parameters.

```{r lrtacid, fig.width = 5, fig.height = 4, results='hide', message=FALSE, warning=FALSE}
acidity.obs <- unlist(acidity)

acidity.dM <- datMix(acidity.obs, dist = "norm", discrete = FALSE, 
                     MLE.function = MLE.norm.list, 
                     theta.bound.list = norm.bound.list)

set.seed(0)
res <- mix.lrt(acidity.dM, B = 50, quantile = 0.95)
plot(res)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/lrt-norm.png">

## Section 6. Non-standard mixtures

In all preceding examples, the families of component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976784.jpg"> belonged to one of the "standard" probability distributions included in the **stats** package, which provides the density/mass function, cumulative distribution function, quantile function and random variate generation for selected distributions. The function names are of the form `dxxx`, `pxxx`, `qxxx` and `rxxx` respectively. With some additional effort, it is possible to use the **mixComp** package on "non-standard" distributions -- the user merely has to provide functions evaluating the density and generating random numbers for the component distribution. In accordance with **R** conventions, the user-generated function `dxxx` has to take `x` and the distribution parameters as input and returns the value of the density function specified by the parameters at the point `x`. Likewise, `rxxx` requires `n` and the distribution parameters as input and returns `n` random numbers based on the distribution specified by the aforementioned parameters.

As an example, consider a sample <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> from a 3-component mixture of normals with means equal to 10, 11 and 13. Assume that the standard deviation of all components is known to be 0.5, yet the number of components and their means are unknown. Then each of the components follows a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650052707.jpg"> distribution, which shall be called `norm0.5`. The first step is always that of creating the `dxxx` and `rxxx` functions, since they will be called by the **mixComp** functions.

The following example creates the `Mix` and `rMix` objects
based on the density of a normal mixture with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650048978.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650052857.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650052889.jpg"> and plots the obtained mixture density and the corresponding random sample. 

```{r, figures-side, fig.show="hold", out.width="50%", results='hide', message=FALSE, warning=FALSE}
dnorm0.5 <- function(x, mean){
  dnorm(x, mean = mean,  sd = 0.5)
}
rnorm0.5 <- function(n, mean){
  rnorm(n, mean = mean,  sd = 0.5)
}
## create objects `Mix` and `rMix`:
set.seed(1)
norm0.5Mix <- Mix("norm0.5", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 11, 13))
norm0.5RMix <- rMix(1000, obj = norm0.5Mix)
## plot the results:
plot(norm0.5Mix)
plot(norm0.5RMix)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/norm0.5Mix.png">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/norm0.5RMix.png">

Below we will estimate of the mixture density using `mix.lrt` given a sample from the 3-component normal mixture with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650048978.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650052857.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650052889.jpg">.

We start by creating all necessary inputs:
```{r}
norm0.5.list <- vector(mode = "list", length = 1)
names(norm0.5.list) <- c("mean")
norm0.5.list$mean <- c(-Inf, Inf)

MLE.norm0.5 <- function(dat) mean(dat)

norm0.5.dM <- RtoDat(norm0.5RMix, theta.bound.list = norm0.5.list,
                     MLE.function = MLE.norm0.5)
```


Finally, the **mixComp** procedures can be used on the `datMix` object as usual. The results can be printed and plotted using `print` and `plot` functions.
```{r, results='hide', message=FALSE, warning=FALSE}
set.seed(1)
res <- mix.lrt(norm0.5.dM, B = 50, quantile = 0.95)
```


```{r, fig.width = 5, fig.height = 4}
print(res)
plot(res)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/images/est-norm0.5.png">

## Section 7. Computational nuance for mixComp functions using the solnp() solver

Several functions in the **mixComp** package (namely, `nonparamHankel`, `paramHankel`, `hellinger.cont`, `hellinger.disc`, `L2.disc`, `mix.lrt`) make use of `solnp()`  function (**Rsolnp** library), which is a solver for general nonlinear programming problems. The above mentioned **mixComp** functions attempt to generate good starting values by clustering the data via `clara` function prior to applying  `solnp()`, which lead to convergence of the algorithm in most of the cases. However when running multiple simulations, `solnp()` might not converge for particular initial values with default control values. This may happen when very few observations are assigned to some of the clusters have, in which case the solver can get "stuck", not even resulting in bad exit status (codes 1 and 2 in returned by `solnp()` convergence value). This issue can be overcome by specifying the control parameters in the functions using the `solnp()` solver (e.g. by defining `L2.boot.disc(geom.dM, j.max = 5, B = 500, ql = 0.025, qu = 0.975, control = list("rho" = 0.1, tol = 0.0000001, trace = 0))`, see description of the `solnp()` function for further details on setting the control parameters) or by setting a time limit for the function execution and going to the next iteration whenever the time limit is exceeded. 


## Section 8. Summary and discussion

In this paper we presented the R package **mixComp**, a collection of routines developed to estimate a mixture's complexity from a data sample <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">. Moreover, it provides the possibility of generating artificial data from specified mixtures as well as proper visualization tools for the complexity estimation result, plotting either the successive determinant values or the final fitted mixture. If estimates of the component weights and parameters are sought, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg"> can be passed to one of the many **R** packages specialized in their calculation, or one of the **mixComp** functions returning weight and parameter estimates can be used. However, it should be noted that the theory on which this package is based concentrates on showing consistency for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg"> -- other estimates obtained in the process are merely by-products. A primary goal of the package was to make it extendible, meaning that it can be utilized on mixtures of less well-known distributions, as long as sufficient information on the density and random variate generation is provided. We hope that this package will be useful for practitioners in the many areas where mixture models are applicable.

## Computational details

All computations and graphics in this paper have been done using **R** version 4.0.0 with the packages **boot** 1.3-24, **cluster** 2.1.0, **expm** 0.999-4, **matrixcalc** 1.0-3, **Rsolnp** 1.16 and **kdensity** 1.0.1. **R** itself
and all packages used are available from the Comprehensive **R** Archive Network (CRAN) at https://CRAN.R-project.org/.


## Appendix 
### Distance and non-parametric estimator definitions

Let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650054486.jpg"> be an i.i.d. sample of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from some continuous distribution with unknown density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027645.jpg">. Its *kernel densty estimator* is defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650054693.jpg">

with kernel <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650054779.jpg"> and bandwidth <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650054831.jpg">.

As an extension, [@adap] defined the *adaptive kernel density estimator* as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055002.jpg">	

with kernel <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650054779.jpg">, bandwidth <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055259.jpg"> and weights (positive and summing to 1) <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055316.jpg">.
	
Let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650054486.jpg"> be an i.i.d. sample of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from some discrete distribution with unknown probability mass function <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027645.jpg">. Its *empirical probability mass function* is defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650054952.jpg">

Let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027658.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027645.jpg"> be two probability mass functions defined on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055651.jpg">.
The *squared* <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973841.jpg"> *distance* between <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027658.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027645.jpg"> is given by

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055147.jpg">	

Let  <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027658.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027645.jpg"> be two density or probability mass functions defined on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055651.jpg">.
The *squared Hellinger distance* between  <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027658.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650027645.jpg"> is given by

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055464.jpg">	

in the discrete case and by

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650055550.jpg">	

in the continuous case.

### Further results on the real-world datasets

The Table below showcases the results of all complexity estimation procedures when applied to the four real-world datasets that were discussed in this paper. The table also includes the results for the Shakespeare dataset borrowed from [@shakespeare]. It was not mentioned in the vignette but is included in the **mixComp** package. The following settings were used to calculate these results (default setting were used unless indicated otherwise):

* `set.seed(1)` was used for all complexity calculations,
* for `hellinger.cont` and `hellinger.boot.cont` the bandwidths suggested by **kdensity** were used,
* for `nonparamHankel`, we used <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025161.jpg"> and `j.max = 6`,
* for the bootstrapped distance and the LRTS methods, we used `B = 50`.

#### Table A: Results for all **mixComp** methods used on real-world data
|  Method                   | Children | Old Faithful | Shakespeare | Acidity |
|:-------------------------:|:--------:|:------------:|:-----------:|:-------:|
| `nonparamHankel`          |    2     |      x       |      2      |    x    |
| `nonparamHankel`(scaled)  |    2     |      x       |      2      |    x    |
| `paramHankel`             |    2     |      x       |      2      |    x    |
| `paramHankel.scaled`      |    2     |      x       |      2      |    x    |
| `L2.disc`                 |    2     |      x       |      3      |    x    |
| `L2.boot.disc`            |    2     |      x       |      4      |    x    |
| `hellinger.disc`          |    2     |      x       |      3      |    x    |
| `hellinger.boot.disc`     |    2     |      x       |      3      |    x    |
| `hellinger.cont`          |    x     |      2       |      x      |    2    |
| `hellinger.boot.cont`     |    x     |      2       |      x      |    2    |
| `mix.lrt`                 |    2     |      2       |      3      |    2    |

# References

<a id="1">[1]</a> 
Azzalini A, Bowman A (1990). “A Look at Some Data on the Old Faithful Geyser.” Journal
of the Royal Statistical Society C, 39(3), 357–365. doi:10.2307/2347385. URL https:
//www.jstor.org/stable/2347385.

<a id="2">[2]</a> 
Balabdaoui F, de Fournas-Labrosse G (2020). “Least Squares Estimation of a Completely
Monotone pm: From Analysis to Statistics.” JSPI, 204, 55–71.

<a id="3">[3]</a> 
Balabdaoui F, Kulagina Y (2020). “Completely Monotone Distributions: Mixing, Approximation 
and Estimation of Number of Species.” Computational Statistics & Data Analysis,
150, 107014, 26. ISSN 0167-9473. doi:10.1016/j.csda.2020.107014.

<a id="4">[4]</a> 
Benaglia T, Chauveau D, Hunter DR, Young D (2009). “mixtools: An R Package for
Analyzing Finite Mixture Models.” Journal of Statistical Software, 32(6), 1–29. URL
http://www.jstatsoft.org/v32/i06/.

<a id="5">[5]</a> 
Chee CS, Wang Y (2016). “Nonparametric Estimation of Species Richness Using Discrete
k-Monotone Distributions.” Computational Statistics & Data Analysis, 93, 107–118. ISSN
0167-9473.

<a id="6">[6]</a> 
Chen J, Kalbfleisch J (1996). “Penalized Minimum-Distance Estimates in Finite Mixture
Models.” Canadian Journal of Statistics, 24(2), 167–175.

<a id="7">[7]</a> 
Crawford SL (1994). “An Application of the Laplace Method to Finite Mixture Distributions.”
Journal of the American Statistical Association, 89(425), 259–267. ISSN 01621459. URL
http://www.jstor.org/stable/2291222.

<a id="8">[8]</a> 
Crawford SL, DeGroot MH, Kadane JB, Small MJ (1992). “Modeling Lake-Chemistry Distri-
butions: Approximate Bayesian Methods for Estimating a Finite-Mixture Model.” Techno-
metrics, 34(4), 441–453. ISSN 00401706. URL http://www.jstor.org/stable/1268943.

<a id="9">[9]</a> 
Cutler A, Cordero-Braña OI (1996). “Minimum Hellinger Distance Estimation for Finite
Mixture Models.” Journal of the American Statistical Association, 91(436), 1716–1723.
ISSN 01621459. doi:10.2307/2291601. URL http://www.jstor.org/stable/2291601.

<a id="10">[10]</a> 
Dacunha-Castelle D, Gassiat E (1997). “The Estimation of the Order of a Mixture Model.”
Bernoulli, 3(3), 279–299. doi:10.2307/3318593. URL https://www.jstor.org/stable/
3318593.

<a id="11">[11]</a> 
Dempster AP, Laird NM, Rubin DB (1977). “Maximum Likelihood from Incomplete Data
via the EM Algorithm.” Journal of the Royal Statistical Society B, 39(1), 1–38. ISSN
0035-9246. With discussion, URL http://www.jstor.org/stable/2984875.

<a id="12">[12]</a> 
Efron B, Thisted R (1976). “Estimating the Number of Unseen Species: How Many Words
Did Shakespeare Know?” Biometrka, 63, 435–447. ISSN 0006-341X.

<a id="13">[13]</a> 
Figueiredo MAT, Jain AK (2002). “Unsupervised Learning of Finite Mixture Models.” IEEE
Transactions on Pattern Analysis and Machine Intelligence, 24(3), 381–396.

<a id="14">[14]</a> 
Ghalanos A, Theussl S (2015). Rsolnp: General Non-linear Optimization Using Augmented
Lagrange Multiplier Method. R package version 1.16.

<a id="15">[15]</a> 
Goulet V, Dutang C, Maechler M, Firth D, Shapira M, Stadelmann M (2019). expm: Matrix
Exponential, Log, ’etc’. R package version 0.999-4, URL https://CRAN.R-project.org/
package=expm.

<a id="16">[16]</a> 
Grün B, Leisch F (2007). “Fitting Finite Mixtures of Generalized Linear Regressions in
R.” Computational Statistics & Data Analysis, 51(11), 5247–5252. doi:10.1016/j.csda.
2006.08.014.

<a id="17">[17]</a> 
Grün B, Leisch F (2008). “FlexMix Version 2: Finite Mixtures with Concomitant Variables
and Varying and Constant Parameters.” Journal of Statistical Software, 28(4), 1–35. doi:
10.18637/jss.v028.i04. URL http://www.jstatsoft.org/v28/i04/.

<a id="18">[18]</a> 
Härdle W (1991). Smoothing Techniques : With Implementation in S. Springer-Verlag New
York, New York, NY, USA. doi:10.1007/978-1-4612-4432-5.

<a id="19">[19]</a> 
Leisch F (2004). “FlexMix: A General Framework for Finite Mixture Models and Latent
Class Regression in R.” Journal of Statistical Software, 11(8), 1–18. doi:10.18637/jss.
v011.i08. URL http://www.jstatsoft.org/v11/i08/.

<a id="20">[20]</a> 
Lindsay BG (1983a). “The Geometry of Mixture Likelihoods: A General Theory.” The Annals
of Statistics, 11(1), 86–94. ISSN 0090-5364. doi:10.1214/aos/1176346059.30

<a id="21">[21]</a> 
Lindsay BG (1983b). “The Geometry of Mixture Likelihoods: The Exponential Family.” The
Annals of Statistics, 11(3), 783–792. ISSN 0090-5364. doi:10.1214/aos/1176346245.

<a id="22">[22]</a> 
Macdonald P, Du J (2018). mixdist: Finite Mixture Distribution Models. R package version
0.5-5, URL https://CRAN.R-project.org/package=mixdist.

<a id="23">[23]</a> 
Maechler M, Rousseeuw P, Struyf A, Hubert M, Hornik K (2019). cluster: Cluster Analysis
Basics and Extensions. R package version 2.1.0 — For new features, see the ’Changelog’
file (in the package source).

<a id="24">[24]</a> 
McLachlan G, Peel D (2000). Finite Mixture Models. Wiley Series in Probability and Statis-
tics: Applied Probability and Statistics. John Wiley & Sons. ISBN 0-471-00626-2. doi:
10.1002/0471721182.

<a id="25">[25]</a> 
Melnykov V, Chen WC, Maitra R (2012). “MixSim: An R Package for Simulating Data
to Study Performance of Clustering Algorithms.” Journal of Statistical Software, 51(12),
1–25. URL http://www.jstatsoft.org/v51/i12/.

<a id="26">[26]</a> 
Miller JW, Harrison MT (2018). “Mixture Models with a Prior on the Number of Compo-
nents.” Journal of the American Statistical Association, 113(521), 340–356.

<a id="27">[27]</a> 
Müller L (2020). Clustering with Mixtures: A Comparative Study. Master’s thesis, ETH.

<a id="28">[28]</a> 
Moss J, Tveten M (2019). kdensity: Kernel Density Estimation with Parametric Starts
and Asymmetric Kernels. R package version 1.0.1, URL https://CRAN.R-project.org/
package=kdensity.

<a id="29">[29]</a> 
R Core Team (2020). R: A Language and Environment for Statistical Computing. R Foun-
dation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

<a id="30">[30]</a> 
Richardson S, Green PJ (1997). “On Bayesian Analysis of Mixtures with an Unknown Number
of Components (with Discussion).” Journal of the Royal Statistical Society B, 59(4), 731–
792. doi:10.1111/1467-9868.00095. URL https://rss.onlinelibrary.wiley.com/
doi/abs/10.1111/1467-9868.00095.

<a id="31">[31]</a> 
Scrucca L, Fop M, Murphy TB, Raftery AE (2016). “mclust 5: Clustering, Classification and
Density Estimation Using Gaussian Finite Mixture Models.” The R Journal, 8(1), 289–317.
URL 10.32614/RJ-2016-021.

<a id="32">[32]</a> 
Spevack M (1968). “A Complete and Systematic Concordance to the Works of Shakespeare.”
Vols. 1–6. George Olms, Hildesheim.

<a id="33">[33]</a> 
Steutel FW (1969). “Note on Completely Monotone Densities.” The Annals of Mathematical
Statistics, 40, 1130–1131. ISSN 0003-4851.

<a id="34">[34]</a> 
Teicher H (1963). “Identifiability of Finite Mixtures.” The Annals of Mathematical Statistics,
34, 1265–1269. ISSN 0003-4851. doi:10.1214/aoms/1177703862.

<a id="35">[35]</a> 
Thisted RA (1988). Elements of Statistical Computing: Numerical Computation. Chapman
& Hall, Ltd., GBR. ISBN 0412013711.

<a id="36">[36]</a> 
Titterington DM, Smith AFM, Makov UE (1985). Statistical Analysis of Finite Mixture
Distributions. Wiley Series in Probability and Mathematical Statistics: Applied Probability
and Statistics. John Wiley & Sons. ISBN 0-471-90763-4.

<a id="37">[37]</a> 
Umashanger T, Sriram T (2009). “L2E Estimation of Mixture Complexity for Count Data.”
Computational Statistics & Data Analysis, 53(12), 4243 – 4254. ISSN 0167-9473. doi:10.
1016/j.csda.2009.05.013. URL http://www.sciencedirect.com/science/article/
pii/S0167947309002023.

<a id="38">[38]</a> 
Wang HX, Luo B, Zhang QB, Wei S (2004). “Estimation for the Number of Components
in a Mixture Model Using Stepwise Split-and-Merge EM Algorithm.” Pattern Recognition
Letters, 25(16), 1799–1809. ISSN 0167-8655.

<a id="39">[39]</a> 
Woo MJ, Sriram T (2007). “Robust Estimation of Mixture Complexity for Count Data.”
Computational Statistics & Data Analysis, 51(9), 4379 – 4392. ISSN 0167-9473. doi:10.
1016/j.csda.2006.06.006. URL http://www.sciencedirect.com/science/article/
pii/S0167947306001964.

<a id="40">[40]</a> 
Woo MJ, Sriram TN (2006). “Robust Estimation of Mixture Complexity.” Journal of the 
American Statistical Association, 101(476), 1475–1486. doi:10.1198/016214506000000555.

<a id="41">[41]</a> 
Xekalaki E, Karlis D (1999). “On Testing for the Number of Components in a Mixed Poisson
Model.” The Annals of the Institute of Statistical Mathematics, 51, 149–162. doi:10.1023/
A:1003839420071.

<a id="42">[42]</a> 
Ye Y (1987). Interior Algorithms for Linear, Quadratic, and Linearly Constrained Non-Linear
Programming. Ph.D. thesis, Department of ESS, Stanford University.


