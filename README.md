---
title: 'Estimating the Complexity of a Finite Mixture with the \proglang{R} Package \pkg{mixComp}'
tags:
  - R
  - mixture complexity estimation
  - mixture models
  - Hankel matrix
  - minimum distance
  - likelihood-rato test
authors:
  - name: Anja Weigel^[author] # note this makes a footnote 
    affiliation: 2  
  - name: Fadoua Balabdaoui^[co-author] # note this makes a footnote 
    affiliation: 1 
  - name: Yulia Kulagina^[co-author, maintainer] # note this makes a footnote 
    affiliation: 1
  - name: Lilian Mueller^[contributor]
    affiliation: 2
  - name: Martin Maechler^[contributor]
    (package 'nor1mix' as model, <https://orcid.org/0000-0002-8685-9910>)
    affiliation: 1  
affiliations:
  - name: ETH Zurich, Seminar for Statistics, Switzerland
   index: 1
  - name: ETH Zurich, Switzerland
   index: 2
date: 17 April 2022
bibliography: refs.bib

---

# Summary

The **mixComp** package provides a number of methods for obtaining a consistent estimate of the complexity of a finite mixture (the focus is made on the univariate case). The considered approaches can be loosely grouped into three categories:
<ul>
  <li> methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; </li>
  <li> methods based on penalized minimum distance between the unknown probability density and a consistent estimator thereof. The distances considered in this survey are the Hellinger and the $L_2$-distances; </li>
  <li> likelihood ratio test (LRT) - based techniques. </li>
</ul>
While not the primary goal, most methods simultaneously estimate the component weights and parameters. In this document, we give a brief overview of the methodology, and demonstrate the package's functionality in both real world examples and synthetically generated data. Moreover, we show how the package can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution.

# Statement of need

Mixture models occur in numerous settings including random
and fixed effects models, clustering, deconvolution, empirical Bayes prob-
lems and many others. In particular, they are often used to model the data
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

Mixture models have been used extensively in statistical applications and therefore have attracted a lot of attention from both theoretical and computational perspectives. Although the list of works on mixture models is too long to make an exhaustive inventory, we can refer to the following important papers and books: [@Teicher63], [@LindsayI], [@LindsayII], [@Titterington] and [@McLachlan].

The popularity of such models stems from the fact that they allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. To account for such heterogeneity, the (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method. This includes the well-known nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were well studied in the seminal work of [@LindsayI], [@LindsayII]. One remarkable property of the NPMLE of the mixing distribution is that it is, under some simple conditions, a discrete distribution function with at most $k$ number of jumps, where $k$ is the number of distinct observations in the random sample drawn from the mixture. This interesting feature is one reason, among others, for considering the smaller class of finite mixture models, i.e., mixture models with a discrete mixing distribution with a finite number of jumps. The model has the following simple interpretation: the population under study is assumed to consist of several homogeneous subpopulations. These subpopulations, typically referred to as the mixture's components, usually have a specific meaning depending on the problem at hand. In some very simple situations, the number of components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric maximum likelihood estimator (MLE) is known to occur at the fast rate $\sqrt{n}$ (under some regularity conditions). Also, the well-known expectation-maximization (EM) algorithm can be used to find the MLE of all the unknown parameters; see for example [@Dempster]. However, in many statistical applications such knowledge is rarely available and the number of components has to be estimated from the data. Although the mixture is still finite and the distribution of each component is assumed to belong to some parametric family, the estimation framework in this case is much harder than in the fully parametric one, where the number of components is known. In this paper, the terms *order*, *complexity*  and *number of components*  will  be used interchangeably to refer to this unknown number. The main goal of the package **mixComp** is to estimate the unknown complexity using several methods known from the statistical literature. These methods, which are discussed below in more detail, all come with theoretical guarantees for consistency as the sample size gets larger. Of course, consistency in this case means that an estimator is able to exactly recover the unknown complexity for large sample sizes. As expected, the performance of the methods varies according to the underlying mixture distribution and the sample size. This will be illustrated below through several synthetic as well as real datasets.

To describe the estimation problem, we start with some formal notation.  A distribution $F$ is called a *finite mixture* if its density (we write density throughout and keep in mind that it may be taken with respect to the Lebesgue or the counting measure) is of the form

$$f(x) = \sum_{i=1}^p w_i g_i(x),$$

where $p \in \mathbb{N}$ is the mixture complexity, $(w_1, \dots w_p : \sum_{i=1}^p w_i = 1$, $w_i \geq 0,$ for $i=1,\dots,p)$ are the component weights and the density $g_i(x)$ is the $i$-th component of the mixture. As the scope of **mixComp** is limited to mixtures where the family of the component distributions is known, we replace $g_i(x)$ by a parametric density $g(x; \theta_i)$ indexed by the (possibly multivariate, say $d$-dimensional) parameter $\theta_i$ in the parameter space $\Theta \subseteq \mathbb{R}^d$.

Given some complexity $j$, the two relevant parameter spaces can therefore be defined as

$$\Theta_j = \{\theta_1 \dots \theta_j: \theta_i \in \Theta \subseteq \mathbb{R}^d, \text{ for } i = 1,\dots,j\}$$

and

$$W_j = \{w_1, \dots, w_j: \sum_{i=1}^j w_i = 1, w_i \geq 0, \text{ for } i = 1,\dots,j\}.$$

Throughout this document, it is assumed that the family of the component densities $\{g(x; \theta):\theta \in \Theta\}$ is known, but the component parameters $ \textbf{\theta}= (\theta_1, \dots, \theta_p) \in \Theta_p$, the component weights $\textbf{w} = (w_1, \dots, w_p) \in W_p$ and the mixture complexity $p \in \mathbb{N}$ are unknown, with $p$ being the parameter of interest. Assume now that $F$ is a finite mixture distribution with density $f(x) = \sum_{i=1}^p w_i g(x; \theta_i)$ and $\textbf{X} = \{X_1, \dots, X_n\}$ is an i.i.d. sample of size $n$ from $F$. The **mixComp** package aims to estimate the smallest such $p$ on the basis of $\textbf{X}$, either on its own or by simultaneously estimating the weights $w_i$ and the component parameters $\theta_i$, $i \in 1, \dots, p$.

In this setup, it seems natural to test for the number of components by comparing two consecutive models. Traditionally, the problem of choosing between nested models may be approached by applying the generalized likelihood ratio test and referring to the $\chi^2_r$ distribution to assess significance, where $r$ is given by the number of constraints imposed on the alternative hypothesis $H_1$ to arrive at the null hypothesis $H_0$. However, in the context of mixture models, there are several issues hindering application of this classical theory. One of them is that there is no unique way of obtaining $H_0$ from $H_1$. As an example, the two null hypotheses $H_0: w_{j+1} = 0$ and $H_0: \theta_{j+1} = \theta_{1}$ both yield the smaller model, showcasing the difficulties of applying the classical asymptotic theory of the likelihood ratio. This problem has been studied extensively in the literature and numerous alternative approaches to mixture complexity estimation have been suggested, laying the theoretical foundation for the subsequently described algorithms.

This document discusses various categories of functions found in the **mixComp** package, ranging from methods based on Hankel matrices, to techniques based upon distances between densities and likelihood ratio tests. The examples provided in the first sections all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be performed using the functions from the **stats** package. The last section illustrates how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and valuating the density thereof.

Two main features distinguish this package from other mixture-related **R** [@R] packages: Firstly, it is focused on the estimation of the complexity rather than the component weights and parameters. While these are often estimated as a by-product, all methods contained in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture of interest. Secondly, it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** [see @mixtools] and **flexmix** [see @flexmix1; @flexmix2; @flexmix3]  should both be mentioned at this point: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm. Notably, it also contains routines for selecting the number of components based on information criteria and parametric bootstrapping of the likelihood ratio test statistic values. However, they are limited to multinomial and (a variety of) normal mixtures as well as mixtures-of-regressions. Second, while **flexmix** was developed to deal with mixtures-of-regression, it sets itself apart from other packages by its extensibility, a design principle that we also aimed for when creating the  **mixComp** package. Other widely used packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], which is used for grouped conditional data. Interested readers can find a comprehensive list of mixture-related packages on the CRAN Task View: Cluster Analysis and Finite Mixture Models website.

Before moving to the description of the different methods implemented in **mixComp** we would like to briefly mention pther theoretical work on the estimation of mixture complexity not currently included in the package. [@chen] propose a method that is reminiscent of the ones described in Section 4. The main difference is that the authors consider distribution functions instead densities, i.e. they consider minimizing a penalized distance between the distribution function of the mixture and the empirical distribution function. The approach of [@figueiredo] is based on a minimum message length-like criterion, however, their method struggles to deal with mixtures with very different weights. [@xian] propose a procedure based on alternating between splitting and merging the components in an EM-algorithm. This algorithm requires selecting two thresholds, the choice of which is somewhat unclear when working with a specific dataset. [@miller] follow a Bayesian approach, taking the usual finite mixture model with Dirichlet weights and putting a prior distribution on the unknown number of components. 

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


## Section 3. Functions using Hankel matrices

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
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_art_1.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_art_2.png" />
</p>

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
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_art_1.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_art_2.png" />
</p>

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
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/images/np_real.png" />
<img src="https://github.com/yuliadm/mixComp/blob/main/images/p_real.png" />
</p>

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
