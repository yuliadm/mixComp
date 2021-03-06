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
affiliations:
 - name: ETH Zurich, Seminar for Statistics, Switzerland
   index: 1
 - name: ETH Zurich, Switzerland
   index: 2
date: 19 April 2022
bibliography: paper/refs.bib
---

# Installation

To install from `CRAN`, use:
``` r
install.packages("mixComp")
```
For installing the development version from GitHub:
``` r
# install.packages("devtools")
devtools::install_github("yuliadm/mixComp")
```

# Community Guidelines

Contributions in the form of feedback, comments, code, and bug reports are welcome. To contribute code or report a bug, please 

  - fork the source code, modify, and issue a pull request or
   
  - file a Github issue: https://github.com/yuliadm/mixComp/issues.

Any other contribution questions and requests for support can be directed to the package maintainer Yulia Kulagina (yulia.kulagina@stat.math.ethz.ch).

# **mixComp** for mixture complexity estimation 

Mixture models have been used extensively in statistical applications and therefore have attracted a lot of attention from both theoretical and computational perspectives. Although the list of works on mixture models is too long to make an exhaustive inventory, we can refer to the following important papers and books: [[21]](#21), [[12]](#12),  [[13]](#13), [[23]](#23) and [[15]](#15).

The popularity of such models stems from the fact that they allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. To account for such heterogeneity, the (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method. This includes the well-known nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were well studied in the seminal work of [[12]](#12), [[13]](#13). One remarkable property of the NPMLE of the mixing distribution is that it is, under some simple conditions, a discrete distribution function with at most *k* number of jumps, where *k* is the number of distinct observations in the random sample drawn from the mixture. This interesting feature is one reason, among others, for considering the smaller class of finite mixture models, i.e., mixture models with a discrete mixing distribution with a finite number of jumps. The model has the following simple interpretation: the population under study is assumed to consist of several homogeneous subpopulations. These subpopulations, typically referred to as the mixture's components, usually have a specific meaning depending on the problem at hand. In some very simple situations, the number of components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric maximum likelihood estimator (MLE) is known. Also, the well-known expectation-maximization (EM) algorithm can be used to find the MLE of all the unknown parameters; see for example [[6]](#6). However, in many statistical applications such knowledge is rarely available and the number of components has to be estimated from the data. Although the mixture is still finite and the distribution of each component is assumed to belong to some parametric family, the estimation framework in this case is much harder than in the fully parametric one, where the number of components is known. In this paper, the terms *order*, *complexity*  and *number of components*  will  be used interchangeably to refer to this unknown number. The main goal of the package **mixComp** is to estimate the unknown complexity using several methods known from the statistical literature. These methods, which are discussed below in more detail, all come with theoretical guarantees for consistency as the sample size gets larger. Of course, consistency in this case means that an estimator is able to exactly recover the unknown complexity for large sample sizes. As expected, the performance of the methods varies according to the underlying mixture distribution and the sample size.

The methods that were included in the package can be roughly devided into three categories: methods based on Hankel matrices, following the theory as described in [[5]](#5) and selected because of the fact that computation of the mixture parameters is not required, a method based on the likelihood ratio test statistic (LRTS) following [[28]](#28) since a likelihood ratio test seems like a natural approach in this setting and methods employing minimum distance calculations based on several works and included as a computationally more efficient alternative to the LRTS method for certain distributions and distances; see [[26]](#26), [[27]](#27), [[24]](#24), [[4]](#4). For example, when the distance is taken to be the Hellinger distance, such an approach is especially fast for discrete distributions. For a more fluid reading, the relevant theory will be portrayed at the beginning of each of the respective sections. The examples depicted in these first chapters all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be done by functions available from the **stats** package. The last chapter showcases how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and evaluating the density thereof.

The mathematical details of each approached can be found in the `documentation.md` file. 

Two main features distinguish this package from other mixture-related **R** [[19]](#19) packages: Firstly, it is focused on the estimation of the complexity rather than the component weights and parameters. While these are often estimated as a by-product, all methods contained in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture of interest. Secondly, it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** (see [[2]](#2)) and **flexmix** (see [[8]](#8), [[9]](#9), [[11]](#11)), should both be mentioned at this point: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm. Notably, it also contains routines for selecting the number of components based on information criteria and parametric bootstrapping of the likelihood ratio test statistic values. However, they are limited to multinomial and (a variety of) normal mixtures as well as mixtures-of-regressions. Second, while **flexmix** was developed to deal with mixtures-of-regression, it sets itself apart from other packages by its extensibility, a design principle that we also aimed for when creating the  **mixComp** package. Other widely used packages dealing with mixture models are **mclust** [[20]](#20), which fits mixtures of Gaussians using the EM algorithm, **MixSim** [[16]](#16), which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [[14]](#14), which is used for grouped conditional data. Interested readers can find a comprehensive list of mixture-related packages on the CRAN Task View: Cluster Analysis and Finite Mixture Models website.

Before moving to the description of the different methods implemented in **mixComp** we would like to briefly mention other theoretical work on the estimation of mixture complexity not currently included in the package. [[3]](#3) propose a method that is reminiscent of the ones described in Section 4. The main difference is that the authors consider distribution functions instead densities, i.e. they consider minimizing a penalized distance between the distribution function of the mixture and the empirical distribution function. The approach of [[7]](#7) is based on a minimum message length-like criterion, however, their method struggles to deal with mixtures with very different weights. [[25]](#25) propose a procedure based on alternating between splitting and merging the components in an EM-algorithm. This algorithm requires selecting two thresholds, the choice of which is somewhat unclear when working with a specific dataset. [[17]](#17) follow a Bayesian approach, taking the usual finite mixture model with Dirichlet weights and putting a prior distribution on the unknown number of components. 

# Objects and functions defined in mixComp

Table 1 depicts five object classes defined in **mixComp**. The first two respectively represent a finite mixture distribution and a random sample drawn from such a distribution. The `Mix` object is printed as a matrix of component weights and parameters and is plotted as the density of the specified mixture, showing the overall as well as the component densities. The `rMix` object prints as the vector of observations and plots as a histogram, showcasing the individual components as well as the full sample. Both objects contain a number of attributes giving additional information, details of which can be found in the corresponding **R** help files. 

#### Table 1: Objects and functions defined in mixComp
|  Object class  | Created via                                  | Description                                     |
|:--------------:|:--------------------------------------------:|:-----------------------------------------------:|
| `Mix`          | `Mix`                                        | Represents a finite mixture                     |
| `rMix`         | `rMix`                                       | Randomly generated data from a finite mixture   |
| `datMix`       | `datMix` or `RtoDat`                         | Observed data from (presumably) a finite mixture|
| `hankDet`      | `nonparamHankel`                             | Vector of estimated Hankel matrix determinants  |
| `paramEst`     | `paramHankel(.scaled)`, `L2(.boot).disc`, `hellinger(.boot).disc`, `hellinger(.boot).cont` or `mix.lrt`  | Complexity estimate , together with estimates of the weights and the component parameters|

The generation of an object of class `Mix` hinges on four central arguments: a string `dist` specifying the name of the family of component densities (or kernels), a boolean `discrete` stating whether the distribution is discrete, a vector `w` giving the weights and a list `theta.list` (the component parameters can also be supplied via the `...` argument) containing the parameters of the component densities. While the creation of `Mix` objects is mostly straightforward, two things should be noted in this regard: First, **mixComp** procedures will search for functions called `rdist` and `ddist` in the accessible namespaces. For most "standard" distributions, these functions are contained in the **stats** package and do not need to be user-written (compare with the Section 6). To make use of these functions, it is essential that the string `dist` is named correctly (e.g. to create a gaussian mixture on the basis of the **stats** package, `dist` has to be specified as `norm` instead of `normal`, `gaussian` etc. for the package to find the functions `dnorm` and `rnorm`). Second, the names of the list elements of `theta.list` (for the names of the `...` arguments) have to match the names of the formal arguments of the functions `ddist` and `rdist` exactly (e.g. for a gaussian mixture, the list elements have to be named `mean` and `sd`, as these are the formal arguments used by `rnorm` and `dnorm` functions of the **stats** package).

The third object class shown in Table 1, called `datMix`, represents the data vector based on which the mixture complexity is supposed to be estimated. These objects are most central to the package, as every procedure estimating the order of a mixture takes a `datMix` object as input. Apart from the vector of observations, it contains other "static" information needed for the estimation procedure (in contrast to "tuning parameters", which can be changed with every function call. An example of such a tuning parameter is the number of bootstrap replicates for a function employing a bootstrap procedure). A brief overview of which "static" attributes need to be supplied for each complexity estimation routine is given in Table 2. 

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

Unlike the above mentioned objects whose creation precedes any type of mixture complexity estimation, objects of the bottom two classes (see Table 1) contain the results of the estimation procedures. Generally, the functions estimating the number of components differ in the types of families of component densities for which they allow and in whether they provide estimates of the weights and the component parameters, the latter determining the object class of the estimation result. These differences are shown in Table 3. The function `nonparamHankel` returns an object of class `hankDet`, which is a vector of determinants (scaled and/or penalized), each entry belonging to a certain complexity estimate. The link between these determinants and the mixture complexity will be discussed in the Section 3. `paramEst` objects arise when using any other function estimating the mixture complexity, all of which additionally return estimates of the component weights and parameters. For both object classes, print and plot methods are available to summarize and visualize the estimation results. 

#### Table 3: Distribution restrictions and output types of different functions contained in mixComp
|     R function                    |Restrictions on the family of the component density |Estimation of the mixture parameters|
|:---------------------------------:|:--------------------------------------------------:|:----------------------------------------------:|
|`nonparamHankel`                   |Compatible with `explicit`, `translation` or `scale`|                                                |
|`nonparamHankel(.scaled)`          |Compatible with `explicit`, `translation` or `scale`|                       x                        |  
|`L2(.boot).disc`                   |Discrete distributions                              |                       x                        | 
|`hellinger(.boot).disc`            |Discrete distributions                              |                       x                        | 
|`hellinger(.boot).cont`            |Continuous distributions                            |                       x                        | 
|`mix.lrt`                          |                                                    |                       x                        | 

# Examples using simulated data

The following example creates two `Mix` objects, a 3-component mixture of normal distributions and a 3-component mixture of Poisson distributions. 

``` r
set.seed(0)
# construct a Nix object:
normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
# plot the mixtures:
plot(normLocMix, main = "3-component normal mixture", cex.main = 0.9)
plot(poisMix, main = "3-component poisson mixture", cex.main = 0.9)
```
<p float="left">
  <img src="https://github.com/yuliadm/mixComp/blob/main/figures/normMix.png" width="500" height="400"/>
  <img src="https://github.com/yuliadm/mixComp/blob/main/figures/poisMix.png" width="500" height="400"/>
</p>

If required, random samples can be generated from these mixtures.
``` r
# generate random samples:
normLocRMix <- rMix(1000, obj = normLocMix)
poisRMix <- rMix(1000, obj = poisMix)
# plot the histograms of the random samples:
plot(normLocRMix, main = "Three component normal mixture", cex.main = 0.9)
plot(poisRMix, main = "Three component poisson mixture", cex.main = 0.9)
```
<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/normRMix.png" width="500" height="400"/>
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/poisRMix.png" width="500" height="400"/>
</p>

We now define the functions computing the moments (for the Gaussian and Poisson mixtures respectively):
``` r
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}

explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
         matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}
```

Define the MLE functions and construct the `datMix` objects:
``` r 
MLE.norm.mean <- function(dat) mean(dat)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)

MLE.pois <- function(dat) mean(dat)

# create datMix objects:
pois.dM <- RtoDat(poisRMix, theta.bound.list = list(lambda = c(0, Inf)), 
                  MLE.function = MLE.pois, Hankel.method = "explicit",
                  Hankel.function = explicit.pois)

normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list, Hankel.method = "translation",
                     Hankel.function = mom.std.norm)
```

In the case of the scaled version of the method, the penalty should be scaled as well. 
``` r
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

``` r
# print the results (for the Poisson mixture)
print(poisdets_sca_pen)
# plot results for both mixtures:
par(mar = c(5, 5, 1, 1))
plot(poisdets_sca_pen, main = "3-component Poisson mixture", cex.main = 0.9)
plot(normdets_sca_pen, main = "3-component Normal mixture", cex.main = 0.9)
```

<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/np_art_1.png" width="500" height="400"/>
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/np_art_2.png" width="500" height="400"/>
</p>

Having created the data ourselves, we know that it comes from a 3-component Poisson mixture and a 3-component Gaussian mixture respectively. The resulting plots indicate that while theoretically sound, the scaled version of the Hankel method can struggle to correctly identify the number of components in practice.

Applying `paramHankel.scaled` to the same Poisson and Normal mixtures results in the correct identification of the mixture complexity in both cases as can be seen in the plot.

``` r
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
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/p_art_1.png" width="500" height="400"/>
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/p_art_2.png" width="500" height="400"/>
</p>


We now apply the minimum Hellinger distance (with AIC-based penalty term) to the artificially created samples from the 3-component Poisson and normal mixtures using the procedures `hellinger.disc` for the former and `hellinger.cont` for the latter. 

``` r
set.seed(0)
h_disc_pois <- hellinger.disc(pois.dM, threshold = "AIC")
h_cont_norm <- hellinger.cont(normLoc.dM, bandwidth = 0.5, sample.n = 5000, 
                      threshold = "AIC")
par(mar = c(5, 5, 1, 1))
plot(h_disc_pois)
plot(h_cont_norm)
```

<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/dist_art_1.png" width="500" height="400"/>
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/dist_art_2.png" width="500" height="400"/>
</p>


# Examples using real-world data

### The Old Faithful dataset

As a simple example of a given dataset to which mixture models have been applied extensively, take the Old Faithful dataset [[19]](#19), [[1]](#1), [[10]](#10). In the context of mixture model estimation, the variable `waiting`, which gives the time in minutes between eruptions of the Old Faithful geyser in the Yellowstone National Park, is often considered to be the variable of interest. To estimate the number of components of the mixture distribution that provides a suitable approximation to the `waiting` data via **mixComp**, the raw data vector of observations has to be converted to a `datMix` object first. For the sake of exposition we specify all arguments of  the `datMix` function. As has often been done in the relevant literature, we assume that the data comes from a normal mixture.

``` r
faithful.obs <- faithful$waiting
norm.dist <- "norm"
norm.discrete <- FALSE
```
Now a named list containing the bounds for the component parameters (the mean and standard deviation) has to be created:

``` r
# define the range for parameter values:
norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
```
Next, the argument `MLE.function` for the mean and for the standard deviation have to be defined:

``` r
# define the MLE functions for the mean and sd: 
MLE.norm.mean <- function(dat) mean(dat)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)
```
The last two arguments, `Hankel.method` and `Hankel.function`, need to be supplied if the mixture complexity is to be estimated based on the Hankel matrix of the moments of the mixing distribution. 

``` r
method <- "translation"
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}
```

Finally, all previously generated objects are combined to a `datMix` object.

``` r
# construct a datMix object that summarizes all the necessary information:
faithful.dM <- datMix(faithful.obs, dist = norm.dist, discrete = norm.discrete,
                      theta.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = method,
                      Hankel.function = mom.std.norm)
```

We will now check how the minimum Hellinger distance method works for these data. Fitting the distance methods to a continuous density requires a choice of bandwidth. While using the adaptive bandwidth is an option, if the user does not want to do so, it is recommended to use the function `kdensity` from the package **kdensity** [[28]](#28) which automatically selects an optimal bandwidth (can be accessed via `kdensity(data)$bw`). 

`hellinger.cont` fits a 2-component mixture to the data, which fits the data well and comprises similar parameter estimates to those found in the literature.

``` r
# estimate the number of components:
library(kdensity)
res <- hellinger.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw,
                      sample.n = 5000, threshold = "AIC")
plot(res)
```
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/hell_cont_norm.png" width="500" height="400">

### The Children dataset

As another a real-world example, we look at the Children dataset whose content was taken from the Annual Report of the pension fund S.P.P. of 1952. The dataset initially appeared in work of [[22]](#22) and was subsequently analysed by many authors. It entails data on 4075 widows who recieved pension from the fund, with their number of children being our variable of interest. For example, there are 3062 widows without children, 587 widows with one child, etc. Many authors have noted that this data is not consistent with being a random sample from a Poisson distribution since the number of zeros found in the data is too large. Thisted approached this by fitting a mixture of two populations, one which is always zero and one which follows a Poisson distribution. **mixComp** includes this data stored as a dataframe. Here, we want to investigate 
how the Hankel matrix methods compare when fitting the data to a mixture of Poissons.

The estimation process starts with defining the MLE function and constructing of the `datMix` object.

``` r
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


First, we define the penalty term and check the nonparametric method. The result suggests that the data comes from a 2-component mixture.

``` r
# define the penalty:
pen <- function(j, n) j * log(n)
# estimate the number of components:
set.seed(0)
(det_sca_pen <- nonparamHankel(children.dM, j.max = 5, scaled = TRUE, 
                              B = 1000, pen.function = pen))
#< 
#< Estimation of the scaled and penalized determinants for a 'pois' mixture model:
#<  Number of components Determinant
#<                     1    21.61041
#<                     2    17.15443
#<                     3    25.60157
#<                     4    33.67663
#<                     5    41.79636
#plot the results:                              
#plot the results:
plot(det_sca_pen, main = "Non-parametric Hankel method for Children dataset",
     cex.main = 0.9)
```

Next, we check the fit of the parametric version. The printed result of `paramHankel.scaled` shows that this method also suggests 2 to be the number of components, with the first component corresponding to a Poisson distribution with the rate of 0.0306. Note that the limit case proposed by [[22]](#22) results in a point mass at 0, and that this fit therefore nicely lines up with the idea of a component accounting for only the zero observations. The plot shows that this method yields a sensible fit overall.

``` r
set.seed(0)
param_sca <- paramHankel.scaled(children.dM, j.max = 5, B = 1000, ql = 0.025, 
                          qu = 0.975)
#< 
#< Parameter estimation for a 1 component 'pois' mixture model:
#< Function value: 3640.3094
#<              w lambda
#< Component 1: 1 0.3995
#< Optimization via user entered MLE-function.
#< ----------------------------------------------------------------------
#< 
#< Parameter estimation for a 2 component 'pois' mixture model:
#< Function value: 3350.9289
#<                    w lambda
#< Component 1: 0.65959 0.0306
#< Component 2: 0.34041 1.1143
#< Converged in 3 iterations.
#< ----------------------------------------------------------------------
#< 
#< The estimated order is 2.
plot(param_sca, breaks = 8, ylim = c(0, 0.8))
```

<p float="left">
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/np_real.png" width="500" height="400"/>
<img src="https://github.com/yuliadm/mixComp/blob/main/figures/p_real.png" width="500" height="400"/>
</p>

# Computational nuance for mixComp functions using the solnp() solver

Several functions in the **mixComp** package (namely, `nonparamHankel`, `paramHankel`, `hellinger.cont`, `hellinger.disc`, `L2.disc`, `mix.lrt`) make use of `solnp()`  function (**Rsolnp** library), which is a solver for general nonlinear programming problems. The above mentioned **mixComp** functions attempt to generate good starting values by clustering the data via `clara` function prior to applying  `solnp()`, which lead to convergence of the algorithm in most of the cases. However when running multiple simulations, `solnp()` might not converge for particular initial values with default control values. This may happen when very few observations are assigned to some of the clusters have, in which case the solver can get "stuck", not even resulting in bad exit status (codes 1 and 2 in returned by `solnp()` convergence value). This issue can be overcome by specifying the control parameters in the functions using the `solnp()` solver (e.g. by defining `L2.boot.disc(geom.dM, j.max = 5, B = 500, ql = 0.025, qu = 0.975, control = list("rho" = 0.1, tol = 0.0000001, trace = 0))`, see description of the `solnp()` function for further details on setting the control parameters) or by setting a time limit for the function execution and going to the next iteration whenever the time limit is exceeded. 


# Further details

For more information on the functions used in **mixComp** and for further examples see the vignette or the documentation.md file.

# References

<a id="1">[1]</a> 
Azzalini A, Bowman A (1990). ???A Look at Some Data on the Old Faithful Geyser.??? Journal
of the Royal Statistical Society C, 39(3), 357???365. doi:10.2307/2347385. URL https:
//www.jstor.org/stable/2347385.

<a id="2">[2]</a> 
Benaglia T, Chauveau D, Hunter DR, Young D (2009). ???mixtools: An R Package for
Analyzing Finite Mixture Models.??? Journal of Statistical Software, 32(6), 1???29. URL
http://www.jstatsoft.org/v32/i06/.

<a id="3">[3]</a> 
Chen J, Kalbfleisch J (1996). ???Penalized Minimum-Distance Estimates in Finite Mixture
Models.??? Canadian Journal of Statistics, 24(2), 167???175.

<a id="4">[4]</a> 
Cutler A, Cordero-Bra??a OI (1996). ???Minimum Hellinger Distance Estimation for Finite
Mixture Models.??? Journal of the American Statistical Association, 91(436), 1716???1723.
ISSN 01621459. doi:10.2307/2291601. URL http://www.jstor.org/stable/2291601.

<a id="5">[5]</a> 
Dacunha-Castelle D, Gassiat E (1997). ???The Estimation of the Order of a Mixture Model.???
Bernoulli, 3(3), 279???299. doi:10.2307/3318593. URL https://www.jstor.org/stable/
3318593.

<a id="6">[6]</a> 
Dempster AP, Laird NM, Rubin DB (1977). ???Maximum Likelihood from Incomplete Data
via the EM Algorithm.??? Journal of the Royal Statistical Society B, 39(1), 1???38. ISSN
0035-9246. With discussion, URL http://www.jstor.org/stable/2984875.

<a id="7">[7]</a> 
Figueiredo MAT, Jain AK (2002). ???Unsupervised Learning of Finite Mixture Models.??? IEEE
Transactions on Pattern Analysis and Machine Intelligence, 24(3), 381???396.

<a id="8">[8]</a> 
Gr??n B, Leisch F (2007). ???Fitting Finite Mixtures of Generalized Linear Regressions in
R.??? Computational Statistics & Data Analysis, 51(11), 5247???5252. doi:10.1016/j.csda.
2006.08.014.

<a id="9">[9]</a> 
Gr??n B, Leisch F (2008). ???FlexMix Version 2: Finite Mixtures with Concomitant Variables
and Varying and Constant Parameters.??? Journal of Statistical Software, 28(4), 1???35. doi:
10.18637/jss.v028.i04. URL http://www.jstatsoft.org/v28/i04/.

<a id="10">[10]</a> 
H??rdle W (1991). Smoothing Techniques : With Implementation in S. Springer-Verlag New
York, New York, NY, USA. doi:10.1007/978-1-4612-4432-5.

<a id="11">[11]</a> 
Leisch F (2004). ???FlexMix: A General Framework for Finite Mixture Models and Latent
Class Regression in R.??? Journal of Statistical Software, 11(8), 1???18. doi:10.18637/jss.
v011.i08. URL http://www.jstatsoft.org/v11/i08/.

<a id="12">[12]</a> 
Lindsay BG (1983a). ???The Geometry of Mixture Likelihoods: A General Theory.??? The Annals
of Statistics, 11(1), 86???94. ISSN 0090-5364. doi:10.1214/aos/1176346059.30

<a id="13">[13]</a> 
Lindsay BG (1983b). ???The Geometry of Mixture Likelihoods: The Exponential Family.??? The
Annals of Statistics, 11(3), 783???792. ISSN 0090-5364. doi:10.1214/aos/1176346245.

<a id="14">[14]</a> 
Macdonald P, Du J (2018). mixdist: Finite Mixture Distribution Models. R package version
0.5-5, URL https://CRAN.R-project.org/package=mixdist.

<a id="15">[15]</a> 
McLachlan G, Peel D (2000). Finite Mixture Models. Wiley Series in Probability and Statis-
tics: Applied Probability and Statistics. John Wiley & Sons. ISBN 0-471-00626-2. doi:
10.1002/0471721182.

<a id="16">[16]</a> 
Melnykov V, Chen WC, Maitra R (2012). ???MixSim: An R Package for Simulating Data
to Study Performance of Clustering Algorithms.??? Journal of Statistical Software, 51(12),
1???25. URL http://www.jstatsoft.org/v51/i12/.

<a id="17">[17]</a> 
Miller JW, Harrison MT (2018). ???Mixture Models with a Prior on the Number of Compo-
nents.??? Journal of the American Statistical Association, 113(521), 340???356.

<a id="18">[18]</a> 
Moss J, Tveten M (2019). kdensity: Kernel Density Estimation with Parametric Starts
and Asymmetric Kernels. R package version 1.0.1, URL https://CRAN.R-project.org/
package=kdensity.

<a id="19">[19]</a> 
R Core Team (2020). R: A Language and Environment for Statistical Computing. R Foun-
dation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

<a id="20">[20]</a> 
Scrucca L, Fop M, Murphy TB, Raftery AE (2016). ???mclust 5: Clustering, Classification and
Density Estimation Using Gaussian Finite Mixture Models.??? The R Journal, 8(1), 289???317.
URL 10.32614/RJ-2016-021.

<a id="21">[21]</a> 
Teicher H (1963). ???Identifiability of Finite Mixtures.??? The Annals of Mathematical Statistics,
34, 1265???1269. ISSN 0003-4851. doi:10.1214/aoms/1177703862.

<a id="22">[22]</a> 
Thisted RA (1988). Elements of Statistical Computing: Numerical Computation. Chapman
& Hall, Ltd., GBR. ISBN 0412013711.

<a id="23">[23]</a> 
Titterington DM, Smith AFM, Makov UE (1985). Statistical Analysis of Finite Mixture
Distributions. Wiley Series in Probability and Mathematical Statistics: Applied Probability
and Statistics. John Wiley & Sons. ISBN 0-471-90763-4.

<a id="24">[24]</a> 
Umashanger T, Sriram T (2009). ???L2E Estimation of Mixture Complexity for Count Data.???
Computational Statistics & Data Analysis, 53(12), 4243 ??? 4254. ISSN 0167-9473. doi:10.
1016/j.csda.2009.05.013. URL http://www.sciencedirect.com/science/article/
pii/S0167947309002023.

<a id="25">[25]</a> 
Wang HX, Luo B, Zhang QB, Wei S (2004). ???Estimation for the Number of Components
in a Mixture Model Using Stepwise Split-and-Merge EM Algorithm.??? Pattern Recognition
Letters, 25(16), 1799???1809. ISSN 0167-8655.

<a id="26">[26]</a> 
Woo MJ, Sriram T (2007). ???Robust Estimation of Mixture Complexity for Count Data.???
Computational Statistics & Data Analysis, 51(9), 4379 ??? 4392. ISSN 0167-9473. doi:10.
1016/j.csda.2006.06.006. URL http://www.sciencedirect.com/science/article/
pii/S0167947306001964.

<a id="27">[27]</a> 
Woo MJ, Sriram TN (2006). ???Robust Estimation of Mixture Complexity.??? Journal of the 
American Statistical Association, 101(476), 1475???1486. doi:10.1198/016214506000000555.

<a id="28">[28]</a> 
Xekalaki E, Karlis D (1999). ???On Testing for the Number of Components in a Mixed Poisson
Model.??? The Annals of the Institute of Statistical Mathematics, 51, 149???162. doi:10.1023/
A:1003839420071.
