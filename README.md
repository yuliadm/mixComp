# mixComp

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

The popularity of such models stems from the fact that they allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. To account for such heterogeneity, the (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method. This includes the well-known nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were well studied in the seminal work of [@LindsayI], [@LindsayII]. One remarkable property of the NPMLE of the mixing distribution is that it is, under some simple conditions, a discrete distribution function with at most $k$ number of jumps, where $k$ is the number of distinct observations in the random sample drawn from the mixture. This interesting feature is one reason, among others, for considering the smaller class of finite mixture models, i.e., mixture models with a discrete mixing distribution with a finite number of jumps. The model has the following simple interpretation: the population under study is assumed to consist of several homogeneous subpopulations. These subpopulations, typically referred to as the mixture's components, usually have a specific meaning depending on the problem at hand. In some very simple situations, the number of components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric maximum likelihood estimator (MLE) is known to occur at the fast rate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973529.jpg"> (under some regularity conditions). Also, the well-known expectation-maximization (EM) algorithm can be used to find the MLE of all the unknown parameters; see for example [@Dempster]. However, in many statistical applications such knowledge is rarely available and the number of components has to be estimated from the data. Although the mixture is still finite and the distribution of each component is assumed to belong to some parametric family, the estimation framework in this case is much harder than in the fully parametric one, where the number of components is known. In this paper, the terms *order*, *complexity*  and *number of components*  will  be used interchangeably to refer to this unknown number. The main goal of the package ***mixComp** is to estimate the unknown complexity using several methods known from the statistical literature. These methods, which are discussed below in more detail, all come with theoretical guarantees for consistency as the sample size gets larger. Of course, consistency in this case means that an estimator is able to exactly recover the unknown complexity for large sample sizes. As expected, the performance of the methods varies according to the underlying mixture distribution and the sample size. This will be illustrated below through several synthetic as well as real datasets.

To describe the estimation problem, we start with some formal notation.  A distribution <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> is called a *finite mixture* if its density (we write density throughout and keep in mind that it may be taken with respect to the Lebesgue or the counting measure) is of the form

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973590.jpg">

where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973785.jpg"> is the mixture complexity, 
<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973708.jpg"> are the component weights and the density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973751.jpg"> is the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975606.jpg">th component of the mixture. As the scope of **mixComp** is limited to mixtures where the family of the component distributions is known, we replace <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973751.jpg"> by a parametric density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976027.jpg"> indexed by the (possibly multivariate, say <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg">-dimensional) parameter <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976173.jpg"> in the parameter space <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922033.jpg">.

Given some complexity <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">, the two relevant parameter spaces can therefore be defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976536.jpg">

and

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976599.jpg">

Throughout this document, it is assumed that the family of the component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976784.jpg"> is known, but the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976825.jpg">, the component weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977021.jpg"> and the mixture complexity <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973785.jpg"> are unknown, with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> being the parameter of interest. Assume now that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> is a finite mixture distribution with density <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649973590.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977574.jpg"> is an i.i.d. sample of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg">. 

The **mixComp** package aims to estimate the smallest such <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> on the basis of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">, either on its own or by simultaneously estimating the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977761.jpg"> and the component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976173.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978071.jpg">.

In this setup, it seems natural to test for the number of components by comparing two consecutive models (with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978303.jpg"> components, say) and regarding the mixture with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> components as the imposition of a null hypothesis on the larger model (intuitively, the null hypothesis 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978460.jpg"> 

comes to mind). Traditionally, the problem of testing between nested models may be approached by applying the generalized likelihood ratio test and referring to the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978511.jpg">  distribution to assess significance, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978628.jpg"> is given by the number of constraints imposed on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978742.jpg"> in order to arrive at <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978727.jpg">. However, in the context of mixture models, there is no unique way of obtaining <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978727.jpg"> from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978742.jpg">. As an example, the two null hypotheses <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978884.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978911.jpg"> both yield the smaller model, demonstrating the difficulties of applying this method. This problem has been studied extensively in the literature and numerous other approaches to complexity estimation have been suggested, laying the theoretical foundation for the algorithms described further.

This document discusses various categories of functions found in the **mixComp** package, ranging from methods based on Hankel matrices, to techniques based upon distances between densities and likelihood ratio tests. The examples provided in the first sections all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be performed using the functions from the **stats** package. The last section illustrates how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and valuating the density thereof.

Two main features distinguish this package from other mixture-related `R` [@R] packages: Firstly, it is focused on the estimation of the complexity rather than the component weights and parameters. While these are often estimated as a by-product, all methods contained in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture of interest. Secondly, it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** [see @mixtools] and **flexmix** [see @flexmix1; @flexmix2; @flexmix3]  should both be mentioned at this point: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm. Notably, it also contains routines for selecting the number of components based on information criteria and parametric bootstrapping of the likelihood ratio test statistic values. However, they are limited to multinomial and (a variety of) normal mixtures as well as mixtures-of-regressions. Second, while **flexmix** was developed to deal with mixtures-of-regression, it sets itself apart from other packages by its extensibility, a design principle that we also aimed for when creating the  **mixComp** package. Other widely used packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], which is used for grouped conditional data. Interested readers can find a comprehensive list of mixture-related packages on the CRAN Task View: Cluster Analysis and Finite Mixture Models website.

We also want to briefly mention further theoretical work on the estimation of mixture complexity not included in the **mixComp** functionality for now. [@chen] propose a method that is reminiscent of the ones described in Section 4. The main difference is that the authors consider distribution functions instead densities, i.e. they consider minimizing a penalized distance between the distribution function of the mixture and the empirical distribution function. The approach of [@figueiredo] is based on a minimum message length-like criterion, however, their method struggles to deal with mixtures with very different weights. [@xian] propose a procedure based on alternating between splitting and merging the components in an EM-algorithm. This algorithm requires selecting two thresholds, the choice of which is somewhat unclear when working with a specific dataset. [@miller] follow a Bayesian approach, taking the usual finite mixture model with Dirichlet weights and putting a prior distribution on the unknown number of components. 

The methods that were included in the package can be roughly devided into three categories: methods based on Hankel matrices, following the theory of [@hankel] and selected because of the fact that computation of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979158.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979228.jpg"> is not required, a method based on the likelihood ratio test statistic (LRTS) following [@lrt] since a likelihood ratio test seems like a natural approach in our setting and methods employing minimum distance calculations based on works from multiple authors -- [@hell; @hellcont; @l2; @adap] and included as a computationally more efficient alternative to the LRTS method for certain distributions and distances (this method is especially fast for discrete distributions when using Hellinger distance). The relevant theory will be portrayed at the beginning of each of the respective chapters. The examples depicted in these first chapters all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be done by functions available from the **stats** package. The last chapter showcases how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and evaluating the density thereof.


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

The generation of an object of class `Mix` hinges on four central arguments: a string `dist` specifying the name of the family of component densities (or kernels) <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014854.jpg">, a boolean`discrete` stating whether the distribution is discrete, a vector `w` giving the weights <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977761.jpg"> and a list `theta.list` (the component parameters can also be supplied via the `...` argument) containing the parameters of the component densities <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015122.jpg">. While the creation of `Mix` objects is mostly straightforward, two things should be noted in this regard: First, the **mixComp** procedures will search for functions called `rdist` and `ddist` in the accessible namespaces. For most "standard" distributions, these functions are contained in the **stats** package and do not need to be user-written (compare with the Section 6). To make use of these functions, it is essential that the string `dist` is named correctly (e.g. to create a gaussian mixture on the basis of the **stats** package, `dist` has to be specified as `norm` instead of `normal`, `gaussian` etc. for the package to find the functions `dnorm` and `rnorm`). Second, the names of the list elements of `theta.list`(for the names of the `...` arguments) have to match the names of the formal arguments of the functions `ddist` and `rdist` exactly (e.g. for a gaussian mixture, the list elements have to be named `mean` and `sd`, as these are the formal arguments used by `rnorm` and `dnorm`).

The following example creates two `Mix` objects, a 3-component mixture of normal distributions and a 3-component mixture of Poisson distributions. 

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#<"
)
options(knitr.duplicate.label = "allow")
```

```{r setup}
library(mixComp)
```

```{r mixobj}
set.seed(0)
# construct a Nix object:
normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
```

The created `Mix` objects can then be plotted.
```{r mixplot, figures-side, fig.show="hold", out.width="50%"}
# plot the mixtures:
par(mar = c(5, 5, 1, 1))
plot(normLocMix, main = "3-component normal mixture", cex.main = 0.9)
plot(poisMix, main = "3-component poisson mixture", cex.main = 0.9)
```

If required, random samples can be generated from these mixtures.
```{r rmix}
# generate random samples:
normLocRMix <- rMix(1000, obj = normLocMix)
poisRMix <- rMix(1000, obj = poisMix)
```

The histograms of the generated random samples can be plotted.
```{r plotrmix, figures-side, fig.show="hold", out.width="50%"}
# plot the histograms of the random samples:
par(mar = c(5, 5, 1, 1))
plot(normLocRMix, main = "Three component normal mixture", cex.main = 0.9)
plot(poisRMix, main = "Three component poisson mixture", cex.main = 0.9)
```

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

In the preceding example, the data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> was taken from an already existing dataset. As seen before, the `rMix` function can be used to generate a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg">-sized sample from a specific mixture. If this synthesized data is to be used in simulations (i.e. passed to one of the functions estimating the mixture complexity) an `rMix` object can be converted to a `datMix` object via the `RtoDat` function. Apart from `dist` and `discrete`, all `datMix` arguments have to be supplied to `RtoDat` likewise. We will introduce the relevant examples in Section 3.

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

In 1997, [@hankel] showed that the order of a mixture <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> is characterized by the smallest integer <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> such that the determinant of the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017103.jpg"> Hankel matrix of the first <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017156.jpg"> moments of the mixing distribution equals zero. Moreover, it can be shown that this determinant is zero for all <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017206.jpg">. Formally, for any vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017273.jpg"> with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017336.jpg"> equal to 1, the *Hankel matrix* of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017386.jpg"> is defined as the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017441.jpg"> matrix given by

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017499.jpg">

Now, let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017555.jpg"> be the vector containing the first <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017605.jpg"> (raw) moments of the mixing distribution. For finite mixture models, this amounts to

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017661.jpg">

Then, for all <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017720.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017786.jpg"> is non-negative and

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017845.jpg">

Making use of this fact, the first approach to estimating the order of a mixture that is implemented in **mixComp** relies on initially finding a consistent estimator of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018010.jpg"> based on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">, say <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg">, to then iteratively calculate the applicable Hankel matrix while increasing the assumed order <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> until a sufficiently small value of 
<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015122.jpg"> is attained. However, since <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650015122.jpg"> should be close to 0 for all <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017206.jpg">, this would lead to choosing <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg"> rather larger than the true value and it seems natural to introduce a penalty term. Therefore [@hankel] define the empirical penalized objective function as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018549.jpg">

with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018991.jpg"> being a positive function converging to 0 as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019043.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019115.jpg"> being positive and strictly increasing. 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019176.jpg">

is then a consistent estimator of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg">.

As an extension to simply adding a penalty term to the determinant, a scaling approach was considered by [@lilian]. Let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019301.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019380.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019564.jpg">. Since the estimated moments <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg"> are asymptotically normal, one can apply the delta method giving

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019736.jpg">

Instead of inspecting the vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019799.jpg">, one could therefore also base the complexity analysis on a vector of scaled determinants, employing a nonparametric bootstrap procedure on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg">.
To this end, let <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019902.jpg"> denote the covariance matrix of the determinants <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021468.jpg"> calculated on the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021576.jpg">th bootstrap sample for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021631.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021686.jpg">. Note that 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650021686.jpg">

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

with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> the true cumulative distribution function. This expression is called <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022588.jpg"> only for notational consisteny; the RHS is simply a closed-form expression for the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">th moment of the mixing distribution of a geometric mixture. Hence one may take

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
 
Example 3.1. in [@hankel, p.284] describes how <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018010.jpg"> can be estimated if the family of component distributions <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023420.jpg"> is given by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023472.jpg">, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is a known probability distribution whose moments can be given explicitly. In this case, a triangular linear system can be solved for the estimated moments of the mixing distribution <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg"> using the empirical moments of the mixture distribution and the known moments of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">. The former can be estimated from the data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> whereas the latter has to be supplied by the user. Thus, `Hankel.function` contains a function of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> returning the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">th (raw) moment of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">.

As an example, consider a mixture of normal distributions with unknown mean and unit variance. Then <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is the standard normal distribution, and its <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">th moment <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024376.jpg"> is defined as

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024444.jpg">

```{r}
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}
```


#### 3. `Hankel.method = "scale"`

Similarly, example 3.2. in [@hankel, p.285] describes how <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018010.jpg"> can be estimated if the family of component distributions <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023420.jpg"> is given by <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024662.jpg">, where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is a known probability distribution whose moments can be given explicitly. Likewise, a triangular linear system can be solved for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650018121.jpg">, using the empirical moments of the mixture distribution and the known moments of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">. `Hankel.function` contains a function of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> returning the <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">th moment of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg">. Note that squares have to be taken everywhere if for some integer <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650024922.jpg"> (compare with [@hankel, p.285]).

As an example, consider a mixture of normal distributions with zero mean and unknown variance. Then <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650023535.jpg"> is again the standard normal distribution, and its <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">th moment is defined as above.

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

The resulting plots indicate that while theoretically sound, the scaled version of the Hankel method can struggle to correctly identify the number of components in practice.

As the preceding example shows, it can be quite difficult to determine the order estimate from the vector of estimated determinants alone. Thus, the package includes another option of estimating <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> based on Hankel matrices, however, using a more "parametric" approach which goes hand in hand with estimating <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979158.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649979228.jpg">. The `paramHankel` procedure initially assumes the mixture to only contain a single component, setting <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650022011.jpg">, and then sequentially tests <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025663.jpg"> versus <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025642.jpg"> for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025687.jpg">, until the algorithm terminates. To do so, it determines the MLE for a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-component mixture <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025877.jpg">, generates `B` parametric bootstrap samples of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from the distribution corresponding to <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030786.jpg"> and calculates `B` determinants of the corresponding <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017103.jpg"> Hankel matrices. The null hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026196.jpg"> is rejected and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> increased by 1 if the determinant value based on the original data vector <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977791.jpg"> lies outside of the interval <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026269.jpg">, a range specified by the `ql` and `qu` empirical quantiles of the bootstrapped determinants. Otherwise, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is returned as the order estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg">. 

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
# convert the data to vetor:
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

The preceding notation was held as broad as possible, since different distance measures $D$ and non-parametric estimators <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg"> can be used. Those relevant to the package are mostly well-known, still, definitions can be found in the Appendix. Three procedures are implemented in **mixComp** based on the foregoing methodology: `L2.disc`, `hellinger.disc` and `hellinger.cont`.

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

(recall that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg"> is the number of component parameters, i.e. <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922033.jpg">). If a customized function is to be used, its arguments have to named `j` and `n` once more, so if the user wants to include the number of component parameters <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649976201.jpg">, it has to be entered explicitly. 

#### 3. `hellinger.cont`
 
Unlike the two preceding functions, this procedure is applicable to *continuous* mixture models and uses a kernel density estimator (KDE) as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg">. Its `bandwidth` can be chosen by the user, or the adaptive KDE found in [@adap, p. 1720, equation (2)] may be used by specifying `bandwidth = "adaptive"`. The calculations are based on the continuous version of the squared Hellinger distance, where the "best" estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650028366.jpg"> for a given <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> corresponds to

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650047870.jpg">

Since the computational burden of optimizing over an integral to find the "best" weights and component parameters is immense, the algorithm approximates the objective function defined in the previous equation by sampling <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650047870.jpg"> `sample.n` observations <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650048345.jpg"> from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg"> and setting

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

`hellinger.cont` fits a 2-component mixture to the data, which fits the data well and comprises similar parameter estimates to those found in the literature.

```{r faithplothel, fig.width = 5, fig.height = 4}
# estimate the number of components:
library(kdensity)
res <- hellinger.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw,
                      sample.n = 5000, threshold = "AIC")
plot(res)
```


At this point, it is worth having a closer look at the thresholds. They each satisfy <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049424.jpg"> as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019043.jpg">, the sole condition the authors require. Now, the consistency proofs for the estimators defined in this Section all rely on the fact that, as <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650019043.jpg">,

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049583.jpg">

and

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049649.jpg">

where <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977354.jpg"> is the true complexity (compare with [@l2, p. 4253, Proof of the Theorem], [@hell, p. 4383, Proof] and [@hellcont, p. 1485, Proof of Theorem 1]. If however <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049424.jpg"> faster than <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049788.jpg"> for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650017206.jpg">, asymptotically, the decision rule outlined above will always lead to <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> being rejected. Therefore, a second condition should be placed on <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030894.jpg">, namely choosing it in accordance with 

<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049788.jpg">

Both the LIC and the AIC as well as the SBC in the continuous case do not satisfy this condition, yet they are still part of the package for two reasons. First, since they were used in the original papers, they are included for the sake of completeness and reproducibility of original results. Second, consistency is an asymptotic property, and while the aforementioned thresholds do not fulfill it, they still perform well (and not rarely better than consistent thresholds) for smaller sample sizes. In the example above, the number of components is correctly identified under the non-consistent AIC threshold. Nonetheless, the user will get a warning when using one of non-consistent predefined thresholds.

The preceding example shows that <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg"> directly depends on the chosen threshold <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030894.jpg">. While some thresholds can be motivated better than others from a theoretical perspective, the choice will ultimately always remain somewhat arbitrary. It would thus be desirable to have versions of the preceding functions which do not suffer from this drawback. `L2.boot.disc`, `hellinger.boot.disc` and `hellinger.boot.cont` all work similarly to their counterparts, with the exception that the difference in distances is not compared to a predefined threshold but a value generated by a bootstrap procedure.  At every iteration (of <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">), the procedure sequentially tests <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025663.jpg"> versus <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025642.jpg"> for <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650025687.jpg">, using a parametric bootstrap to generate `B` samples of size <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649977745.jpg"> from a <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg">-component mixture given the previously calculated "best" parameter values <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030786.jpg">. For each of the bootstrap samples, again the "best" estimates corresponding to densities with <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649978303.jpg"> components are calculated, as well as their difference in distances from <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650030450.jpg">. The null hypothesis <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026196.jpg"> is rejected and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> increased by 1 if the original difference <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049788.jpg"> lies outside of the interval <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650026269.jpg">, specified by the `ql` and `qu` empirical quantiles of the bootstrapped differences. Otherwise, <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649975639.jpg"> is returned as the order estimate <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650014579.jpg">. 

As an example, consider the poisson `rMix` object generated in the Section 2. It was a sample from the 3-component poisson mixture corresponding to <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049280.jpg"> and <img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1650049340.jpg">.
The corresponding `datMix` object can be generated to use the function  `hellinger.boot.disc` for identifying the number of components as follows: 

```{r poishelex, results='hide', message=FALSE, warning=FALSE, fig.width = 5, fig.height = 4}
poisList <- vector(mode = "list", length = 1)
names(poisList) <- "lambda"
poisList$lambda <- c(0, Inf)
MLE.pois <- function(dat) mean(dat)
pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, 
                  MLE.function = MLE.pois)
set.seed(1)
res <- hellinger.boot.disc(pois.dM, B = 50, ql = 0.025, qu = 0.975)

plot(res)
```
