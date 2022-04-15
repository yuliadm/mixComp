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
| `paramEst`     | `paramHankel(.scaled)`, `L2(.boot).disc`, `hellinger(.boot).disc`, `hellinger(.boot).cont` or `mix.lrt`  | Complexity estimate $\hat{p}$, together with estimates of the weights $\hat{\mathbf{w}}$ and the component parameters $\hat{\mathbf{\theta}}$|

The generation of an object of class `Mix` hinges on four central arguments: a string `dist` specifying the name of the family of component densities (or kernels) $\{g(x;\theta):\theta \in \Theta \}$, a boolean`discrete` stating whether the distribution is discrete, a vector `w` giving the weights $w_i$ and a list `theta.list` (the component parameters can also be supplied via the `...` argument) containing the parameters of the component densities $\theta_i, i \in 1, \dots, p$. While the creation of `Mix` objects is mostly straightforward, two things should be noted in this regard: First, the **mixComp** procedures will search for functions called `rdist` and `ddist` in the accessible namespaces. For most "standard" distributions, these functions are contained in the **stats** package and do not need to be user-written (compare with the Section 6). To make use of these functions, it is essential that the string `dist` is named correctly (e.g. to create a gaussian mixture on the basis of the **stats** package, `dist` has to be specified as `norm` instead of `normal`, `gaussian` etc. for the package to find the functions `dnorm` and `rnorm`). Second, the names of the list elements of `theta.list`(for the names of the `...` arguments) have to match the names of the formal arguments of the functions `ddist` and `rdist` exactly (e.g. for a gaussian mixture, the list elements have to be named `mean` and `sd`, as these are the formal arguments used by `rnorm` and `dnorm`).

The following example creates two `Mix` objects, a $3$-component mixture of normal distributions and a $3$-component mixture of Poisson distributions. 

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

The third object class shown in Table 1, called `datMix`, represents the data vector $\textbf{X}$ based on which the mixture complexity is supposed to be estimated. These objects are most central to the package, as every procedure estimating the order of a mixture takes a `datMix` object as input. Apart from $\mathbf{X}$, it contains other "static" information needed for the estimation procedure (in contrast to "tuning parameters", which can be changed with every function call. An example of such a tuning parameter is the number of bootstrap replicates for a function employing a bootstrap procedure). A brief overview of which "static" attributes need to be supplied for each complexity estimation routine is given in Table 2. 

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
As a simple example of a given dataset to which mixture models have been applied extensively, take the Old Faithful dataset [@R; @faithful1; @faithful2]. In the context of mixture model estimation, the variable `waiting`, which gives the time in minutes between eruptions of the Old Faithful geyser in the Yellowstone National Park, is often considered to be the variable of interest. To estimate the number of components of the mixture distribution that provides a suitable approximation to the `waiting` data via **mixComp**, the raw data vector of observations has to be converted to a `datMix` object first. For the sake of exposition we specify all arguments of  the `datMix` function, starting with the vector of observations $\textbf{X}$ and the string `dist`, specifying $\{g(x;\theta):\theta \in \Theta \}$ and the boolean `discrete`. As has often been done in the relevant literature, we assume that the data comes from a normal mixture.

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

Next, the argument `MLE.function` contains a single function if $d = 1$ or a list of functions of length $d$ otherwise, specifying the $d$ functions needed to estimate the MLE's of $\theta$ based on $\textbf{X}$ if $p$ were equal to 1 (i.e. the MLE's of the component distribution). If this argument is supplied and the `datMix` object is handed over to a complexity estimation procedure relying on optimizing over a likelihood function, the `MLE.function` attribute will be used for the single component case. In case the objective function is either not a likelihood or corresponds to a mixture with more than $1$ component, numerical optimization will be used based on **Rsolnp**'s function `solnp` [@solnp, @Rsolnp]. The initial values (for the parameters of a $j$-component mixture, say) supplied to the solver are then calculated as follows: the data $\textbf{X}$ is clustered into $j$ groups by the function `clara` (of the **cluster** package by [@cluster] and the data corresponding to each group is given to `MLE.function`. The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates. Specifying `MLE.function` is optional and if it is not, for example because the MLE solution does not exists in closed form, numerical optimization is used to find the relevant MLE's.

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

In the preceding example, the data vector $\textbf{X}$ was taken from an already existing dataset. As seen before, the `rMix` function can be used to generate a $n$-sized sample from a specific mixture. If this synthesized data is to be used in simulations (i.e. passed to one of the functions estimating the mixture complexity) an `rMix` object can be converted to a `datMix` object via the `RtoDat` function. Apart from `dist` and `discrete`, all `datMix` arguments have to be supplied to `RtoDat` likewise. We will introduce the relevant examples in Section 3.

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
