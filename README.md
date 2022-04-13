# mixComp

The **mixComp** package provides a number of methods for estimating the complexity of a finite mixture.
The considered approaches can be loosely grouped into three categories:
<ul>
  <li> methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; </li>
  <li> methods based on penalized minimum distance between the unknown probability density and a consistent estimator thereof. The distances considered in this survey are the Hellinger and the L_2-distances; </li>
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


Please cite this software as:
???

For a BibTeX entry, use the output from citation(package = "mixComp").

## Installation

To install from 'CRAN', use:
```
install.packages("mixComp")
```

You can install the development version from GitHub with:
```
devtools::install_github("yuliamd/mixComp")
```

## Section 1. Introduction to finite mixture models and mixComp


Mixture models have been used extensively in statistical applications and therefore have attracted a lot of attention from both theoretical and computational perspectives. Although the list of works on mixture models is too long to make an exhaustive inventory, we can refer to the following important papers and books: [@Teicher63], [@LindsayI], [@LindsayII], [@Titterington] and [@McLachlan].

The popularity of such models stems from the fact that they allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. To account for such heterogeneity, the (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method. This includes the well-known nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were well studied in the seminal work of [@LindsayI], [@LindsayII]. One remarkable property of the NPMLE of the mixing distribution is that it is, under some simple conditions, a discrete distribution function with at most $k$ number of jumps, where $k$ is the number of distinct observations in the random sample drawn from the mixture. This interesting feature is one reason, among others, for considering the smaller class of finite mixture models, i.e., mixture models with a discrete mixing distribution with a finite number of jumps. The model has the following simple interpretation: the population under study is assumed to consist of several homogeneous subpopulations. These subpopulations, typically referred to as the mixture's components, usually have a specific meaning depending on the problem at hand. In some very simple situations, the number of components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric maximum likelihood estimator (MLE) is known to occur at the fast rate $\sqrt{n}$ (under some regularity conditions). Also, the well-known expectation-maximization (EM) algorithm can be used to find the MLE of all the unknown parameters; see for example [@Dempster]. However, in many statistical applications such knowledge is rarely available and the number of components has to be estimated from the data. Although the mixture is still finite and the distribution of each component is assumed to belong to some parametric family, the estimation framework in this case is much harder than in the fully parametric one, where the number of components is known. In this paper, the terms *order*, *complexity*  and *number of components*  will  be used interchangeably to refer to this unknown number. The main goal of the package ***mixComp** is to estimate the unknown complexity using several methods known from the statistical literature. These methods, which are discussed below in more detail, all come with theoretical guarantees for consistency as the sample size gets larger. Of course, consistency in this case means that an estimator is able to exactly recover the unknown complexity for large sample sizes. As expected, the performance of the methods varies according to the underlying mixture distribution and the sample size. This will be illustrated below through several synthetic as well as real datasets.
