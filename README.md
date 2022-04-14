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

To describe the estimation problem, we start with some formal notation.  A distribution<img src="https://github.com/yuliadm/mixComp/blob/main/misc/Tex2Img_1649922138.jpg"> is called a *finite mixture* if its density (we write density throughout and keep in mind that it may be taken with respect to the Lebesgue or the counting measure) is of the form

<img src="http://www.sciweavers.org/tex2img.php?eq=f%28x%29%20%3D%20%5Csum_%7Bi%3D1%7D%5Ep%20w_i%20g_i%28x%29%2C&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="f(x) = \sum_{i=1}^p w_i g_i(x)," width="150" height="51" />

where <img src="http://www.sciweavers.org/tex2img.php?eq=p%20%5Cin%20%5Cmathbb%7BN%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="p \in \mathbb{N}" width="207" height="18" /> is the mixture complexity, <img src="http://www.sciweavers.org/tex2img.php?eq=%28w_1%2C%20%5Cdots%20w_p%20%3A%20%5Csum_%7Bi%3D1%7D%5Ep%20w_i%20%3D%201&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="(w_1, \dots w_p : \sum_{i=1}^p w_i = 1" width="183" height="51" />, <img src="http://www.sciweavers.org/tex2img.php?eq=w_i%20%5Cgeq%200%2C%24%20for%20%24i%3D1%2C%5Cdots%2Cp%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="w_i \geq 0,$ for $i=1,\dots,p)" width="190" height="19" /> are the component weights and the density <img src="http://www.sciweavers.org/tex2img.php?eq=g_i%28x%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="g_i(x)" width="44" height="19" /> is the $i$th component of the mixture. As the scope of **mixComp** is limited to mixtures where the family of the component distributions is known, we replace <img src="http://www.sciweavers.org/tex2img.php?eq=g_i%28x%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="g_i(x)" width="44" height="19" /> by a parametric density <img src="http://www.sciweavers.org/tex2img.php?eq=g%28x%3B%20%5Ctheta_i%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="g(x; \theta_i)" width="64" height="19" /> indexed by the (possibly multivariate, say <img src="http://www.sciweavers.org/tex2img.php?eq=d&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="d" width="17" height="15" />-dimensional) parameter <img src="http://www.sciweavers.org/tex2img.php?eq=%5Ctheta_i&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\theta_i" width="18" height="18" /> in the parameter space <img src="http://www.sciweavers.org/tex2img.php?eq=%5CTheta%20%5Csubseteq%20%5Cmathbb%7BR%7D%5Ed&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\Theta \subseteq \mathbb{R}^d" width="215" height="21" />.

Given some complexity $j$, the two relevant parameter spaces can therefore be defined as
$$\Theta_j = \{\theta_1 \dots \theta_j: \theta_i \in \Theta \subseteq \mathbb{R}^d, \text{ for } i = 1,\dots,j\}$$
and
$$W_j = \{w_1, \dots, w_j: \sum_{i=1}^j w_i = 1, w_i \geq 0, \text{ for } i = 1,\dots,j\}.$$

Throughout this document, it is assumed that the family of the component densities $\{g(x; \theta):\theta \in \Theta\}$ is known, but the component parameters $\mathbf{\theta}  = (\theta_1, \dots, \theta_p) \in \Theta_p$, the component weights $\textbf{w} = (w_1, \dots, w_p) \in W_p$ and the mixture complexity $p \in \mathbb{N}$ are unknown, with $p$ being the parameter of interest. Assume now that $F$ is a finite mixture distribution with density $$f(x) = \sum_{i=1}^p w_i g(x; \theta_i)$$ and $\textbf{X} = \{X_1, \dots, X_n\}$ is an i.i.d. sample of size $n$ from $F$. The **mixComp** package aims to estimate the smallest such $p$ on the basis of $\textbf{X}$, either on its own or by simultaneously estimating the weights $w_i$ and the component parameters $\theta_i$, $i \in 1, \dots, p$.

In this setup, it seems natural to test for the number of components by comparing two consecutive models (with $j$ and $j+1$ components, say) and regarding the mixture with $j$ components as the imposition of a null hypothesis on the larger model (intuitively, the null hypothesis 
$$H_0: \left\{ \exists i \in \{1, 2, \dots, j+1\}: w_{i} = 0 \right\}$$ 
comes to mind). Traditionally, the problem of testing between nested models may be approached by applying the generalized likelihood ratio test and referring to the $\chi^2_r$ distribution to assess significance, where $r$ is given by the number of constraints imposed on $H_1$ in order to arrive at $H_0$. However, in the context of mixture models, there is no unique way of obtaining $H_0$ from $H_1$. As an example, the two null hypotheses $H_0: w_{j+1} = 0$ and $H_0: \theta_{j+1} = \theta_{1}$ both yield the smaller model, demonstrating the difficulties of applying this method. This problem has been studied extensively in the literature and numerous other approaches to complexity estimation have been suggested, laying the theoretical foundation for the algorithms described further.

This document discusses various categories of functions found in the **mixComp** package, ranging from methods based on Hankel matrices, to techniques based upon distances between densities and likelihood ratio tests. The examples provided in the first sections all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be performed using the functions from the **stats** package. The last section illustrates how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and valuating the density thereof.

Two main features distinguish this package from other mixture-related `R` [@R] packages: Firstly, it is focused on the estimation of the complexity rather than the component weights and parameters. While these are often estimated as a by-product, all methods contained in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture of interest. Secondly, it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** [see @mixtools] and **flexmix** [see @flexmix1; @flexmix2; @flexmix3]  should both be mentioned at this point: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm. Notably, it also contains routines for selecting the number of components based on information criteria and parametric bootstrapping of the likelihood ratio test statistic values. However, they are limited to multinomial and (a variety of) normal mixtures as well as mixtures-of-regressions. Second, while **flexmix** was developed to deal with mixtures-of-regression, it sets itself apart from other packages by its extensibility, a design principle that we also aimed for when creating the  **mixComp** package. Other widely used packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], which is used for grouped conditional data. Interested readers can find a comprehensive list of mixture-related packages on the CRAN Task View: Cluster Analysis and Finite Mixture Models website.

We also want to briefly mention further theoretical work on the estimation of mixture complexity not included in the **mixComp** functionality for now. [@chen] propose a method that is reminiscent of the ones described in Section 4. The main difference is that the authors consider distribution functions instead densities, i.e. they consider minimizing a penalized distance between the distribution function of the mixture and the empirical distribution function. The approach of [@figueiredo] is based on a minimum message length-like criterion, however, their method struggles to deal with mixtures with very different weights. [@xian] propose a procedure based on alternating between splitting and merging the components in an EM-algorithm. This algorithm requires selecting two thresholds, the choice of which is somewhat unclear when working with a specific dataset. [@miller] follow a Bayesian approach, taking the usual finite mixture model with Dirichlet weights and putting a prior distribution on the unknown number of components. 

The methods that were included in the package can be roughly devided into three categories: methods based on Hankel matrices, following the theory of [@hankel] and selected because of the fact that computation of $\mathbf{w}$ and $\mathbf{\theta}$ is not required, a method based on the likelihood ratio test statistic (LRTS) following [@lrt] since a likelihood ratio test seems like a natural approach in our setting and methods employing minimum distance calculations based on works from multiple authors -- [@hell; @hellcont; @l2; @adap] and included as a computationally more efficient alternative to the LRTS method for certain distributions and distances (this method is especially fast for discrete distributions when using Hellinger distance). The relevant theory will be portrayed at the beginning of each of the respective chapters. The examples depicted in these first chapters all contain mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well as random variate generation may be done by functions available from the **stats** package. The last chapter showcases how the **mixComp** package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the component distribution and evaluating the density thereof.
