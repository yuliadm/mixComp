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
bibliography: refs.bib

---

# Summary

Mixture models, studied extensively in [@Teicher63; @LindsayI; @LindsayII; @Titterington; @McLachlan] etc., allow for modeling heterogeneous data. The (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed as a random variable drawn from some unknown mixing distribution. The number of mixture components could be known in advance, in which case the model is *fully parametric* and the well-known expectation-maximization (EM) algorithm [@Dempster] can be used for finding the maximum likelihood estimates (MLE) of the unknown parameters. However, in many applications the number of components is unknown and has to be estimated from the data. 

**mixComp** provides three categories of methods for estimating the unknown complexity of a (univariate) finite mixture. All methods come with theoretical guarantees for consistency. Their performance varies according to the underlying mixture distribution and the sample size.

# Statement of need

Two main features distinguish **mixComp** from other mixture-related **R** [@R] packages: 

- while mixture component weights and parameters are often estimated as a by-product, all **mixComp** methods are based on theory specifically developed to consistently estimate mixture complexity;

- it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

Aside from **mixtools**'s [@mixtools] focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm, also containing routines for selecting the number of components based on information criteria and parametric bootstrapping of the LRT statistic values. However, they are limited to multinomial and normal mixtures as well as mixtures-of-regressions. While **flexmix** [@flexmix1; @flexmix2; @flexmix3] was developed to deal with mixtures-of-regression, it stands out from other packages by its extensibility, a design principle that we also aimed for when creating **mixComp**. Other packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], used for grouped conditional data. 

**mixComp** can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution. The package is aimed at practitioners studying phenomena that can be effectively modelled using mixture distributions. 

# Methods

A distribution $F$ is called a *finite mixture* if its probability density/mass is of the form
$$f(x) = \sum_{i=1}^p w_i g_i(x),$$
$p \in \mathbb{N}$ being the mixture complexity, $(w_1, \dots w_p : \sum_{i=1}^p w_i = 1$, $w_i \geq 0,$ for $i=1,\dots,p)$ - the component weights and $g_i(x)$ - the $i$-th component density. 

Given some complexity $j$, the relevant parameter spaces are
$$\Theta_j = \{\theta_1, \dots, \theta_j: \theta_i \in \Theta \subseteq \mathbb{R}^d, \quad d \in \mathbb{N}, \quad \text{ for } i = 1,\dots,j\}, \text { and }$$
$$W_j = \{w_1, \dots, w_j: \sum_{i=1}^j w_i = 1, w_i \geq 0, \text{ for } i = 1,\dots,j\}.$$

Assume the family of the component densities $\{g(x; \theta)\}$ is known, while $\theta=(\theta_1, \dots, \theta_p) \in \Theta_p$, $\textbf{w} = (w_1, \dots, w_p) \in W_p$ and $p \in \mathbb{N}$ are unknown. 

### 1. Functions using Hankel matrices of moments of the mixing distribution

The basic approach [@hankel] estimates the mixture order (based on $\textbf{X} = \{X_1, \dots, X_n\}$, an i.i.d. $n$-sample from $F$) as
$$\hat{p} := \text{argmin}_{j \in \mathbb{N}} J_n(j),$$
where 
$$J_n(j) := \lvert \det H(\hat{\textbf{c}}_{2j+1}) \rvert + A(j)l(n),$$
with $l(n)$ being a positive function converging to $0$ as $n\to\infty$, $A(j)$ - positive and strictly increasing function,
$H(\mathbf{c_{2j+1}})$ - the Hankel matrix built on $\hat{\textbf{c}}_{2j+1}$, the consistent estimator of the first $2j+1$ moments of the mixing distribution. 

**mixComp** offers several methods for calculating $\hat{\textbf{c}}^{2j+1}$ and provides extensions of the basic approach.


# References
