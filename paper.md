---
title: 'mixComp: An R package for estimating complexity of a mixture mixture'
tags:
  - R
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Adrian M. Price-Whelan^[Co-first author] # note this makes a footnote saying 'Co-first author'
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID^[Co-first author] # note this makes a footnote saying 'Co-first author'
    affiliation: 2
  - name: Author with no affiliation^[Corresponding author]
    affiliation: 2
affiliations:
 - name: ETH Zurich, Seminar for Statistics, Switzerland
   index: 1
 - name: ETH Zurich, Switzerland
   index: 2
date: 13 August 2017
bibliography: refs.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The **mixComp** package provides a number of methods for obtaining a consistent estimate of the complexity of a 
finite mixture (the focus is made on the univariate case). The considered approaches can be loosely grouped into three categories:
<ul>
  <li> methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; </li>
  <li> methods based on penalized minimum distance between the unknown probability density and a consistent estimator thereof. 
    The distances considered in this survey are the Hellinger and the $L_2$-distances; </li>
  <li> likelihood ratio test (LRT) - based techniques. </li>
</ul>
While not the primary goal, most methods simultaneously estimate the component weights and parameters. 
In this document, we give a brief overview of the methodology, and demonstrate the package's functionality in 
both real world examples and synthetically generated data. Moreover, we show how the package can be used on virtually 
any parametric mixture as long as functions generating random variates and evaluating the density are provided 
for the component distribution.

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
It is applicable to parametric mixtures well beyond those whose component distributions are 
included in the **R** package **stats**, making it more customizable than most packages for model-based clustering. 
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

Mixture models have been used extensively in statistical applications and therefore have attracted a lot of attention from both theoretical 
and computational perspectives. Although the list of works on mixture models is too long to make an exhaustive inventory, we can 
refer to the following important papers and books: [@Teicher63], [@LindsayI], [@LindsayII], [@Titterington] and [@McLachlan].
The popularity of such models stems from the fact that they allow for modeling heterogeneous data whose distribution cannot be 
captured by a single parametric distribution. To account for such heterogeneity, the (unknown) distribution is assumed to result 
from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from 
some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution 
functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method. 
This includes the well-known nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were well studied in the seminal 
work of [@LindsayI], [@LindsayII]. One remarkable property of the NPMLE of the mixing distribution is that it is, under some simple conditions, 
a discrete distribution function with at most $k$ number of jumps, where $k$ is the number of distinct observations in the random sample drawn 
from the mixture. This interesting feature is one reason, among others, for considering the smaller class of finite mixture models, i.e., 
mixture models with a discrete mixing distribution with a finite number of jumps. The model has the following simple interpretation: 
the population under study is assumed to consist of several homogeneous subpopulations. These subpopulations, typically referred to as 
the mixture's components, usually have a specific meaning depending on the problem at hand. In some very simple situations, the number 
of components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the 
parametric maximum likelihood estimator (MLE) is known to occur at the fast rate $\sqrt{n}$ (under some regularity conditions). Also, 
the well-known expectation-maximization (EM) algorithm can be used to find the MLE of all the unknown parameters; see for example [@Dempster]. 
However, in many statistical applications such knowledge is rarely available and the number of components has to be estimated from the data. 
Although the mixture is still finite and the distribution of each component is assumed to belong to some parametric family, the estimation 
framework in this case is much harder than in the fully parametric one, where the number of components is known. In this paper, the terms *order*, 
*complexity*  and *number of components*  will  be used interchangeably to refer to this unknown number. The main goal of the package **mixComp** 
is to estimate the unknown complexity using several methods known from the statistical literature. These methods, which are discussed below 
in more detail, all come with theoretical guarantees for consistency as the sample size gets larger. Of course, consistency in this case 
means that an estimator is able to exactly recover the unknown complexity for large sample sizes. As expected, the performance of the 
methods varies according to the underlying mixture distribution and the sample size. This will be illustrated below through several 
synthetic as well as real datasets.

To describe the estimation problem, we start with some formal notation.  A distribution $F$ is called a *finite mixture* if its density 
(we write density throughout and keep in mind that it may be taken with respect to the Lebesgue or the counting measure) is of the form

$$f(x) = \sum_{i=1}^p w_i g_i(x),$$

where $p \in \mathbb{N}$ is the mixture complexity, $(w_1, \dots w_p : \sum_{i=1}^p w_i = 1$, $w_i \geq 0,$ for $i=1,\dots,p)$ are the 
component weights and the density $g_i(x)$ is the $i$-th component of the mixture. As the scope of **mixComp** is limited to mixtures 
where the family of the component distributions is known, we replace $g_i(x)$ by a parametric density $g(x; \theta_i)$ indexed by the 
(possibly multivariate, say $d$-dimensional) parameter $\theta_i$ in the parameter space $\Theta \subseteq \mathbb{R}^d$.
Given some complexity $j$, the two relevant parameter spaces can therefore be defined as
$$\Theta_j = \{\theta_1 \dots \theta_j: \theta_i \in \Theta \subseteq \mathbb{R}^d, \text{ for } i = 1,\dots,j\}$$

and

$$W_j = \{w_1, \dots, w_j: \sum_{i=1}^j w_i = 1, w_i \geq 0, \text{ for } i = 1,\dots,j\}.$$

Throughout this document, it is assumed that the family of the component densities $\{g(x; \theta):\theta \in \Theta\}$ is known, but 
the component parameters $ \textbf{\theta}= (\theta_1, \dots, \theta_p) \in \Theta_p$, the component weights 
$\textbf{w} = (w_1, \dots, w_p) \in W_p$ and the mixture complexity $p \in \mathbb{N}$ are unknown, with $p$ being the parameter of interest. 
Assume now that $F$ is a finite mixture distribution with density $f(x) = \sum_{i=1}^p w_i g(x; \theta_i)$ and 
$\textbf{X} = \{X_1, \dots, X_n\}$ is an i.i.d. sample of size $n$ from $F$. The **mixComp** package aims to 
estimate the smallest such $p$ on the basis of $\textbf{X}$, either on its own or by simultaneously estimating the 
weights $w_i$ and the component parameters $\theta_i$, $i \in 1, \dots, p$.
In this setup, it seems natural to test for the number of components by comparing two consecutive models. 
Traditionally, the problem of choosing between nested models may be approached by applying the generalized likelihood 
ratio test and referring to the $\chi^2_r$ distribution to assess significance, where $r$ is given by the number of constraints 
imposed on the alternative hypothesis $H_1$ to arrive at the null hypothesis $H_0$. However, in the context of mixture models, 
there are several issues hindering application of this classical theory. One of them is that there is no unique way of obtaining 
$H_0$ from $H_1$. As an example, the two null hypotheses $H_0: w_{j+1} = 0$ and $H_0: \theta_{j+1} = \theta_{1}$ both yield the smaller model, 
showcasing the difficulties of applying the classical asymptotic theory of the likelihood ratio. This problem has been studied extensively 
in the literature and numerous alternative approaches to mixture complexity estimation have been suggested, laying the theoretical 
foundation for the subsequently described algorithms.

This document discusses various categories of functions found in the **mixComp** package, ranging from methods based on Hankel matrices, 
to techniques based upon distances between densities and likelihood ratio tests. The examples provided in the first sections all contain 
mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well 
as random variate generation may be performed using the functions from the **stats** package. The last section illustrates how the **mixComp** 
package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from 
the component distribution and valuating the density thereof.
Two main features distinguish this package from other mixture-related **R** [@R] packages: Firstly, it is focused on the estimation of 
the complexity rather than the component weights and parameters. While these are often estimated as a by-product, all methods contained 
in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture of interest. 
Secondly, it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, 
making it more customizable than most packages for model-based clustering. 
The packages **mixtools** [see @mixtools] and **flexmix** [see @flexmix1; @flexmix2; @flexmix3]  should both be mentioned at this point: 
aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, 
it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm. Notably, it also contains 
routines for selecting the number of components based on information criteria and parametric bootstrapping of the likelihood ratio 
test statistic values. However, they are limited to multinomial and (a variety of) normal mixtures as well as mixtures-of-regressions. 
Second, while **flexmix** was developed to deal with mixtures-of-regression, it sets itself apart from other packages by its extensibility, 
a design principle that we also aimed for when creating the  **mixComp** package. Other widely used packages dealing with mixture models 
are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation 
from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], which is used for grouped conditional 
data. Interested readers can find a comprehensive list of mixture-related packages on the CRAN Task View: Cluster Analysis and 
Finite Mixture Models website.

Before moving to the description of the different methods implemented in **mixComp** we would like to briefly mention pther theoretical 
work on the estimation of mixture complexity not currently included in the package. [@chen] propose a method that is reminiscent of the 
ones described in Section 4. The main difference is that the authors consider distribution functions instead densities, i.e. they consider 
minimizing a penalized distance between the distribution function of the mixture and the empirical distribution function. The approach 
of [@figueiredo] is based on a minimum message length-like criterion, however, their method struggles to deal with mixtures with very 
different weights. [@xian] propose a procedure based on alternating between splitting and merging the components in an EM-algorithm. 
This algorithm requires selecting two thresholds, the choice of which is somewhat unclear when working with a specific dataset. [@miller] 
follow a Bayesian approach, taking the usual finite mixture model with Dirichlet weights and putting a prior distribution on the unknown 
number of components. 

The methods that were included in the package can be roughly devided into three categories: methods based on Hankel matrices, following 
the theory as described in [@hankel] and selected because of the fact that computation of $\textbf{w}$ and $\textbf{\theta}$ is not required, 
a method based on the likelihood ratio test statistic (LRTS) following [@lrt] since a likelihood ratio test seems like a natural approach 
in this setting and methods employing minimum distance calculations based on several works and included as a computationally more efficient 
alternative to the LRTS method for certain distributions and distances; see [@hell; @hellcont; @l2; @adap]. For example, when the distance 
is taken to be the Hellinger distance, such an approach is especially fast for discrete distributions. For a more fluid reading, the relevant 
theory will be portrayed at the beginning of each of the respective sections. The examples depicted in these first chapters all contain 
mixtures of "standard" distributions for which evaluation of the density, cumulative distribution function and quantile function as well 
as random variate generation may be done by functions available from the **stats** package. The last chapter showcases how the **mixComp** 
package can be used to estimate the complexity of any mixture as long as the user provides functions generating random variates from the 
component distribution and evaluating the density thereof.

# References
