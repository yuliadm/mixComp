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

Mixture models [see @LindsayI; @LindsayII; @McLachlan; @Teicher63; @Titterington] are a popular statistical tool for modeling heterogeneous data. The number of mixture components could be known in advance, in which case they can be easily estimated ( e.g. the maximum likelihood estimates (MLE) of the unknown parameters can be estimated via the well-known expectation-maximization (EM) algorithm [@Dempster]). However, in many applications the number of components is unknown and has to be estimated from the data. 

**mixComp** provides three categories of methods for estimating the unknown complexity of a (univariate) finite mixture:

  - methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; 
  
  - methods based on penalized minimum distance between the unknown probability density and its consistent estimator;
   
  - likelihood ratio test (LRT) - based techniques. 

All methods come with theoretical guarantees for consistency. 

# Statement of need

Two main features distinguish **mixComp** from other mixture-related **R** [@R] packages: 

- while mixture component weights and parameters are often estimated as a by-product, all **mixComp** methods are based on theory specifically developed to consistently estimate mixture complexity;

- **mixComp** is applicable to parametric mixtures beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

Aside from **mixtools**'s [@mixtools] focus on mixture-of-regressions and non-parametric mixtures, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm, also containing routines for selecting the number of components based on information criteria and parametric bootstrapping of the LRT statistic values. However, they are limited to multinomial and normal mixtures as well as mixtures-of-regressions. While **flexmix** [@flexmix1; @flexmix2; @flexmix3] was developed to deal with mixtures-of-regression, it stands out due to its extensibility, a design principle that we also aimed for. Other packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], used for grouped conditional data. 

**mixComp** can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution. The package is aimed at practitioners studying phenomena that can be effectively modelled using mixture distributions. 

# Methods

A distribution $F$ is called a *finite mixture* if its probability density/mass is of the form
$$f(x) = \sum_{i=1}^p w_i g_i(x, \theta_i),$$
$p \in \mathbb{N}$ being the mixture complexity, $(w_1, \dots w_p : \sum_{i=1}^p w_i = 1$, $w_i \geq 0,$ for $i=1,\dots,p)$ - component weights and $g_i(x, \theta_i)$ - $i$-th component density. 

Given some complexity $j$, the relevant parameter spaces are
$$\Theta_j = \{\theta_1, \dots, \theta_j: \theta_i \in \Theta \subseteq \mathbb{R}^d, \quad d \in \mathbb{N}, \quad \text{ for } i = 1,\dots,j\}, \text { and }$$
$$W_j = \{w_1, \dots, w_j: \sum_{i=1}^j w_i = 1, w_i \geq 0, \text{ for } i = 1,\dots,j\}.$$

Assume the family of the component densities $\{g(x; \theta)\}$ is known, while $\textbf{\theta} = (\theta_1, \dots, \theta_p) \in \Theta_p$, $\textbf{w} = (w_1, \dots, w_p) \in W_p$ and $p \in \mathbb{N}$ are unknown. 

### 1. Functions using Hankel matrices

The basic approach [@hankel] estimates (based on $\textbf{X} = \{X_1, \dots, X_n\}$, an i.i.d. $n$-sample from $F$)
$$\hat{p} := \text{argmin}_{j \in \mathbb{N}} \Big\{ \lvert \det H(\hat{\textbf{c}}_{2j+1}) \rvert + A(j)l(n) \Big\},$$
with positive function $l(n) \to 0 \text{ as } n \to \infty$; positive, strictly increasing function $A(j)$; 
$H(\hat{\textbf{c}}_{2j+1})$ - Hankel matrix built on $\hat{\textbf{c}}_{2j+1}$, the consistent estimator of the first $2j+1$ moments of the mixing distribution. 

**mixComp** offers several methods for calculating $\hat{\textbf{c}}^{2j+1}$ and provides extensions of the basic approach.

### 2. Functions using distances

Consider the parametric family $\mathcal{F}_j = \{ f_{j, \mathbf{w},\mathbf{\theta}} : (\mathbf{w}, \mathbf{\theta}) \in W_j \times \Theta_j \},$
$f_{j,\mathbf{w},\mathbf{\theta}}(x) = \sum_{i = 1}^j w_i g(x; \theta_i), \quad \{g(x;\theta): \theta \in \Theta \},$ 
$\mathcal{F}_j \subseteq \mathcal{F}_{j+1}, \forall j = 1,2, \dots$.

Find the 'best' estimate (e.g. MLE) $(\hat{\mathbf{w}}^j, \hat{\mathbf{\theta}}^j) \in W_j \times \Theta_j$ for a given $j$ and thereby specified probability density/mass function $\hat{f}_j(x) = f_{j, \hat{\mathbf{w}}^j, \hat{\mathbf{\theta}}^j}(x),$
and the non-parametric density/mass estimate $\tilde{f}_n(x)$.  
$$\hat{p} = \min_j \big\{D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n) \leq t(j,n) \big\},$$ 
where $D$ denotes the distance between the distributions, $t(j,n)$ - the penalty term.

**mixComp** offers several distance-based procedures described in [@l2; @hell; @hellcont]. 

### 3. Functions using LRTS

Find the MLE for the mixture density with $j$ and $j+1$ components ($j=1,2,\dots$), yielding $(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j}) \in W_j \times \Theta_j$ and $(\hat{\mathbf{w}}^{j+1}, \hat{\mathbf{\theta}}^{j+1}) \in W_{j+1} \times \Theta_{j+1}$,
$$\text{LRTS}= -2\ln\left(\frac{L_{\textbf{X}}(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j})}{L_{\textbf{X}}(\hat{\mathbf{w}}^{j+1}, \hat{\mathbf{\theta}}^{j+1})}\right) \quad \text{, with}$$

$L_{\textbf{X}}$ being the likelihood function given ${\textbf{X}}$.

Use a parametric bootstrap to generate `B` $n$-samples from a $j$-component mixture given $(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j})$. For each bootstrap sample, compute the MLEs and LRTS corresponding to the mixture densities with $j$ and $j+1$ components. Reject $H_0: p = j$, set $j \leftarrow j+1$ if the LRTS is larger than the specified quantile of its bootstrapped counterparts; set $\hat{p} = j$ otherwise. 

# Example

`paramHankel.scaled` function applied to a Poisson mixture: $\mathbf{w}=(0.45,0.45,0.1), \textrm{ and } \mathbf{\lambda}=(1,5,10)$. 

``` r
set.seed(0)
# construct a Mix object:
poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), 
		lambda = c(1, 5, 10))
# plot the density:
plot(poisMix, main = "3-component poisson mixture", cex.main = 0.9)
```

![3-component poisson mixture](figures/poisMix.png) 

``` r
# generate a random sample:
poisRMix <- rMix(1000, obj = poisMix)
# plot the histogram:
plot(poisRMix, main = "3-component poisson mixture", cex.main = 0.9)
```

![3-component poisson mixture](figures/poisRMix.png) 

Use that for $Y \sim Pois(\lambda)$,
$$\lambda^j = \mathbb{E}[Y(Y-1)\dots(Y-j+1)] \quad \textrm{and}$$
$$\hat{c}^{2j+1}_j = \frac{1}{n} \sum_{i=1}^n X_i(X_i-1)\dots(X_i-j+1).$$

``` r
# define the function for computing the moments:
explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
         matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}
# define the MLE function:
MLE.pois <- function(dat) mean(dat)
# convert to datMix object:
pois.dM <- RtoDat(poisRMix, theta.bound.list = list(lambda = c(0, Inf)), 
                  MLE.function = MLE.pois, Hankel.method = "explicit",
                  Hankel.function = explicit.pois)
# define the penalty function:
pen <- function(j, n){
  j * log(n)
}
# apply the nonparamHankel function to the datMix objects:
set.seed(1)
pois_sca_pen <- paramHankel.scaled(pois.dM)
# plot the results (estimated component densities & estimated mixture):
plot(pois_sca_pen)
```
![Scaled Hankel determinants for a poisson mixture](figures/p_art_1.png)

# References
