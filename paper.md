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

Mixture models have been used broadly in statistical applications and studied extensively in [@Teicher63; @LindsayI; @LindsayII; @Titterington; @McLachlan] etc. They allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. The (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed as a random variable drawn from some unknown mixing distribution. In simple situations, the number of mixture components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric MLE occurs at the rate $\sqrt{n}$ (given certain conditions). Also, the well-known expectation-maximization (EM) algorithm [@Dempster] can be used to find the MLE of the unknown parameters. However, in many applications such knowledge is rarely available and the number of components has to be estimated from the data. 

The **mixComp** package provides several methods for estimating the unknown complexity of a (univariate) finite mixture that can be arranged into 3 categories:

  - methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; 
  
  - methods based on penalized minimum distance between the unknown probability density and its consistent estimator;
   
  - likelihood ratio test (LRT) - based techniques. 

All methods come with theoretical guarantees for consistency. Their performance varies according to the underlying mixture distribution and the sample size.

# Statement of need

Two main features distinguish **mixComp** from other mixture-related **R** [@R] packages: 

- while mixture component weights and parameters are often estimated as a by-product, all **mixComp** methods are based on theory specifically developed to consistently estimate mixture complexity;

- it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** [@mixtools] and **flexmix** [@flexmix1; @flexmix2; @flexmix3] should be mentioned: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm, also containing routines for selecting the number of components based on information criteria and parametric bootstrapping of the LRT statistic values. However, they are limited to multinomial and normal mixtures as well as mixtures-of-regressions. While **flexmix** was developed to deal with mixtures-of-regression, it stands out from other packages by its extensibility, a design principle that we aimed for when creating **mixComp**. Other packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], used for grouped conditional data. 

**mixComp** can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution. The estimation results can be printed and plotted for further analysis. The package is aimed at practitioners studying phenomena that can be effectively modelled using mixture distributions. 

# General Framework

A distribution $F$ is called a *finite mixture* if its probability density/mass is of the form
$$f(x) = \sum_{i=1}^p w_i g_i(x),$$
$p \in \mathbb{N}$ being the mixture complexity, $(w_1, \dots w_p : \sum_{i=1}^p w_i = 1$, $w_i \geq 0,$ for $i=1,\dots,p)$ - the component weights and $g_i(x)$ - the $i$-th component density of the mixture. If the family of the component distributions is known, one can replace $g_i(x)$ by a parametric density/mass $g(x; \theta_i)$ indexed by $\theta_i \in \Theta \subseteq \mathbb{R}^d, d \in \mathbb{N}$, $d$ #dimensions.

Given some complexity $j$, the relevant parameter spaces are
$$\Theta_j = \{\theta_1 \dots \theta_j: \theta_i \in \Theta \subseteq \mathbb{R}^d, \text{ for } i = 1,\dots,j\}, \text { and }$$
$$W_j = \{w_1, \dots, w_j: \sum_{i=1}^j w_i = 1, w_i \geq 0, \text{ for } i = 1,\dots,j\}.$$

Assume the family of the component densities $\{g(x; \theta):\theta \in \Theta\}$ is known, the component parameters $\textbf{\theta}= (\theta_1, \dots, \theta_p) \in \Theta_p$, weights $\textbf{w} = (w_1, \dots, w_p) \in W_p$ and mixture complexity $p \in \mathbb{N}$ are unknown. **mixComp** selects the smallest $p$ yielding the 'best' fit (in one of the discussed below senses) to $\textbf{X} = \{X_1, \dots, X_n\}$, an i.i.d. $n$-sample from $F$.

### 1. Functions using Hankel matrices

The basic approach [@hankel] estimates the mixture order as
$$\hat{p} := \text{argmin}_{j \in \mathbb{N}} J_n(j),$$
where 
$$J_n(j) := \lvert \det H(\hat{\textbf{c}}^{2j+1}) \rvert + A(j)l(n),$$
with $l(n)$ being a positive function converging to $0$ as $n\to\infty$, $A(j)$ - positive and strictly increasing function; 
$\hat{\textbf{c}}^{2j+1}$ - a consistent estimator of $\textbf{c}^{2j+1}$, the first $2j+1$ moments of the mixing distribution,
and $$H(\mathbf{c})$$ - the Hankel matrix of $\mathbf{c} = (c_0=1, c_1 \dots, c_{2k}) \in \mathbb{R}^{2k+1}$. 

**mixComp** offers several methods for calculating the consistent moment estimates $\hat{\textbf{c}}^{2j+1}$ and provides extensions of the basic approach.

### 2. Functions using distances

Consider the parametric family
$$\mathcal{F}_j = \{ f_{j, \mathbf{w},\mathbf{\theta}} : (\mathbf{w}, \mbox{\boldmath$\theta$}) \in W_j \times \Theta_j \},$$ 

with $\{g(x;\theta): \theta \in \Theta \}$ and
$$f_{j,\mathbf{w},\mathbf{\theta}}(x) = \sum_{i = 1}^j w_i g(x; \theta_i).$$

The estimation procedure is based on finding the 'best' estimate $(\hat{\mathbf{w}}^j, \hat{\mathbf{\theta}}^j) \in W_j \times \Theta_j$ for a given $j$, to compare the thereby specified probability density/mass function 
$$\hat{f}_j(x) = f_{j, \hat{\mathbf{w}}^j, \hat{\mathbf{\theta}}^j}(x),$$
with a non-parametric probability density/mass estimate $\tilde{f}_n(x)$. As $\mathcal{F}_j \subseteq \mathcal{F}_{j+1}$, a distance measure $D(\hat{f}_j, \tilde{f}_n)$ is monotonically non-increasing with $j$, a penalty $t(j,n)$ is required: 
\begin{equation}\label{eq:distances}
\hat{p} = min_j \big\{D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n) \leq t(j,n) \big\}.
\end{equation}
**mixComp** offers 3 distance-based procedures: `L2.disc`, `hellinger.disc` and `hellinger.cont`. Consistency of these estimators has been shown in [@l2; @hell; @hellcont]. 

### 3. Functions using LRTS

By iteratively increasing $j$, find the MLE for the density of a mixture with $j$ and $j+1$ components, giving $(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j}) \in W_j \times \Theta_j$ and $(\hat{\mathbf{w}}^{j+1}, \hat{\mathbf{\theta}}^{j+1}) \in W_{j+1} \times \Theta_{j+1}$, calculate 
$$\text{LRTS}= -2\ln\left(\frac{L_0}{L_1}\right) \quad \text{, with}$$

$$L_0 = L_{\textbf{X}}(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j}) \quad\quad \text{and} \quad\quad L_1 = L_{\textbf{X}}(\hat{\mathbf{w}}^{j+1}, \hat{\mathbf{\theta}}^{j+1})\text{,}$$
$L_{\textbf{X}}$ being the likelihood function given ${\textbf{X}}$.

Then use a parametric bootstrap to generate `B` $n$-samples from a $j$-component mixture given the calculated MLE $(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j})$. For each bootstrap sample, the MLEs and LRTS corresponding to the mixture densities with $j$ and $j+1$ components are computed. $H_0: p = j$ is rejected and $j$ increased by $1$ if the LRTS is larger than the specified quantile of its bootstrapped counterparts. Otherwise, $\hat{p} = j$. 

# Example

Consider a 3-component Poisson mixture with $\mathbf{w}=(0.45,0.45,0.1), \textrm{ and } \mathbf{\lambda}=(1,5,10)$ and apply the `paramHankel` function with scaling. 

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
# plot the results
plot(pois_sca_pen)
```
![Scaled Hankel determinants for a poisson mixture](figures/p_art_1.png)

The reader is referred to **mixComp** documentation for more technical details and examples. 

# References
