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

Mixture models have been used broadly in statistical applications, and have been studied extensively in [@Teicher63], [@LindsayI], [@LindsayII], [@Titterington], [@McLachlan], and other works. They allow for modeling heterogeneous data whose distribution cannot be captured by a single parametric distribution. The (unknown) distribution is assumed to result from mixing over some latent parameter in the following sense: the latent parameter is viewed itself as a random variable drawn from some unknown mixing distribution. When this mixing distribution is only assumed to belong to the ensemble of all possible distribution functions, the mixture model is called *nonparametric* and estimation of the mixing distribution requires using some nonparametric estimation method, such as the nonparametric maximum likelihood estimator (NPMLE) whose fundamental properties were scrutinized in [@LindsayI], [@LindsayII]. One remarkable property of the NPMLE of the mixing distribution is that, it is, under some conditions, a discrete distribution function with at most $k$ number of jumps, where $k$ is the number of distinct observations in the random sample drawn from the mixture. This allows for considering the smaller class of finite mixture models (characterized by a discrete mixing distribution with a finite number of jumps). In some very simple situations, the number of mixture components could be known in advance, in which case the model is *fully parametric* and convergence of classical estimators such as the parametric maximum likelihood estimator (MLE) occurs at the rate $\sqrt{n}$ (given certain conditions). Also, the well-known expectation-maximization (EM) algorithm [@Dempster] can be used to find the MLE of the unknown parameters. However, in many applications such knowledge is rarely available and the number of components has to be estimated from the data. The goal of **mixComp** is to estimate the unknown complexity using several methods. 

The **mixComp** package provides several methods for estimating the unknown complexity of a (univariate) finite mixture that can be loosely grouped into 3 categories:

  - methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; 
  
  - methods based on penalized minimum distance between the unknown probability density and a consistent estimator thereof;
   
  - likelihood ratio test (LRT) - based techniques. 

These methods all come with theoretical guarantees for consistency. The performance of the methods varies according to the underlying mixture distribution and the sample size. While not the primary goal, most methods simultaneously estimate the component weights and parameters. 

# Statement of need

Two main features distinguish **mixComp** from other mixture-related **R** [@R] packages: 

- while mixture component weights and parameters are often estimated as a by-product, all methods contained in **mixComp** are based on theory specifically developed to consistently estimate the number of components in the mixture; 

- it is applicable to parametric mixtures well beyond those whose component distributions are included in the **stats** package, making it more customizable than most packages for model-based clustering. 

The packages **mixtools** [@mixtools] and **flexmix** [@flexmix1; @flexmix2; @flexmix3] should both be mentioned at this point: aside from **mixtools**'s focus on mixture-of-regressions and non-parametric mixtures which are less relevant to this package, it is widely used to fit (multivariate) normal, multinomial or gamma mixtures with the EM algorithm, also containing routines for selecting the number of components based on information criteria and parametric bootstrapping of the LRT statistic values. However, they are limited to multinomial and normal mixtures as well as mixtures-of-regressions. While **flexmix** was developed to deal with mixtures-of-regression, it stands out from other packages by its extensibility, a design principle that we aimed for when creating **mixComp**. Other packages dealing with mixture models are **mclust** [@mclust], which fits mixtures of Gaussians using the EM algorithm, **MixSim** [@mixsim], which allows for simulation from mixtures and comparing the performance of clustering algorithms, and **mixdist** [@mixdist], used for grouped conditional data. 

**mixComp** can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution. The estimation results can be printed and plotted for further analysis. The package is aimed at practitioners studying phenomena that can be effectively modelled using mixture distributions. In particular, it can be used to identify settings and conditions, under which a certain method provides more accurate estimates than the others.

# General Framework

A distribution $F$ is called a *finite mixture* if its density (we write density throughout and keep in mind that it may be taken with respect to the Lebesgue or the counting measure) is of the form
$$f(x) = \sum_{i=1}^p w_i g_i(x),$$
where $p \in \mathbb{N}$ is the mixture complexity, $(w_1, \dots w_p : \sum_{i=1}^p w_i = 1$, $w_i \geq 0,$ for $i=1,\dots,p)$ are the component weights and the density $g_i(x)$ is the $i$-th component of the mixture. As the scope of **mixComp** is limited to mixtures where the family of the component distributions is known, we replace $g_i(x)$ by a parametric density $g(x; \theta_i)$ indexed by the (possibly multivariate, say $d$-dimensional) parameter $\theta_i$ in the parameter space $\Theta \subseteq \mathbb{R}^d$.

Given some complexity $j$, the two relevant parameter spaces can therefore be defined as
$$\Theta_j = \{\theta_1 \dots \theta_j: \theta_i \in \Theta \subseteq \mathbb{R}^d, \text{ for } i = 1,\dots,j\}$$
and
$$W_j = \{w_1, \dots, w_j: \sum_{i=1}^j w_i = 1, w_i \geq 0, \text{ for } i = 1,\dots,j\}.$$

Assume that the family of the component densities $\{g(x; \theta):\theta \in \Theta\}$ is known, the component parameters $\textbf{\theta}= (\theta_1, \dots, \theta_p) \in \Theta_p$, weights $\textbf{w} = (w_1, \dots, w_p) \in W_p$ and the mixture complexity $p \in \mathbb{N}$ are unknown. Suppose that $F$ is a finite mixture distribution with density $f(x) = \sum_{i=1}^p w_i g(x; \theta_i)$ and $\textbf{X} = \{X_1, \dots, X_n\}$ is an i.i.d. sample of size $n$ from $F$. **mixComp** selects the smallest such $p$, either on its own or by simultaneously estimating the weights $w_i$ and the component parameters $\theta_i$, $i \in 1, \dots, p$.

# 1. Functions using Hankel matrices

The basic approach (due to [@hankel]) estimates the mixture order as
$$\hat{p} := \text{argmin}_{j \in \mathbb{N}} J_n(j),$$
where 
$$J_n(j) := \lvert \det H(\hat{\textbf{c}}^{2j+1}) \rvert + A(j)l(n),$$
with $l(n)$ being a positive function converging to $0$ as $n\to\infty$, $A(j)$ positive and strictly increasing; 
$\hat{\textbf{c}}^{2j+1}$ is a consistent estimator  of $\textbf{c}^{2j+1}$, the first $2j+1$ moments of the mixing distribution,
and for any vector $\mathbf{c} = (c_0, \dots, c_{2k}) \in \mathbb{R}^{2k+1}$ with $c_0 = 1$, the Hankel matrix of $\mathbf{c}$ is defined as the $(k+1)\times(k+1)$ matrix:
$$H(\mathbf{c})_{i,j} = c_{i+j-2}, \quad \quad 1 \leq i,j \leq k+1.$$

**mixComp** offers 3 methods (`explicit`, `translation` and `scale`) for calculating the consistent moment estimates $\hat{\textbf{c}}^{2j+1}$ and also provides several extensions of the basic approach.

# 2. Functions using distances

Consider the parametric family of mixture densities
$$\mathcal{F}_j = \{ f_{j, \mathbf{w},\mathbf{\theta}} : (\mathbf{w}, \mbox{\boldmath$\theta$}) \in W_j \times \Theta_j \}.$$ 

With $\{g(x;\theta): \theta \in \Theta \}$ set in advance, elements of $\mathcal{F}_j$ can be written as
$$f_{j,\mathbf{w},\mathbf{\theta}}(x) = \sum_{i = 1}^j w_i g(x; \theta_i).$$

The support of $f$ will depend on the support of $g$ and $\mathcal{F}_j \subseteq \mathcal{F}_{j+1}$ (by setting $w_{j+1} = 0$) for all $j$. Now take a specific mixture $f_0 = f_{p_0, \mathbf{w}_0,\mathbf{\theta}_0}$, where $(\mathbf{w}_0,\mathbf{\theta}_0) \in W_{p_0} \times \Theta_{p_0}$. Clearly, the mixture complexity is defined as
$$p_0 = \min\{j:f_0 \in \mathcal{F}_j\}.$$

This suggests an estimation procedure based on initially finding the 'best' possible estimate (in a sense to be determined) $(\hat{\mathbf{w}}^j, \hat{\mathbf{\theta}}^j) \in W_j \times \Theta_j$ for a given value of $j$, in order to compare the thereby specified probability density/mass function 
$$\hat{f}_j(x) = f_{j, \hat{\mathbf{w}}^j, \hat{\mathbf{\theta}}^j}(x),$$
with a non-parametric density/probability mass estimate $\tilde{f}_n(x)$. As the classes $\mathcal{F}_j$ and $\mathcal{F}_{j+1}$ are nested, the distance $D$ between $\hat{f}_j$ and $\tilde{f}_n$ will only decrease with $j$. Thus, it makes sense to add some penalty term $t(j,n)$ (increasing in $j$) to $D(\hat{f}_j, \tilde{f}_n)$ and find the first value of $j$ where the penalized distance for $j$ is smaller than that for $j+1$. The algorithm starts at $j = 1$, if $j$ is the first integer satisfying
\begin{equation}\label{eq:distances}
D(\hat{f}_j, \tilde{f}_n) - D(\hat{f}_{j+1}, \tilde{f}_n) \leq t(j,n),
\end{equation}
$j$ is taken as the estimate $\hat{p}$, otherwise $j$ is increased by $1$ and the procedure is repeated. Consistency of these estimators has been shown in [@l2; @hell; @hellcont]. **mixComp** includes 3 procedures based on the foregoing methodology: `L2.disc`, `hellinger.disc` and `hellinger.cont`.

# 3. Functions using LRTS

The function `mix.lrt` iteratively increases the assumed order $j$ and finds the MLE for both, the density of a mixture with $j$ and $j+1$ components, giving $(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j}) \in W_j \times \Theta_j$ and $(\hat{\mathbf{w}}^{j+1}, \hat{\mathbf{\theta}}^{j+1}) \in W_{j+1} \times \Theta_{j+1}$. It then calculates the corresponding LRTS, defined as
$$\text{LRTS}= -2\ln\left(\frac{L_0}{L_1}\right) \quad \text{, with}$$

$$L_0 = L_{\textbf{X}}(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j}) \quad\quad \text{and} \quad\quad L_1 = L_{\textbf{X}}(\hat{\mathbf{w}}^{j+1}, \hat{\mathbf{\theta}}^{j+1})\text{,}$$
$L_{\textbf{X}}$ being the likelihood function given the data ${\textbf{X}}$.

Next, a parametric bootstrap is used to generate `B` $n$-samples from a $j$-component mixture given the previously calculated MLE $(\hat{\mathbf{w}}^{j}, \hat{\mathbf{\theta}}^{j})$. For each bootstrap sample, the MLEs and LRTS corresponding to the mixture densities with $j$ and $j+1$ components are calculated. $H_0: p = j$ is rejected and $j$ increased by $1$ if the LRTS based on the original data vector $\textbf{X}$ is larger than the chosen quantile of its bootstrapped counterparts. Otherwise, $j$ is returned as the order estimate $\hat{p}$. 

# Examples

The following example creates a `Mix` objects for a 3-component Poisson mixture and plots the mixture density with $\mathbf{w}=(0.45,0.45,0.1), \mathbf{\lambda}=(1,5,10)$ and generates a random sample. 

``` r
set.seed(0)
# construct a Mix object:
poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), 
		lambda = c(1, 5, 10))
# plot the mixtures:
plot(poisMix, main = "3-component poisson mixture", cex.main = 0.9)
```

![3-component poisson mixture](figures/poisMix.png) 

``` r
# generate a random sample:
poisRMix <- rMix(1000, obj = poisMix)
# plot the histogram of the random sample:
plot(poisRMix, main = "3-component poisson mixture", cex.main = 0.9)
```

![3-component poisson mixture](figures/poisRMix.png) 

We now apply the Hankel matrix method to these simulated data.

If $Y \sim Pois(\lambda)$, it is known that
$$\lambda^j = \mathbb{E}[Y(Y-1)\dots(Y-j+1)].$$
Thus
$$\hat{c}^{2j+1}_j = \frac{1}{n} \sum_{i=1}^n X_i(X_i-1)\dots(X_i-j+1).$$

``` r
# define the function for computing the moments:
explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
         matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}
```

Convert the previously obtained samples yielding the objects of class `rMix` to objects of class `datMix`.

``` r
MLE.pois <- function(dat) mean(dat)
# create datMix objects:
pois.dM <- RtoDat(poisRMix, theta.bound.list = list(lambda = c(0, Inf)), 
                  MLE.function = MLE.pois, Hankel.method = "explicit",
                  Hankel.function = explicit.pois)

```

Define the penalty $A(j)l(n) = j \cdot log(n)$ (for the scaled version, the original penalty should be multiplied by $\sqrt{n}$), apply the paramHankel function with scaling, print and plot the results: 

``` r
# define the penalty function:
pen <- function(j, n){
  j * log(n)
}
# apply the nonparamHankel function to the datMix objects:
set.seed(1)
pois_sca_pen <- paramHankel.scaled(pois.dM)
#< 
#< Parameter estimation for a 1 component 'pois' mixture model:
#< Function value: 3008.6529
#<              w lambda
#< Component 1: 1  3.661
#< Optimization via user entered MLE-function.
#< ----------------------------------------------------------------------
#< 
#< Parameter estimation for a 2 component 'pois' mixture model:
#< Function value: 2439.6581
#<                    w lambda
#< Component 1: 0.53139 1.1141
#< Component 2: 0.46861 6.5491
#< Converged in 3 iterations.
#< ----------------------------------------------------------------------
#< 
#< Parameter estimation for a 3 component 'pois' mixture model:
#< Function value: 2401.5429
#<                     w  lambda
#< Component 1: 0.455849  0.8742
#< Component 2: 0.464642  5.1340
#< Component 3: 0.079509 11.0305
#< Converged in 3 iterations.
#< ----------------------------------------------------------------------
#< 
#< The estimated order is 3.
# plot the results
plot(pois_sca_pen)
```
![Scaled Hankel determinants for a poisson mixture](figures/p_art_1.png)

The reader is referred to consult the package documentation for more examples. 

# References
