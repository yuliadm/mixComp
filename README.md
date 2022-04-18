---
title: 'Estimating the Complexity of a Finite Mixture with the \proglang{R} Package \pkg{mixComp}'
tags:
  - R
  - mixture complexity estimation
  - mixture models
  - Hankel matrix
  - minimum distance
  - likelihood-rato test
authors:
  - name: Anja Weigel^[author] # note this makes a footnote 
    affiliation: 2  
  - name: Fadoua Balabdaoui^[co-author] # note this makes a footnote 
    affiliation: 1 
  - name: Yulia Kulagina^[co-author, maintainer] # note this makes a footnote 
    affiliation: 1
  - name: Lilian Mueller^[contributor]
    affiliation: 2
  - name: Martin Maechler^[contributor]
    (package 'nor1mix' as model, <https://orcid.org/0000-0002-8685-9910>)
    affiliation: 1  
affiliations:
  - name: ETH Zurich, Seminar for Statistics, Switzerland
   index: 1
  - name: ETH Zurich, Switzerland
   index: 2
date: 17 April 2022
bibliography: refs.bib

---

# Summary

The **mixComp** package provides a number of methods for estimating the complexity of a finite mixture.
The considered approaches can be loosely grouped into three categories:
<ul>
  <li> methods built upon the determinants of the Hankel matrix of moments of the mixing distribution; </li>
  <li> methods based on penalized minimum distance between the unknown probability density and a consistent estimator thereof. The distances considered in this survey are the Hellinger and the $L_2$-distances; </li>
  <li> likelihood ratio test (LRT) - based techniques. </li>
</ul>
While not the primary goal, most methods simultaneously estimate the component weights and parameters. In this document, we give a brief overview of the methodology, and demonstrate the package's functionality in both real world examples and synthetically generated data. Moreover, we show how the package can be used on virtually any parametric mixture as long as functions generating random variates and evaluating the density are provided for the component distribution.

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
The estimation results can be printed out and plotted for further analysis and study.
The documentation contains multiple examples based on simulated as well as real data, several real-world 
datasets have been built-in to the package for convenience.

The use of the **mixComp** package might be of interest to researchers and practitioners who 
are studying phenomena that can be effectively modelled using mixture distributions. 
Among other things it can be used to identify settings and conditions, under which a certain 
method provides more accurate estimates than the others.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
