## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#<"
)
options(knitr.duplicate.label = "allow")

## ----setup--------------------------------------------------------------------
library(mixComp)

## ----mixobj-------------------------------------------------------------------
set.seed(0)
# construct a Nix object:
normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))
poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))

## ----mixplot, figures-side, fig.show="hold", out.width="50%"------------------
# plot the mixtures:
par(mar = c(5, 5, 1, 1))
plot(normLocMix, main = "3-component normal mixture", cex.main = 0.9)
plot(poisMix, main = "3-component poisson mixture", cex.main = 0.9)

## ----rmix---------------------------------------------------------------------
# generate random samples:
normLocRMix <- rMix(1000, obj = normLocMix)
poisRMix <- rMix(1000, obj = poisMix)

## ----plotrmix, figures-side, fig.show="hold", out.width="50%"-----------------
# plot the histograms of the random samples:
par(mar = c(5, 5, 1, 1))
plot(normLocRMix, main = "Three component normal mixture", cex.main = 0.9)
plot(poisRMix, main = "Three component poisson mixture", cex.main = 0.9)

## ----faithopts----------------------------------------------------------------
faithful.obs <- faithful$waiting
norm.dist <- "norm"
norm.discrete <- FALSE

## ----normlist-----------------------------------------------------------------
# define the range for parameter values:
norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))

## ----normfun------------------------------------------------------------------
# define the MLE functions for the mean and sd: 
MLE.norm.mean <- function(dat) mean(dat)
MLE.norm.sd <- function(dat){
sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)

## ----normmom------------------------------------------------------------------
method <- "translation"
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}

## ----faithdatmix--------------------------------------------------------------
# construct a datMix object that summarizes all the necessary information:
faithful.dM <- datMix(faithful.obs, dist = norm.dist, discrete = norm.discrete,
                      theta.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = method,
                      Hankel.function = mom.std.norm)

## ----geommom------------------------------------------------------------------
# define the function for computing the moments:
explicit.geom <- function(dat, j){
  1 - ecdf(dat)(j - 1)
}

## ----geompois-----------------------------------------------------------------
# define the function for computing the moments:
explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
         matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}

## -----------------------------------------------------------------------------
# define the function for computing the moments:
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}

## ----rtodat-------------------------------------------------------------------
MLE.pois <- function(dat) mean(dat)

# create datMix objects:
pois.dM <- RtoDat(poisRMix, theta.bound.list = list(lambda = c(0, Inf)), 
                  MLE.function = MLE.pois, Hankel.method = "explicit",
                  Hankel.function = explicit.pois)


normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list, Hankel.method = "translation",
                     Hankel.function = mom.std.norm)

## ----nonph--------------------------------------------------------------------
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

## ----plotnonph, figures-side, fig.show="hold", out.width="50%"----------------
# print the results (for the Poisson mixture)
print(poisdets_sca_pen)
# plot results for both mixtures:
par(mar = c(5, 5, 1, 1))
plot(poisdets_sca_pen, main = "3-component Poisson mixture", cex.main = 0.9)
plot(normdets_sca_pen, main = "3-component Normal mixture", cex.main = 0.9)

## ----plotph, figures-side, fig.show="hold", out.width="50%"-------------------
# apply papamHankel.scaled to datMix objects:
set.seed(1)
pois_sca_pen <- paramHankel.scaled(pois.dM)
norm_sca_pen <- paramHankel.scaled(normLoc.dM)
# plot the results for both mixtures:
par(mar=c(5, 5, 1, 1))
plot(pois_sca_pen,)
plot(norm_sca_pen)

## ----geomex-------------------------------------------------------------------
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

## ----childex------------------------------------------------------------------
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

## ----childplotnph, fig.width = 5, fig.height = 4------------------------------
# define the penalty:
pen <- function(j, n) j * log(n)
# estimate the number of components:
set.seed(0)
(det_sca_pen <- nonparamHankel(children.dM, j.max = 5, scaled = TRUE, 
                              B = 1000, pen.function = pen))
#plot the results:
plot(det_sca_pen, main = "Non-parametric Hankel method for Children dataset",
     cex.main = 0.9)

## ----childplotph, fig.width = 5, fig.height = 4-------------------------------
set.seed(0)
param_sca <- paramHankel.scaled(children.dM, j.max = 5, B = 1000, ql = 0.025, 
                          qu = 0.975)
plot(param_sca, breaks = 8, ylim = c(0, 0.8))

## ----plothel, figures-side, fig.show="hold", out.width="50%"------------------
set.seed(0)
h_disc_pois <- hellinger.disc(pois.dM, threshold = "AIC")
h_cont_norm <- hellinger.cont(normLoc.dM, bandwidth = 0.5, sample.n = 5000, 
                      threshold = "AIC")
par(mar = c(5, 5, 1, 1))
plot(h_disc_pois)
plot(h_cont_norm)

## ----faithplothel, fig.width = 5, fig.height = 4------------------------------
# estimate the number of components:
library(kdensity)
res <- hellinger.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw,
                      sample.n = 5000, threshold = "AIC")
plot(res)

## ----poishelex, results='hide', message=FALSE, warning=FALSE, fig.width = 5, fig.height = 4----
poisList <- vector(mode = "list", length = 1)
names(poisList) <- "lambda"
poisList$lambda <- c(0, Inf)
MLE.pois <- function(dat) mean(dat)
pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, 
                  MLE.function = MLE.pois)
set.seed(1)
res <- hellinger.boot.disc(pois.dM, B = 50, ql = 0.025, qu = 0.975)

plot(res)

## ---- fig.width = 5, fig.height = 4, results='hide', message=FALSE, warning=FALSE----
set.seed(1)
res <- mix.lrt(faithful.dM, B = 50, quantile = 0.95)
print(res)
plot(res)

## ----lrtacid, fig.width = 5, fig.height = 4, results='hide', message=FALSE, warning=FALSE----
acidity.obs <- unlist(acidity)

acidity.dM <- datMix(acidity.obs, dist = "norm", discrete = FALSE, 
                     MLE.function = MLE.norm.list, 
                     theta.bound.list = norm.bound.list)

set.seed(0)
res <- mix.lrt(acidity.dM, B = 50, quantile = 0.95)
plot(res)

## ---- figures-side, fig.show="hold", out.width="50%", results='hide', message=FALSE, warning=FALSE----
dnorm0.5 <- function(x, mean){
  dnorm(x, mean = mean,  sd = 0.5)
}
rnorm0.5 <- function(n, mean){
  rnorm(n, mean = mean,  sd = 0.5)
}
## create objects `Mix` and `rMix`:
set.seed(1)
norm0.5Mix <- Mix("norm0.5", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 11, 13))
norm0.5RMix <- rMix(1000, obj = norm0.5Mix)
## plot the results:
plot(norm0.5Mix)
plot(norm0.5RMix)

## -----------------------------------------------------------------------------
norm0.5.list <- vector(mode = "list", length = 1)
names(norm0.5.list) <- c("mean")
norm0.5.list$mean <- c(-Inf, Inf)

MLE.norm0.5 <- function(dat) mean(dat)

norm0.5.dM <- RtoDat(norm0.5RMix, theta.bound.list = norm0.5.list,
                     MLE.function = MLE.norm0.5)

## ---- results='hide', message=FALSE, warning=FALSE----------------------------
set.seed(1)
res <- mix.lrt(norm0.5.dM, B = 50, quantile = 0.95)

## ---- fig.width = 5, fig.height = 4-------------------------------------------
print(res)
plot(res)

