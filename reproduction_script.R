#----------------------------------------------------------------------
# Preliminaries, load packages
#----------------------------------------------------------------------

library("mixComp")
library("kdensity")


#----------------------------------------------------------------------
# SECTION 2, FIGURE 1
# Create two 'Mix' objects representing a Poisson and a Gaussian
# mixture and plot them
#----------------------------------------------------------------------

poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1),
               lambda = c(1, 5, 10))
normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3),
                  mean = c(10, 13, 17), sd = c(1, 1, 1))

plot(poisMix)
plot(normLocMix)


#----------------------------------------------------------------------
# SECTION 2, FIGURE 2
# Create two 'rMix' objects representing randomly drawn data from the 
# distributions specified in 'poisMix' and 'normLocMix'
#----------------------------------------------------------------------

set.seed(0)
normLocRMix <- rMix(1000, obj = normLocMix)
poisRMix <- rMix(1000, obj = poisMix)

plot(poisRMix)
plot(normLocRMix)


#----------------------------------------------------------------------
# SECTION 2
# Create 'datMix' object from the Old Faithful 'waiting' variable
#----------------------------------------------------------------------

faithful.obs <- faithful$waiting
norm.dist <- "norm"
norm.discrete <- FALSE
norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))

MLE.norm.mean <- function(dat) mean(dat)
MLE.norm.sd <- function(dat){
  sqrt((length(dat) - 1) / length(dat)) * sd(dat)
} 
MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
                      "MLE.norm.sd" = MLE.norm.sd)

method <- "translation"
mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}

faithful.dM <- datMix(faithful.obs, dist = norm.dist, 
                      discrete = norm.discrete, theta.bound.list = norm.bound.list,
                      MLE.function = MLE.norm.list, Hankel.method = method,
                      Hankel.function = mom.std.norm)


#----------------------------------------------------------------------
# SECTION 3
# Create the relevant 'Hankel.function' for each a geometric, a 
# Poisson and a Gaussian (sd = 1) mixture
#----------------------------------------------------------------------

explicit.geom <- function(dat, j){
  1 - ecdf(dat)(j - 1)
}

explicit.pois <- function(dat, j){
  mat <- matrix(dat, nrow = length(dat), ncol = j) - 
  matrix(0:(j-1), nrow = length(dat), ncol = j, byrow = TRUE)
  return(mean(apply(mat, 1, prod)))
}

mom.std.norm <- function(j){
  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
}


#----------------------------------------------------------------------
# SECTION 3, FIGURE 3
# Example: use 'nonparamHankel' on the randomly generated Poisson and 
# Gaussian mixture data
#----------------------------------------------------------------------

# create 'datMix' objects
MLE.pois <- function(dat) mean(dat)

pois.dM <- RtoDat(poisRMix, theta.bound.list = list(lambda = c(0, Inf)), 
                  MLE.function = MLE.pois, Hankel.method = "explicit",
                  Hankel.function = explicit.pois)

normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                     MLE.function = MLE.norm.list, Hankel.method = "translation",
                     Hankel.function = mom.std.norm)


# define penalty function and estimate complexity
pen <- function(j, n) {j * log(n)}

set.seed(0)
res1 <- nonparamHankel(pois.dM, j.max = 5, pen.function = pen,
                       scaled = TRUE)

res2 <- nonparamHankel(normLoc.dM, j.max = 5, pen.function = pen,
                       scaled = TRUE)

plot(res1, main = "3-component Poisson mixture")
plot(res2, main = "3-component Gaussian mixture")


#----------------------------------------------------------------------
# SECTION 3, FIGURE 4
# Example: use 'paramHankel.scaled' on the randomly generated Poisson
# and Gaussian mixture data
#----------------------------------------------------------------------

set.seed(0)

res1 <- paramHankel.scaled(pois.dM)
res2 <- paramHankel.scaled(normLoc.dM)

plot(res1)
plot(res2)


#----------------------------------------------------------------------
# SECTION 3, FIGURE 5
# Example: use 'nonparamHankel' and paramHankel.scaled' on the 
# Children dataset
#----------------------------------------------------------------------

# create 'datMix' object
children.obs <- unlist(children)

MLE.pois <- function(dat) mean(dat)

children.dM <- datMix(children.obs, dist = "pois", discrete = TRUE, 
                      Hankel.method = "explicit", Hankel.function = explicit.pois,
                      theta.bound.list = list(lambda = c(0, Inf)), MLE.function = MLE.pois)

# define penalty and estimate complexity with 'nonparamHankel'
pen <- function(j, n) j * log(n)

set.seed(0)
det_sca_pen <- nonparamHankel(children.dM, j.max = 5, scaled = TRUE, 
                              B = 1000, pen.function = pen)
det_sca_pen

plot(det_sca_pen, main = "Non-parametric Hankel method for Children dataset",
     cex.main = 0.9)


# estimate complexity with 'paramHankel.scaled'
set.seed(0)
res <- paramHankel.scaled(children.dM, j.max = 5, B = 1000, ql = 0.025, 
                          qu = 0.975)
res
plot(res, breaks = 8, ylim = c(0, 0.8))

#----------------------------------------------------------------------
# SECTION 4, FIGURE 6
# Example: use 'hellinger.disc' and 'hellinger.cont' on the randomly 
# generated Poisson and Gaussian mixture data
#----------------------------------------------------------------------

set.seed(0)
res1 <- hellinger.disc(pois.dM, threshold = "AIC")

set.seed(0)
res2 <- hellinger.cont(normLoc.dM, bandwidth = kdensity(normLocRMix)$bw,
                       sample.n = 5000, threshold = "AIC")

plot(res1)
plot(res2)


#----------------------------------------------------------------------
# SECTION 4, FIGURE 7
# Show the effect of different bandwith values on the KDE of the Old 
# Faithful 'waiting' variable using 'kdensity'
#----------------------------------------------------------------------

plot(kdensity(faithful$waiting, bw = 1), main = "Bandwidth = 1", xlab = "")
plot(kdensity(faithful$waiting, bw = 4), main = "Bandwidth = 4", xlab = "")
plot(kdensity(faithful$waiting, bw = 8), main = "Bandwidth = 8", xlab = "")


#----------------------------------------------------------------------
# SECTION 4, FIGURE 8
# Example: use 'hellinger.cont' on the 'waiting' variable
#----------------------------------------------------------------------

set.seed(0)
res <- hellinger.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw,
                      sample.n = 5000, threshold = "AIC")
res
plot(res)


#----------------------------------------------------------------------
# SECTION 4, FIGURE 10
# Example: use 'hellinger.boot.disc' on the Shakespeare dataset
#----------------------------------------------------------------------

# create 'datMix' object
shakespeare.obs <- unlist(shakespeare) - 1

MLE.geom <- function(dat) 1 / (mean(dat) + 1)

Shakespeare.dM <- datMix(shakespeare.obs, dist = "geom", discrete = TRUE, 
                         MLE.function = MLE.geom, theta.bound.list = list(prob = c(0, 1)))


# use 'hellinger.boot.disc' for complexity estimation
set.seed(0)
res <- hellinger.boot.disc(Shakespeare.dM, B = 50, ql = 0.025, qu = 0.975)
res
plot(res, breaks = 100, xlim = c(0, 20))


#----------------------------------------------------------------------
# SECTION 5, FIGURE 11 
# Example: use 'mix.lrt' on the Acidity dataset
#----------------------------------------------------------------------

# create 'datMix' object
acidity.obs <- unlist(acidity)

acidity.dM <- datMix(acidity.obs, dist = "norm", discrete = FALSE, 
                     MLE.function = MLE.norm.list, theta.bound.list = norm.bound.list)


# use 'mix.lrt' for complexity estimation
set.seed(0)
res <- mix.lrt(acidity.dM, B = 50, quantile = 0.95)
res
plot(res)


#----------------------------------------------------------------------
# SECTION 6, FIGURE 12
# Create 'Mix' and 'rMix' objects representing a Gaussian mixture with 
# known standard deviation (= 0.5)
#----------------------------------------------------------------------

# create functions evaulating the density and generating random numbers
# for the component distribution
dnorm0.5 <- function(x, mean){
  dnorm(x, mean = mean, sd = 0.5)
}

rnorm0.5 <- function(n, mean){
  rnorm(n, mean = mean, sd = 0.5)
}


# create 'Mix' and 'rMix' object
set.seed(0)
norm0.5Mix <- Mix("norm0.5", discrete = FALSE, w = c(0.3, 0.4, 0.3), 
                  mean = c(10, 11, 13))
norm0.5RMix <- rMix(1000, obj = norm0.5Mix)
plot(norm0.5Mix)
plot(norm0.5RMix)


#----------------------------------------------------------------------
# SECTION 6
# Create corresponding 'datMix' object
#----------------------------------------------------------------------

norm0.5.list <- list("mean" = c(-Inf, Inf))

MLE.norm0.5 <- function(dat) mean(dat)

norm0.5.dM <- RtoDat(norm0.5RMix, theta.bound.list = norm0.5.list,
                     MLE.function = MLE.norm0.5)


#----------------------------------------------------------------------
# SECTION 6, FIGURE 13
# Example: use 'mix.lrt' on the 'norm0.5.dM' object
#----------------------------------------------------------------------

set.seed(0)
res <- mix.lrt(norm0.5.dM, B = 50, quantile = 0.95)
res
plot(res)



#----------------------------------------------------------------------
# APPENDIX B, TABLE 4
# Applying all applicable methods to the four real-world examples:
# Children, Old Faithful, Shakespeare and Acidity
#----------------------------------------------------------------------

# default values used for the estimation:
# penalty function for 'nonparamHankel' when 'scaled = FALSE'
pen <- function(j, n){
  (j * log(n) / sqrt(n))
}
# penalty function for 'nonparamHankel' when 'scaled = TRUE'
pen_sca <- function(j, n){
  (j * log(n))
}


# Children dataset, poisson components so only use methods for discrete
# distributions

# nonparamHankel
set.seed(1)
res <- nonparamHankel(children.dM, j.max = 6, pen.function = pen)
res

# nonparamHankel (scaled)
set.seed(1)
res <- nonparamHankel(children.dM, j.max = 6, scaled = TRUE, pen.function = pen_sca)
res

# paramHankel
set.seed(1)
res <- paramHankel(children.dM)
res

# paramHankel.scaled
set.seed(1)
res <- paramHankel.scaled(children.dM)
res

# L2.disc
set.seed(1)
res <- L2.disc(children.dM)
res

# L2.boot.disc
set.seed(1)
res <- L2.boot.disc(children.dM, B = 50)
res

# hellinger.disc
set.seed(1)
res <- hellinger.disc(children.dM)
res

# hellinger.boot.disc
set.seed(1)
res <- hellinger.boot.disc(children.dM, B = 50)
res

# mix.lrt
set.seed(1)
res <- mix.lrt(children.dM, B = 50)
res      


# Old Faithful dataset (waiting variable), Gaussian components so only use methods for
# continuous distributions; no Hankel methods because neither mean = 0 nor sd = 1
# can be assumed

# hellinger.cont
set.seed(1)
res <- hellinger.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw)
res

# hellinger.boot.cont
set.seed(1)
res <- hellinger.boot.cont(faithful.dM, bandwidth = kdensity(faithful.obs)$bw, B = 50)
res

# mix.lrt
set.seed(1)
res <- mix.lrt(faithful.dM, B = 50)
res


# Shakespeare data, geometric components so only use methods for discrete
# distributions

# create 'datMix' object now adding the information for the Hankel methods
Shakespeare.dM <- datMix(shakespeare.obs, dist = "geom", discrete = TRUE, 
                         Hankel.method = "explicit", Hankel.function = explicit.geom, 
                         MLE.function = MLE.geom, theta.bound.list = list(prob = c(0, 1)))

# nonparamHankel
set.seed(1)
res <- nonparamHankel(Shakespeare.dM, j.max = 6, scaled = FALSE, pen.function = pen)
res

# nonparamHankel (scaled)
set.seed(1)
res <- nonparamHankel(Shakespeare.dM, j.max = 6, scaled = TRUE, pen.function = pen_sca)
res

# paramHankel
set.seed(1)
res <- paramHankel(Shakespeare.dM)
res

# paramHankel.scaled
set.seed(1)
res <- paramHankel.scaled(Shakespeare.dM)
res

# L2.disc
set.seed(1)
res <- L2.disc(Shakespeare.dM)
res

# L2.boot.disc
set.seed(1)
res <- L2.boot.disc(Shakespeare.dM, B = 50)
res

# hellinger.disc
set.seed(1)
res <- hellinger.disc(Shakespeare.dM)
res

# hellinger.boot.disc
set.seed(1)
res <- hellinger.boot.disc(Shakespeare.dM, B = 50)
res

# mix.lrt
set.seed(1)
res <- mix.lrt(Shakespeare.dM, B = 50)
res


# Acidity dataset, Gaussian components so only use methods for continuous distributions;
# no Hankel methods because neither mean = 0 nor sd = 1 can be assumed

# hellinger.cont
set.seed(1)
res <- hellinger.cont(acidity.dM, bandwidth = kdensity(acidity.obs)$bw)
res

# hellinger.boot.cont
set.seed(1)
res <- hellinger.boot.cont(acidity.dM, bandwidth = kdensity(acidity.obs)$bw, B = 50)
res

# mix.lrt
set.seed(1)
res <- mix.lrt(acidity.dM, B = 50)
res


#----------------------------------------------------------------------
# Call to 'sessionInfo'
#----------------------------------------------------------------------

sessionInfo() 