
## Purpose: estimating the adaptive Kernel density estimator found in Cutler &
##          Cordero-Brana and sampling from it

ADAPkde <- function(dat, ndistparams, n, j, init, dist, formals.dist, dist_call,
                    sample.n, sample.plot, bs_iter = NULL){

  if(j == 1) w <- 1
  else w <- c(init[1:(j - 1)], 1 - sum(init[1:(j - 1)]))

  theta.list.long <- vector(mode = "list", length = ndistparams)
  names(theta.list.long) <- formals.dist
  for(i in 1:ndistparams){
    theta.list.long[[i]] <- matrix(init[(i*j):((1 + i)*j - 1)], nrow = n, ncol = j, byrow = TRUE)
  }
  theta.list.long$x <- dat

  # get matrix of a's given an estimate of the weights w and an estimate of the component
  # parameters theta.list.long (both supplied via 'init') at points dat

  w.array <- matrix(w, nrow = n, ncol = j, byrow = TRUE)
  f.array <- array(suppressWarnings(do.call(dist_call, args = theta.list.long)), dim = c(n, j))
  a.array <- w.array * f.array / rowSums(array(w.array * f.array, dim = c(n, j)))
  # a.array will be infinite if f.array is infinite for any point (can happen when
  # solnp does not converge).
  a.array[which(is.na(a.array))] <- max(a.array, na.rm = TRUE)

  # calculate \hat{\sigma}_i, i in 1,...j i.e. for every component as the empirical
  # standard deviation of those observations whose weighted density is highest in
  # component i

  ind <- apply(f.array, 1, function(x) which(x == max(x)))
  # if two components have the same value (should not really happen) choose the
  # component with the lower index
  if(is.list(ind)) ind <- sapply(ind, min)

  sigma <- mapply(function(i) sd(dat[ind == i]), 1:j)
  # default value in case no observation has its highest density in component i
  sigma[is.na(sigma)] <- 1
  bw <- matrix(rep(2.283 * sigma / n^(0.283), each = n), nrow = n, ncol = j, byrow = FALSE)

  theta.list.long$mean <- matrix(dat, nrow = n, ncol = j, byrow = FALSE)
  theta.list.long$sd <- bw

  # adaptive kernel density estimate
  kde <- function(y, tll = theta.list.long, a = a.array, N = n){
    return(1/N*sum(a*dnorm(x = y, mean = tll$mean, sd = tll$sd)))
  }
  kde <- Vectorize(kde)

  # sampling from the adaptive kde
  mean.ind <- sample(1:n, size = sample.n, replace = TRUE)
  means <- dat[mean.ind]
  weights <- array(a.array[mean.ind, ]/bw[1, ], dim = c(sample.n, j))
  sd.sample <- function(i) sample(1:j, size = 1, prob = weights[i,])
  sd.sample <- Vectorize(sd.sample)
  sd.ind <- sd.sample(1:sample.n)
  sds <- bw[1, sd.ind]
  sample <- rnorm(sample.n, means, sds)

  if(sample.plot == TRUE){

    if(is.null(bs_iter)){ # not in bootstrap loop yet

      txt <- "Sample from KDE based on data"
      hist(sample, freq = FALSE, col = "light grey", main = txt, breaks = 100,
           xlab = "Sample")
      vals <- seq(min(sample), max(sample), length.out = 100)
      kdevals <- kde(vals)
      lines(vals, kdevals)

    } else if(bs_iter != 0){
      # bs_iter equal 0 just recomputes statistic based on original data

      txt <- paste("Sample from Bootstrap KDE: Iteration ", bs_iter, sep = "")
      hist(sample, freq = FALSE, col = "light grey", breaks = 100,
           main = txt, xlab = "Sample")
      vals <- seq(min(sample), max(sample), length.out = 100)
      lines(vals, kde(vals))

    }
  }

  return(list(kde = kde, sample = sample))

}



## Purpose: returns the approximated Hellinger distance function corresponding to
##          parameters x

.get.fmin.hellinger.c <- function(kde, dat, formals.dist, ndistparams, dist,
                                  sample, kdevals, dist_call){

   fmin <- function(x){

    j <- (length(x) + 1)/(ndistparams + 1)
    w <- c(x[1:(j - 1)], 1 - sum(x[1:(j - 1)]))
    if(any(w < 0)) return(0)

    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- matrix(x[(i*j):((1 + i)*j - 1)], nrow = length(sample), ncol = j,
                                     byrow = TRUE)
    }
    theta.list.long$x <- sample

    # NAs or warnings can happen as solnp sometimes uses intermediate
    # parameter values outside of the box constraints (to speed up convergence
    # and avoid numerical ill conditioning)
    mat <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    fvals <- as.matrix(mat[, ]) %*% w
    # Additional constraint of not allowing all weight being put on a single observation;
    # this may happen because we approximate the integral with a sum and this sum can be
    # maximized by increasing the value of a single summand sufficiently while letting
    # the others go to 0
    if(any(is.na(fvals)) || any(fvals > 1000)) return(0)

    return (- sum(sqrt(fvals / kdevals)))
  }

}



## Purpose: returns the approximated Hellinger distance function corresponding to
##          parameters x when the mixture consists of only a single component

.get.fmin.hellinger.c.0 <- function(kde, dat, formals.dist, ndistparams, dist, sample,
                                    kdevals, dist_call){

  fmin <- function(x){

    theta <- setNames(x, formals.dist)
    theta.list <- split(unname(theta), names(theta))
    theta.list$x <- sample

    # NAs or warnings can happen as solnp sometimes uses intermediate
    # parameter values outside of the box constraints (to speed up convergence
    # and avoid numerical ill conditioning)
    fvals <- suppressWarnings(do.call(dist_call, args = theta.list))
    # Additional constraint of not allowing all weight being put on a single observation;
    # this may happen because we approximate the integral with a sum and this sum can be
    # maximized by increasing the value of a single summand sufficiently while letting
    # the others go to 0
    if(any(is.na(fvals)) || any(fvals > 1000)) return(0)

    return (- sum(sqrt(fvals / kdevals)))
  }

}



## Purpose: returns the true Hellinger distance between the kde and the mixture
#           corresponding to parameters x

.get.hellingerD <- function(x, j, ndistparams, formals.dist, kde, dist, discrete, values){

  if(length(x) == ndistparams) w <- 1
  else w <- c(x[1:(j - 1)], 1 - sum(x[1:(j - 1)]))
  theta <- x[j:length(x)]

  theta <- setNames(as.vector(theta), rep(formals.dist, each = j))
  theta.list <- split(unname(theta), names(theta))
  Mix.obj <- Mix(dist, discrete = discrete, w = w, theta.list = theta.list)

  objective <- try(integrate(function(x){sqrt(suppressWarnings(dMix(x, Mix.obj)) * kde(x))}, -Inf, Inf,
                             subdivisions = 1000L)[[1]], silent = TRUE)
  if(inherits(objective, "try-error")){
    warning(" Error while calculating the value of the true objective function. \n Returning the value of the approximated objective function instead. \n")
    objective <- values[length(values)]
  }

  return(2 - 2 * objective)
}



## Purpose: Hellinger distance based method of estimating the mixture complexity of a
##          continuous mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object

#' @title Estimation of a Continuous Mixture Complexity Based on Hellinger Distance
#'
#' @description Estimation of a continuous mixture complexity as well as its component weights and parameters by minimizing the squared Hellinger distance to a kernel density estimate.
#'
#' @details Define the \eqn{complexity} of a finite continuous mixture \eqn{F} as the smallest integer \eqn{p}, such that its probability density function (pdf) \eqn{f} can be written as
#' \deqn{f(x) = w_1*g(x;\theta_1) + \dots + w_p*g(x;\theta_p).}
#' Further, let \eqn{g, f} be two probability density functions. The squared Hellinger distance between \eqn{g} and \eqn{f} is given by
#' \deqn{H^2(g,f) = \int (\sqrt{g(x)}-\sqrt{f(x)})^2 = 2 - 2\int\sqrt{f(x)}\sqrt{g(x)},}
#' where \eqn{\sqrt{g(x)}}, respectively \eqn{\sqrt{f(x)}} denotes the square root of the probability density functions at point x.
#' To estimate \eqn{p}, \code{hellinger.cont} iteratively increases the assumed complexity \eqn{j} and finds the "best" estimate for both, the pdf of a mixture with \eqn{j} and \eqn{j+1} components, ideally by calculating the parameters that minimize the sum of squared Hellinger distances to a kernel density estimate evaluated at each point.  Since the computational burden of optimizing over an integral to find the "best" component weights and parameters is immense, the algorithm approximates the objective function by sampling \code{sample.n} observations \eqn{Y_i} from the kernel density estimate and using
#' \deqn{2 - 2\sum \sqrt{f(Y_i)} / \sqrt{g(Y_i)},}
#' instead, with \eqn{f} being the mixture density and \eqn{g} being the kernel density estimate. Once the "best" parameters have been obtained, the difference in squared distances is compared to a predefined \code{threshold}. If this difference is smaller than the threshold, the algorithm terminates and the true complexity is estimated as \eqn{j}, otherwise \eqn{j} is increased by 1 and the procedure is started over. The predefined thresholds are the \code{"AIC"} given by 
#' \deqn{(d+1)/n}
#' and the \code{"SBC"} given by 
#' \deqn{(d+1)log(n)/(2n),} 
#' \eqn{n} being the sample size and \eqn{d} the number of component parameters, i.e. \eqn{\theta} is in \eqn{R^d}. Note that, if a customized function is to be used, it may only take the arguments \code{j} and \code{n}, so if the user wants to include the number of component parameters \eqn{d}, it has to be entered explicitly.
#' \code{hellinger.boot.cont} works similarly to \code{hellinger.cont} with the exception that the difference in squared distances is not compared to a predefined threshold but a value generated by a bootstrap procedure. At every iteration (\eqn{j}), the function sequentially tests \eqn{p = j} versus \eqn{p = j+1} for \eqn{j = 1,2, \dots}, using a parametric bootstrap to generate \code{B} samples of size \eqn{n} from a \eqn{j}-component mixture given the previously calculated "best" parameter values. For each of the bootstrap samples, again the "best" estimates corresponding to pdfs with \eqn{j} and \eqn{j+1} components are calculated, as well as their difference in approximated squared Hellinger distances from the kernel density estimate. The null hypothesis \eqn{H_0: p = j} is rejected and \eqn{j} increased by 1 if the original difference in squared distances lies outside of the interval \eqn{[ql, qu]}, specified by the \code{ql} and \code{qu} empirical quantiles of the bootstrapped differences. Otherwise, \eqn{j} is returned as the complexity estimate.
#' To calculate the minimum of the Hellinger distance (and the corresponding parameter values), the solver \code{\link[Rsolnp]{solnp}} is used. The initial values supplied to the solver are calculated as follows: the data is clustered into \eqn{j} groups by the function \code{\link[cluster]{clara}} and the data corresponding to each group is given to \code{MLE.function} (if supplied to the \code{\link{datMix}} object \code{obj}, otherwise numerical optimization is used here as well). The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates.
#' @aliases hellinger.cont hellinger.boot.cont
#' @usage
#' hellinger.cont(obj, bandwidth, j.max = 10, threshold = "SBC", sample.n = 5000,
#'                sample.plot = FALSE, control = c(trace = 0))
#'
#' hellinger.boot.cont(obj, bandwidth, j.max = 10, B = 100, ql = 0.025,
#'                     qu = 0.975, sample.n = 3000, sample.plot = FALSE,
#'                     control = c(trace = 0), ...)
#' @param obj object of class \code{\link{datMix}}.
#' @param bandwidth numeric, indicating the bandwidth to be used. Can also be set to "adaptive" if the adaptive kernel density estimator as defined by Cutler & Cordero-Brana (1996, page 1720, Equation 2) should be employed.
#' @param j.max integer stating the maximal number of components to be considered.
#' @param threshold function or character string in \code{c("AIC", "SBC")} specifying which threshold should be used to compare two mixture estimates of complexities \eqn{j} and \eqn{j+1}. If the difference in minimized squared distances is smaller than the relevant threshold, \eqn{j} will be returned as complexity estimate.
#' @param sample.n integer, specifying the sample size to be used for approximation of the objective function (see details).
#' @param sample.plot logical, indicating whether the histogram of the sample drawn to approximate the objective function should be plotted.
#' @param B integer, specifying the number of bootstrap replicates.
#' @param ql numeric between \eqn{0} and \eqn{1}, specifying the lower quantile to which the observed difference in minimized squared distances will be compared.
#' @param qu numeric between \eqn{0} and \eqn{1}, specifying the upper quantile to which the observed difference in minimized squared distances will be compared.
#' @param control control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.
#' @param \dots further arguments passed to the \code{\link[boot]{boot}} function.
#' @return Object of class \code{paramEst} with the following attributes:
#'      \item{dat}{data based on which the complexity is estimated.}
#'      \item{dist}{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density/ mass function and \code{rdist} generates random variates.}
#'      \item{ndistparams}{integer specifying the number of parameters identifying the component distribution, i.e. if \eqn{\theta} is in \eqn{R^d} then \code{ndistparams}\eqn{ = d}.}
#'      \item{formals.dist}{string vector specifying the names of the formal arguments identifying the distribution \code{dist} and used in \code{ddist} and \code{rdist}, e.g. for a gaussian mixture (\code{dist = norm}) amounts to 
#'      \code{mean} and \code{sd}, as these are the formal arguments used by \code{dnorm} and \code{rnorm}.}
#'      \item{discrete}{logical indicating whether the underlying mixture distribution is discrete. Will always be \code{FALSE} in this case.}
#'      \item{mle.fct}{attribute \code{MLE.function} of \code{obj}.}
#'      \item{pars}{say the complexity estimate is equal to some \eqn{j}. Then \code{pars} is a numeric vector of size \eqn{(d+1)*j-1} specifying the component weight and parameter estimates, given as 
#' \deqn{(w_1, ... w_{j-1}, \theta 1_1, ... \theta 1_j, \theta 2_1, ... \theta d_j).}}
#'      \item{values}{numeric vector of function values gone through during optimization at iteration \eqn{j}, the last entry being the value at the optimum.}
#'      \item{convergence}{integer, indicating whether the solver has converged (0) or not (1 or 2) at iteration \eqn{j}.}
#' @references Details can be found in 
#' \enumerate{
#' \item M.-J. Woo and T. Sriram, "Robust Estimation of Mixture Complexity", Journal of the American Statistical Association, Vol. 101, No. 476, 1475-1486, Dec. 2006.
#' \item A. Cutler, O.I. Cordero-Brana, "Minimum Hellinger Distance Estimation for Finite Mixture Models." Journal of the American Statistical Association, Vol. 91, No. 436, 1716-1723,  Dec. 1996.}
#' @seealso \code{\link{hellinger.disc}} for the same estimation method for discrete mixtures, \code{\link[Rsolnp]{solnp}} for the solver, \code{\link{datMix}} for the creation of the \code{datMix} object.
#' @keywords cluster
#' @examples
#'
#' ### generating 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#'
#' ### generating 'rMix' from 'Mix' object (with 1000 observations)
#' set.seed(1)
#' normLocRMix <- rMix(10000, normLocMix)
#'
#' ### generating 'datMix' from 'R' object
#'
#' ## generate list of parameter bounds
#'
#' norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
#'
#' ## generate MLE functions
#'
#' # for "mean"
#' MLE.norm.mean <- function(dat) mean(dat)
#' # for "sd" (the sd function uses (n-1) as denominator)
#' MLE.norm.sd <- function(dat){
#' sqrt((length(dat) - 1) / length(dat)) * sd(dat)
#' }
#' # combining the functions to a list
#' MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
#'                       "MLE.norm.sd" = MLE.norm.sd)
#'
#' ## generating 'datMix' object
#' normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
#'                      MLE.function = MLE.norm.list)
#'
#'
#' ### complexity and parameter estimation
#' \dontrun{
#' set.seed(0)
#' res <- hellinger.cont(normLoc.dM, bandwidth = 0.5, sample.n = 5000)
#' plot(res)
#' }
#'
#' @export hellinger.cont
hellinger.cont <- function(obj, bandwidth, j.max = 10, threshold = "SBC",
                           sample.n = 5000, sample.plot = FALSE, control = c(trace = 0)){

  # check relevant inputs
  .input.checks.functions(obj, bandwidth = bandwidth, j.max = j.max, thrshHel = threshold,
                          sample.n = sample.n, sample.plot = sample.plot,
                          assert_continuous = TRUE, Hankel = FALSE, param = TRUE)
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())

  if(is.character(threshold)){
    # otherwise it is a function and will be calculated further down
    if(threshold == "AIC") thresh <- (ndistparams + 1)/N
    if(threshold == "SBC") thresh <- ((ndistparams + 1) * log(N))/(2 * N)
  }

  j0 <- 0

  repeat{

    j0 <- j0 + 1  # current complexity estimate
    j1 <- j0 + 1

    if(is.function(threshold)){
      thresh <- threshold(n = N, j = j0)
    }

    if(j0 != 1){ # pass on constraints for bootstrap

      ineq.j0 <- ineq.j1
      lx.j0 <- lx.j1
      ux.j0 <- ux.j1

    } else { # need to estimate j0 so need initial values and restrictions

      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)
      restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                           upper = upper)
      lx.j0 <- restrictions.j0$lx
      ux.j0 <- restrictions.j0$ux
    }

    if(bandwidth == "adaptive"){ # use the adaptive Kernel density estimate (kde)

      # will not really be used in the computation of the adaptive kde anyways (for j0 == 1);
      # just supplying something so the same function (ADAPkde) can be used
      if(j0 == 1) theta.j1 <- initial.j0

      # compute the kde and draw a sample from it (used to calculate the approximate
      # hellinger distance); has to be recomputed for every j0
      kde.list <- ADAPkde(dat, ndistparams, N, j0, theta.j1, dist, formals.dist, dist_call,
                          sample.n, sample.plot)
      kde <- kde.list$kde
      sample <- kde.list$sample
      kdevals <- kde(sample)

      # calculate optimal parameters for j0; cannot reuse them with the adaptive kde
      # since the kde and kdevals will change in every iteration of j0
      if(j0 > 1){ # need to include weight restrictions in optimization

        initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                       formals.dist)
        fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)
        ineq.j0 <- restrictions.j0$ineq

        opt <- solnp(initial.j0, fun = fmin, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                     LB = lx.j0, UB = ux.j0, control = control)
        # if we estimate multiple components check that all weights satisfy the constraints
        theta.j0 <- opt$pars <- .augment.pars(opt$pars, j0)

      } else { # already know w = 1 (single component mixture)

        fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                        kdevals, dist_call)

        opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)
        theta.j0 <- opt$pars
      }

      Hellinger.j0 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j0, j0,  ndistparams, formals.dist,
                                                                        kde, dist, discrete, opt$values[length(opt$values)])
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
      .printresults(opt, j0, dist, formals.dist, ndistparams)

    } else if(j0 != 1) { # if bandwidth != "adaptive"
      # for the standard kde we can reuse the values calculated for j1 since the kde
      # does not change with j0

      theta.j0 <- theta.j1
      Hellinger.j0 <- Hellinger.j1
      conv.j0 <- conv.j1
      values.j0 <- values.j1

    } else { # if bandwidth != "adaptive" and j0 == 1

      # compute the kde and draw a sample from it (used to calculate the approximate
      # hellinger distance) in the first iteration (i.e. j0 == 1)
      kde <- kdensity(dat, bw = bandwidth, kernel = "gaussian")
      rkernel <- function(n) rnorm(n, sd = bandwidth)
      sample <- sample(dat, size = sample.n, replace = TRUE) + rkernel(n = sample.n)
      kdevals <- kde(sample)

      if(sample.plot == TRUE){
        hist(sample, freq = FALSE, breaks = 100)
        lines(seq(min(sample), max(sample), length.out = 100), kde(seq(min(sample), max(sample), length.out = 100)))
      }

      # calculate optimal parameters for j0 = 1

      fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)

      opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)


      theta.j0 <- opt$pars
      Hellinger.j0 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j0, j0,  ndistparams, formals.dist,
                                                                        kde, dist, discrete, opt$values[length(opt$values)])
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
      .printresults(opt, j0, dist, formals.dist, ndistparams)
    }


    # calculate optimal parameters for j1 (always need weight restrictions since j1
    # starts from 2)

    fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                  kdevals, dist_call)

    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq.j1 <- restrictions.j1$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)

    opt <- solnp(initial.j1, fun = fmin, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1,
                 LB = lx.j1, UB = ux.j1, control = control)
    theta.j1 <- opt$pars <- .augment.pars(opt$pars, j1)
    Hellinger.j1 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j1, j1,  ndistparams, formals.dist,
                                                                      kde, dist, discrete, opt$values[length(opt$values)])
    conv.j1 <- opt$convergence
    values.j1 <- opt$values

    .printresults(opt, j1, dist, formals.dist, ndistparams)

    if((Hellinger.j0 - Hellinger.j1) < thresh){
      # so that the printed result reflects that the order j.max was actually estimated
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if(j0 == j.max){
      break
    }

  }

  .return.paramEst(j0, j.max, dat, theta.j0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete = discrete, MLE.function)
}



## Purpose: Hellinger distance based method of estimating the mixture complexity of a
##          continuous mixture (as well as the weights and component parameters) returning
##          a 'paramEst' object (using bootstrap)

#' @export hellinger.boot.cont
hellinger.boot.cont <- function(obj, bandwidth, j.max = 10, B = 100, ql = 0.025,
                                qu = 0.975, sample.n = 3000, sample.plot = FALSE,
                                control = c(trace = 0), ...){

  # check relevant inputs
  .input.checks.functions(obj, j.max = j.max,  B = B, ql = ql, qu = qu,
                          assert_continuous = TRUE, Hankel = FALSE, param = TRUE)
  
  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())

  j0 <- 0

  repeat{

    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1

    if(j0 != 1){ # pass on constraints for bootstrap

      ineq.j0 <- ineq.j1
      lx.j0 <- lx.j1
      ux.j0 <- ux.j1

    } else { # need to estimate j0 so need initial values and restrictions

      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)
      restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                           upper = upper)
      lx.j0 <- restrictions.j0$lx
      ux.j0 <- restrictions.j0$ux
    }

    if(bandwidth == "adaptive"){ # use the adaptive Kernel density estimate (kde)

      # will not really be used in the computation of the adaptive kde anyways (for j0 == 1);
      # just supplying something so the same function (ADAPkde) can be used
      if(j0 == 1) theta.adap <- initial.j0

      # compute the kde and draw a sample from it (used to calculate the approximate
      # hellinger distance); has to be recomputed for every j0
      kde.list <- ADAPkde(dat, ndistparams, N, j0, theta.adap, dist, formals.dist, dist_call,
                          sample.n, sample.plot)
      kde <- kde.list$kde
      sample <- kde.list$sample
      kdevals <- kde(sample)

      # calculate optimal parameters for j0; cannot reuse them with the adaptive kde
      # since the kde and kdevals will change in every iteration of j0
      if(j0 > 1){ # need to include weight restrictions in optimization

        initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                       formals.dist)
        fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)

        opt <- solnp(initial.j0, fun = fmin, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                     LB = lx.j0, UB = ux.j0, control = control)
        # if we estimate multiple components check that all weights satisfy the constraints
        theta.j0 <- opt$pars <- .augment.pars(opt$pars, j0)

      } else { # already know w = 1 (single component mixture)

        fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                        kdevals, dist_call)

        opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)
        theta.j0 <- opt$pars
      }

      Hellinger.j0 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j0, j0,  ndistparams, formals.dist,
                                                                        kde, dist, discrete, opt$values[length(opt$values)])
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
      .printresults(opt, j0, dist, formals.dist, ndistparams)

    } else if(j0 != 1) { # and if bandwidth != "adaptive"
      # for the standard kde we can reuse the values calculated for j1 since the kde
      # does not change with j0

      theta.j0 <- theta.j1
      Hellinger.j0 <- Hellinger.j1
      conv.j0 <- conv.j1
      values.j0 <- values.j1

    } else { # if bandwidth != "adaptive" and j0 == 1

      # compute the kde and draw a sample from it (used to calculate the approximate
      # hellinger distance) in the first iteration (i.e. j0 == 1)
      kde <- kdensity(dat, bw = bandwidth, kernel = "gaussian")
      rkernel <- function(n) rnorm(n, sd = bandwidth)
      sample <- sample(dat, size = sample.n, replace = TRUE) + rkernel(n = sample.n)
      kdevals <- kde(sample)

      if(sample.plot == TRUE){
        hist(sample, freq = FALSE, breaks = 100)
        lines(seq(min(sample), max(sample), length.out = 100), kde(seq(min(sample), max(sample), length.out = 100)))
      }

      # calculate optimal parameters for j0 = 1

      fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)

      opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)


      theta.j0 <- opt$pars
      Hellinger.j0 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j0, j0,  ndistparams, formals.dist,
                                                                        kde, dist, discrete, opt$values[length(opt$values)])
      conv.j0 <- opt$convergence
      values.j0 <- opt$values
      .printresults(opt, j0, dist, formals.dist, ndistparams)
    }

    # calculate optimal parameters for j1 (always need weight restrictions since j1
    # starts from 2)

    fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                  kdevals, dist_call)

    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq.j1 <- restrictions.j1$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)

    opt <- solnp(initial.j1, fun = fmin, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1,
                 LB = lx.j1, UB = ux.j1, control = control)
    theta.j1 <- opt$pars <- .augment.pars(opt$pars, j1)
    Hellinger.j1 <- opt$values[length(opt$values)] <- .get.hellingerD(theta.j1, j1,  ndistparams, formals.dist,
                                                                      kde, dist, discrete, opt$values[length(opt$values)])
    conv.j1 <- opt$convergence
    values.j1 <- opt$values

    .printresults(opt, j1, dist, formals.dist, ndistparams)

    diff.0 <- Hellinger.j0 - Hellinger.j1

    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams,
                                            mle.est = theta.j0, j = j0)
    Mix.boot <- Mix(dist = dist, discrete = discrete, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                    name = "Mix.boot")

    ran.gen <- function(dat, mle){
      rMix(n = length(dat), obj = mle)
    }

    # counting bootstrap iterations to print progression
    bs_iter <- - 1

    stat <- function(dat){

      assign("bs_iter", bs_iter + 1, inherits = TRUE)
      if(bs_iter != 0){

        # don't include first iteration as this just uses the original data
        # to calculate t0
        message(paste("Running bootstrap iteration ", bs_iter, " testing for ", j0,
                  " components.\n", sep = ""))

      } else message(paste("\n"))

      # calculate optimal parameters for j0

      initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)


      if(bandwidth == "adaptive"){ # use the adaptive Kernel density estimate (kde)

        # compute the kde and draw a sample from it
        kde.list <- ADAPkde(dat, ndistparams, N, j0, theta.adap, dist, formals.dist, dist_call,
                            sample.n, sample.plot, bs_iter)
        kde <- kde.list$kde
        sample <- kde.list$sample

      } else { # use the standard gaussian Kernel estimate with the supplied bandwidth

        # compute the kde and draw a sample from it
        kde <- kdensity(dat, bw = bandwidth, kernel = "gaussian")
        rkernel <- function(n) rnorm(n, sd = bandwidth)
        sample <- sample(dat, size = sample.n, replace = TRUE) + rkernel(n = sample.n)

        if(sample.plot == TRUE && bs_iter != 0){
          txt <- paste("Sample from Bootstrap KDE: Iteration ", bs_iter, sep = "")
          hist(sample, freq = FALSE, breaks = 100, col = "light grey",
               main = txt, xlab = "Sample")
          lines(seq(min(sample), max(sample), length.out = 100), kde(seq(min(sample), max(sample), length.out = 100)))
        }

      }

      kdevals <- kde(sample)

      if(j0 > 1){ # need to include weight restrictions in optimization

        fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                      kdevals, dist_call)
        opt <- solnp(initial.j0, fun = fmin, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                     LB = lx.j0, UB = ux.j0, control = control)
        # if we estimate multiple components check that all weights satisfy the constraints
        theta.boot0 <- .augment.pars(opt$pars, j0)

      } else { # already know w = 1 (single component mixture)

        fmin <- .get.fmin.hellinger.c.0(kde, dat, formals.dist, ndistparams, dist, sample,
                                        kdevals, dist_call)
        opt <- solnp(initial.j0, fun = fmin, LB = lx.j0, UB = ux.j0, control = control)
        theta.boot0 <- opt$pars
      }

      Hellinger.boot0 <- .get.hellingerD(theta.boot0, j0,  ndistparams, formals.dist, kde, dist, discrete, opt$values[length(opt$values)])

      # calculate optimal parameters for j1 (always need weight restrictions since j1
      # starts from 2)

      fmin <- .get.fmin.hellinger.c(kde, dat, formals.dist, ndistparams, dist, sample,
                                    kdevals, dist_call)

      initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                     formals.dist)

      opt <- solnp(initial.j1, fun = fmin, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1,
                   LB = lx.j1, UB = ux.j1, control = control)

      theta.boot1 <- .augment.pars(opt$pars, j1)
      Hellinger.boot1 <- .get.hellingerD(theta.boot1, j1,  ndistparams, formals.dist, kde, dist, discrete, opt$values[length(opt$values)])

      return(Hellinger.boot0 - Hellinger.boot1)

    }

    bt <- boot(dat, statistic = stat, R = B, sim = "parametric", ran.gen = ran.gen,
               mle = Mix.boot, ...)
    diff.boot <- bt$t

    q_lower <- quantile(diff.boot, probs = ql)
    q_upper <- quantile(diff.boot, probs = qu)

    theta.adap <- theta.j1 # adaptive KDE is based on previous theta.j1 estimate

    if(diff.0 >= q_lower && diff.0 <= q_upper){
      # so that the printed result reflects that the order j.max was actually estimated
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if (j0 == j.max){
      break
    }

  }

  .return.paramEst(j0, j.max, dat, theta.j0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete = discrete, MLE.function)
}
