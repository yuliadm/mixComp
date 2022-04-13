
## Purpose: returns the negative log likelihood function of parameters x
##          (length of x determines the mixture complexity)

.get.negloglik.dist <- function(dat, dist, formals.dist, ndistparams, dist_call){

  function(x){

    j <- (length(x) + 1)/(ndistparams + 1) # current complexity estimate
    n <- length(dat)
    w <- x[1:(j-1)]

    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- matrix(x[(i*j):((1 + i)*j-1)], nrow = n, ncol = j, byrow = TRUE)
    }
    theta.list.long$x <- dat

    # # NAs or warnings can happen as solnp sometimes uses intermediate
    # # parameter values outside of the box constraints (to speed up convergence
    # # and avoid numerical ill conditioning)
    mat <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    w <- c(x[1:(j - 1)], 1 - sum(x[1:(j-1)]))
    if(any(w < 0)) return(0)
    res <- as.matrix(mat) %*% w
    if(any(is.na(res))) return(0)

    # avoid log(0)
    res[res == 0] <- ifelse(suppressWarnings(min(res[res != 0])) < .Machine$double.xmin,
                            min(res[res != 0]), .Machine$double.xmin)

    return(-sum(log(res)))
  }
}



## Purpose: returns the negative log likelihood function of parameters x
##          when the mixture consists of only a single component

.get.negloglik.dist.0 <- function(dat, dist, formals.dist, ndistparams, dist_call){

  function(x){

    n <- length(dat)
    dist_call <- get(paste("d", dist, sep = ""))

    theta.list.long <- vector(mode = "list", length = ndistparams)
    names(theta.list.long) <- formals.dist
    for(i in 1:ndistparams){
      theta.list.long[[i]] <- rep(x[i], n)
    }
    theta.list.long$x <- dat

    # NAs or warnings can happen as solnp sometimes uses intermediate
    # parameter values outside of the box constraints (to speed up convergence
    # and avoid numerical ill conditioning)
    res <- suppressWarnings(do.call(dist_call, args = theta.list.long))
    if(any(is.na(res))) return(0)

    # avoid log(0)
    res[res == 0] <- ifelse(suppressWarnings(min(res[res != 0])) < .Machine$double.xmin,
                            min(res[res != 0]), .Machine$double.xmin)

    return(-sum(log(res)))
  }
}


## Purpose: returns the MLE function contained in a 'datMix' object as a list

.get.MLE.function <- function(obj){

  if(!is.null(attr(obj, "MLE.function"))){
    if(is.list(attr(obj, "MLE.function"))) MLE.function <- attr(obj, "MLE.function")
    else MLE.function <- list(attr(obj, "MLE.function"))
  }
  else MLE.function <- NULL
  return(MLE.function)

}


## Purpose: returns the restrictions set when optimizing the negative log likelihood
##          function

.get.restrictions <- function(j, ndistparams, lower, upper){

  # first j-1 weights need to sum up to less than 1
  ineq <- function(x){
    return(sum(x[1:(j - 1)]))
  }

  # lower bound for all parameters
  lx <- rep(0, j - 1) # lower bound for weights
  for(i in 1:ndistparams){lx <- c(lx, rep(lower[i], j))} # lower bound for component parameters

  # upper bound for all parameters
  ux <- rep(1, j - 1) # upper bound for weights
  for(i in 1:ndistparams){ux <- c(ux, rep(upper[i], j))} # upper bound for component parameters

  return(list("ineq" = ineq, "lx" = lx, "ux" = ux))
}


## Purpose: returns the parameters (weights and component parameters)
##          used in parametric bootstrap

.get.bootstrapparams <- function(formals.dist, ndistparams, mle.est, j)  {

  theta.list <- vector(mode = "list", length = ndistparams)
  names(theta.list) <- formals.dist
  for (i in 1:ndistparams){
    theta.list[[i]] <-  mle.est[(i*j):((1 + i)*j - 1)]
  }

  if(length(mle.est) != ndistparams){
    w <- c(mle.est[1:(j - 1)], 1 - sum(mle.est[1:(j - 1)]))
  } else w <- 1

  return(list("theta.list" = theta.list, "w" = w))
}



## Purpose: returns initial values for the calculation of the initial values
##          passed to solnp (in case numerical optimization has to be used
##          to determine the initial values)

.get.initialvals.init <- function(j, ndistparams, lower, upper){

  lower[lower == -Inf] <- -100
  upper[upper == Inf] <- 100
  lower.new <- lower + (upper-lower)/4 # such that we are in the midde 50% of the range
  upper.new <- upper - (upper-lower)/4 # such that we are in the midde 50% of the range

  if(j != 1) initial <- rep(1/j, (j - 1)) # uniform weights
  else initial <- NULL

  initial <- c(initial, runif(j*ndistparams, rep(lower.new, each = j),
                              rep(upper.new, each = j)))
  initial
}



## Purpose: returns initial values to be passed to solnp

.get.initialvals <- function(dat, j, ndistparams, MLE.function = NULL, lower, upper, dist,
                             formals.dist, dist_call){

  # cluster data into j groups
  cluster.vec <- clara(dat, j, rngR = TRUE, pamLike = TRUE, medoids.x = FALSE)$clustering
  groups <- mapply(function(i){sum(abs(cluster.vec - i) <= .Machine$double.eps)}, 1:j)

  # check that each group has at least 2 observations
  if(any(groups == 1)){
    ind <- which(groups == 1)
    for(i in ind){
      distance <- abs(dat[cluster.vec == i] - dat)
      distance[distance == 0] <- max(distance)
      closest <- which(distance == min(distance))
      cluster.vec[closest] <- i
    }
  }

  # initial weights proportional to number of observations in each cluster
  initial <- (mapply(function(i){sum(abs(cluster.vec - i) <= .Machine$double.eps)}, 1:j)/length(cluster.vec))[-j]

  if(!is.null(MLE.function)){

    # apply MLE function to each of the clusters
    init.mat <- mapply(function(i){sapply(MLE.function,
                                          function(fun) fun(dat[cluster.vec == i]))}, 1:j)

  } else {

    # use numerical optimization to find the MLE for each of the clusters
    init <- .get.initialvals.init(j = 1, ndistparams = ndistparams, lower = lower,
                                  upper = upper)
    likelihood0list <- mapply(function(i){.get.negloglik.dist.0(dat[cluster.vec == i],
                                                                dist, formals.dist, ndistparams,
                                                                dist_call)}, 1:j)

    # solnp at times has a hard time converging if some groups have very few observations:
    # stopping criteria
    setTimeLimit(cpu = 10, elapsed = 10, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
    env <- environment()
    tryCatch(
      init.mat <- mapply(function(i){solnp(init, likelihood0list[[i]], LB = lower, UB = upper,
                                           control = c(trace = 0))$pars}, 1:j),
      error = function(e) assign("init.mat", rep(init, j), envir = env))
  }

  initial.theta <- as.vector(t(init.mat))

  # don't want starting values on the boarder of the feasible region
  initial.theta[initial.theta == lower] <- 1.05*initial.theta[initial.theta == lower]
  initial.theta[initial.theta == lower] <- 0.05 # if lower == 0
  initial.theta[initial.theta == upper] <- 0.95*initial.theta[initial.theta == upper]
  initial.theta[initial.theta == upper] <- -0.05 # if upper == 0

  initial <- c(initial, initial.theta)
  initial
}



## Purpose: brings weights back into feasible region in case solnp goes outside of it
##          (this may happen for convergence and conditioning reasons)

.augment.pars <- function(pars, j){

  if(any(pars[1:(j-1)] < 0)) pars[1:(j - 1)][pars[1:(j - 1)] < 0] <- 1e-08
  if(sum(pars[1:(j-1)]) >= 1) pars[1:(j - 1)] <- pars[1:(j - 1)]/(sum(pars[1:(j - 1)]) + 1e-08)
  pars

}



## Purpose: returns an object of class 'paramEst', which contains the result
##          of the function estimating the mixture complexity (alongside the
##          weights and component parameters)

.return.paramEst <- function(j, j.max, dat, pars, values, conv, dist, ndistparams, formals.dist,
                             discrete, MLE.function){

  if(j < j.max) message(paste("\nThe estimated order is ", j, ".", sep = ""))
  else message(paste("\nThe estimation procedure reached the value 'j.max' equal to ", j.max,
                 ".\nThe result showcases the last parameter estimates.", sep = ""))

  paramEst <- j
  attr(paramEst, "dat") <- as.vector(dat)
  attr(paramEst, "pars") <- pars
  attr(paramEst, "values") <- values
  attr(paramEst, "convergence") <- conv
  attr(paramEst, "dist") <- dist
  attr(paramEst, "ndistparams") <- ndistparams
  attr(paramEst, "formals.dist") <- formals.dist
  attr(paramEst, "discrete") <- discrete
  attr(paramEst, "mle.fct") <- MLE.function
  class(paramEst) <- "paramEst"
  invisible(paramEst)

}


## Purpose: prints the results of the numerical optimization procedure

.printresults <- function(res, j, dist, formals.dist, ndistparams){

  niter <- length(res$values)
  funv <- res$values[niter] # function value
  conv <- res$convergence
  if(conv == 0) conv <- TRUE
  else conv <- FALSE

  pars <- res$pars
  pars_complete <- numeric(j*(ndistparams + 1))
  if(j != 1){
    # add last weight completely determined by the others as they sum to 1
    pars_complete[1:(j - 1)] <- pars[1:(j - 1)]
    pars_complete[j] <- 1-sum(pars[1:(j - 1)])
    pars_complete[(j + 1):length(pars_complete)] <- pars[j:length(pars)]
  } else {
    pars_complete[1] <- 1
    pars_complete[2:length(pars_complete)] <- pars[1:length(pars)]
  }

  cmat <- matrix(pars_complete, nrow = j, ncol = ndistparams + 1)
  colnames(cmat) <- c("w", formals.dist)
  rownames(cmat) <- paste("Component ", 1:j, ":", sep = "")

  message(paste("\nParameter estimation for a ", j, " component '", dist, "' mixture model:",
            "\nFunction value: ",
            format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
            "\n", "\n", sep = ""))
  printCoefmat(cmat)

  if(conv == TRUE){
    message(paste( "Converged in ", niter, " iterations.\n", strrep("-", 70), sep = ""))
  }
  else message("Optimization algorithm did not converge.\n", strrep("-", 70), sep = "")
}



## Purpose: prints the results obtained by applying the MLE function

.printresultsMLE <- function(res, dist, formals.dist, ndistparams, likelihood0){

  pars <- c(1, res) # parameters
  funv <- likelihood0(res) # function value

  cmat <- matrix(pars, nrow = 1, ncol = ndistparams + 1)
  colnames(cmat) <- c("w", formals.dist)
  rownames(cmat) <- paste("Component ", 1, ":", sep = "")

  message(paste("\nParameter estimation for a 1 component '", dist, "' mixture model:",
            "\nFunction value: ",
            format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
            "\n", "\n", sep = ""))
  printCoefmat(cmat)
  message("Optimization via user entered MLE-function.\n", strrep("-", 70), sep = "")

}


## Purpose: returs a list of commonly needed variables; called at the beginning of
##          most function estimating the mixture complexity

.get.list <- function(obj){

  theta.bound.list <- attr(obj, "theta.bound.list")
  formals.dist <- names(theta.bound.list)
  ndistparams <- length(formals.dist)
  bounds <- as.vector(matrix(unlist(theta.bound.list, use.names = FALSE),
                             nrow = ndistparams, byrow = TRUE))
  dat <- as.numeric(obj)
  dist <-  attr(obj, "dist")
  dist_call <- get(paste("d", dist, sep = ""))


  return(list(

    dist = dist,
    theta.bound.list = theta.bound.list,
    formals.dist = formals.dist,
    ndistparams = ndistparams,
    dist_call = dist_call,

    bounds = bounds,
    lower = bounds[1:ndistparams],
    upper = bounds[(ndistparams + 1):(2*ndistparams)],

    dat = dat,
    N = length(dat),
    n.max = max(dat),
    discrete = attr(obj, "discrete"),
    continuous = !attr(obj, "discrete"),

    Hankel.method = attr(obj, "Hankel.method"),
    Hankel.function = attr(obj, "Hankel.function"),
    MLE.function = .get.MLE.function(obj)))

}



## Purpose: method of estimating the mixture complexity (as well as the weights
##          and component parameters) returning a 'paramEst' object

#' @title Estimation of Mixture Complexity (and Component Weights/Parameters) Based on Hankel Matrix Approach
#'
#' @description Estimation method of mixture complexity as well as component weights and parameters based on estimating the determinant of the Hankel matrix of the moments of the mixing distribution and comparing it to determinant values generated by a parametric bootstrap.
#'
#' @details Define \eqn{complexity} of a finite mixture \eqn{F} as the smallest integer \eqn{p}, such that its pdf/pmf \eqn{f} can be written as
#' \deqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p).}
#' The \code{paramHankel} procedure initially assumes that the mixture only contains one component, setting \eqn{j = 1}, then sequentially tests \eqn{p = j} versus \eqn{p = j+1} for \eqn{j = 1,2, \dots}. It determines the MLE for a \eqn{j}-component mixture, generates \code{B} parametric bootstrap samples of size \eqn{n} from the distribution the MLE corresponds to and calculates \code{B} determinants of the corresponding \eqn{(j+1) x (j+1)} Hankel matrices of the first \eqn{2j} raw moments of the mixing distribution (for details see \code{\link{nonparamHankel}}). The null hypothesis \eqn{H_0: p = j} is rejected and \eqn{j} increased by 1 if the determinant value based on the original data lies outside of the interval \eqn{[ql, qu]}, a range specified by \code{ql} and \code{qu}, empirical quantiles of the bootstrapped determinants. Otherwise, \eqn{j} is returned as the complexity estimate.
#' \code{paramHankel.scaled} functions similarly to \code{paramHankel} with the exception that the bootstrapped determinants are scaled by the empirical standard deviation of the bootstrap sample. To scale the original determinant, \code{B} nonparametric bootstrap samples of size \eqn{n} are generated from the data, the corresponding determinants are calculated and their empirical standard deviation is used.
#' The MLEs are calculated via the \code{MLE.function} attribute of the \code{datMix} object \code{obj} for \eqn{j = 1}, if it is supplied. For all other \eqn{j} (and also for \eqn{j = 1} in case \code{MLE.function = NULL}) the solver \code{\link[Rsolnp]{solnp}} is used to calculate the minimum of the negative log-likelihood. The initial values supplied to the solver are calculated as follows: the data is clustered into \eqn{j} groups by the function \code{\link[cluster]{clara}} and the data corresponding to each group is supplied to \code{MLE.function} (if supplied to the \code{datMix} object, otherwise numerical optimization is used). The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates.
#'
#' @aliases paramHankel paramHankel.scaled print.paramEst plot.paramEst
#' @usage
#' paramHankel(obj, j.max = 10, B = 1000, ql = 0.025, qu = 0.975, 
#'             control = c(trace = 0), \dots)
#'
#' paramHankel.scaled(obj, j.max = 10, B = 100, ql = 0.025, qu = 0.975, 
#'                    control = c(trace = 0), \dots)
#' @param obj object of class \code{\link{datMix}}.
#' @param j.max integer stating the maximal number of components to be considered.
#' @param B integer specifying the number of bootstrap replicates.
#' @param ql numeric between \eqn{0} and \eqn{1} specifying the lower bootstrap quantile to which the observed determinant value will be compared.
#' @param qu numeric between \eqn{0} and \eqn{1} specifying the upper bootstrap quantile to which the observed determinant value will be compared.
#' @param control control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.
#' @param x object of class \code{paramEst}.
#' @param mixture logical flag, indicating whether the estimated mixture density should be plotted, set to \code{TRUE} by default.
#' @param components logical flag, indicating whether the individual mixture components should be plotted, set to \code{TRUE} by default.
#' @param ylim range of y values to use; if not specified (or containing \code{NA}), the function tries to construct reasonable default values.
#' @param cex.main magnification to be used for main titles relative to the current setting of \code{cex}, see \code{\link[graphics]{par}}.
#' @param \dots
#'    \describe{
#'       \item{in \code{paramHankel()} and \code{paramHankel.scaled()}:}{further arguments passed to the \code{\link[boot]{boot}} function.}
#'       \item{in \code{plot.hankDet()}:}{further arguments passed to the \code{\link[graphics]{hist}} function plotting the data.}
#'       \item{in \code{print.hankDet()}:}{further arguments passed to the \code{\link[stats]{printCoefmat}} function.}
#'     }
#' @return Object of class \code{paramEst} with the following attributes:
#'       \item{dat}{data based on which the complexity is estimated.}
#'       \item{dist}{character string giving the abbreviated name of the component distribution, such that the function \code{ddist} evaluates its density/mass and \code{rdist} generates random variates.}
#'       \item{ndistparams}{integer specifying the number of parameters identifying the component distribution, i.e. if \eqn{\theta} is in \eqn{R^d} then \code{ndistparams}\eqn{ = d}.}
#'       \item{formals.dist}{string vector specifying the names of the formal arguments identifying the distribution \code{dist} and used in \code{ddist} and \code{rdist}, e.g. for a gaussian mixture (\code{dist = norm}) amounts to \code{mean} and \code{sd}, as these are the formal arguments used by \code{dnorm} and \code{rnorm}.}
#'       \item{discrete}{logicalflag, indicating whether the underlying mixture distribution is discrete.}
#'       \item{mle.fct}{attribute \code{MLE.function} of \code{obj}.}
#'       \item{pars}{Say the complexity estimate is equal to some \eqn{j}. Then \code{pars} is a numeric vector of size \eqn{(d+1)*j-1} specifying the component weight and parameter estimates, given as 
#' \deqn{(w_1, ... w_{j-1}, \theta 1_1, ... \theta 1_j, \theta 2_1, ... \theta d_j).}}
#'       \item{values}{numeric vector of function values gone through during optimization at iteration \eqn{j}, the last entry being the value at the optimum.}
#'       \item{convergence}{indicates whether the solver has converged (0) or not (1 or 2) at iteration \eqn{j}.}
#'
#' @seealso \code{\link{nonparamHankel}} for estimation of the mixture complexity based on the Hankel matrix without parameter estimation, \code{\link[Rsolnp]{solnp}} for the solver, \code{\link{datMix}} for creation of the \code{datMix} object.
#' @examples
#'
#' ## create 'Mix' object
#' poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
#'
#' ## create random data based on 'Mix' object (gives back 'rMix' object)
#' set.seed(1)
#' poisRMix <- rMix(1000, obj = poisMix)
#'
#' ## create 'datMix' object for estimation
#' # generate list of parameter bounds
#' poisList <- list("lambda" = c(0, Inf))
#'
#' # generate MLE function
#' MLE.pois <- function(dat){
#'   mean(dat)
#' }
#'
#' # generate function needed for estimating the j^th moment of the
#' # mixing distribution via Hankel.method "explicit"
#'
#' explicit.pois <- function(dat, j){
#'   res <- 1
#'   for (i in 0:(j-1)){
#'     res <- res*(dat-i)
#'   }
#'   return(mean(res))
#' }
#'
#' # generating 'datMix' object
#' pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, MLE.function = MLE.pois,
#'                   Hankel.method = "explicit", Hankel.function = explicit.pois)
#'
#'
#' ## complexity and parameter estimation
#' \donttest{
#' set.seed(1)
#' res <- paramHankel(pois.dM)
#' plot(res)
#' }
#'
#' @rdname paramHankel
#' @export paramHankel
paramHankel <- function(obj, j.max = 10, B = 1000, ql = 0.025,  qu = 0.975,
                        control = c(trace = 0), ...){
  
  # check relevant inputs
  .input.checks.functions(obj, B = B, j.max = j.max, ql = ql, qu = qu,
                          Hankel = TRUE, param = TRUE)

  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())

  fun <- .moments.map(Hankel.method = Hankel.method)
  likelihood <- .get.negloglik.dist(dat, dist, formals.dist, ndistparams, dist_call)
  likelihood0 <- .get.negloglik.dist.0(dat, dist, formals.dist, ndistparams, dist_call)
  j <- 0

  repeat{

    j <- j + 1 # current complexity estimate

    if (j == 1 && !is.null(MLE.function)){ # Calculate MLE via the MLE function
      # only possible for j equal to 1 and if a MLE function was supplied

      mle.est <- sapply(MLE.function, function(fun) fun(dat))
      .printresultsMLE(mle.est, dist, formals.dist, ndistparams, likelihood0)
      values <- likelihood0(mle.est)
      conv <- NULL

    } else { # Calculate MLE of a j component mixture numerically

      restrictions <- .get.restrictions(j = j, ndistparams = ndistparams, lower = lower,
                                        upper = upper)
      ineq <- restrictions$ineq
      lx <- restrictions$lx
      ux <- restrictions$ux

      initial <- .get.initialvals(dat, j, ndistparams, MLE.function, lower, upper,
                                  dist, formals.dist, dist_call)

      if(j != 1){ # need to include weight restrictions in optimization

        opt <- solnp(initial, fun = likelihood, ineqfun = ineq, ineqLB = 0, ineqUB = 1,
                     LB = lx, UB = ux, control = control)
        opt$pars <- .augment.pars(opt$pars, j)
        opt$values[length(opt$values)] <- likelihood(opt$pars)

        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)

      } else { # already know w = 1 (single component mixture)

        opt <- solnp(initial, fun = likelihood0, LB = lx, UB = ux, control = control)
        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)

      }
    }

    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams,
                                            mle.est = mle.est, j = j)
    Mix.mle <- Mix(dist = dist, discrete = discrete, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                   name = "Mix-boot")

    ran.gen <- function(dat, mle){
      rMix(n = length(dat), obj = mle)
    }

    bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
               j = j, sim = "parametric", ran.gen = ran.gen, mle = Mix.mle, ...)
    det0 <- bt$t0 # original determinant value
    det.boot <- bt$t # bootstrapped determinant values

    q_lower <- quantile(det.boot, probs = ql)
    q_upper <- quantile(det.boot, probs = qu)

    if(det0 >= q_lower && det0 <= q_upper){
      # so that the printed result reflects that the order j.max was actually estimated
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if(j == j.max){
      break
    }
  }

  .return.paramEst(j, j.max, dat, mle.est, values, conv, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}


## Purpose: method of estimating the mixture complexity (as well as the weights
##          and component parameters) returning a 'paramEst' object
#' @export paramHankel.scaled
paramHankel.scaled <- function(obj, j.max = 10, B = 100, ql = 0.025, qu = 0.975,
                               control = c(trace = 0), ...){

  # check relevant inputs
  .input.checks.functions(obj, B = B, j.max = j.max, ql = ql, qu = qu, Hankel = TRUE,
                          param = TRUE)

  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())

  fun <- .moments.map(Hankel.method = Hankel.method)
  likelihood <- .get.negloglik.dist(dat, dist, formals.dist, ndistparams, dist_call)
  likelihood0 <- .get.negloglik.dist.0(dat, dist, formals.dist, ndistparams, dist_call)
  j <- 0

  repeat{

    j <- j + 1 # current complexity estimate

    if (j == 1 && !is.null(MLE.function)){ # Calculate MLE via the MLE function
      # only possible for j equal to 1 and if a MLE function was supplied

      mle.est <- sapply(MLE.function, function(fun) fun(dat))
      .printresultsMLE(mle.est, dist, formals.dist, ndistparams, likelihood0)
      values <- likelihood0(mle.est)
      conv <- NULL

    } else { # Calculate MLE of a j component mixture numerically

      restrictions <- .get.restrictions(j = j, ndistparams = ndistparams, lower = lower,
                                        upper = upper)
      ineq <- restrictions$ineq
      lx <- restrictions$lx
      ux <- restrictions$ux
      initial <- .get.initialvals(dat, j, ndistparams, MLE.function, lower, upper,
                                  dist, formals.dist, dist_call)

      if(j != 1){ # need to include weight restrictions in optimization

        opt <- solnp(initial, fun = likelihood, ineqfun = ineq, ineqLB = 0, ineqUB = 1,
                     LB = lx, UB = ux, control = control)
        opt$pars <- .augment.pars(opt$pars, j)
        opt$values[length(opt$values)] <- likelihood(opt$pars)

        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)

      } else { # already know w = 1 (single component mixture)

        opt <- solnp(initial, fun = likelihood0, LB = lx, UB = ux, control = control)
        mle.est <- opt$pars
        values <- opt$values
        conv <- opt$convergence
        .printresults(opt, j, dist, formals.dist, ndistparams)

      }
    }

    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams,
                                            mle.est = mle.est, j = j)
    Mix.mle <- Mix(dist = dist, discrete = discrete, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                   name = "Mix.boot")

    ran.gen <- function(dat, mle){
      rMix(n = length(dat), obj = mle)
    }

    # parametric bootstrap for bootstrapped scaled determinants
    bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
               j = j, sim = "parametric", ran.gen = ran.gen, mle = Mix.mle, ...)
    det0 <- bt$t0 # original determinant value
    det.boot.param <- bt$t # bootstrapped determinants
    det.boot.scaled <- det.boot.param/(sd(det.boot.param))

    q_lower <- quantile(det.boot.scaled, probs = ql)
    q_upper <- quantile(det.boot.scaled, probs = qu)

    # non parametric bootstrap for scaling of det0
    bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
               j = j, ...)
    det.boot.np <- bt$t
    det0.scaled <- det0/(sd(det.boot.np))

    if(det0.scaled >= q_lower && det0.scaled <= q_upper){
      # so that the printed result reflects that the order j.max was actually estimated
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if(j == j.max){
      break
    }

  }

  .return.paramEst(j, j.max, dat, mle.est, values, conv, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}



# Purpose: plot method for 'paramEst' objects
#' @rdname paramHankel
#' @method plot paramEst
#' @export
plot.paramEst <- function(x, mixture = TRUE, components = TRUE, ylim = NULL, cex.main = 0.9,
                          ...){

  obj <- x
  dat <- attr(obj, "dat")
  class(dat) <- "rMix"
  n <- length(dat)
  pars <- attr(obj, "pars")
  dist <- attr(obj, "dist")
  dist_call <- get(paste("d", dist, sep = ""))
  ndistparams <- attr(obj, "ndistparams")
  formals.dist <- attr(obj, "formals.dist")
  discrete <- attr(obj, "discrete")

  # construct reasonable x values
  if(discrete == TRUE){
    x <- seq(min(dat), max(dat), by = 1)
  } else {
    x <- seq(min(dat), max(dat), by = 0.01)
  }

  theta.list.long <- vector(mode = "list", length = ndistparams)
  names(theta.list.long) <- formals.dist
  for(i in 1:ndistparams){
    theta.list.long[[i]] <- matrix(pars[(i*obj):((1 + i)*obj - 1)], nrow = length(x),
                                   ncol = obj, byrow = TRUE)
  }
  theta.list.long$x <- x

  if(length(pars) == ndistparams){ # i.e. the estimated mixture complexity equals 1
    w <- 1
    y <- do.call(dist_call, theta.list.long)
  } else {
    w <- c(pars[1:(obj - 1)], 1-sum(pars[1:(obj-1)]))
    y <- do.call(dist_call, theta.list.long) %*% w
  }

  if(is.null(ylim) || anyNA(ylim)){ # construct reasonable ordinate values
    plt.inv <- suppressWarnings(plot(dat, plot = FALSE, right = FALSE))
    mx <- max(plt.inv$density, y)
    ylim <- c(0, mx)
  }

  cols <- .get.colors(alpha = 0.9)[c(6, 7)]
  main <- paste("Estimated ", obj, " component '", dist, "' mixture model", sep = "")
  plot(dat, freq = FALSE, ylim = ylim, main = main, xlab = "Data", border = "darkgrey",
       components = FALSE, cex.main = cex.main, ...)

  if(components == TRUE){ # plot components

    y.comp <- matrix(w, nrow = length(x), ncol = obj,
                byrow = TRUE) * do.call(dist_call, theta.list.long)

    mapply(function(i) lines(x = x, y = y.comp[, i], col = cols[1], lwd = 1, lty = 2),
           i = 1:obj)
  }

  if(mixture == TRUE){ # plot mixture

    if(discrete == TRUE) lines(theta.list.long$x, y, col = cols[2], lwd = 1, type = "b",
                               pch = 20, cex = 0.70)
    else lines(theta.list.long$x, y, col = cols[2], lwd = 1)

  }

}



# Purpose: print method for 'paramEst' objects
#' @rdname paramHankel
#' @method print paramEst
#' @export
print.paramEst <- function(x, ...){

  obj <- x
  j <- as.numeric(obj)
  pars <- attr(obj, "pars")
  values <- attr(obj, "values")
  conv <- attr(obj, "conv")
  dist <- attr(obj, "dist")
  ndistparams <- attr(obj, "ndistparams")
  formals.dist <- attr(obj, "formals.dist")
  mle.fct <- attr(obj, "mle.fct")
  niter <- length(values)
  funv <- values[niter]

  if(!is.null(mle.fct) & j == 1){ # parameters estimated via MLE.function

    pars <- c(1, pars)

    cmat <- matrix(pars, nrow = 1, ncol = ndistparams + 1)
    colnames(cmat) <- c("w", formals.dist)
    rownames(cmat) <- paste("Component ", 1, ":", sep = "")

    message(paste("\nParameter estimation for a 1 component '", dist, "' mixture model:",
              "\nFunction value: ",
              format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
              "\n", "\n", sep = ""))
    printCoefmat(cmat, ...)
    message("Optimization via user entered MLE-function.\n", strrep("-", 70), sep = "")

  } else { # parameters estimated numerically

    if(conv == 0) conv <- TRUE
    else conv <- FALSE

    # add last weight completely determined by the others as they sum to 1
    pars_complete <- numeric(j*(ndistparams + 1))
    if(j != 1){
      pars_complete[1:(j-1)] <- pars[1:(j-1)]
      pars_complete[j] <- 1-sum(pars[1:(j-1)])
      pars_complete[(j+1):length(pars_complete)] <- pars[j:length(pars)]
    } else {
      pars_complete[1] <- 1
      pars_complete[2:length(pars_complete)] <- pars[1:length(pars)]
    }

    cmat <- matrix(pars_complete, nrow = j, ncol = ndistparams +1)
    colnames(cmat) <- c("w", formals.dist)
    rownames(cmat) <- paste("Component ", 1:j, ":", sep = "")

    message(paste("\nParameter estimation for a ", j, " component '", dist, "' mixture model:",
              "\nFunction value: ",
              format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE),
              "\n", "\n", sep = ""))
    printCoefmat(cmat, ...)

    if(conv == TRUE){
      message(paste( "Converged in ", niter, " iterations.\n", sep = ""))
    }
    else message("Optimization algorithm did not converge.\n", sep = "")
  }
  message(paste("\nThe estimated order is ", j, ".", sep = ""))

}
