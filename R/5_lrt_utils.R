
## Purpose: LRT based method of estimating the mixture complexity (as well as
##          the weights and component parameters) returning a 'paramEst' object
##          (using bootstrap)

#' @title Estimation of a Mixture Complexity Based on Likelihood Ratio Test Statistics
#'
#' @description Estimation of a mixture complexity as well as its component weights and parameters based on comparing the likelihood ratio test statistic (LRTS) to a bootstrapped quantile.
#'
#' @details Define the \eqn{complexity} of a finite mixture \eqn{F} as the smallest integer \eqn{p}, such that its pdf/pmf \eqn{f} can be written as
#' \deqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p).}
#' To estimate \eqn{p}, \code{mix.lrt} sequentially tests \eqn{p = j} versus \eqn{p = j+1} for \eqn{j = 1,2, \dots}, by finding the maximum likelihood estimator (MLE) for the density of a mixture with \eqn{j} and \eqn{j+1} components and calculating the corresponding likelihood ratio test statistic (LRTS). Next, a parametric bootstrap procedure is used to generate \code{B} samples of size \eqn{n} from a \eqn{j}-component mixture given the previously calculated MLE. For each of the bootstrap samples, the MLEs corresponding to densities of mixtures with \eqn{j} and \eqn{j+1} components are calculated, as well as the LRTS. The null hypothesis \eqn{H_0: p = j} is rejected and \eqn{j} increased by 1 if the LRTS based on the original data is larger than the chosen \code{quantile} of its bootstrapped counterparts. Otherwise, \eqn{j} is returned as the complexity estimate.
#' The MLEs are calculated via the \code{MLE.function} attribute (of the \code{\link{datMix}} object \code{obj}) for \eqn{j = 1}, if it is supplied. For all other \eqn{j} (and also for \eqn{j = 1} in case \code{MLE.function = NULL}) the solver \code{\link[Rsolnp]{solnp}} is used to calculate the minimum of the negative log-likelihood. The initial values supplied to the solver are calculated as follows: the data is clustered into \eqn{j} groups by the function \code{\link[cluster]{clara}} and the data corresponding to each group is given to \code{MLE.function} (if supplied to the \code{datMix} object, otherwise numerical optimization is used here as well). The size of the groups is taken as initial component weights and the MLE's are taken as initial component parameter estimates.
#' @usage
#' mix.lrt(obj, j.max = 10, B = 100, quantile = 0.95, control = c(trace = 0), ...)
#' @param obj object of class \code{\link{datMix}}.
#' @param j.max integer, giving the maximal complexity to be considered.
#' @param B integer, specifying the number of bootstrap replicates.
#' @param quantile numeric between \eqn{0} and \eqn{1} specifying the bootstrap quantile to which the observed LRTS will be compared.
#' @param control control list of optimization parameters, see \code{\link[Rsolnp]{solnp}}.
#' @param \dots further arguments passed to the \code{\link[boot]{boot}} function.
#' @return Object of class \code{paramEst} with the following attributes:
#'      \item{dat}{data based on which the complexity is estimated.}
#'      \item{dist}{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density/ mass function and \code{rdist} generates random variates.}
#'      \item{ndistparams}{integer specifying the number of parameters identifying the component distribution, i.e. if \eqn{\theta} is in \eqn{R^d} then \code{ndistparams}\eqn{ = d}.}
#'      \item{formals.dist}{string vector, specifying the names of the formal arguments identifying the distribution \code{dist} and used in \code{ddist} and \code{rdist}, e.g. for a gaussian mixture (\code{dist = norm}) amounts to \code{mean} and \code{sd}, as these are the formal arguments used by \code{dnorm} and \code{rnorm}.}
#'      \item{discrete}{logical indicating whether the underlying mixture distribution is discrete.}
#'      \item{mle.fct}{attribute \code{MLE.function} of \code{obj}.}
#'      \item{pars}{Say the complexity estimate is equal to some \eqn{j}. Then \code{pars} is a numeric vector of size \eqn{(d+1)*j-1} specifying the component weight and parameter estimates, given as
#' \deqn{(w_1, ... w_{j-1}, \theta 1_1, ... \theta 1_j, \theta 2_1, ... \theta d_j).}}
#'      \item{values}{numeric vector of function values gone through during optimization at iteration \eqn{j}, the last entry being the value at the optimum.}
#'      \item{convergence}{integer indicating whether the solver has converged (0) or not (1 or 2) at iteration \eqn{j}.}
#' @seealso \code{\link[Rsolnp]{solnp}} for the solver, \code{\link{datMix}} for the creation of the \code{datMix} object.
#' @keywords cluster
#' @examples
#'
#' ### generating 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#'
#' ### generating 'rMix' from 'Mix' object (with 1000 observations)
#' set.seed(0)
#' normLocRMix <- rMix(1000, normLocMix)
#'
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
#' ### complexity and parameter estimation
#' \donttest{
#' set.seed(0)
#' res <- mix.lrt(normLoc.dM, B = 30)
#' plot(res)
#' }
#'
#' @export mix.lrt
mix.lrt <- function(obj, j.max = 10, B = 100, quantile = 0.95, control = c(trace = 0), ...){

  # check relevant inputs
  .input.checks.functions(obj, j.max = j.max, quantile = quantile, Hankel = FALSE,
                          param = TRUE)

  # get standard variables
  variable_list <- .get.list(obj)
  list2env(variable_list, envir = environment())


  likelihood0 <- .get.negloglik.dist.0(dat, dist, formals.dist, ndistparams, dist_call)
  likelihood <- .get.negloglik.dist(dat = dat, dist = dist, formals.dist = formals.dist,
                                    ndistparams = ndistparams, dist_call)
  j0 <- 0

  repeat{

    j0 <- j0 + 1 # current complexity estimate
    j1 <- j0 + 1

    if(j0 > 1){ # if j1 was calculated in the last interation, pass it over to j0...

      mle.est0 <- mle.est1
      L0 <- L1
      conv.j0 <- conv.j1
      values.j0 <- values.j1

      # also need to pass over the restrictions as they will be used in the bootstrap
      ineq.j0 <- ineq.j1
      lx.j0 <- lx.j1
      ux.j0 <- ux.j1

    } else {  # ... or calculate j0 directly if j0 = 1 (j1 has not been calculated yet)
      # in this case we already know w = 1 (single component mixture)

      restrictions.j0 <- .get.restrictions(j = j0, ndistparams = ndistparams, lower = lower,
                                           upper = upper)
      ineq.j0 <- restrictions.j0$ineq
      lx.j0 <- restrictions.j0$lx
      ux.j0 <- restrictions.j0$ux

      if (!is.null(MLE.function)){ # Calculate MLE via the MLE function

        mle.est0 <- sapply(MLE.function, function(fun) fun(dat))
        L0 <- likelihood0(mle.est0)
        .printresultsMLE(mle.est0, dist, formals.dist, ndistparams, likelihood0)
        conv.j0 <- NULL
        values.j0 <- L0

      } else { # Calculate MLE of a j component mixture numerically

        initial.j0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                       formals.dist)

        opt <- solnp(initial.j0, fun = likelihood0, ineqfun = ineq.j0, ineqLB = 0, ineqUB = 1,
                     LB = lx.j0, UB = ux.j0, control = control)

        .printresults(opt, j0, dist, formals.dist, ndistparams)
        mle.est0 <- opt$pars
        L0 <- likelihood0(mle.est0)
        conv.j0 <- opt$convergence
        values.j0 <- opt$values

      }
    }

    # optimization for j1. Starts from j1 = 2 so we always need to include weight
    # restrictions in optimization

    restrictions.j1 <- .get.restrictions(j = j1, ndistparams = ndistparams, lower = lower,
                                         upper = upper)
    ineq.j1 <- restrictions.j1$ineq
    lx.j1 <- restrictions.j1$lx
    ux.j1 <- restrictions.j1$ux
    initial.j1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper, dist,
                                   formals.dist)

    opt <- solnp(initial.j1, fun = likelihood, ineqfun = ineq.j1, ineqLB = 0, ineqUB = 1,
                 LB = lx.j1, UB = ux.j1, control = control)
    mle.est1 <- opt$pars <- .augment.pars(opt$pars, j1)
    L1 <- opt$values[length(opt$values)] <- likelihood(opt$pars)
    conv.j1 <- opt$convergence
    values.j1 <- opt$values

    .printresults(opt, j1, dist, formals.dist, ndistparams)

    LRT <- 2*(L0 - L1)

    # parameters used for parametric bootstrap and corresponding 'Mix' object
    param.list.boot <- .get.bootstrapparams(formals.dist = formals.dist, ndistparams = ndistparams,
                                            mle.est = mle.est0, j = j0)
    Mix.mle <- Mix(dist = dist, discrete = discrete, w = param.list.boot$w, theta.list = param.list.boot$theta.list,
                   name = "Mix.boot")

    ran.gen <- function(dat, mle){
      rMix(n = length(dat), obj = mle)
    }

    # counting bootstrap iterations to print progression
    bs_iter <- -1

    stat <- function(dat){

      assign("bs_iter", bs_iter + 1, inherits = TRUE)
      if(bs_iter != 0){

        # don't include first iteration as this just uses the original data
        # to calculate t0
        message(paste("Running bootstrap iteration ", bs_iter, " testing for ", j0,
                  " components.\n", sep = ""))

      } else message(paste("\n"))

      # in the bootstrap we have to calculate the values for j0 and j1 as the bootstrap
      # data changes in every iteration (cannot reuse last j1 values as j0)

      likelihood_boot <- .get.negloglik.dist(dat = dat, dist = dist, formals.dist = formals.dist,
                                             ndistparams = ndistparams, dist_call)

      # calculate optimal parameters for j0

      if (j0 == 1){  # already know w = 1 (single component mixture)

        likelihood_boot0 <- .get.negloglik.dist.0(dat = dat, dist = dist, formals.dist = formals.dist,
                                                  ndistparams = ndistparams, dist_call)

        if(!is.null(MLE.function)){ # Calculate MLE via the MLE function

          opt.boot0 <- sapply(MLE.function, function(fun) fun(dat))
          L0.boot <- likelihood_boot0(opt.boot0)

        } else {  # Calculate MLE of a j0 = 1 component mixture numerically

          initial.boot0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                            formals.dist)

          opt.boot0 <- solnp(initial.boot0, fun = likelihood_boot0, LB = lx.j0, UB = ux.j0,
                             control = control)$values
          L0.boot <- opt.boot0[length(opt.boot0)]

        }

      } else { # need to include weight restrictions in optimization for j0 != 1

        initial.boot0 <- .get.initialvals(dat, j0, ndistparams, MLE.function, lower, upper, dist,
                                          formals.dist)
        opt.boot0 <- solnp(initial.boot0, fun = likelihood_boot, ineqfun = ineq.j0, ineqLB = 0,
                           ineqUB = 1, LB = lx.j0, UB = ux.j0, control = control)
        opt.boot0$pars <- .augment.pars(opt.boot0$pars, j0)
        L0.boot <- likelihood_boot(opt.boot0$pars)

      }

      # calculate optimal parameters for j1 (always need restrictions since j1
      # starts from 2)

      initial.boot1 <- .get.initialvals(dat, j1, ndistparams, MLE.function, lower, upper,
                                        dist, formals.dist)

      opt.boot1 <- solnp(initial.boot1, fun = likelihood_boot, ineqfun = ineq.j1, ineqLB = 0,
                         ineqUB = 1, LB = lx.j1, UB = ux.j1, control = control)
      opt.boot1$pars <- .augment.pars(opt.boot1$pars, j1)
      L1.boot <- likelihood_boot(opt.boot1$pars)

      return(2*(L0.boot - L1.boot))
    }


    bt <- boot(dat, statistic = stat, R = B, sim = "parametric", ran.gen = ran.gen,
               mle = Mix.mle, ...)
    LRT.boot <- bt$t

    qu <- quantile(LRT.boot, quantile)

    if(LRT <= qu){
      # so that the printed result reflects that the order j.max was actually estimated
      # rather than just returned as the default
      j.max <- j.max + 1
      break
    } else if (j0 == j.max){
      break
    }
  }

  .return.paramEst(j0, j.max, dat, mle.est0, values.j0, conv.j0, dist, ndistparams, formals.dist,
                   discrete, MLE.function)
}
