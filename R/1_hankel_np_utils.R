

## Purpose: calculate the determinant of the Hankel matrix of the moments of the
##          mixing distribution via the method "explicit"

.deterHankel.explicit <- function(dat, inds = 1:length(dat), Hankel.function, j.max = NULL,
                                  j = NULL){ # wierd argument order for bootstrap
  # need j.max for nonparamHankel, j for paramHankel

  dat <- dat[inds] # for bootstrap

  if(is.null(j.max)){ # only calculate determinant for a single complexity estimate
    mn <- j
    mx <- j
  } else { # calculate determinant for all complexity estimates up to j.max
    det.vec <- numeric(j.max)
    mn <- 1
    mx <- j.max
  }

  Hankel.functionV <- Vectorize(Hankel.function, vectorize.args = "j")
  cp_moments <- Hankel.functionV(dat, 1:(2*mx))
  cp_add1 <- c(1, cp_moments)

  for(i in mn:mx){

    H <- hankel.matrix((i + 1), cp_add1[1:((2*i)+1)])
    if(is.null(j.max)){ # return single determinant
      det.vec <- det(H)
    } else det.vec[i] <- det(H) # return vector of determinants

  }

  return(det.vec)
}


## Purpose: calculate the determinant of the Hankel matrix of the moments of the
##          mixing distribution via the method "translation"

.deterHankel.translation <- function(dat, inds = 1:length(dat), Hankel.function, j.max = NULL,
                                     j = NULL){ # wierd argument order for bootstrap
  # need j.max for nonparamHankel, j for paramHankel

  dat <- dat[inds] # for bootstrap
  n <- length(dat)

  if(is.null(j.max)){ # only calculate determinant for a single complexity estimate
    mn <- j
    mx <- j
  } else { # calculate determinant for all complexity estimates up to j.max
    det.vec <- numeric(j.max)
    mn <- 1
    mx <- j.max
  }

  # construct elements of triagular linear system
  # empirical moments of mixture distribution
  X_hat <- (1/n) * mapply(function(x, y){sum(x^y)}, rep(list(dat), 2*mx), 1:(2*mx))
  Hankel.functionV <- Vectorize(Hankel.function)
  # theoretical moments of G
  EY <- Hankel.functionV(1:(2*mx))

  b <- X_hat - EY # vector b
  EY_A <- c(1, EY)

  ## generate matrix A
  A <- matrix(0, nrow = 2*mx, ncol = 2*mx)
  for (i_row in 1:(2*mx)){
    for (i_col in 1:i_row){
      A[i_row, i_col] <- choose(i_row, i_row-i_col)*EY_A[i_row-i_col+1]
    }
  }

  # solve system
  cp_moments <- solve(a = A, b = b)

  cp_add1 <- c(1, cp_moments)
  for(i in mn:mx){

    H <- hankel.matrix((i + 1), cp_add1)
    if(is.null(j.max)){ # return single determinant
      det.vec <- det(H)
    } else det.vec[i] <- det(H) # return vector of determinants

  }

  return(det.vec)
}



## Purpose: calculate the determinant of the Hankel matrix of the moments of the
##          mixing distribution via the method "scale"

.deterHankel.scale <- function(dat, Hankel.function, inds = 1:length(dat), j.max = NULL, j = NULL,
                               message = FALSE){ # wierd argument order for bootstrap
  # need j.max for nonparamHankel, j for paramHankel

  dat <- dat[inds] # for bootstrap
  n <- length(dat)

  if(is.null(j.max)){ # only calculate determinant for a single complexity estimate
    mn <- j
    mx <- j
  } else { # calculate determinant for all complexity estimates up to j.max
    det.vec <- numeric(j.max)
    mn <- 1
    mx <- j.max
  }

  # construct elements of triagular linear system
  # empirical moments of mixture distribution
  X_hat <- (1/n) * mapply(function(x, y){sum(x^y)}, rep(list(dat), 2*mx), 1:(2*mx))
  Hankel.functionV <- Vectorize(Hankel.function)
  # theoretical moments of G
  EY <- Hankel.functionV(1:(2*mx))

  if (sum(EY == 0) != 0){  # cannot devide by zero, in that case take squares of everything

    if(message == TRUE) message("Moments of mixing distribution M^k are 0 for some k. Calculating M^(2k) instead.")

    X_hat <- (1/n) * mapply(function(x,y){sum(x^y)}, rep(list(dat),2*mx),
                            seq(from = 2, to = (4*mx), by = 2))
    EY <- Hankel.functionV(seq(from = 2, to = (4*mx), by = 2))

  }

  cp_moments <- X_hat/EY
  cp_add1 <- c(1, cp_moments)

  for(i in mn:mx){

    H <- hankel.matrix((i + 1), cp_add1)
    if(is.null(j.max)){ # return single determinant
      det.vec <- det(H)
    } else det.vec[i] <- det(H) # return vector of determinants

  }

  return(det.vec)
}


## Purpose: return the correct function for calculating the determinant based on
##          "Hankel.method"

.moments.map <- function(Hankel.method) {

  if (Hankel.method == "explicit"){
    fun <- .deterHankel.explicit
  }
  else if (Hankel.method == "translation"){
    fun <- .deterHankel.translation
  }
  else if (Hankel.method == "scale"){
    fun <- .deterHankel.scale
  }

  fun
}



## Purpose: check whether the inputs to any of the functions estimating mixtuer complexity
##          are of the correct form

.input.checks.functions <- function(obj, B, j.max, ql, qu, quantile, n.inf, thrshL2,
                                    thrshHel, pen.function, scaled, Hankel, param,
                                    bandwidth, sample.n, sample.plot, assert_discrete,
                                    assert_continuous){

  if(Hankel == TRUE){
    # if the function uses the Hankel matrix of the moments of the mixing distribution

    if(is.null(attr(obj, "Hankel.method")))
      stop("'Hankel.method' has to be defined as a datMix object attribute!")
    if(!attr(obj, "Hankel.method") %in% c("explicit", "translation", "scale")){
      stop("Hankel.method is not one of the implemented methods.
           Please enter \"explicit\", \"translation\" or \"scale\"")
    }

    if(is.null(attr(obj, "Hankel.function")))
      stop("'Hankel.function' has to be defined as a datMix object attribute!")
    if(!is.function(attr(obj, "Hankel.function")))
      stop("'Hankel.function' has to be a function!")
  }

  if(param == TRUE){
    # if the functions estimates the weights as well as the component parameters
    if(is.null(attr(obj, "theta.bound.list")))
      stop("'theta.bound.list' has to be specified as a datMix object attribute!")
    
    if(is.null(attr(obj, "discrete"))) stop("'discrete' has to be specified as a datMix object attribute!")
    
    if(missing(assert_discrete)){} # allowed be missing, needed for discrete minimum distance
    else if(attr(obj, "discrete") != TRUE)
      stop("The 'discrete' attribute of the Rdat object is not set to TRUE, however
         this function only works for discrete data.")
    
    if(missing(assert_continuous)){} # allowed be missing, needed for continuous minimum distance
    else if(attr(obj, "discrete") != FALSE)
      stop("The 'discrete' attribute of the Rdat object is not set to FALSE, however
         this function only works for continuous data.")
  }

  if(is.null(obj) || !is.datMix(obj)) stop("'obj' has to be a 'datMix' object!")

  if(is.null(j.max) || !is.numeric(j.max) || !(as.integer(j.max) == j.max))
    stop("'j.max' has to be an integer!")

  if(missing(pen.function)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(!is.null(pen.function) && !is.function(pen.function))
    stop("'pen.function' has to be NULL or a function!")
  else if(!is.null(pen.function) && (!all(names(formals(pen.function)) %in% c("n", "j"))
                                     || !all(c("n", "j") %in% names(formals(pen.function)))))
    stop("'pen.function' must contain the arguments \"n\" and \"j\".")

  if(missing(scaled)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(scaled) || !(is.logical(scaled))) stop("'scaled' has to be logical!")

  if(missing(B)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if (is.null(B) || !is.numeric(B) || !(as.integer(B) == B)) stop("'B' has to be an integer!")

  if(missing(ql) & missing(qu)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(ql) ||  is.null(qu) || !is.numeric(ql) || !is.numeric(qu) ||
          !(ql >= 0 & ql <= 1) || !(qu >= 0 & qu <= 1) || ql > qu)
    stop("'ql' and 'qu' have to be numerics between 0 and 1 with qu > ql!")

  if(missing(quantile)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(!is.numeric(quantile) || !(quantile >= 0 & quantile <= 1))
    stop("'quantile' has to be a numeric between 0 and 1!")

  if(missing(n.inf)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(n.inf) || !is.numeric(n.inf) || !(as.integer(n.inf) == n.inf))
    stop("'n.inf' has to be an integer!")

  if(missing(thrshL2)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(thrshL2) || !(is.function(thrshL2)|| thrshL2 %in% c("LIC", "SBC")))
    stop("'threshold' has to be a user-entered function or an element of c(\"LIC\", \"SBC\")")
  else if(is.function(thrshL2) && (!all(names(formals(thrshL2)) %in% c("n", "j"))
                                   || !all(c("n", "j") %in% names(formals(thrshL2)))))
    stop("The function 'threshold' needs to have arguments \"n\" and \"j\".")
  else if(!is.function(thrshL2) && thrshL2 == "LIC")
    warning("While being used in Umashanger & Sriram's original paper,  asymptotically, the threshold 'LIC' does not go to 0 slower than the difference in squared L2 distances once the correct order p is reached and the estimator is therefore not consistent.",
            call. = FALSE)

  if(missing(thrshHel)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(thrshHel) || !(is.function(thrshHel)|| thrshHel %in% c("AIC", "SBC")))
    stop("'threshold' has to be a user-entered function or an element of c(\"AIC\", \"SBC\")")
  else if(is.function(thrshHel) && (!all(names(formals(thrshHel)) %in% c("n", "j"))
                                    || !all(c("n", "j") %in% names(formals(thrshHel)))))
    stop("The function 'threshold' needs to have arguments \"n\" and \"j\".")
  else if(!is.function(thrshHel) && thrshHel == "AIC")
    warning("While being used in Woo & Sriram's original paper, asymptotically, the threshold 'AIC' does not go to 0 slower than the difference in squared Hellinger distances once the correct order p is reached and the estimator is therefore not consistent.",
            call. = FALSE)

  if(missing(bandwidth)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(bandwidth) ||
          (bandwidth != "adaptive" && (!is.numeric(bandwidth) || length(bandwidth) != 1)))
    stop("'bandwidth' has to be either \"adaptive\" or a single numeric!")

  if(missing(sample.n)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(sample.n) || !is.numeric(sample.n) || !(as.integer(sample.n) == sample.n))
    stop("'sample.n' has to be an integer!")

  if(missing(sample.plot)){} # allowed be missing
  # but if it is not missing needs to fulfil requirements
  else if(is.null(sample.plot) || !(is.logical(sample.plot)))
    stop("'sample.plot' has to be logical!")

}


## Purpose: calculate the scaled vector of determinants of the Hankel matrix of the
##          moments of the mixing distribution

.deterHankel.scaled <- function(dat, Hankel.method, Hankel.function, j.max = 10, B = 1000, ...){

  fun <- .moments.map(Hankel.method = Hankel.method)

  if(Hankel.method == "scale") # get warning if squared moments have to be calculated
    D_hat <- fun(dat = dat, Hankel.function = Hankel.function, j.max = j.max, message = TRUE)
  else D_hat <- fun(dat = dat, Hankel.function = Hankel.function, j.max = j.max)

  # bootstrapped vector of determinants
  bt <- boot(dat, statistic = fun, R = B, Hankel.function = Hankel.function,
             j.max = j.max, ...)
  D_boot <- bt$t

  cov_boot <- cov(D_boot)
  cov_boot_sqrt <- sqrtm(cov_boot)
  cov_boot_sqrt_inv <- try(solve(cov_boot_sqrt), silent = TRUE)
  if(inherits(cov_boot_sqrt_inv, "try-error")) {
    stop("The matrix squareroot of the estimated covariance matrix could not be inverted: singular system.",
         call. = FALSE)
  }
  Y_boot <- cov_boot_sqrt_inv %*% D_hat

  det.vec <- abs(Y_boot)
  return(as.vector(det.vec))
}


## Purpose: function returning the vector of estimated determinants for orders up to
##          j.max; returns an object of class "hankDet"

#' @title Estimation of Mixture Complexity Based on Hankel Matrix
#'
#' @description Estimation of mixture complexity based on estimating the determinant of the Hankel matrix of the moments of the mixing distribution. The estimated determinants can be scaled and/or penalized.
#'
#' @details Define the \eqn{complexity} of a finite mixture \eqn{F} as the smallest integer \eqn{p}, such that its pdf/pmf \eqn{f} can be written as
#' \deqn{f(x) = w_1*g(x;\theta _1) + \dots + w_p*g(x;\theta _p).}
#' \code{nonparamHankel} estimates \eqn{p} by iteratively increasing the assumed complexity \eqn{j} and calculating the determinant of the \eqn{(j+1) x (j+1)} Hankel matrix made up of the first \eqn{2j} raw moments of the mixing distribution. As shown by Dacunha-Castelle & Gassiat (1997), once the correct complexity is reached (i.e. for all \eqn{j >= p}), this determinant is zero.
#' This suggests an estimation procedure for \eqn{p} based on initially finding a consistent estimator of the moments of the mixing distribution and then choosing the estimator \eqn{estim_p} as the value of \eqn{j} which yields a sufficiently small value of the determinant. Since the estimated determinant is close to 0 for all \eqn{j >= p}, this could lead to choosing \eqn{estim_p} rather larger than the true value. The function therefore returns all estimated determinant values corresponding to complexities up to \code{j.max}, so that the user can pick the lowest \eqn{j} generating a sufficiently small determinant. In addition, the function allows the inclusion of a penalty term as a function of the sample size \code{n} and the currently assumed complexity \code{j} which will be added to the determinant value (by supplying \code{pen.function}), and/or scaling of the determinants (by setting \code{scaled  = TRUE}). For scaling, a nonparametric bootstrap is used to calculate the covariance of the estimated determinants, with \code{B} being the size of the bootstrap sample. The inverse of the square root of this covariance matrix (i.e. the matrix \eqn{S^{(-1)}} such that \eqn{A = SS} (see \code{\link[expm]{sqrtm}}), where A is the covariance matrix) is then multiplied with the estimated determinant vector to get the scaled determinant vector. Note that in the case of the scaled version the penalty function chosen should be multiplied by \eqn{\sqrt{n}} before it is entered as \code{pen.function}: let \eqn{S*} denote a \eqn{j_m x j_m} covariance matrix of the determinants calculated for the \eqn{b}th bootstrap sample (\eqn{b=1,...,B} and j=1,...,j_m). Then \eqn{S*} goes to \eqn{S/n} as \eqn{B,n} go to infinity. Write \deqn{S*^{-1/2}=\sqrt{n}*\hat{S}^{-1/2}.} Define the rescaled vector \deqn{(y_1,...,y_{j_m})^T = \sqrt{n}*\hat{S}^{-1/2}(\hat{d}_1,...,\hat{d}_{j_m})^T.} Then the creterion to be minimized becomes  \deqn{|y_j|  + pen.function*\sqrt{n}.} See further sections for examples.
#' For a thorough discussion of the methods that can be used for the estimation of the moments see the details section of \code{\link{datMix}}.
#'
#' @aliases nonparamHankel print.hankDet plot.hankDet
#' @usage
#' nonparamHankel(obj, j.max = 10, pen.function = NULL, scaled = FALSE, B = 1000, ...)
#' @param obj object of class \code{\link{datMix}}.
#' @param j.max integer specifying the maximal number of components to be considered.
#' @param pen.function function with arguments \code{j} and \code{n} specifying the penalty added to the determinant value in the objective function, given sample size \eqn{n} and the assumed complexity at current iteration \eqn{j}. If left empty, no penalty will be added. If non-empty and \code{scaled} is \code{TRUE}, the penalty function will be added after the determinants are scaled.
#' @param scaled logical flag specifying whether the vector of estimated determinants should be scaled.
#' @param B integer specifying the number of bootstrap replicates used for scaling of the determinants. Ignored if \code{scaled} is \code{FALSE}.
#' @param x object of class \code{hankDet}.
#' @param type character denoting type of plot, see, e.g. \code{\link[graphics]{lines}}. Defaults to \code{"b"}.
#' @param xlab,ylab labels for the x and y axis with defaults (the default for \code{ylab} is created within the function, if no value is supplied).
#' @param mar numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot, see \code{\link[graphics]{par}}.
#' @param ylim range of y values to use.
#' @param \dots
#'    \describe{
#'      \item{in \code{nonparamHankel()}:}{further arguments passed to the \code{\link[boot]{boot}} function if \code{scaled} is \code{TRUE}.}
#'      \item{in \code{plot.hankDet()}:}{further arguments passed to \code{\link[base]{plot}}.}
#'      \item{in \code{print.hankDet()}:}{further arguments passed to \code{\link[base]{print}}.}
#'   }
#' @return Vector of estimated determinants (optionally scaled and/or penalized) as an object of class \code{hankDet} with the following attributes:
#'       \item{scaled}{logical flag indicating whether the determinants are scaled.}
#'       \item{pen}{logical flag indicating whether a penalty was added to the determinants.}
#'       \item{dist}{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers.}
#' @references D. Dacunha-Castelle and E. Gassiat, "The estimation of the order of a mixture model", Bernoulli, Volume 3, Number 3, 279-299, 1997.
#' @seealso \code{\link{paramHankel}} for a similar approach which additionally estimates the component weights and parameters, \code{\link{datMix}} for construction of a \code{\link{datMix}} object.
#' @keywords cluster
#' @examples
#' \donttest{
#' ## create 'Mix' object
#' geomMix <- Mix("geom", discrete = TRUE, w = c(0.1, 0.6, 0.3), prob = c(0.8, 0.2, 0.4))
#'
#' ## create random data based on 'Mix' object (gives back 'rMix' object)
#' set.seed(1)
#' geomRMix <- rMix(1000, obj = geomMix)
#'
#' ## create 'datMix' object for estimation
#'
#' # explicit function giving the estimate for the j^th moment of the
#' # mixing distribution, needed for Hankel.method "explicit"
#'
#' explicit.fct.geom <- function(dat, j){
#'  1 - ecdf(dat)(j - 1)
#' }
#'
#' ## generating 'datMix' object
#' geom.dM <- RtoDat(geomRMix, Hankel.method = "explicit",
#'                  Hankel.function = explicit.fct.geom)
#'
#' ## function for penalization w/o scaling
#' pen <- function(j, n){
#'   (j*log(n))/(sqrt(n))
#' }
#'
#' ## estimate determinants w/o scaling
#' set.seed(1)
#' geomdets_pen <- nonparamHankel(geom.dM, pen.function = pen, j.max = 5,
#'                                scaled = FALSE)
#' plot(geomdets_pen, main = "Three component geometric mixture")
#'
#'
#' ## function for penalization with scaling
#' pen <- function(j, n){
#'   j*log(n)
#' }
#'
#' ## estimate determinants using the same penalty with scaling
#' geomdets_pen <- nonparamHankel(geom.dM, pen.function = pen, j.max = 5,
#'                                scaled = TRUE)
#' plot(geomdets_pen, main = "Three component geometric mixture")
#' }
#' @rdname nonparamHankel
#' @export nonparamHankel
nonparamHankel <- function(obj, j.max = 10, pen.function = NULL, scaled = FALSE, B = 1000,
                           ...){

  .input.checks.functions(obj, j.max = j.max, pen.function = pen.function, scaled = scaled,
                          B = B, Hankel = TRUE, param = FALSE)

  Hankel.method <- attr(obj, "Hankel.method")
  Hankel.function <- attr(obj, "Hankel.function")

  dat <- as.numeric(obj)

  if(!is.null(pen.function)){
    pen <- TRUE
    n <- length(dat)
  }
  else pen <- FALSE

  if (scaled == TRUE){
    res <- .deterHankel.scaled(dat, Hankel.function = Hankel.function,
                               Hankel.method = Hankel.method, j.max = j.max, B = B, ...)

  } else {

    fun <- .moments.map(Hankel.method = Hankel.method)
    res <- abs(fun(dat = dat, j.max = j.max, Hankel.function = Hankel.function))
  }

  if (pen == TRUE) res <- res + pen.function(j = 1:j.max, n = n)

  class(res) <- "hankDet"
  attr(res, "scaled") <- scaled
  attr(res, "pen") <- pen
  attr(res, "dist") <- attr(obj, "dist")
  return(res)
}


## Purpose: print method for "hankDet" objects
#' @rdname nonparamHankel
#' @method print hankDet
#' @export
print.hankDet <- function(x, ...){

  obj <- x
  scaled <- attr(obj, "scaled")
  pen <- attr(obj, "pen")
  dist <- attr(obj, "dist")

  if(scaled == TRUE & pen == TRUE){
    header <- paste("\nEstimation of the scaled and penalized determinants for a '", dist,
                    "' mixture model:\n", "\n", sep = "")
  } else if(scaled == TRUE & pen == FALSE){
    header <- paste("\nEstimation of the scaled determinants for a '", dist,
                    "' mixture model:\n", "\n", sep = "")
  } else if(scaled == FALSE & pen == TRUE){
    header <- paste("\nEstimation of the penalized determinants for a '", dist,
                    "' mixture model:\n", "\n", sep = "")
  } else {
    header <- paste("\nEstimation of the determinants for a '", dist, "' mixture model:\n",
                    "\n", sep = "")
  }

  message(header)
  cmat <- matrix(c(1:length(obj), obj), nrow = length(obj), ncol = 2)
  colnames(cmat) <- c("Number of components", "Determinant")
  rownames(cmat) <- rep("", length(obj))
  print(cmat, ...)

}



## Purpose: plot method for "hankDet" objects
#' @rdname nonparamHankel
#' @method plot hankDet
#' @export
plot.hankDet <- function(x, type = "b", xlab = "j",  ylab = NULL, mar = NULL,
                         ylim = c(min(0, min(obj)), max(obj)), ...){

  obj <- x
  scaled <- attr(obj, "scaled")
  pen <- attr(obj, "pen")

  if(is.null(ylab)){
    if(scaled == TRUE & pen == TRUE){
      ylab <- bquote(atop(NA, atop(textstyle("Determinant value"),
                                   textstyle("(scaled and penalized)"))))
      if(is.null(mar)) mar <-  c(5,6,4,2) + 0.1
    }
    else if (scaled == TRUE){
      ylab <- bquote(atop(NA, atop(textstyle("Determinant value"),
                                   textstyle("(scaled)"))))
      if(is.null(mar)) mar <-  c(5,6,4,2) + 0.1
    }
    else if (pen == TRUE){
      ylab <- bquote(atop(NA, atop(textstyle("Determinant value"),
                                   textstyle("(penalized)"))))
      if(is.null(mar)) mar <-  c(5,6,4,2) + 0.1
    } else {
      ylab <- "Determinant value"
      if(is.null(mar)) mar <- c(5,5,4,2) + 0.1 # ylab needs less space on the left
    }
  }

  opar <- par(mar = mar)
  on.exit(par(opar))
  plot(x = 1:length(obj), y = obj, type = type, ylab = ylab, xlab = xlab, ylim = ylim, ...)

}
