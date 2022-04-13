## Purpose: variable that are defined globally within this package

globalVariables(c("dist", "theta.bound.list", "formals.dist", "ndistparams", "dist_call",
                  "bounds", "lower", "upper",
                  "dat", "N", "n.max", "discrete", "continuous",
                  "Hankel.method", "Hankel.function", "MLE.function"))


## Purpose: Gives colors for plots

.get.colors <- function(alpha){

  cols <- list(c(0, 0, 1), # blue
               c(0, 1, 0), # green
               c(0, 1, 1), # turqoise
               c(1, 0.6, 0), # orange
               c(0.57, 0.1, 0.33), # lilac
               c(0.35, 0, 0.22), # dark lilac
               c(0, 0.3, 0.3), # dark turqoise
               c(0, 0.5, 0.5)) # medium turqoise
  cols <- sapply(cols, function(x) rgb(x[1], x[2], x[3], alpha = alpha))

}


## Purpose: Constructor for 'Mix' mixture objects

#' @title Mixtures of Univariate Distributions
#'
#' @description Function constructing objects of class \code{Mix} that represent finite mixtures of any univariate distribution. Additionally methods for printing and plotting are provided.
#'
#' @aliases Mix print.Mix is.Mix
#' @usage
#' Mix(dist, discrete, w = NULL, theta.list = NULL, name = NULL, \dots)
#'
#' is.Mix(x)
#' @param dist character string providing the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers. The function sources functions for the density/mass estimation and random variate generation from distributions in \code{\link[stats]{distributions}}, so the abbreviations should be specified accordingly. Thus to create a gaussian mixture, set \code{dist = "norm"}, for a poisson mixture, set \code{dist = "pois"}. The \code{Mix} function will find the functions \code{dnorm}, \code{rnorm} and \code{dpois}, \code{rpois} respectively.
#' @param discrete logical flag, should be set to TRUE if the mixture distribution is discrete and to FALSE if continuous.
#' @param w numeric vector of length \eqn{p}, specifying the mixture weights \eqn{w[i]} of the components, \eqn{i = 1,\dots,p}. If the weights do not add up to 1, they will be scaled accordingly. Equal weights for all components are used by default.
#' @param theta.list named list specifying the component parameters. The names of the list elements have to match the names of the formal arguments of the functions \code{ddist} and \code{rdist} exactly. For a gaussian mixture, the list elements would have to be named \code{mean} and \code{sd}, as these are the formal arguments used by \code{rnorm} and \code{dnorm} functions from \code{\link[stats]{distributions}}. Alternatively, the component parameters can be supplied directly as named vectors of length \eqn{p} via \dots
#' @param name optional name tag of the result (used for printing and plotting).
#' @param x
#'    \describe{
#'      \item{in \code{is.Mix()}:}{returns TRUE if the argument is a \code{datMix} object and FALSE otherwise.}
#'      \item{in \code{print.Mix()}:}{object of class \code{Mix}.}
#' }
#' @param \dots
#'    \describe{
#'      \item{in \code{Mix()}:}{alternative way of supplying the component parameters (instead of using \code{theta.list}).}
#'      \item{in \code{print.Mix()}:}{further arguments passed to the print method.}
#' }
#' @return An object of class \code{Mix} (implemented as a matrix) with the following attributes:
#'     \item{dim}{dimensions of the matrix.}
#'     \item{dimnames}{a \code{\link[base]{dimnames}} attribute for the matrix.}
#'     \item{name}{optional name tag for the result passed on to printing and plotting methods.}
#'     \item{dist}{character string giving the abbreviated name of the component distribution, such that the function \code{ddist} evaluates its density/mass and \code{rdist} generates random variates.}
#'     \item{discrete}{logical flag indicating whether the mixture distribution is discrete.}
#'     \item{theta.list}{named list specifying component parameters.}
#' @seealso \code{\link{dMix}} for the density, \code{\link{rMix}} for random numbers (and construction of an \code{rMix} object) and \code{\link{plot.Mix}} for the plot method.
#' @keywords cluster
#' @examples
#'
#' # define 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#' poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
#'
#' # plot 'Mix' object
#' plot(normLocMix)
#' plot(poisMix)
#'
#' @export Mix
Mix <- function(dist, discrete, w = NULL, theta.list = NULL, name = NULL, ...){

  if(!is.character(dist)) stop("'dist' needs to be a character string specifying the
                               distribution of the mixture components!")

  if(!is.logical(discrete)) stop("'discrete' needs to be a boolean specifying whether the
                                   distribution is discrete!")

  if(!is.null(theta.list))
    if(!is.list(theta.list))
      stop("Input to theta.list has to be of class 'list'!")
    else theta.list <- theta.list
  # component parameters can also be entered via ...
  else theta.list <- list(...)

  if(is.null(unlist(theta.list)))
    stop("The component parameters have to be entered as 'theta.list' or as input to ...!")

  formals.dist <- names(theta.list)
  ndistparams <- length(formals.dist)

  # check if the vectors of component parameters (e.g. mean and sd for normal distribution)
  # are of equal length
  equal.length <- function(x){
    len <- sapply(x, length)
    diff(range(len)) < .Machine$double.eps
  }

  if(!equal.length(theta.list))
    stop("The elements of theta.list (or the inputs to ...) must be of equal length!")

  if(!is.numeric(unlist(theta.list)))
    stop("The elements of theta.list (or the inputs to ...) must all be numeric!")

  # true mixture complexity
  p <- length(theta.list[[1]])

  if(is.null(w)){ # construct default value for the weights
    w <- rep.int(1/p, p)
  } else {

    if(length(w) != p || !is.numeric(w) || any(w < 0))
      stop("'w' must be a numeric >= 0 with same length as the elements of theta.list
           (or the inputs to ...)")
    s <- sum(w)
    if(abs(s - 1) > .Machine$double.eps) w <- w/s

  }

  if(is.null(name)) { # construct default name

    sformat <- function(v) sapply(v, format, digits = 1)
    pPar <- function(pp) {
      pp <- if(p >= 10) c(sformat(pp[1:9]), "....") else sformat(pp)
      paste(pp, collapse = "'")
    }
    name <- paste0(dist, "Mix", format(p))
    for(i in 1:ndistparams){
      name <- paste0(name, "_", pPar(theta.list[[i]]))
    }

  }
  if(!is.character(name)) name <- as.character(name)

  # construct matrix to be printed with weights as first column and other parameters
  # as further columns
  dat <- matrix(NA, nrow = p, ncol = ndistparams + 1)
  colnames(dat) <- as.character(1:(ndistparams + 1))
  for(i in 1:ndistparams){
    dat[, i+1] <- theta.list[[i]]
    colnames(dat)[i+1] <- formals.dist[i]
  }
  dat[, 1] <- w
  colnames(dat)[1] <- "w"

  # check whether there exists a function generating random numbers and returning the
  # probability denisty/mass for the distribution specified as 'dist' and if their
  # arguments are compatible with the component parameter names specified in 'theta.list'
  # (or via ...)
  dist_call_r <- try(get(paste("r", dist, sep = "")), silent = TRUE)
  dist_call_d <- try(get(paste("d", dist, sep = "")), silent = TRUE)
  if(inherits(dist_call_r, "try-error") || inherits(dist_call_d, "try-error"))
    stop("combining the string \"dist\" with \"d\" or \"r\" has to yield an existing function name!")

  if(!all(formals.dist %in% names(formals(dist_call_r))))
    stop(paste("The names of theta.list do not match the names of the formal arguments
               of the function r", dist, sep = ""))
  if(!all(formals.dist %in% names(formals(dist_call_d))))
    stop(paste("The names of theta.list do not match the names of the formal arguments
               of the function d", dist, sep = ""))

  structure(name = name, dist = dist, discrete = discrete, theta.list = theta.list, class = "Mix",
            .Data = dat)
}



## Purpose: check whether 'obj' is a "Mix" object (auxiliary function)

is.Mix <- function(x){

  obj <- x
  dist_call_d <- try(get(paste("d", attr(obj, "dist"), sep = "")), silent = TRUE)
  dist_call_r <- try(get(paste("r", attr(obj, "dist"), sep = "")), silent = TRUE)
  if(inherits(dist_call_d, "try-error") || inherits(dist_call_r, "try-error")) return(FALSE)
  theta.list <- attr(obj, "theta.list")
  formals.dist <- names(theta.list)

  inherits(obj, "Mix") && is.matrix(obj) &&
    (!is.null(w <- obj[, "w"])) && dim(obj)[1] == length(w) &&
    dim(obj)[2] == (length(theta.list) + 1) &&
    is.numeric(w) && all(w >= 0) && abs(sum(w) - 1) < 1000*.Machine$double.eps &&
    length(w) == length(theta.list[[1]]) &&
    length(unique(sapply(theta.list, function(x) length(x)))) == 1 &&
    all(formals.dist %in% names(formals(dist_call_d))) &&
    all(formals.dist %in% names(formals(dist_call_r))) &&
    is.logical(attr(obj, "discrete"))

}


## Purpose: plot method for "Mix" objects (mixtures)
#' @title \code{plot} Method for \code{\link{Mix}} Objects
#'
#' @description \code{plot} method for \code{\link{Mix}} objects visualizing the mixture density, with an option of showing the component densities.
#'
#' @param x object of class \code{Mix}.
#' @param ylim range of y values to use, if not specified (or containing \code{NA}), the function tries to construct reasonable default values.
#' @param xlim range of x values to use, particularly important if \code{xout} is not specified. If not specified, the function tries to construct reasonable default values.
#' @param xout numeric or \code{NULL} giving the abscissae at which to draw the density.
#' @param n number of points to generate if \code{xout} is unspecified (for continuous distributions).
#' @param type character denoting the type of plot, see e.g. \code{\link[graphics]{lines}}. Defaults to \code{"l"} if the mixture distribution is continuous and to \code{"h"} if discrete.
#' @param xlab,ylab labels for the x and y axis with defaults.
#' @param main main title of plot, defaulting to the \code{\link{Mix}} object name.
#' @param lwd line width for plotting, a positive number.
#' @param log logical flag, if \code{TRUE}, probabilities/densities \eqn{f} are plotted as \eqn{log(f)}. Only works if \code{components} is set to \code{FALSE}.
#' @param h0 logical flag indicating whether the line \eqn{y = 0} should be drawn.
#' @param components logical flag indicating whether the individual mixture components should be plotted, set to \code{TRUE} by default.
#' @param parH0 graphical parameters for drawing the line \eqn{y = 0} if \code{h0} is set to \code{TRUE}.
#' @param parComp graphical parameters for drawing the individual components if \code{components} is set to \code{TRUE}.
#' @param  \dots further arguments passed to the function for plotting the mixture density.
#' @seealso \code{\link{Mix}} for the construction of \code{Mix} objects, \code{\link{dMix}} for the density/mass of a mixture.
#' @keywords cluster
#' @examples
#'
#' # define 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#' poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
#'
#' # plot 'Mix' object
#' plot(normLocMix)
#' plot(poisMix)
#'
#' @method plot Mix
#' @export
plot.Mix <- function(x, ylim, xlim = NULL, xout = NULL, n = 511, type = NULL,
                     xlab = "x", ylab = "f(x)", main = attr(obj, "name"), lwd = 1.4,
                     log = FALSE, components = TRUE, h0 = FALSE,
                     parComp = list(col = NULL, lty = 3, lwd = 1),
                     parH0 = list(col = NULL, lty = 3, lwd = 1), ...){

  obj <- x
  if(!is.numeric(n)) stop("'n' has to be an integer!")
  if(!is.logical(components)) stop("'components' has to be logical!")
  if(!is.list(parH0)) stop("'parH0' has to be a list!")
  if(!is.list(parComp)) stop("'parComp' has to be a list!")

  if(components == TRUE && is.null(parComp$col)){ # default colors
    parComp$col <- .get.colors(alpha = 0.9)[8]
  }

  if(h0 == TRUE && is.null(parH0$col)){ # default colors
    parH0$col <- .get.colors(alpha = 0.9)[6]
  }

  if(is.null(type)){ # default type
    if (attr(obj, "discrete") == FALSE) type = "l"
    else type = "h"
  }

  if(is.null(xlim) && is.null(xout)){  # construct "reasonable" abscissa values

    set.seed(1)
    rand <- rMix(1000, obj = obj)
    d.o <- dMix(obj, x = rand)
    if (attr(obj, "discrete") == FALSE)
      rand_trunc <- rand[d.o >= max(d.o)/200] # discard values with too low density
    else rand_trunc <- rand
    xlim <- c(min(rand_trunc), max(rand_trunc))

  } # here we definitely have xlim values

  if(!is.null(xlim) && is.null(xout)){ # construct xout from xlim

    if (attr(obj, "discrete") == FALSE)
      xout <- seq(xlim[1], xlim[2], length = n)
    else xout <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 1)

  } # here we definitely have xout values


  if(log == TRUE & components == TRUE){
    warning("Only either 'log' or 'components' can take the value true. Setting components = FALSE")
    components <- FALSE
  }

  d.o <- dMix(obj, x = xout, log = log)

  if(missing(ylim) || anyNA(ylim)){ # construct "reasonable" ordinate values
    if(log == FALSE)
      ylim <- c(0, max(d.o))
    else ylim <- c(min(d.o), max(d.o))
  }

  # plot mixture density
  plot(y = d.o, x = xout, type = type, xlim = xlim, ylim = ylim,
       main = main, xlab = xlab, ylab = ylab, lwd = lwd, ...)

  if(components == TRUE) { # plot individual components

    w <- obj[ ,"w"]
    p <- length(w)
    dist_call <- get(paste("d", attr(obj, "dist"), sep = ""))
    theta.list <- attr(obj, "theta.list")
    theta.list <- lapply(theta.list, function(y) matrix(y, nrow = length(xout), ncol = p,
                                                        byrow = TRUE))
    theta.list$x <- xout
    y <- matrix(w, nrow = length(xout), ncol = p,
                byrow = TRUE) * do.call(dist_call, theta.list)
    mapply(function(i) do.call(lines, c(list(x = xout, y = y[, i]), parComp)), i = 1:p)

  }

  if(h0 == TRUE) do.call(abline, c(list(h = 0), parH0))
}


## Purpose: print method for "Mix" objects (mixtures)
#' @rdname Mix
#' @method print Mix
#' @export
print.Mix <- function(x, ...){

  obj <- x
  if(!is.Mix(obj)) stop("obj is not a 'Mix' object!")

  message(paste("'", attr(obj, "dist"), sep = ""), "Mixture' object",
      paste("\t ``", attr(obj, "name"), "''", sep=''), "\n")

  att <- attributes(obj);
  att <- att[names(att) != "dist" & names(att) != "theta.list" & names(att) != "discrete" & names(att) != "name"]
  attributes(obj) <- if(length(att) > 0) att

  class(obj) <- character(0)
  print(obj, ...)
  invisible(x)

}


## Purpose: density evaluation for "Mix" objects (mixtures)
#' @title Mixture density
#'
#' @description Evaluation of the (log) density function of a mixture specified as a \code{Mix} object.
#'
#' @usage
#' dMix(x, obj, log = FALSE)
#' @param x vector of quantiles.
#' @param obj object of class \code{\link{Mix}}.
#' @param log logical flag, if \code{TRUE}, probabilities/densities \eqn{f} are returned as \eqn{log(f)}.
#' @return \code{dMix(x)} returns a numeric vector of probability values \eqn{f(x)} and logarithm thereof if \code{log} is \code{TRUE}.
#' @seealso \code{\link{Mix}} for the construction of \code{Mix} objects, \code{\link{rMix}} for random number generation (and construction of \code{rMix} objects) and \code{\link{plot.Mix}} for plotting the densities computed using \code{\link{dMix}}.
#' @keywords cluster
#' @examples
#'
#' # define 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#'
#' # evaluate density at points x
#' x <- seq(7, 20, length = 501)
#' dens <- dMix(x, normLocMix)
#' plot(x, dens, type = "l")
#'
#' # compare to plot.Mix
#' plot(normLocMix)
#'
#' @export dMix
dMix <- function(x, obj, log = FALSE){

  if(!is.numeric(x)) stop("'x' has to be numeric!")
  if(!is.Mix(obj)) stop("'obj' has to be a 'Mix' object!")
  if(!is.logical(log)) stop("'log' has to be logical!")

  w <- obj[,"w"]
  p <- length(w) # number of components
  theta.list <- attr(obj, "theta.list")
  theta.list <- lapply(theta.list, function(y) matrix(y, nrow = length(x), ncol = p,
                                                      byrow = TRUE))
  theta.list$x <- x

  dist_call <- get(paste("d", attr(obj, "dist"), sep = ""))
  y <- rowSums(matrix(w, nrow = length(x),
                      ncol = p, byrow = TRUE) * do.call(dist_call, theta.list))
  if(log) log(y) else y

}



## Purpose: Generate random numbers according to "Mix" object;
##          simultaneously creates an "rMix" object

#' @title Random Variate Generation from a Mixture Distribution
#'
#' @description Function for generating a random sample of size \code{n}, distributed according to a mixture specified as \code{\link{Mix}} object. Returns an object of class \code{\link{rMix}}.
#'
#' @details For a mixture of \eqn{p} components, generates the number of observations in each component as multinomial, and then use an implemented random variate generation function for each component. The integer (multinomial) numbers are generated via \code{\link[base]{sample}}.
#'
#' @aliases rMix is.rMix print.rMix
#' @usage
#' rMix(n, obj)
#'
#' is.rMix(x)
#' @param n integer specifying the number of observations.
#' @param obj object of class \code{\link{Mix}}.
#' @param x
#'    \describe{
#'      \item{in \code{is.rMix()}:}{R object.}
#'      \item{in \code{print.rMix()}:}{object of class \code{\link{rMix}}.}
#' }
#' @param \dots further arguments passed to the print method.
#' @return An object of class \code{\link{rMix}} with the following attributes (for further explanations see \code{\link{Mix}}):
#'     \item{name}{name of the \code{Mix} object that was given as input.}
#'     \item{dist}{character string stating the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers.}
#'     \item{discrete}{logical flag indicating whether the underlying mixture distribution is discrete.}
#'     \item{theta.list}{named list specifying the parameter values of the \eqn{p} components.}
#'     \item{w}{numeric vector of length \eqn{p} specifying the mixture weights \eqn{w[i]} of the components, \eqn{i = 1,\dots,p}.}
#'     \item{indices}{numeric vector of length \code{n} containing integers between \eqn{1} and \eqn{p} specifying which mixture component each observation belongs to.}
#' @seealso \code{\link{dMix}} for the density, \code{\link{Mix}} for the construction of \code{Mix} objects and \code{\link{plot.rMix}} for the plot method.
#' @keywords cluster
#' @examples
#'
#' # define 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#'
#' # generate n random samples
#' set.seed(1)
#' x <- rMix(1000, normLocMix)
#' hist(x)
#'
#' @export rMix
rMix <- function(n, obj){

  if(!is.Mix(obj)) stop("'obj' has to be a 'Mix' object!")
  if(!is.numeric(n) || !(as.integer(n) == n) || length(n) != 1)
    stop("'n' has to be a single integer!")

  w <- obj[, "w"]
  p <- length(w)
  theta.list <- attr(obj, "theta.list")

  ind <- sample(1:p, prob = w, size = n, replace = TRUE)
  theta.list.expanded <- lapply(theta.list, function(x) x[ind])

  dist_call <- get(paste("r", attr(obj, "dist"), sep = ""))
  dat <- do.call(dist_call, args = c(n = n, theta.list.expanded))

  att <- attributes(obj)
  attributes(dat) <- att[names(att) != "class" & names(att) != "dimnames" & names(att) != "dim"]
  attr(dat, "w") <- w
  attr(dat, "indices") <- ind
  class(dat) <- "rMix"

  return(dat)
}


## Purpose: check whether 'obj' is an "rMix" object (auxiliary function)
is.rMix <- function(x){

  obj <- x
  dist_call_d <- try(get(paste("d", attr(obj, "dist"), sep = "")), silent = TRUE)
  dist_call_r <- try(get(paste("r", attr(obj, "dist"), sep = "")), silent = TRUE)
  if(inherits(dist_call_d, "try-error") || inherits(dist_call_r, "try-error")) return(FALSE)

  theta.list <- attr(obj, "theta.list")
  formals.dist <- names(theta.list)

  inherits(obj, "rMix") && is.numeric(obj) &&
    (!is.null(w <- attr(obj, "w"))) &&
    is.numeric(w) && all(w >= 0) && abs(sum(w)-1) < 1000*.Machine$double.eps &&
    length(w) == length(theta.list[[1]]) &&
    length(unique(sapply(theta.list, function(x) length(x)))) == 1 &&
    all(formals.dist %in% names(formals(dist_call_d))) &&
    all(formals.dist %in% names(formals(dist_call_r))) &&
    is.logical(attr(obj, "discrete")) && all(unique(attr(obj, "indices")) %in% 1:length(w))
}



## Purpose: plot method for "rMix" objects (random numbers generiated via a mixture)

#' @title \code{plot} Method for \code{\link{rMix}} Objects
#'
#' @description plot method for \code{rMix} objects, plotting the histogram of the random sample, with the option of additionally plotting the components (stacked or plotted over one another).
#'
#' @param x object of class \code{rMix}.
#' @param xlab label for the x axis with default.
#' @param ylim range of y values to use, if not specified (or containing \code{NA}), default values are used.
#' @param main main title of the plot, defaulting to the \code{\link{rMix}} object name.
#' @param breaks see \code{\link[graphics]{hist}}. If left unspecified the function tries to construct reasonable default values.
#' @param col colour to be used to fill the bars of the histogram evaluated on the whole data.
#' @param components logical flag indicating whether the plot should show to which component the observations belong (either by plotting individual histograms or by overlaying a stacked barplot), defaulting to \code{TRUE}. Ignored if \code{plot} is \code{FALSE}.
#' @param stacked logical flag indicating whether the component plots should be stacked or plotted one over another, defaulting to \code{FALSE}. Ignored if \code{components} is \code{FALSE} or ignored itself.
#' @param component.colors colors for the component plots. If unspecified, default colors are used.
#' @param freq logical flag, if \code{TRUE}, the histogram graphic is a representation of frequencies, if \code{FALSE}, probability densities. See \code{\link[graphics]{hist}}.
#' @param plot logical flag, if \code{TRUE} (default), a histogram is plotted, otherwise a list of breaks and counts is returned. See \code{\link[graphics]{hist}}.
#' @param \dots further arguments passed to the histogram function evaluated on the whole data as well as the component data (if \code{components} is \code{TRUE} and \code{stacked} is \code{FALSE}).
#' @seealso \code{\link{rMix}} for the construction of \code{rMix} objects.
#' @keywords cluster
#' @examples
#'
#' # define 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#' # generate n random samples
#' set.seed(1)
#' x <- rMix(1000, normLocMix)
#' plot(x)
#'
#' @method plot rMix
#' @export
plot.rMix <- function(x, xlab = attr(obj, "name"), ylim = NULL,
                      main = paste("Histogram of", attr(obj, "name")), breaks = NULL,
                      col = "grey", components = TRUE, stacked = FALSE, component.colors = NULL,
                      freq = TRUE, plot = TRUE, ...){

  obj <- x
  if(!is.logical(components)) stop("'components' has to be logical!")
  if(!is.logical(stacked)) stop("'stacked' has to be logical!")

  # from which component did an observation arise?
  ind <- attr(obj, "indices")

  if(is.null(breaks)){ # create default breaks

    hist <- hist(obj, plot = FALSE, right = FALSE)
    breaks <- hist$breaks
    diffbreaks <- unique(diff(breaks)*0.5)
    if(diff(range(diffbreaks)) < .Machine$double.eps^0.5) diffbreaks <- diffbreaks[1]
    nbreaks <- length(breaks)
    breaks <- rep(breaks, each = 2) + rep(c(0, diffbreaks), nbreaks)
    # making breaks twice as fine so that components can also be shown nicely
  }


  if(is.null(component.colors) && components == TRUE){ # default colors
    if(stacked == FALSE)
      component.colors <- .get.colors(alpha = 0.4)
    else component.colors <- .get.colors(alpha = 0.3)
  }

  if(freq == FALSE && components == TRUE){ # if we plot probabilities and individual components
    if(stacked == TRUE){
      warning("Stacked components cannot be shown when 'freq' is FALSE. Setting freq = TRUE.")
      freq <- TRUE
    } else if(is.null(ylim) || anyNA(ylim)){

      # ordinate values such that component probabilities surely fit in the plot
      plt.inv <- hist(obj, breaks = breaks, plot = FALSE, right = FALSE, ...)
      breaks <- plt.inv$breaks
      ylim <- c(0, 1/min(diff(breaks)))
    }
  }

  if(plot == FALSE){
    plt.inv <- hist(obj, breaks = breaks, plot = plot, ...)
    return(plt.inv)
  } else {
    hist(obj, col = col, breaks = breaks, main = main, xlab = xlab, freq = freq,
         ylim = ylim, plot = plot, right = FALSE, ...)
    # invisible plot to reuse breaks later for components

    plt.inv <- suppressWarnings(hist(obj, breaks = breaks, plot = FALSE, right = FALSE,
                                     ...))
  }

  if(components == TRUE){

    breaks <- plt.inv$breaks
    p <- length(attr(obj, "theta.list")[[1]]) # number of components

    if(stacked == FALSE){ # just plot the individual histograms
      for (i in 1:p){
        hist(obj[ind == i], col = component.colors[i], breaks = breaks,
             add = TRUE, freq = freq, plot = plot, right = FALSE, ...)
      }

    } else { # plot stacked rectangles

      nbreaks <- length(breaks)
      mat <- matrix(0, nrow = p + 1, ncol = nbreaks - 1)
      # cumulative version of mat
      matcum <- matrix(0, nrow = p + 1, ncol = nbreaks - 1)

      for(i in 2:(p + 1)){ # for every component i...

        # How many observations lie within the bins
        mat[i, 1:(ncol(mat)-1)] <- mapply(function(j) length(obj[obj < breaks[j + 1] &
                                                          obj >= breaks[j] &
                                                          ind == (i - 1)]), 1:(ncol(mat)-1))
        matcum[i, 1:(ncol(mat)-1)] <- mapply(function(j) matcum[i - 1, j] + mat[i, j], 1:(ncol(mat)-1))

        # ... how many observations lie within last bin
        # (seperate because both breaks value need to be included)
        mat[i, ncol(mat)] <- length(obj[obj <= breaks[ncol(mat) + 1] & obj >= breaks[ncol(mat)] & ind == (i - 1)])
        matcum[i, ncol(mat)] <- matcum[i - 1, ncol(mat)] + mat[i, ncol(mat)]

        rect(xleft = breaks[1:(nbreaks - 1)], ybottom = matcum[i - 1, ], xright = breaks[2:(nbreaks)],
             ytop = matcum[i, ], col = component.colors[i - 1])
      }

    }
  }
}


## Purpose: print method for "rMix" objects (random numbers generated via a mixture
#' @rdname rMix
#' @method print rMix
#' @export
print.rMix <- function(x, ...){
  print(as.vector(x), ...)
  return(x)
}


# Purpose: Constructor for 'datMix' objects, to be passed to functions estimating
#          mixture complexity, contains all "static" information about the data

#' @title Constructor for Objects for Which to Estimate the Mixture Complexity
#'
#' @description Function to generate objects of class \code{datMix} to be passed to other \code{mixComp} functions used for estimating mixture complexity.
#'
#' @details If the \code{datMix} object is supposed to be passed to a function that calculates the Hankel matrix of the moments of the mixing distribution (i.e. \code{\link{nonparamHankel}}, \code{\link{paramHankel}} or \code{\link{paramHankel.scaled}}), the arguments \code{Hankel.method} and \code{Hankel.function} have to be specified. The \code{Hankel.method}s that can be used to generate the estimate of the (raw) moments of the mixing distribution and the corresponding \code{Hankel.function}s are the following, where \eqn{j} specifies an estimate of the number of components:
#'  \describe{
#'     \item{\code{"explicit"}}{For this method, \code{Hankel.function} contains a function with arguments called \code{dat} and \code{j}, explicitly estimating the moments of the mixing distribution from the data and assumed mixture complexity at current iteration. Note that what Dacunha-Castelle & Gassiat (1997) called the "natural" estimator in their paper is equivalent to using \code{"explicit"} with \code{Hankel.function}
#' \deqn{f_j((1/n) * \sum_i(\psi_j(X_i))).}}
#'     \item{\code{"translation"}}{This method corresponds to Dacunha-Castelle & Gassiat's (1997) example 3.1. It is applicable if the family of component distributions \eqn{(G_\theta)} is given by
#' \deqn{dG_\theta(x) = dG(x-\theta),}
#' where \eqn{G} is a known probability distribution, such that its moments can be expressed explicitly. \code{Hankel.function} contains a function of \eqn{j} returning the \eqn{j}th (raw) moment of \eqn{G}.}
#'     \item{\code{"scale"}}{This method corresponds to Dacunha-Castelle & Gassiat's (1997) example 3.2. It is applicable if the family of component distributions \eqn{(G_\theta)} is given by
#' \deqn{dG_\theta(x) = dG(x / \theta),}
#' where \eqn{G} is a known probability distribution, such that its moments can be expressed explicitly. \code{Hankel.function} contains a function of \eqn{j} returning the \eqn{j}th (raw) moment of \eqn{G}.}
#'  }
#' If the \code{datMix} object is supposed to be passed to a function that estimates the component weights and parameters (i.e. all but \code{\link{nonparamHankel}}), the arguments \code{discrete} and \code{theta.bound.list} have to be specified, and \code{MLE.function} will be used in the estimation process if it is supplied (otherwise the MLE is found numerically).
#' @aliases datMix is.datMix print.datMix
#' @usage
#' datMix(dat, dist, discrete = NULL, theta.bound.list = NULL,
#'        MLE.function = NULL, Hankel.method = NULL, Hankel.function = NULL)
#'
#' is.datMix(x)
#' @param dat numeric vector containing observations from the mixture model.
#' @param dist character string providing the (abbreviated) name of the component distribution, such that the function \code{ddist} evaluates its density function and \code{rdist} generates random numbers. The function sources functions for the density/mass estimation and random variate generation from distributions in \code{\link[stats]{distributions}}, so the abbreviations should be specified accordingly. Thus to create a gaussian mixture, set \code{dist = "norm"}, for a poisson mixture, set \code{dist = "pois"}. The \code{MixComp} functions will find the functions \code{dnorm}, \code{rnorm} and \code{dpois}, \code{rpois} respectively.
#' @param discrete logical flag indicating whether the mixture distribution is discrete, required for methods that estimate component weights and parameters.
#' @param theta.bound.list named list specifying the upper and lower bounds for the component parameters. The names of the list elements have to match the names of the formal arguments of the functions \code{ddist} and \code{rdist} exactly as specified in the distributions in \code{\link[stats]{distributions}}. For a gaussian mixture, the list elements would have to be named \code{mean} and \code{sd}, as these are the formal arguments used by \code{rnorm} and \code{dnorm}. Has to be supplied if a method that estimates the component weights and parameters is to be used.
#' @param MLE.function function (or a list of functions) which takes the data as input and outputs the maximum likelihood estimator for the parameter(s) the component distribution \code{dist}. If the component distribution has more than one parameter, a list of functions has to be supplied and the order of the MLE functions has to match the order of the component parameters in \code{theta.bound.list} (e.g. for a normal mixture, if the first entry of \code{theta.bound.list} is the bounds of the mean, then then first entry of \code{MLE.function} has to be the MLE of the mean). If this argument is supplied and the \code{datMix} object is handed over to a complexity estimation procedure relying on optimizing over a likelihood function, the \code{MLE.function} attribute will be used for the single component case. In case the objective function is neither a likelihood nor corresponds to a mixture with more than 1 component, numerical optimization will be used based on \code{\link{Rsolnp}}'s function \code{\link[Rsolnp]{solnp}}, but \code{MLE.function} will be used to calculate the initial values passed to \code{solnp}. Specifying \code{MLE.function} is optional. If not supplied, for example because the MLE solution does not exist in a closed form, numerical optimization is used to find the relevant MLE.
#' @param Hankel.method character string in \code{c("explicit", "translation", "scale")},  specifying the method of estimating the moments of the mixing distribution used to calculate the relevant Hankel matrix. Has to be specified when using \code{\link{nonparamHankel}}, \code{\link{paramHankel}} or \code{\link{paramHankel.scaled}}. For further details see below.
#' @param Hankel.function function required for the moment estimation via \code{Hankel.method}. This normally depends on \code{Hankel.method} as well as \code{dist}. For further details see below.
#' @param x
#'   \describe{
#'      \item{in \code{is.datMix()}:}{returns TRUE if the argument is a \code{datMix} object and FALSE otherwise.}
#'      \item{in \code{print.datMix()}:}{object of class \code{datMix}.}
#'  }
#' @param \dots further arguments passed to the print method.
#' @return Object of class \code{datMix} with the following attributes (for further explanations see above):
#'      \item{dist}{character string giving the abbreviated name of the component distribution, such that the function \code{ddist} evaluates its density/mass and \code{rdist} generates random variates.}
#'      \item{discrete}{logical flag indicating whether the mixture distribution is discrete.}
#'      \item{theta.bound.list}{named list specifying the upper and lower bounds for the component parameters.}
#'      \item{MLE.function}{function which computes the MLE of the component distribution \code{dist}.}
#'      \item{Hankel.method}{character string taking on values \code{"explicit"}, \code{"translation"}, or \code{"scale"}, specifying the method of estimating the moments of the mixing distribution to compute the corresponding Hankel matrix.}
#'      \item{Hankel.function}{function required for the moment estimation via \code{Hankel.method}. See details for more information.}
#' @seealso \code{\link{RtoDat}} for conversion of \code{\link{rMix}} to \code{datMix} objects.
#' @keywords cluster
#' @examples
#'
#' ## observations from a (presumed) mixture model
#' obs <- faithful$waiting
#'
#' ## generate list of parameter bounds (assuming gaussian components)
#' norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
#'
#' ## generate MLE functions
#' # for "mean"
#' MLE.norm.mean <- function(dat) mean(dat)
#' # for "sd" (the sd function uses (n-1) as denominator)
#' MLE.norm.sd <- function(dat){
#'   sqrt((length(dat) - 1) / length(dat)) * sd(dat)
#' }
#' # combining the functions to a list
#' MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
#'                       "MLE.norm.sd" = MLE.norm.sd)
#'
#' ## function giving the j^th raw moment of the standard normal distribution,
#' ## needed for calculation of the Hankel matrix via the "translation" method
#' ## (assuming gaussian components with variance 1)
#'
#' mom.std.norm <- function(j){
#'   ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
#' }
#'
#' ## generate 'datMix' object
#' faithful.dM <- datMix(obs, dist = "norm", discrete = FALSE,
#'                       theta.bound.list = norm.bound.list, MLE.function = MLE.norm.list,
#'                       Hankel.method = "translation", Hankel.function = mom.std.norm)
#'
#' ## using 'datMix' object to estimate the mixture complexity
#' set.seed(1)
#' res <- paramHankel.scaled(faithful.dM)
#' plot(res)
#'
#' @export datMix
datMix <- function(dat, dist, discrete = NULL, theta.bound.list = NULL, MLE.function = NULL,
                   Hankel.method = NULL, Hankel.function = NULL){


  # return 'datMix' object with relevant attributes
  structure(class = "datMix", dist = dist,
            theta.bound.list = theta.bound.list, discrete = discrete,
            MLE.function = MLE.function, Hankel.method = Hankel.method,
            Hankel.function = Hankel.function, .Data = dat)
}


## Purpose: is 'obj' a "datMix" object?
#' @export is.datMix
is.datMix <- function(x){

  obj <- x
  dist_call <- try(get(paste("d", attr(obj, "dist"), sep = "")), silent = TRUE)
  if(inherits(dist_call, "try-error")) return(FALSE)

  inherits(obj, "datMix") && is.numeric(obj)
}


## Purpose: print method for "datMix" objects (don't print attributes)
#' @rdname datMix
#' @method print datMix
#' @export
print.datMix <- function(x, ...){
  print(as.vector(x), ...)
}


## Purpose: Convert "RMix" object to "datMix" object

#' @title Converting \code{rMix} to \code{datMix} Objects
#'
#' @description Function for converting objects of class \code{\link{rMix}} to objects of class \code{\link{datMix}}, so that they could be passed to functions estimating the mixture complexity.
#'
#' @usage
#' RtoDat(obj, theta.bound.list = NULL, MLE.function = NULL, Hankel.method = NULL,
#'        Hankel.function = NULL)
#' @param obj object of class \code{rMix}.
#' @param theta.bound.list named list specifying the upper and lower bounds for the component parameters. The names of the list elements have to match the names of the formal arguments of the functions \code{ddist} and \code{rdist} exactly as specified in the distributions in \code{\link[stats]{distributions}}. For a gaussian mixture, the list elements would have to be named \code{mean} and \code{sd}, as these are the formal arguments used by \code{rnorm} and \code{dnorm}. Has to be supplied if a method that estimates the component weights and parameters is to be used.
#' @param MLE.function function (or a list of functions) which takes the data as input and outputs the maximum likelihood estimator for the parameter(s) the component distribution \code{dist}. If the component distribution has more than one parameter, a list of functions has to be supplied and the order of the MLE functions has to match the order of the component parameters in \code{theta.bound.list} (e.g. for a normal mixture, if the first entry of \code{theta.bound.list} is the bounds of the mean, then then first entry of \code{MLE.function} has to be the MLE of the mean). If this argument is supplied and the \code{datMix} object is handed over to a complexity estimation procedure relying on optimizing over a likelihood function, the \code{MLE.function} attribute will be used for the single component case. In case the objective function is neither a likelihood nor corresponds to a mixture with more than 1 component, numerical optimization will be used based on \code{\link{Rsolnp}}'s function \code{\link[Rsolnp]{solnp}}, but \code{MLE.function} will be used to calculate the initial values passed to \code{solnp}. Specifying \code{MLE.function} is optional. If not supplied, for example because the MLE solution does not exist in a closed form, numerical optimization is used to find the relevant MLE.
#' @param Hankel.method character string in \code{c("explicit", "translation", "scale")}, specifying the method of estimating the moments of the mixing distribution used to calculate the relevant Hankel matrix. Has to be specified when using \code{\link{nonparamHankel}}, \code{\link{paramHankel}} or \code{\link{paramHankel.scaled}}. For further details see the details section of \code{\link{datMix}}.
#' @param Hankel.function function required for the moment estimation via \code{Hankel.method}. This normally depends on \code{Hankel.method} as well as \code{dist}. For further details see the \code{\link{datMix}} details section.
#' @return Object of class \code{datMix} with the following attributes (for further explanations see above):
#'      \item{dist}{character string giving the abbreviated name of the component distribution, such that the function \code{ddist} evaluates its density/mass and \code{rdist} generates random variates.}
#'      \item{discrete}{logical flag indicating whether the mixture distribution is discrete.}
#'      \item{theta.bound.list}{named list specifying the upper and lower bounds for the component parameters.}
#'      \item{MLE.function}{function which computes the MLE of the component distribution \code{dist}.}
#'      \item{Hankel.method}{character string taking on values \code{"explicit"}, \code{"translation"}, or \code{"scale"}, specifying the method of estimating the moments of the mixing distribution to compute the corresponding Hankel matrix.}
#'      \item{Hankel.function}{function required for the moment estimation via \code{Hankel.method}. See details for more information.}
#' @seealso \code{\link{datMix}} for direct generation of a \code{datMix} object from a vector of observations.
#' @keywords cluster
#' @examples
#' ### generating 'Mix' object
#' normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
#'                   sd = c(1, 1, 1))
#'
#' ### generating 'rMix' from 'Mix' object (with 1000 observations)
#' set.seed(1)
#' normLocRMix <- rMix(1000, normLocMix)
#'
#' ### generating 'datMix' from 'R' object
#'
#' ## generate list of parameter bounds
#'
#' norm.bound.list <- vector(mode = "list", length = 2)
#' names(norm.bound.list) <- c("mean", "sd")
#' norm.bound.list$mean <- c(-Inf, Inf)
#' norm.bound.list$sd <- c(0, Inf)
#'
#' ## generate MLE functions
#'
#' # for "mean"
#' MLE.norm.mean <- function(dat) mean(dat)
#' # for "sd" (the sd function uses (n-1) as denominator)
#' MLE.norm.sd <- function(dat){
#'   sqrt((length(dat) - 1) / length(dat)) * sd(dat)
#' }
#' # combining the functions to a list
#' MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
#'                       "MLE.norm.sd" = MLE.norm.sd)
#'
#' ## function giving the j^th raw moment of the standard normal distribution,
#' ## needed for calculation of the Hankel matrix via the "translation" method
#' ## (assuming gaussian components with variance 1)
#'
#' mom.std.norm <- function(j){
#'  ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
#' }
#'
#' normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
#'                      MLE.function = MLE.norm.list, Hankel.method = "translation",
#'                      Hankel.function = mom.std.norm)
#'
#' ### using 'datMix' object to estimate the mixture
#'
#' set.seed(0)
#' res <- paramHankel.scaled(normLoc.dM)
#' plot(res)
#'
#' @export RtoDat
RtoDat <- function(obj, theta.bound.list = NULL, MLE.function = NULL, Hankel.method = NULL,
                   Hankel.function = NULL){

  new <- obj
  att <- attributes(obj)
  # discard other information like component parameter values as these are
  # supposed to be estimated from the "datMix" object
  attributes(new) <- att[names(att) == "dist" | names(att) == "discrete"]
  class(new) <- "datMix"
  attr(new, "theta.bound.list") <- theta.bound.list
  attr(new, "MLE.function") <- MLE.function
  attr(new, "Hankel.method") <- Hankel.method
  attr(new, "Hankel.function") <- Hankel.function
  return(new)

}
