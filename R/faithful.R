#' Faithful Dataset
#'
#' Waiting time between eruptions and the duration of the eruption 
#' for the Old Faithful geyser in Yellowstone National Park, 
#' Wyoming, USA from \code{\link{datasets}}.
#'
#' @name faithful
#' @docType data
#' @usage data(faithful)
#' @format A data frame with 272 observations on 2 variables:
#'      \describe{
#'          \item{eruptions}{numeric, eruption time in mins}
#'          \item{waiting}{numeric, waiting time to next eruption (in mins)}
#' }
#' @keywords datasets
#' @source 
#' \enumerate{
#' \item
#' Azzalini, A. and Bowman, A. W. (1990). 
#' A look at some data on the Old Faithful geyser. Applied Statistics, 39, 357--365.
#' 
#'
#' \item \code{\link{datasets}}
#' }
#' @examples
#' data(faithful)
#' 
#' faithful.obs <- faithful$waiting
#' 
#' # function giving the j^th raw moment of the standard normal distribution,
#' # needed for calculation of the Hankel matrix via the "translation" method
#' # (assuming gaussian components with variance 1)
#' mom.std.norm <- function(j){
#'   ifelse(j %% 2 == 0, prod(seq(1, j - 1, by = 2)), 0)
#' }
#' 
#' # generate list of parameter bounds
#' norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
#' 
#' # define the MLE functions for the mean and sd: 
#' MLE.norm.mean <- function(dat) mean(dat)
#' MLE.norm.sd <- function(dat){
#' sqrt((length(dat) - 1) / length(dat)) * sd(dat)
#' } 
#' MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)
#' 
#' # construct a 'datMix' object that summarizes all the necessary information:
#' faithful.dM <- datMix(faithful.obs, dist = "norm", discrete = FALSE,
#'                       theta.bound.list = norm.bound.list,
#'                       MLE.function = MLE.norm.list, Hankel.method = "translation",
#'                       Hankel.function = mom.std.norm)
#' 
#' # estimate the number of components and plot the results:
#' res <- hellinger.cont(faithful.dM, bandwidth = 4,
#'                       sample.n = 5000, threshold = "AIC")
#' plot(res)
NULL
