#' Acidity Dataset
#'
#' @description Data for the Acidity index on the log-scale for 155 lakes 
#' in North-Central Wisconsin from the Eastern Lake Survey.
#' The measurements are acid neutralizing capacity (ANC) on the log scale; 
#' specifically, log(ANC + 50).
#' The Acidity index describes the capability of a lake to absorb acid; 
#' low ANC values can lead to a loss of biological resources, see Crawford (1994).
#'
#' @name acidity
#' @docType data
#' @usage data(acidity)
#' @format A data frame with 155 observations on 1 variable.
#' @keywords datasets
#' @source Crawford et al. (1992) Modeling Lake-Chemistry Distributions: 
#' Approximate Bayesian Methods for Estimating a Finite-Mixture Model, Technometrics, 34:4, 441-453
#' @examples
#' data(acidity)
#' 
#' acidity.obs <- unlist(acidity)
#' 
#' # define the MLE functions for the mean and sd: 
#' MLE.norm.mean <- function(dat) mean(dat)
#' MLE.norm.sd <- function(dat){
#'   sqrt((length(dat) - 1) / length(dat)) * sd(dat)
#' } 
#' MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean, "MLE.norm.sd" = MLE.norm.sd)
#'
#' # define the range for parameter values:
#' norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
#' 
#' # create 'datMix' object:
#' acidity.dM <- datMix(acidity.obs, dist = "norm", discrete = FALSE, 
#'                      MLE.function = MLE.norm.list, 
#'                      theta.bound.list = norm.bound.list)
#'                      
#' \donttest{                      
#' set.seed(0)
#' res <- mix.lrt(acidity.dM, B = 50, quantile = 0.95)
#' plot(res)
#' }
NULL
