#' Shakespeare Dataset
#'
#' Shakespeare's word type frequencies data from Efron and Thisted (1976).
#'
#' @name shakespeare
#' @docType data
#' @usage data(shakespeare)
#' @format A data frame with 30792 observations on 1 variable. Replicates are generated
#' to reflect the frequencies of word types (words used exactly n times n = 1, 2, ..., 100).
#' As there are 14376 word types that were used once, 1 appears 14376 times in the data,
#' as there are 4343 word types that were used twice, 2 appears 4343 times in the data, etc.
#' @keywords datasets
#' @source Efron, B. and Thisted, R. (1976). Estimating the number of unseen species: 
#' how many words did Shakespeare know? Biometrka 63 435-447.
#' @examples
#' data(shakespeare)
#' 
#' shakespeare.obs <- unlist(shakespeare) - 1
#' 
#' # define the MLE function:
#' MLE.geom <- function(dat) 1 / (mean(dat) + 1)
#' 
#' Shakespeare.dM <- datMix(shakespeare.obs, dist = "geom", discrete = TRUE, 
#'                          MLE.function = MLE.geom,
#'                          theta.bound.list = list(prob = c(0, 1)))
#'
#' # estimate the number of components and plot the results:
#' \donttest{
#' set.seed(0)
#' res <- hellinger.boot.disc(Shakespeare.dM, B = 50, ql = 0.025, qu = 0.975)
#' plot(res, breaks = 100, xlim = c(0, 20))
#' }
NULL
