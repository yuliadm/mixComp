
## tests for mix.lrt function:    

testthat::test_that("errors", {
  ### generating 'Mix' object
  normLocMix <- Mix("norm", discrete = FALSE, w = c(0.25, 0.25, 0.25, 0.25), mean = c(5, 10, 13, 17),
                    sd = c(1, 1, 1, 1))
  ### generating 'rMix' from 'Mix' object (with 1000 observations)
  set.seed(0)
  normLocRMix <- rMix(1000, normLocMix)
  ## specify list of parameter bounds
  norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
  ## generate MLE functions
  MLE.norm.mean <- function(dat) mean(dat)
  MLE.norm.sd <- function(dat){
    sqrt((length(dat) - 1) / length(dat)) * sd(dat)
  }
  # combining the functions to a list
  MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
                        "MLE.norm.sd" = MLE.norm.sd)
  
  ## generating 'datMix' object
  normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                       MLE.function = MLE.norm.list)
  
  set.seed(0)
  ## missing mandatory argument:  
  testthat::expect_error(
    mix.lrt(normLoc.dM, B = NULL))  
  
  set.seed(0)
  ## input of class other than datMix:  
  testthat::expect_error(
    mix.lrt(normLocRMix, B = 50))  
})   


testthat::test_that("obects checks", {    
  ### generate a 'Mix' object:
  normLocMix <- Mix("norm", discrete = FALSE, w = c(0.25, 0.25, 0.25, 0.25), mean = c(5, 10, 13, 17),
                    sd = c(1, 1, 1, 1))
  ### generate 'rMix' from 'Mix' object (with 1000 observations):
  set.seed(0)
  normLocRMix <- rMix(1000, normLocMix)
  ## specify list of parameter bounds:
  norm.bound.list <- list("mean" = c(-Inf, Inf), "sd" = c(0, Inf))
  ## define MLE functions:
  MLE.norm.mean <- function(dat) mean(dat)
  MLE.norm.sd <- function(dat){
    sqrt((length(dat) - 1) / length(dat)) * sd(dat)
  }
  # combine the functions to a list:
  MLE.norm.list <- list("MLE.norm.mean" = MLE.norm.mean,
                        "MLE.norm.sd" = MLE.norm.sd)
  ## create a 'datMix' object:
  normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list,
                       MLE.function = MLE.norm.list)
  
  set.seed(0)
  res <- mix.lrt(normLoc.dM, B = 100)
  
  ## check the class of the output:   
  testthat::expect_equal(
    class(res), "paramEst")
  
  set.seed(0)
  ## check the output:   
  testthat::expect_equal(
    res[1], 4)
})   


