## tests for hellinger.disc function:    

testthat::test_that("errors", {
  ## Gaussian mixture (continuous distribution)
  ### generate a 'Mix' object:
  normLocMix <- Mix("norm", discrete = FALSE, w = c(0.1, 0.2, 0.4, 0.3), mean = c(-5, -2, 0, 3),
                    sd = c(1, 1, 1, 1))
  ### generate 'rMix' from 'Mix' object (with 1000 observations):
  set.seed(0)
  normLocRMix <- rMix(500, normLocMix)
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
  normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list, MLE.function = MLE.norm.list)
  
  ## Poisson mixture (discrete distribution)
  ## create 'Mix' object
  poisMix <- Mix("pois", discrete = TRUE, w = c(0.2, 0.3, 0.5), lambda = c(1, 3, 6))
  ## generate random sample:
  set.seed(0)
  poisRMix <- rMix(1000, obj = poisMix)
  # specify list of parameter bounds:
  poisList <- list("lambda" = c(0, Inf))
  # define the MLE function:
  MLE.pois <- function(dat){
    mean(dat)
  }
  # generate the 'datMix' object:
  pois.dM <- RtoDat(poisRMix) 
  
  ## missing mandatory argument in datMix::  
  testthat::expect_error(
    hellinger.disc(pois.dM), "'theta.bound.list' has to be specified as a datMix object attribute!")  
  
  # update 'datMix':
  pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, MLE.function = MLE.pois) 
   
  ## input of class other than datMix:  
  testthat::expect_error(
    hellinger.disc(poisRMix))
  
  ## input from a continuous distribution:
  testthat::expect_error(
    hellinger.disc(normLoc.dM), "The 'discrete' attribute of the Rdat object is not set to TRUE, however
         this function only works for discrete data", fixed = TRUE)
  
  # expect warning when using AIC threshold:
  testthat::expect_warning(
    hellinger.disc(pois.dM, threshold="AIC"))
  
})   


testthat::test_that("function checks", {    
  ## create 'Mix' object
  poisMix <- Mix("pois", discrete = TRUE, w = c(0.2, 0.3, 0.5), lambda = c(1, 3, 6))
  ## generate random sample:
  set.seed(0)
  poisRMix <- rMix(1000, obj = poisMix)
  # specify list of parameter bounds:
  poisList <- list("lambda" = c(0, Inf))
  # define the MLE function:
  MLE.pois <- function(dat){
    mean(dat)
  }
  # generate the 'datMix' object:
  pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, MLE.function = MLE.pois) 
  
  # overwrite some default parameters for speeding up the computations:
  set.seed(0)
  res <- hellinger.disc(pois.dM)
  
  ## check the class of the output:   
  testthat::expect_equal(
    class(res), "paramEst")
  
  ## check the output:  
  testthat::expect_equal(
    res[1], 2)
})   


