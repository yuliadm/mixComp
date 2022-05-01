
## tests for nonparamHankel function:    
testthat::test_that("errors", {
  
  poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
  ## create random data based on 'Mix' object (gives back 'rMix' object)
  set.seed(0)
  poisRMix <- rMix(1000, obj = poisMix)
  # specify the list of parameter bounds:
  poisList <- list("lambda" = c(0, Inf))
  # define the MLE function:
  MLE.pois <- function(dat){
    mean(dat)
  }
  # define function for estimating the j^th moment of the mixing distribution:
  explicit.pois <- function(dat, j){
    res <- 1
    for (i in 0:(j-1)){
      res <- res*(dat-i)
    }
    return(mean(res))
  }
  # create a 'datMix' object:
  pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, MLE.function = MLE.pois,
                    Hankel.method = "explicit", Hankel.function = NULL)
  set.seed(0)
  
  ## missing mandatory argument in the input datMix object:  
  testthat::expect_error(
    nonparamHankel(pois.dM),
    "'Hankel.function' has to be defined as a datMix object attribute!", fixed = TRUE) 
  
  # update the 'datMix' object:
  pois.dM <- RtoDat(poisRMix, MLE.function = MLE.pois, Hankel.function = explicit.pois)
  
  set.seed(0)
  ## missing mandatory argument in the input datMix object:  
  testthat::expect_error(
    nonparamHankel(pois.dM),
    "'Hankel.method' has to be defined as a datMix object attribute!", fixed = TRUE) 
  
  set.seed(0)
  ## input is not of class 'datMix':   
  testthat::expect_error(
    nonparamHankel(1))
})

testthat::test_that("function checks", {
  
  poisMix <- Mix("pois", discrete = TRUE, w = c(0.45, 0.45, 0.1), lambda = c(1, 5, 10))
  ## create random data based on 'Mix' object (gives back 'rMix' object)
  set.seed(0)
  poisRMix <- rMix(1000, obj = poisMix)
  # specify the list of parameter bounds:
  poisList <- list("lambda" = c(0, Inf))
  # define the MLE function:
  MLE.pois <- function(dat){
    mean(dat)
  }
  # define function for estimating the j^th moment of the mixing distribution:
  explicit.pois <- function(dat, j){
    res <- 1
    for (i in 0:(j-1)){
      res <- res*(dat-i)
    }
    return(mean(res))
  }
  # create a 'datMix' object:
  pois.dM <- RtoDat(poisRMix, theta.bound.list = poisList, MLE.function = MLE.pois,
                    Hankel.method = "explicit", Hankel.function = explicit.pois)
  
  ## check class of the output:  
  testthat::expect_equal( 
    class(nonparamHankel(pois.dM, j.max = 5)), "hankDet")
  
  set.seed(0)    
  res <- nonparamHankel(pois.dM, j.max = 5, scaled = TRUE, B = 100)
  ## check output:    
  testthat::expect_equal(
    round(res[1],3), 12.638)
  
  ## check output:    
  testthat::expect_equal(
    round(res[4],3), 0.468)
})    

