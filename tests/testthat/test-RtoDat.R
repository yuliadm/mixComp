## tests for RtoDat function:    

testthat::test_that("errors", {
  ## missing mandatory argument:  
  testthat::expect_error(
    RtoDat(theta.bound.list = NULL, MLE.function = NULL, Hankel.method = NULL,
           Hankel.function = NULL),
    'argument "obj" is missing, with no default')  
  
})   


testthat::test_that("obects checks", {    
  normLocMix <- Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17),
                    sd = c(1, 1, 1))
  ### generate 'rMix' from 'Mix' object (1000 obs)
  set.seed(0)
  normLocRMix <- rMix(1000, normLocMix)
  
  ## generate list of parameter bounds
  norm.bound.list <- vector(mode = "list", length = 2)
  names(norm.bound.list) <- c("mean", "sd")
  norm.bound.list$mean <- c(-Inf, Inf)
  norm.bound.list$sd <- c(0, Inf)
  
  normLoc.dM <- RtoDat(normLocRMix, theta.bound.list = norm.bound.list)
  
  ## check class of the output:   
  testthat::expect_equal(
    class(normLoc.dM), "datMix")
})   


