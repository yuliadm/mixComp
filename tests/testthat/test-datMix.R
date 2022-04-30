## tests for datMix function: 
   
testthat::test_that("errors", {
   ## missing mandatory argument:  
    testthat::expect_error(
      datMix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1)),
        'unused arguments (w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))', fixed = TRUE)  
})

testthat::test_that("object checks", {
    
    children.obs <- unlist(children)
    # define the MLE function:
    MLE.pois <- function(dat) mean(dat)
    # construct a datMix object:
    children.dM <- datMix(children.obs, dist = "pois", discrete = TRUE, 
                        Hankel.method = NULL, 
                        Hankel.function = NULL,
                        theta.bound.list = list(lambda = c(0, Inf)), 
                        MLE.function = MLE.pois)
  
    ## class for a Mix object:  
    testthat::expect_equal( 
      class(children.dM), "datMix")
    
     ## is.datMix for a datMix object:    
     testthat::expect_equal(
       is.datMix(children.dM), TRUE)
      
     ## input is not a datMix object:
     testthat::expect_equal(
       is.datMix(children.obs), FALSE)  
})    
