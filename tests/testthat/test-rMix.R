## tests for rMix function:    

testthat::test_that("errors", {
  ## argument obj is not of class 'Mix':  
  testthat::expect_error(
    rMix(1000, c(1, 2, 3)),
    "obj' has to be a 'Mix' object!")  
  
  ## number of samples is not a natural number:   
  testthat::expect_error(
    rMix(0.56, Mix("geom", discrete = TRUE, w = c(0.1, 0.9), prob = c(0.2, 0.5))),
    "'n' has to be a single integer!")
})   

    
testthat::test_that("obects checks", {    
  set.seed(0) 
  testthat::expect_equal(
    as.vector(rMix(5, Mix("geom", discrete = TRUE, w = c(1), prob = c(0.3)))), c(7, 0, 3, 0, 1) )
  
  ## class for an rMix object:  
  testthat::expect_equal( 
    class(rMix(1000, Mix("geom", discrete = TRUE, w = c(0.1, 0.9), prob = c(0.2, 0.5)))), "rMix")
  
  ## input is not an Mix object:
  testthat::expect_equal(
    is.rMix(c(1, 2, 3)),
    FALSE)  
})   

