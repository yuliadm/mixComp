## tests for dMix function:    

testthat::test_that("errors", {
  ## argument obj is not of class 'Mix':  
  testthat::expect_error(
    dMix(c(1,3,5), c(1, 2, 3)),
    "obj' has to be a 'Mix' object!")  
  
  ## missing argument:  
testthat::expect_error(
    dMix(100),
    'argument "obj" is missing, with no default')  
  
  ## the point is not in the 'dgeom' function support:   
  testthat::expect_warning(
    dMix(0.56, Mix("geom", discrete = TRUE, w = c(0.1, 0.9), prob = c(0.2, 0.5))))
})   


testthat::test_that("obects checks", {    
  set.seed(0) 
  testthat::expect_equal(
    dMix(c(1, 2, 3), Mix("geom", discrete = TRUE, w = c(0.1, 0.3, 0.6), prob = c(0.3, 0.4,0.8))) , c(0.18900, 0.07710, 0.04005) )
  
  ## check class of the output:  
  testthat::expect_equal( 
    class(dMix(10, Mix("geom", discrete = TRUE, w = c(0.1, 0.9), prob = c(0.2, 0.5)))), "numeric")
})   
 


