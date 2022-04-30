## tests for Mix function: 
   
testthat::test_that("errors", {
  ## missing mandatory argument:  
  testthat::expect_error(
    Mix("geom", w = c(0.2, 0.8), prob = c(0.2, 0.7)),
    'argument "discrete" is missing, with no default')  
  
  ## weight and parameter vectors of different lengths:   
  testthat::expect_error(
    Mix("geom", discrete = TRUE, w = c(0.2, 0.8), prob = c(0.2, 0.5, 0.7)),
    "'w' must be a numeric >= 0 with same length as the elements of theta.list
           (or the inputs to ...)", fixed = TRUE)
})   

   
testthat::test_that("obects checks", {    
  ## class for a Mix object:  
  testthat::expect_equal( 
    class(Mix("geom", discrete = TRUE, w = c(0.2, 0.8), prob = c(0.2, 0.7))), "Mix")
  ## is.Mix for a Mix object:    
  
  testthat::expect_equal(
    is.Mix(Mix("norm", discrete = FALSE, w = c(0.3, 0.4, 0.3), mean = c(10, 13, 17), sd = c(1, 1, 1))),
    TRUE)
  
  ## input is not a Mix object:
  testthat::expect_equal(
    is.Mix(NULL),
    FALSE)  
})   
