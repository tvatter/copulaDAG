context("Function 'pairwise_direction'")

test_that("function works", {
  
  set.seed(0)
  n <- 1e3
  NX <- rnorm(n)
  NY <- rnorm(n)
  X <- NX
  Y <- X^2 + NY
  
  expect_true(pairwise_direction(X, Y))
  expect_true(pairwise_direction(X, Y, bicop(cbind(rank(X)/(n+1), 
                                                   rank(Y)/(n+1)), 
                                             family_set = "tll", 
                                             nonpar_method = "constant")))
  expect_true(pairwise_direction(X, Y, 
                                 family_set = "tll", 
                                 nonpar_method = "constant"))
})

test_that("sanity checks work", {
  
  expect_error(pairwise_direction(hist(1:3, plot = FALSE)))
  
  expect_error(pairwise_direction(matrix(1,3,3)))
  expect_error(pairwise_direction(matrix(1,3,2), 1:2))
  expect_error(pairwise_direction(data.frame(a="a",b=1)))
  
  expect_error(pairwise_direction(1:3))
  expect_error(pairwise_direction(1:3, 1:2))
  expect_error(pairwise_direction(rep("a", 3), 1:3))
  expect_error(pairwise_direction(1:3, rep("a", 3)))
  
  expect_error(pairwise_direction(1:3, 1:3, 12))
})
