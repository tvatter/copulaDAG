context("Function 'copula_dag'")

test_that("function works", {
  
  
  set.seed(123)
  n <- 5e2
  X <- rnorm(n)
  Y <- X ^ 2 + rnorm(n)
  Z <- Y / 2 + 1.3 * rnorm(n)
  W <- X + rnorm(n) ^ 3
  
  dag_data <- cbind(X, Y, Z, W)
  trueDAG <- cbind(c(0, 0, 0, 0), c(1, 0, 0, 0), c(0, 1, 0, 0), c(1, 0, 0, 0))
  
  estDAG <- copula_dag(dag_data, alpha = 0.01)
  
  expect_true(
    is.matrix(trueDAG) &&
      is.matrix(estDAG$Adj) &&
      dim(trueDAG) == dim(estDAG$Adj) && all(trueDAG == estDAG$Adj)
  )
})

test_that("sanity checks work", {
  
  set.seed(123)
  n <- 5e2
  X <- rnorm(n)
  Y <- X ^ 2 + rnorm(n)
  Z <- Y / 2 + 1.3 * rnorm(n)
  W <- X + rnorm(n) ^ 3
  
  # should throw an error if there are identical variables (? or to handle differently)
  dag_data <- cbind(X, Y, X, W)
  expect_error(copula_dag(dag_data, alpha = 0.01))
  
  dag_data <- cbind(X, Y)
  expect_error(copula_dag(dag_data, alpha = 0.01))
  expect_error(copula_dag(dag_data, alpha = "a"))
  expect_error(copula_dag(data.frame(rep("a",n), rep("1",n), rep("b",n)), alpha = "a"))
  
  
})