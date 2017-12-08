#' Copula-based direction for a single pair
#'
#' @param x a numeric vector, bivariate matrix or data frame.
#' @param y `NULL` (default) or a vector with compatible dimensions to `x`.
#' @param cop `NULL` (default) or a 
#' \code{\link[rvinecopulib:bicop_dist]{bicop_dist}} object.
#' @param ... additional arguments passed to 
#' \code{\link[rvinecopulib:bicop]{bicop}}
#'
#' @details
#' Assuming that either (a) X := N_X and Y = f(X, N_Y) or (b) Y := N_Y and 
#' X = f(Y, N_X), where N_X and N_Y are independent noise processes, 
#' this copula-based method recovers the direction, 
#' namely (a) X -> Y or (b) Y -> X.
#' 
#' @return TRUE if X -> Y and FALSE if Y -> X
#'
#' @examples
#' # Y = X^2 + N_Y with Gaussian noises
#' n <- 1e3
#' NX <- rnorm(n)
#' NY <- rnorm(n)
#' X <- NX
#' Y <- X^2 + NY
#'
#' pairwise_direction(X, Y)
#' @export
#' @importFrom assertthat assert_that
#' @importFrom Hmisc hoeffd
#' @importFrom stats predict
#' @importFrom rvinecopulib bicop
#' @importFrom utils modifyList
#' @importFrom energy eqdist.etest
#' @importFrom statmod gauss.quad.prob
pairwise_direction <- function(x, y = NULL, cop = NULL, ...) {
  
  # basic sanity checks
  assert_that(is.vector(x) || is.matrix(x) || is.data.frame(x),
              msg = "x should be a vector, matrix or data.frame")
  if (is.matrix(x) || is.data.frame(x)) {
    assert_that(ncol(x) == 2, 
                msg = "x is a matrix or data.frame but ncol(x) is not equal to 2")
    assert_that(is.null(y), 
                msg = "x is a matrix or data.frame but y is not NULL")
    assert_that(all(apply(x,2,is.numeric)), 
                msg = "all the elements of x should be numeric")
    y <- x[,2]
    x <- x[,1]
  } else {
    assert_that(!is.null(y),
                msg = "x is a vector but y is NULL")
    assert_that(length(x) == length(y))
    assert_that(is.numeric(x))
    assert_that(is.numeric(y))
  }
  if (!is.null(cop)) {
    assert_that(any(class(cop) == "bicop"), 
                msg = "cop should be a bicop object")
  }
  
  # get pseudo-observations (i.e. F(X) and F(Y))
  n <- length(x)
  u1 <- rank(x)/(n+1)
  u2 <- rank(y)/(n+1)
  
  # if copula not provided, estimate it
  if (is.null(cop)) {
    # basic arguments
    pars <- list(data = cbind(u1,u2),
                 family_set = "tll",
                 nonpar_method = "constant")
    # additional arguments
    pars <- modifyList(pars, list(...))
    # estimation
    cop <- do.call(bicop, pars)
  }
  
  uw <- gauss.quad.prob(5)
  h <- sapply(uw$nodes, function(uu) {
    u1p <- predict(object = cop, newdata = cbind(uu, u2), what = "hinv2")
    u2p <- predict(object = cop, newdata = cbind(u1, uu), what = "hinv1")
    
    xp <- quantile(x, u1p)
    yp <- quantile(y, u2p)
    h1 <- abs(cor(x, xp))
    h2 <- abs(cor(y, yp))
    
    ifelse(h1 != h2, h1 < h2, NA)
  })
  sum(uw$weights[!is.na(h)]*h[!is.na(h)])
}