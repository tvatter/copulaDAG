#' Copula-based Directed Acyclic Graph
#'
#' @param x a matrix or data frame with more than two columns 
#' (see \code{\link{pairwise_direction}} for the bivariate case)
#' @param alpha a scalar number in (0,1) indicating the level of the test.
#' @param skelecton a logical indicating whether the skeleton only should be 
#' returned.
#' @param ... additional arguments passed to 
#' \code{\link[rvinecopulib:bicop]{bicop}}
#'
#' @details
#' TODO
#' 
#' @return TODO
#'
#' @examples
#' TODO 
#' @export
#' @importFrom utils combn
#' @importFrom igraph graph_from_edgelist
#' @importFrom assertthat is.number is.flag
copula_dag <- function(x, alpha = 0.1, skeleton = FALSE, ...) {
  
  # basic sanity checks
  assert_that(is.matrix(x) || is.data.frame(x),
              msg = "x should be a vector, matrix or data.frame")
  assert_that(ncol(x) > 2)
  assert_that(all(apply(x,2,is.numeric)), 
              msg = "all the elements of x should be numeric")
  assert_that(is.number(alpha) && alpha > 0 && alpha < 1,
              msg = "alpha should be a scalar in (0,1)")
  assert_that(is.flag(skeleton))
  
  # get pseudo-observations
  n <- nrow(x)
  d <- ncol(x)
  udata <- apply(x, 2, rank)/(n+1)
  
  # basic arguments
  pars <- list(data = matrix(),
               family_set = "tll", 
               nonpar_method = "constant")
  # additional arguments
  pars <- modifyList(pars, list(...))
  
  # get all combinations of two variables
  all_var <- 1:d
  all_comb <- t(combn(all_var, 2))
  
  # fit a models and returns pdf required for the residuals for every combination
  system.time(copdata <- lapply(1:nrow(all_comb), function(j) {
    # add data
    pars$data <- udata[,all_comb[j,]]
    # fit a bivariate copula
    cop <- do.call(bicop, pars)
    # get pdf required for the residuals
    outer(pars$data[,1], pars$data[,2], 
          function(u1, u2) predict(object = cop, 
                                   newdata = cbind(u1, u2), 
                                   what = "pdf"))
  }))
  
  system.time(p_ures <- lapply(1:nrow(all_comb), function(k) {
    ures <- get_pred(all_comb[k,1], all_comb[k,2], all_comb, udata, copdata)
    list(p = hoeffd(ures[,1], ures[,2])$P[2],
         ures = ures)
  }))

  # select pairs by p-values
  sel_comb <- sapply(p_ures, function(x) x$p < alpha)
  
  if (skeleton) {
    sel_comb <- all_comb[sel_comb, ]
  } else {
    # get direction
    # dir_comb <- sapply(seq_along(sel_comb), 
    #                    function(comb) {
    #                      if (sel_comb[comb] == TRUE) {
    #                        ures <- p_ures[[comb]]$ures
    #                        ifelse(rep(pairwise_direction(ures), 2), 
    #                               all_comb[comb,], rev(all_comb[comb,]))
    #                      } else {
    #                        NULL
    #                      }
    #                    })
    # dir_comb <- matrix(unlist(dir_comb), ncol = 2, byrow = TRUE)
    sel_comb <- apply(all_comb[sel_comb,], 1, function(comb)
      ifelse(rep(pairwise_direction(x[,comb]), 2), rev(comb), comb))
    sel_comb <- t(sel_comb)
  }

  # create and return graph
  return(graph_from_edgelist(sel_comb, directed = !skeleton))
}

get_ij_combs <- function(i, j, all_comb) {
  sapply(1:nrow(all_comb), function(k) {
    ifelse((all_comb[k, 1] == i && all_comb[k, 2] != j) || 
             (all_comb[k, 2] == i && all_comb[k, 1] != j), 
           1, ifelse((all_comb[k, 1] != i && all_comb[k, 2] == j) || 
                       (all_comb[k, 2] != i && all_comb[k, 1] == j), 
                     2, 0))
  }) 
}

if_vec_to_matrix <- function(u) {
  if (NCOL(u) == 1)
    u <- matrix(u, 1, length(u))
  if (!is.matrix(u))
    u <- as.matrix(u)
  
  u
}

get_compl <- function(i, icombs) {
  apply(if_vec_to_matrix(icombs), 1, function(k) ifelse(i == k[1], k[2], k[1]))
}

get_i_data <- function(i, sel_combs, all_comb, copdata) {
  i_combs <- all_comb[sel_combs,]
  i_compl <- get_compl(i, i_combs)
  i_data <- copdata[sel_combs]
  i_data <- lapply(1:length(i_compl), function(k) {
    if (i_compl[k] < i) {
      return(t(i_data[[k]]))
    } else {
      return(i_data[[k]])
    }
  })
  if (length(i_compl) > 1) {
    return(Reduce("*", i_data))
  } else {
    return(i_data[[1]])
  }
}

get_pred <- function(i, j, all_comb, udata, copdata) {
  ij_combs <- get_ij_combs(i, j, all_comb)
  i_data <- get_i_data(i, ij_combs == 1, all_comb, copdata)
  j_data <- get_i_data(j, ij_combs == 2, all_comb, copdata)
  
  ui <- udata[,i]
  uj <- udata[,j]
  predi <- sapply(1:ncol(i_data), function(k) 
    sum(i_data[,k]*ifelse(ui <= ui[k], 1, 0))/sum(i_data[,k]))
  predj <- sapply(1:ncol(j_data), function(k) 
    sum(j_data[,k]*ifelse(uj <= uj[k], 1, 0))/sum(j_data[,k]))
  cbind(predi, predj)
}
