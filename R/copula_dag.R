#' Copula-based Directed Acyclic Graph
#'
#' @param x a matrix or data frame with more than two columns 
#' (see \code{\link{pairwise_direction}} for the bivariate case)
#' @param alpha a scalar number in (0,1) indicating the level of the test.
#' @param pc a logical indicating whether the PC algorithm should be used 
#' to find the skeleton.
#' @param skeleton a logical indicating whether the skeleton only should be 
#' returned.
#' @param fast a logical indicating whether the memory greedy (but faster) 
#' method should be used.
#' @param ... additional arguments passed to 
#' \code{\link[rvinecopulib:bicop]{bicop}}
#'
#' @details
#' TODO
#' 
#' @return estimated causal DAG and its adjacency matrix
#'
#' @examples
#' 
#' set.seed(123)
#' n <- 500
#' X <- rnorm(n)
#' Y <- X ^ 2 + rnorm(n)
#' Z <- Y / 2 + 1.3 * rnorm(n)
#' W <- X + rnorm(n) ^ 3
#' dag_data <- cbind(X, Y, Z, W)
#' trueDAG <- cbind(c(0, 0, 0, 0), c(1, 0, 0, 0), c(0, 1, 0, 0), c(1, 0, 0, 0))
#'
#' estDAG <- copula_dag(dag_data, alpha = 0.01)
#'
#' @export
#' @importFrom utils combn
#' @importFrom igraph graph_from_edgelist graph_from_adjacency_matrix as_edgelist as_adjacency_matrix
#' @importFrom assertthat is.number is.flag
#' @importFrom pcalg skeleton
copula_dag <- function(x, alpha = 0.1, 
                       pc = TRUE, 
                       skeleton = FALSE, 
                       fast = TRUE, ...) {
  
  # basic sanity checks
  assert_that(is.matrix(x) || is.data.frame(x),
              msg = "x should be a vector, matrix or data.frame")
  assert_that(ncol(x) > 2)
  assert_that(all(apply(x,2,is.numeric)), 
              msg = "all the elements of x should be numeric")
  assert_that(any(duplicated(x, MARGIN = 2)) == FALSE,
              msg = "variables in x cannot be identical")
  assert_that(is.number(alpha) && alpha > 0 && alpha < 1,
              msg = "alpha should be a scalar in (0,1)")
  assert_that(is.flag(pc))
  assert_that(is.flag(skeleton))
  assert_that(is.flag(fast))
  
  # get pseudo-observations
  n <- nrow(x)
  if (fast && n > 1e3) {
    msg <- paste0("fast = TRUE is not available for sample sizes larger ",
                  "than 1e3, switching to fast = FALSE.", collapse = "")
    message(msg)
    fast <- FALSE
  }
    
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
  
  # fit a model for every combination
  copdata <- lapply(1:nrow(all_comb), function(j) {
    # add data
    pars$data <- udata[,all_comb[j,]]
    # fit a bivariate copula
    cop <- do.call(bicop, pars)
  })
  
  # for small sample sizes, we can pre-compute the require pdfs
  if (fast) {
    # return pdf required for the residuals for every combination
    copdata <- lapply(1:nrow(all_comb), function(j) {
      # get pdf required for the residuals
      outer(udata[,all_comb[j,1]], udata[,all_comb[j,2]], 
            function(u1, u2) predict(object = copdata[[j]], 
                                     newdata = cbind(u1, u2), 
                                     what = "pdf"))
    })
  }
  
  
  if (pc == TRUE) {
    suffStat <- list(all_comb = all_comb, 
                     udata = udata, 
                     copdata = copdata)
    g_nel <- skeleton(suffStat, hoeffCItest, alpha = alpha, p = d)
    g_igraph <- graph_from_adjacency_matrix(as(g_nel, "matrix"), "undirected")
    sel_comb <- apply(as_edgelist(g_igraph), 2, as.numeric)
  } else {
    p_ures <- lapply(1:nrow(all_comb), function(k) {
      ures <- get_pred(all_comb[k,1], all_comb[k,2], setdiff(all_var, 
                                                             all_comb[k,]), 
                       all_comb, udata, copdata)
      list(p = hoeffd(ures[,1], ures[,2])$P[2],
           ures = ures)
    })
    
    # select pairs by p-values
    sel_comb <- all_comb[sapply(p_ures, function(x) x$p < alpha), ]
  }
  
  if (!skeleton) {
    sel_comb <- apply(sel_comb, 1, function(comb)
      # if TRUE, rev order 
      # pairwise_direction(x[,comb]) on raw data, not on residuals
      ifelse(rep(pairwise_direction(x[,comb]), 2), rev(comb), comb))
    sel_comb <- t(sel_comb)
  }
  
  # create and return graph and its adjacency matrix
  estimated_graph <- graph_from_edgelist(sel_comb, directed = !skeleton)
  return(list(graph <- estimated_graph, Adj = as_adjacency_matrix(estimated_graph, sparse = FALSE)))
}

get_ijk_combs <- function(i, j, k, all_comb) {
  sapply(1:nrow(all_comb), function(l) {
    ifelse((all_comb[l, 1] == i && all_comb[l, 2] %in% k) || 
             (all_comb[l, 2] == i && all_comb[l, 1] %in% k), 
           1, ifelse((all_comb[l, 1] %in% k && all_comb[l, 2] == j) || 
                       (all_comb[l, 2] %in% k && all_comb[l, 1] == j), 
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
  is_bicop <- !is.matrix(i_data[[1]])
  i_data <- lapply(1:length(i_compl), function(k) {
    # compute pdf if not already done
    if (is_bicop) { 
      res <- outer(i_data[[k]]$data[,1], i_data[[k]]$data[,2],
                   function(u1, u2) predict(object = i_data[[k]],
                                            newdata = cbind(u1, u2),
                                            what = "pdf"))
    } else {
      res <- i_data[[k]]
    }
    # transpose when necessary
    if (i_compl[k] < i) {
      return(t(res))
    } else {
      return(res)
    }
  })
  if (length(i_compl) > 1) {
    return(Reduce("*", i_data))
  } else {
    return(i_data[[1]])
  }
}

get_pred <- function(i, j, k, all_comb, udata, copdata) {
  if (is.integer(k) && length(k) == 0L) {
    res <- udata[,c(i,j)]
  } else {
    ij_combs <- get_ijk_combs(i, j, k, all_comb)
    i_data <- get_i_data(i, ij_combs == 1, all_comb, copdata)
    j_data <- get_i_data(j, ij_combs == 2, all_comb, copdata)
    
    ui <- udata[,i]
    uj <- udata[,j]
    predi <- sapply(1:nrow(udata), function(l) 
      sum(i_data[,l]*ifelse(ui <= ui[l], 1, 0))/sum(i_data[,l]))
    predj <- sapply(1:nrow(udata), function(l) 
      sum(j_data[,l]*ifelse(uj <= uj[l], 1, 0))/sum(j_data[,l]))
    res <- cbind(predi, predj)
  }
  return(res)
}

hoeffCItest <- function(x, y, S, suffStat) {
  ures <- get_pred(x, y, S, suffStat$all_comb, suffStat$udata, suffStat$copdata)
  hoeffd(ures[,1], ures[,2])$P[2]
}

