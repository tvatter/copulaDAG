# factory <- function(fun)
#   function(...) {
#     warn <- err <- NULL
#     res <- withCallingHandlers(
#       tryCatch(fun(...), error=function(e) {
#         err <<- conditionMessage(e)
#         NULL
#       }), warning=function(w) {
#         warn <<- append(warn, conditionMessage(w))
#         invokeRestart("muffleWarning")
#       })
#     list(res, warn=warn, err=err)
#   }
# 
# optim_better <- factory(optim)
# 
# computeCausOrder <- function(G)
#   # Copyright (c) 2013  Jonas Peters  [peters@stat.math.ethz.ch]
#   # All rights reserved.
# {
#   p <- dim(G)[2]
#   remaining <- 1:p
#   causOrder <- rep(NA,p)
#   for(i in 1:(p-1))
#   {
#     root <- min(which(colSums(G) == 0))
#     causOrder[i] <- remaining[root]
#     remaining <- remaining[-root]
#     G <- G[-root,-root]
#   }
#   causOrder[p] <- remaining[1]    
#   return(causOrder)
# }