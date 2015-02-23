
#' @title transfer an incidence matrix to an adjacency matrix
#' @param inc, an incidence matrix
#' @return adj, an adiacency matrix
inc_to_adj <- function(inc){
  p <- dim(inc)[1]  # number of Plants
  a <- dim(inc)[2]  # number of Animals
  s <- p + a  # number of all Species
  adj <- matrix(0, s, s)  # initialize the adjacency matrix as a zero-matrix
  adj[1:p, (p + 1):s] <- inc  # the upper right sub-matrix is the incidence matrix
  adj <- adj + t(adj)  # the lower left sub-matrix is transpose of the incidence matrix
  return(adj)
}