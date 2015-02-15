###############################################################################
#' @title Swapping links Algorithm for null model of bipartite networks, that generates random network (ensembles) which keep the node degree distribution of a real network.
#' @param bigraph, incidence matrix of a bipartite network, rows and cols represent two groups of nodes/species
#' @param ntry, the possible maximum times of swapping links to try
#' @return an incidence matrix of bipartite network whose links being randomly swapped.
#' @examples 
#' swaplinks(bigraph)
swaplinks <- function(bigraph, ntry = 5000) {
  count1 <- 0
  B <- bigraph
  NumP <- dim(B)[1]
  NumA <- dim(B)[2]
  while (count1 < ntry){
    count1 <- count1 + 1
    ## pick two rows and two columns
    row1 <- sample(1:NumP, 1)
    row2 <- sample(1:NumP, 1)
    col1 <- sample(1:NumA, 1)
    col2 <- sample(1:NumA, 1)
    ## check swappable
    if (B[row1, col1] == 0.0 && B[row1, col2] > 0.0 && B[row2, col1] > 0.0 && B[row2, col2] == 0.0){
      ## swap
      B[row1, col1] <- B[row1, col2]
      B[row1, col2] <- 0.0
      B[row2, col2] <- B[row2, col1]
      B[row2, col1] <- 0.0
    }
    else{
      if (B[row1, col1] > 0.0 && B[row1, col2] == 0.0 && B[row2, col1] == 0.0 && B[row2, col2] > 0.0){
        ## swap
        B[row1, col2] <- B[row1, col1]
        B[row1, col1] <- 0.0
        B[row2, col1] <- B[row2, col2]
        B[row2, col2] <- 0.0
      }
    }
  }
  return(B)
}

