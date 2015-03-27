
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


###############################################################################
#' @title Generate a connected graph using package [igraph]
#'
#' @param s, size of network. 
#' if graph type is bipartite, s[1], s[2] represent size of two groups; else s is size of network
#' @param k, average degree for the network.
#' 1 < k < s for unipartite network, 1 < k < s[1]*s[2]/(s[1]+s[2]) for bipartite network.
#' @param gtype, Graph type generated: 'bipartite', 'sf', 'er', 'regular'.
#' @param maxtried, the maximum number of tried times. 
#' If have tried [maxtried] times, the function will return no matter whether the connected graph is generated.
#' @param expower exponent coefficient of Scale-Free network
#' @param ... the parms conform to the implementation functions of [igraph]
#' @return the connected graph
#' @details .  
#' @import igraph
graph_connected <- function(s, k, gtype, maxtried = 100, expower = 2.5, ...) {
  #library(igraph)
  if (gtype == 'bipartite' && is.na(s[2])) {  # the bipartite graph need size of two groups of nodes
    warning('sizes of TWO groups of nodes should be designated. 
            we have assumed the size of second group equal to the size of first group.')
    s[2] = s[1]  # if missed second size, we assume it equal to the first size.
  }
  count = 0
  repeat {  # generate a connected graph
    if (gtype == 'bipartite') {
      G = bipartite.random.game(s[1], s[2], type = 'gnm', m = ceiling(k * (s[1] + s[2])))
    } else if (gtype == 'sf') {
      G = static.power.law.game(s, k * s, exponent.out = expower)
    }
    else if (gtype == 'er') {
      G = erdos.renyi.game(s, p.or.m = k * s, type = 'gnm')
    }
    else if (gtype == 'regular') {
      G = k.regular.game(s, k)
    }
    else if (gtype == 'complete') {
      G = graph.full(s)
    }
    if (igraph::is.connected(G)) break  # until a connected graph is generated
    count = count + 1
    if (count == maxtried) {
      warning(paste('Tried', maxtried, 'times, But connected graph still cannot be generated.'))
      break
    }
  }
  G
}
#plot(G, layout = layout.bipartite)

gen_foodweb <- function(s, k, type = 'full') {
  if (type == 'full') {
    g = graph.full(s)
    foodweb = as.matrix(get.adjacency(g))
    for (i in 2:s) {
      for (j in 1:(i-1)) {
        if (foodweb[i, j] == 1) { # if an undirected edge exists between i and j
          if (runif(1) < 0.5)
            foodweb[i, j] = -1
          else
            foodweb[j, i] = -1
        }
      }
    }
  }
  foodweb
}



#' @title another form of uniform distribution between [mean - sd, mean + sd]
runif2 <- function(n, mean, sd) {
  runif(n) * 2 * sd + (mean - sd)
}

display_output <- function(sim.out, coeff, fragility) {
  ode.nstars = laply(sim.out, function(one) {
    one$nstar
  })
  # fragility = fragility(sim.out)
  
  fig_title <- paste('h.mu = ', coeff[['h.mu']], ', gamma.mu = ', coeff[['gamma.mu']], ', delta = ', coeff[['delta']], ', entropy = ', sprintf("%.2f", fragility['entropy']), ', variance = ', fragility['variance'], ', resist = ', fragility['resistance'], sep = '')
  matplot(ode.nstars, main = fig_title, type = 'l', lwd = 1.5, xlab = 'steps of gradual presses', ylab = 'Abundances of species')  
}