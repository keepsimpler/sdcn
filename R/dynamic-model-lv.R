
#' @title Lotka-Volterra (LV) Equations of Holling type II for a community mixed by Competition and Antagonistic interactions
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{A}{a matrix of antagonistic interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#'   \item{g}{the conversion efficiency of antagonistic interactions}
#' }
#' @return the derivation
#' @import deSolve
model_lv2_ca <- function(time, init, parms, ...) {
  r = parms[[1]]  # intrinsic growth rates
  C = parms[[2]]  # the competitive matrix
  A = parms[[3]]  # the antagonistic matrix
  h = parms[[4]]  # handling time
  g = parms[[5]]  # conversion efficiency
  N = init  # initial state
  # seperate the antagonistic matrix to two sub-matrice :
  # First sub-matrix describes the positive interactions among species
  # Second sub-matrix describes the nagative interactions among species
  AP = A[A > 0]
  AN = abs(A[A < 0])
  H = matrix(h, nrow = length(h), ncol = length(h), byrow = T)
  # For the positive sub-matrix, 
  dN <- N * ( r - C %*% N + g * (AP %*% N) / (1 + h * AP %*% N) - (AN %*% N) / (1 + (H * AN) %*% N) )
  list(c(dN))
}

#' @title generate the antagonistic interaction matrix according to the antagonistic network and coefficients
#' @param graph the antagonistic interaction topology of communities, which is a bi-directional network
#' @param gamma.mu
#' @param gamma.sd coefficients that determin a uniform distribution of antagonistic interaction strengths
parms_antago_interactions <- function(graph, gamma.mu, gamma.sd) {
  stopifnot(dim(graph)[1] == dim(graph)[2]) # if not a adjacency matrix, stop
  stopifnot(gamma.mu >= 0, gamma.sd >= 0, gamma.mu >= gamma.sd)
  diag(graph) <- 0  # ensure the diagonal elements of antagonistic matrix equal 0
  edges = sum(graph > 0)  # !leak! number of interactions, number of positive and nagative interactions are same with each other
  foodweb[lower.tri(foodweb) & foodweb != 0] = foodweb[lower.tri(foodweb) & foodweb != 0] * runif(edges, min = gamma.mu - gamma.sd, max = gamma.mu + gamma.sd) # lower tringle of matrix and not equal zero are assigned random variable values
  foodweb[upper.tri(foodweb)] = - t(foodweb)[upper.tri(t(foodweb))] # upper tringle of matrix are nagative of lower tringle
}

#' @title parmaters for antagonism LV2 model according to the network and the coefficients
#' @param antago_graph the antagonistic interaction topology of communities, which is the adjacency matrix of a network
#' @param competitive_graph the competitive interaction topology of communities
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{alpha.mu, alpha.sd}{coefficients of the intrinsic growth rates of species}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta.mu - beta.sd, beta.mu + beta.sd]}
#'    \item{beta1.mu, beta1.sd}{the inter-species competition coefficients}
#'    \item{gamma.mu, gamma.sd}{the inter-species antagonism coefficients}
#'    \item{g.mu, g.sd}{conversion efficiency of antagonistic interactions}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#' }
#' @return a list of parameters for ode model:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{A}{a matrix of antagonistic interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#'   \item{g}{the conversion efficiency of antagonistic interactions}
#' }
parms_lv2_ca <- function(antago_graph, competitive_graph, coeff) {
  stopifnot(dim(antago_graph)[1] == dim(antago_graph)[2], dim(competitive_graph)[1] == dim(competitive_graph)[2], dim(antago_graph)[1] == dim(competitive_graph)[1])
  s = dim(antago_graph)[1]
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- parms_competitive_interactions(competitive_graph, beta0.mu, beta0.sd, beta1.mu, beta1.sd) # generate a competitive interaction matrix
    A <- parms_antago_interactions(antago_graph, gamma.mu, gamma.sd)
    h = runif2(s, h.mu, h.sd)
    g = runif2(s, g.mu, g.sd)
    list(r = r, C = C, A = A, h = h, g = g)     
  })
}


#' @title Lotka-Volterra (LV) Equations of Holling type II for a community mixed by Competition and Mutualism interactions
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{M}{a matrix of mutualism interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
#' @return the derivation
#' @import deSolve
model_lv2_cm <- function(time, init, parms, ...) {
  r = parms[[1]]  # intrinsic growth rates
  C = parms[[2]]  # the competitive matrix
  M = parms[[3]]  # the mutualistic matrix
  h = parms[[4]]  # handling time
  N = init  # initial state
  dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
  list(c(dN))
}

#' @title generate the mutualistic interaction matrix according to the mutualistic network and coefficients
#' @param graph the mutualistic interaction topology of communities, which is the adjacency matrix of a network
#' @param gamma.mu
#' @param gamma.sd coefficients that determin a uniform distribution of mutualistic interaction strengths
#' @param delta coefficient that determin the trade-off between the interaction strength and width(node degree) of species
parms_mutual_interactions <- function(graph, gamma.mu, gamma.sd, delta) {
  stopifnot(dim(graph)[1] == dim(graph)[2]) # if not a adjacency matrix, stop
  stopifnot(gamma.mu >= 0, gamma.sd >= 0, gamma.mu >= gamma.sd, delta >= 0)
  diag(graph) <- 0  # ensure the diagonal elements of mutualistic matrix equal 0
  edges = sum(graph > 0)  # number of all interactions
  M = graph
  degrees = rowSums(M)  # the width (node degree) of species
  M[M > 0] = runif2(edges, gamma.mu, gamma.sd)  # assign inter-species mutualistic interaction strengths
  old_total_strength = sum(M)
  ## !leak!
  M = M / degrees^delta  # trade-off of mutualistic strengths
  new_total_strength = sum(M) # ensure the total strength constant before and after trade-off
  M = M * old_total_strength / new_total_strength
  M = t(M)  ## !leak! ??
  M
}

#' @title generate the competitive interaction matrix according to the competitive network and coefficients
#' @param graph the competitive interaction topology of communities, which is the adjacency matrix of a network
#' @param beta0.mu
#' @param beta0.sd coefficients that determin a uniform distribution of intra-species interaction strengths
#' @param beta1.mu
#' @param beta1.sd coefficients that determin a uniform distribution of inter-species interaction strengths
parms_competitive_interactions <- function(graph, beta0.mu, beta0.sd, beta1.mu, beta1.sd) {
  stopifnot(dim(graph)[1] == dim(graph)[2]) # if not a adjacency matrix, stop
  stopifnot(beta0.mu > 0, beta0.sd >= 0, beta0.mu >= beta0.sd, beta1.mu >= 0, beta1.sd >= 0, beta1.mu >= beta1.sd)
  diag(graph) <- 0
  edges = sum(graph > 0)  # number of all interactions
  s <- dim(graph)[1] # number of total Species
  C <- graph
  C[C > 0] = runif2(edges, beta1.mu, beta1.sd) # assign inter-species competitive interaction strengths
  diag(C) <- runif2(s, beta0.mu, beta0.sd) # assign intra-species competitive interaction strengths
  C
}

#' @title parmaters for mutualism LV2 model according to the network and the coefficients
#' @param mutual_graph the mutualistic interaction topology of communities, which is the adjacency matrix of a network
#' @param competitive_graph the competitive interaction topology of communities
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{alpha.mu, alpha.sd}{coefficients of the intrinsic growth rates of species}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta.mu - beta.sd, beta.mu + beta.sd]}
#'    \item{beta1.mu, beta1.sd}{the inter-species competition coefficients}
#'    \item{gamma.mu, gamma.sd}{the inter-species mutualism coefficients}
#'    \item{delta}{trade-off coefficients of mutualistic interaction strengths}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#' }
#' @return a list of parameters for ode model:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a competitive interaction matrix }
#'   \item{M}{a mutualistic interaction matrix among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
parms_lv2_cm <- function(mutual_graph, competitive_graph, coeff) {
  stopifnot(dim(mutual_graph)[1] == dim(mutual_graph)[2], dim(competitive_graph)[1] == dim(competitive_graph)[2], dim(mutual_graph)[1] == dim(competitive_graph)[1])
  s = dim(mutual_graph)[1]
  with(as.list(coeff), {
    r <- runif2(s, alpha.mu, alpha.sd)
    C <- parms_competitive_interactions(competitive_graph, beta0.mu, beta0.sd, beta1.mu, beta1.sd) # generate a competitive interaction matrix
    M <- parms_mutual_interactions(mutual_graph, gamma.mu, gamma.sd, delta)
    h = runif2(s, h.mu, h.sd)
    list(r = r, C = C, M = M, h = h)     
  })
}

#' @title initial values of state variables, i.e., abundances of species
#' @description Assign initial values according to two criteria: 1. using the equilibrium values of LV1 model as initial values. 2. If any of the initial values is less than 0, using the intrinsic growth rates as initial values.
#' @param parms, the parameters assigned to LV2 model
init_lv2_cm <- function(parms) {
  init = solve(parms$C - parms$M) %*% parms$r
  if (any(init < 0)) {
    warning('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    init = parms$r
  }
  init  
}


