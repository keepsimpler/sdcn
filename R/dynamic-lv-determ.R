
#' @title Lotka-Volterra (LV) Equations of Holling type II for mutualistic communities provided by Bastolla et al.
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param parms, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{the competitive matrix inside plants and animals}
#'   \item{M}{the mutualistic matrix between plants and animals}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
#' @return the derivation
#' @import deSolve
model_lv2 <- function(time, init, parms, ...) {
  r = parms[[1]]  # intrinsic growth rates
  C = parms[[2]]  # the competitive matrix
  M = parms[[3]]  # the mutualistic matrix
  h = parms[[4]]  # handling time
  N = init  # initial state
  dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
  list(c(dN))
}

#' @title parmaters for mutualistic LV2 model according to the network and the coefficients
#' @param graph, the interaction topology of mutualistic communities, which is the incidence matrix of a bipartite network
#' @param coeff, a list of coefficients:
#' \describe{
#'    \item{alpha.mu, alpha.sd}{coefficients of the intrinsic growth rates of species}
#'    \item{beta0.mu, beta0.sd}{the intra-species competition coefficients which determin a uniform distribution in [beta0.mu - beta0.sd, beta0.mu + beta0.sd]}
#'    \item{beta1.mu, beta1.sd}{the inter-species competition coefficients}
#'    \item{gamma.mu, gamma.sd}{the inter-species mutualism coefficients}
#'    \item{delta}{trade-off coefficients of mutualistic interaction strengths}
#'    \item{h.mu, h.sd}{coefficients of the handling time of species}
#' }
#' @return a list of parameters for ode model:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{the competitive matrix inside plants and animals}
#'   \item{M}{the mutualistic matrix between plants and animals}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
parms_lv2 <- function(graph, coeff) {
  p <- dim(graph)[1]  # number of Plants
  a <- dim(graph)[2]  # number of Animals
  s <- p + a          # number of all Species  
  with(as.list(coeff),{
    r = runif(s) * 2 * alpha.sd + (alpha.mu - alpha.sd)
    C = matrix(0, nrow = s, ncol = s) # initialize a zero matrix
    C[1:p, 1:p] = runif(p * p) * 2 * beta1.sd + (beta1.mu - beta1.sd) # assign competitive interactions among plants
    C[(p + 1):s, (p + 1):s] = runif(a * a) * 2 * beta1.sd + (beta1.mu - beta1.sd) # assign competitive interactions among animals
    diag(C) = runif(s) * 2 * beta0.sd + (beta0.mu - beta0.sd) # assign intra-species competitive interactions
    
    edges = sum(graph > 0)  # number of interactions (edges)
    M = inc_to_adj(graph)  # transfer to adjacency matrix
    degrees = rowSums(M)  # 
    M[M > 0] = runif(2 * edges) * 2 * gamma.sd + (gamma.mu - gamma.sd)  # assign inter-species mutualistic interactions
    ## !leak!
    M = M / degrees^delta  # trade-off of mutualistic strengths
    
    h = runif(s) * 2 * h.sd + (h.mu - h.sd)
    list(r = r, C = C, M = M, h = h)     
  })
}

#' @title initial values of state variables, i.e., abundances of species
#' @description Assign initial values according to two criteria: 1. using the equilibrium values of LV1 model as initial values. 2. If any of the initial values is less than 0, using the intrinsic growth rates as initial values.
#' @param parms, the parameters assigned to LV2 model
init_lv2 <- function(parms) {
  init = solve(parms$C - parms$M) %*% parms$r
  if (any(init < 0)) {
    warning('Initial state values is less than 0 !!')
    #stop('Initial state values is less than 0 !!', init(LV2))
    init = parms$r
  }
  init  
}

#' @title Simulate ODE dynamics of autonomous systems.The dynamic starts at initialized state variables, and ends in equilibrium (or error where some values of state variables approach infinity?) 
#' @param model model of ODE dynamics
#' @param parms parameters assigned to the model
#' @param init initial values of the model according to the parameters
#' @param steps steps of simulation
#' @param stepwise step length
#' @param extinct_threshold abundance threshold, species with abundance less than that is considered to be exintct 
#' @return a list of:
#' 
sim_ode_auto <- function(model, parms, init, steps = 1000, stepwise = 1, extinct_threshold) {
  times = seq(from = 0, to = steps * stepwise, by = stepwise)
  ode.out = ode(init, times, model, parms)
  nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)])
  nstar[nstar < extinct_threshold] = 0  # species with abundance less than the threshold is considered to be extinct
  extinct = length(nstar) - sum(nstar > 0) 
  survived = sum(nstar > 0)
  Phi = jacobian.full(y = nstar, func = model, parms = parms)
  ret = list(out = ode.out, nstar = nstar, Phi = Phi, model = model, parms = parms, extinct = extinct, survived = survived)
  ret
}

#' @title Simulate ODE dynamics of non-autonomous systems.A example is ecosystems under "press" perturbations. The dynamic is iteration of successive ODE dynamics of automous sytems (\code{\link{sim_ode_auto}}), while at each iterating step, the parameters and/or state values of systems are changed to reflect "press" perturbations.
#' @inheritParams sim_ode_auto
#' @param perturb a function that change the parameters and state values after each iteration step
#' @param iter_steps iteration steps
#' @param isout if output the transiting trajectory of each ODE iterate step
#' @param ... any arguments which are transfered to perturbation function
sim_ode_press <- function(model, parms, init, steps = 1000, stepwise = 1, extinct_threshold, perturb, iter_steps = 10, isout = TRUE, ...) {
  times = seq(from = 0, to = steps * stepwise, by = stepwise)
  ode.outs = list()
  for(i in 1:iter_steps) {
    print(i)
    ode.out = ode(init, times, model, parms) 
    nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) # species biomass at equilibrium
    nstar[nstar < extinct_threshold] = 0  # species with biomass less than extinct threshold is considered to be extinct
    extinct.species = which(nstar == 0)  # extinct species
    
    Phi = jacobian.full(y = nstar, func = model, parms = parms) # community matrix, Jacobian matrix at equilibrium
    if (isout) {
      ret = list(out = ode.out, nstar = nstar, Phi = Phi, params = parms, extinct.species = extinct.species)
    }
    else {
      ret = list(nstar = nstar, Phi = Phi, params = parms, extinct.species = extinct.species)
    }
    ode.outs[[length(ode.outs) + 1]] = ret
    
    #     if (length(extinct.species) > 0) {
    #       ret = remove.species(parms, nstar, extinct.species)
    #       parms = ret$parms
    #       nstar = ret$nstar
    #     }
    #     if (length(nstar) == 0) break  # if all species are extinct, then stop and quit
    
    perturb.res = perturb(parms, nstar, ...)
    parms = perturb.res$parms
    init = perturb.res$nstar
  }
  ode.outs
}

#' @title perturbation that effect on species by increasing/decreasing the intrinsic growth rates of species
#' @param parms parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param r.delta difference of intrinsic growth rates at each iterating step
perturb_growthrate <- function(parms, nstar, r.delta = 0.01) {
  parms$r = parms$r - r.delta
  list(parms = parms, nstar = nstar)
}



