
#' @title Simulate ODE dynamics of autonomous systems.The dynamic starts at initialized state variables, and ends in equilibrium (or error where some values of state variables approach infinity?) 
#' @param model model of ODE dynamics
#' @param parms parameters assigned to the model
#' @param init initial values of the model according to the parameters
#' @param steps steps of simulation
#' @param stepwise step length
#' @param extinct_threshold abundance threshold, species with abundance less than that is considered to be exintct 
#' @return a list of:
#' \describe{
#'   \item{out}{output of one ODE simulation, including the trajectory of values of state variables}
#'   \item{nstar}{the values of state variables in equilibrium}
#'   \item{Phi}{the Jacobian matrix in equilibrium}
#'   \item{model}{model of ODE dynamics}
#'   \item{parms}{parameters assigned to the model}
#'   \item{extinct}{number of extinct species}
#'   \item{survived}{number of survived species}
#' }
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
#' @param iter_steps possiblely maximum iteration steps
#' @param isout if output the transiting trajectory of each ODE iterate step
#' @param ... any arguments which are transfered to perturbation function
#' @return a list of lists :
#' \describe{
#'   \item{out}{output of one ODE simulation, including the trajectory of values of state variables}
#'   \item{nstar}{the values of state variables in equilibrium}
#'   \item{Phi}{the Jacobian matrix in equilibrium}
#'   \item{parms}{parameters assigned to the model}
#'   \item{extinct.species}{a vector of species that extincted}
#' }
sim_ode_press <- function(model, parms, init, steps = 1000, stepwise = 1, extinct_threshold, perturb, iter_steps = 500, isout = TRUE, ...) {
  times = seq(from = 0, to = steps * stepwise, by = stepwise)
  ode.outs = list()
  for(i in 1:iter_steps) {
    print(i)
    ode.out = ode(init, times, model, parms) 
    nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) # species biomass at equilibrium
    nstar[nstar < extinct_threshold] = 0  # species with biomass less than extinct threshold is considered to be extinct
    extinct.species = which(nstar == 0)  # extinct species

    flag = 0
    # if all species are extinct, will end the simulation
    if (length(nstar) == length(extinct.species)) flag = 1
    # if any species' abundance is NaN, that means the ODE dynamic is unstable, the simulation will also be ended
    if (any(is.nan(nstar))) flag = 2

    Phi = jacobian.full(y = nstar, func = model, parms = parms) # community matrix, Jacobian matrix at equilibrium
    if (isout) {
      ret = list(out = ode.out, nstar = nstar, Phi = Phi, params = parms, extinct.species = extinct.species, flag = flag)
    }
    else {
      ret = list(nstar = nstar, Phi = Phi, params = parms, extinct.species = extinct.species, flag = flag)
    }
    ode.outs[[length(ode.outs) + 1]] = ret
    # if all species are extinct, end the simulation
    if (flag == 1 || flag == 2)
      break;
    
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

#' @title perturbations that effect on species by increasing/decreasing the intrinsic growth rates of all species
#' @param parms parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param r.delta deviation of intrinsic growth rates at each iterating step
perturb_growthrate <- function(parms, nstar, r.delta.mu = 0.01, r.delta.sd = 0.01) {
  #set.seed(1)
  parms$r = parms$r - runif(length(parms$r), min = r.delta.mu - r.delta.sd, max = r.delta.mu + r.delta.sd)
  list(parms = parms, nstar = nstar)
}

#' @title perturbations that effect on species by increasing/decreasing the intrinsic growth rates of a part of species
#' @param parms parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param r.delta deviation of intrinsic growth rates at each iterating step
#' @param perturbed_species the index of perturbed species
perturb_growthrate_part <- function(parms, nstar, r.delta.mu = 0.01, r.delta.sd = 0.01,  perturbed_species) {
  
  parms$r[perturbed_species] = parms$r[perturbed_species] - runif(length(perturbed_species), min = r.delta.mu - r.delta.sd, max = r.delta.mu + r.delta.sd)
  list(parms = parms, nstar = nstar)
}

#' @title perturbations that effect on mutualistic interactions by increasing/decreasing strengths of them
#' @param parms parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param gamma.delta deviation of mutualistic interaction strengths at each iterating step
perturb_mutualistic_strength <- function(parms, nstar, gamma.delta.mu = 0.01, gamma.delta.sd = 0.01) {
  edges = length(parms$M[parms$M > 0])
  parms$M[parms$M > 0] = parms$M[parms$M > 0] - runif(edges, min = gamma.delta.mu - gamma.delta.sd, max = gamma.delta.mu + gamma.delta.sd)
  list(parms = parms, nstar = nstar)
}

#' @title perturbations that remove one species
#' @param parms parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param extinct_species the removed species
perturb_primary_extinct <- function(parms, nstar, extinct_species) {
  nstar = nstar[- extinct_species]  # primary extinction
  parms$r = parms$r[- extinct_species]
  parms$C = parms$C[- extinct_species, - extinct_species]
  parms$M = parms$M[- extinct_species, - extinct_species]
  parms$h = parms$h[- extinct_species]
  list(parms = parms, nstar = nstar)
}

#' @title compute fragility of mutualistic communities in gradual pressed context
#' @param sim.out output of simulation under gradual pressed conditions (\code{\link{sim_ode_press}})
#' @return resistance measured by the length of community trajectory
#' @return fragility measured by the variance of community trajectory
fragility <- function(sim.out) {
  trajectory = laply(sim.out, function(one) {
    length(one$extinct.species)
  })
  resistance2 = sum(trajectory) # the complement area of the trajectory
  trajectory = trajectory[-1] - trajectory[-length(trajectory)]
  fragility.variance = sum(trajectory^2)
  trajectory.positive = trajectory[trajectory > 0]
  fragility.entropy = sum(trajectory.positive * log(trajectory.positive))
  list(trajectory = trajectory, variance = fragility.variance, entropy = fragility.entropy, resistance = length(sim.out), resistance2 = resistance2)
}

fragility.abund <- function(sim.out) {
  trajectory.abund <- laply(sim.out, function(one) {
    sum(one$nstar)
  })
  resistance = sum(trajectory.abund)
  list(resistance = resistance)
}
#' @title compute resistance of mutualistic communities in gradual pressed context
#' @param sim.out output of simulation under gradual pressed conditions (\code{\link{sim_ode_press}})
#' @return resistance measured by the length of community trajectory
resistance <- function(sim.out) {
  length(sim.out)
}


# random simulation times

