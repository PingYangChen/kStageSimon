library(globpso)
library(here)
setwd(here::here())
source("util.R")
source("kStageP2A_Objective.R")

genKStageDesign <- function(designType = "optimal", nStage = 3, 
                            designRequirement = list(p0 = 0.2, p1 = 0.4, alpha = 0.1, beta = 0.1),
                            OmegaL = 70, OmegaU = 30, m1 = 10, mk = 1, qk = 0,
                            psoSetting = list(nSwarm = 32, maxIter = 100), 
                            seed = NULL, verbose = TRUE) {
  
  ### The range of the total sample size
  nMaxRange <- c(OmegaL, OmegaU)
  ### Generate the constraint vector of minimal sample sizes at each stage
  ### This is the input of the objective function
  nMinEachInterim <- c(m1, rep(mk, nStage - 1)) # (do not change!!!)
  
  ### Set the lower and upper bounds for PSO search
  ###  particle = (nMax, nPolarized, rProportionEachInterim) 
  ###  with sizes (1, nStage - 1, nStage)
  upper <- c(nMaxRange[2], rep(0.5*pi, nStage - 1), rep(1, nStage))
  lower <- c(nMaxRange[1], rep(0.0*pi, nStage - 1), rep(0, nStage))
  
  ### Set the PSO configuration
  algSetting <- getPSOInfo(
    nSwarm = psoSetting$nSwarm,     # swarm size 
    maxIter = psoSetting$maxIter,    # number of iterations
    psoType = "basic" # PSO type (one can use "basic" or "quantum" or "cso")
  )
  
  
  if (designType == "optimal") {
    
    ### Run PSO for finding optimal design
    res <- globpso(objFunc = kStageOptimObj, PSO_INFO = algSetting, 
                   lower = lower, upper = upper, seed = seed, verbose = verbose,
                   nMin = nMinEachInterim, cliRequirement = designRequirement)  

  } else if (designType == "minimax") {
    
    ### Run PSO for finding minimax design
    res <- globpso(objFunc = kStageMinMaxObj, PSO_INFO = algSetting, 
                   lower = lower, upper = upper, seed = seed, verbose = verbose,
                   nMin = nMinEachInterim, cliRequirement = designRequirement)
  
  }
  
  if (verbose)
  ### View the optimal design search results
  res$val     # Objective function value
  res$cputime # computing time
  
  ### Transform the PSO outcome into the readable optimal design
  resDesign <- kStageFreqCrit(
    nPolarized = res$par[2:nStage],  
    rProportion = res$par[(nStage + 1):length(res$par)], 
    nMax = res$par[1], nMin = nMinEachInterim, rMin = 0, cliRequirement = cliRequirement)
  
  ### The resulting optimal design
  resDesign$nseq # Sample sizes at each stage (n_1, ..., n_K)
  resDesign$rseq # Stopping cutoff sizes at each stage (r_1, ..., r_K)
  
  ### Properties of the resulting optimal design
  resDesign$t1e # Type I error
  resDesign$t2e # Type II error
  resDesign$en  # Expected sample size under null hypothesis
  resDesign$pet_seq # The probabilities of early termination of the trial at each stage
  
  
}