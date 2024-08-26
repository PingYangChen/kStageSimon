### ======================================================== ###
### Example Code
### ------------
### Using PSO to generate K-stage clinical trial design
### under Simon's frequentist framework
### ======================================================== ###

# Install and load the "globpso" package by following the 
#  instruction in https://github.com/PingYangChen/globpso
library(globpso)
# Import two local R files "util.R" and "kStageP2A_Objective.R"
source("util.R")
source("kStageP2A_Objective.R")


# Set Requirements of the Clinical Trial
cliRequirement <- list(
  p0 = 0.2,    # response rate in the null hypothesis
  p1 = 0.4,    # response rate in the alternative hypothesis
  alpha = 0.1, # upper bound of type I error
  beta = 0.1   # upper bound of type II error
)

# Set the required number of stages
nStage <- 3

# Set Constraints of the Clinical Trial Design
n1Min <- 10 # minimal sample size at the first stage
nrMin <- 1  # minimal sample size at stages after the 1st one
nMaxRange <- c(30, 70) # the range of the total sample size
# Generate the constraint vector of minimal sample sizes at each stage
# This is the input of the objective function
nMinEachInterim <- c(n1Min, rep(nrMin, nStage - 1)) # (do not change)

# Set the PSO configuration
algSetting <- getPSOInfo(
  nSwarm = 128,     # swarm size 
  maxIter = 200,    # number of iterations
  psoType = "basic" # PSO type (one can use "basic" or "quantum" or "cso")
)

# Set Random seed for reproducibility
pso_seed <- 1

# Set the lower and upper bounds for PSO search
#  particle = (nMax, nPolarized, rProportionEachInterim) 
#  with sizes (1, nStage - 1, nStage)
upper <- c(nMaxRange[2], rep(0.5*pi, nStage - 1), rep(1, nStage))
lower <- c(nMaxRange[1], rep(0.0*pi, nStage - 1), rep(0, nStage))

### -------------------------------------------------------- ###
### Find Optimal Design 
### -------------------------------------------------------- ###
# Run PSO for optimal design
optimRes  <- globpso(objFunc = kStageOptimObj, PSO_INFO = algSetting, 
                     lower = lower, upper = upper, 
                     seed = pso_seed, verbose = TRUE,
                     nMin = nMinEachInterim, cliRequirement = cliRequirement)

# View the optimal design search results
optimRes$val     # Objective function value
optimRes$cputime # computing time
optimRes$par



nEachInterim <- c(5, 10, 20, 35, 50)
rseq <- as.integer(nEachInterim/2)

system.time({
  for (i in 1:100) {
    pet_recursive(1, rep(0, 1), nEachInterim[1:1], rseq[1:1], nEachInterim[2], rseq[2], 0.2, 2)  
  }
})[3]
system.time({
  for (i in 1:100) {
    pet_recursive(1, rep(0, 2), nEachInterim[1:2], rseq[1:2], nEachInterim[3], rseq[3], 0.2, 3)
  }
})[3]
system.time({
  for (i in 1:100) {
    pet_recursive(1, rep(0, 3), nEachInterim[1:3], rseq[1:3], nEachInterim[4], rseq[4], 0.2, 4)
  }
})[3]
system.time({
  for (i in 1:100) {
    pet_recursive(1, rep(0, 4), nEachInterim[1:4], rseq[1:4], nEachInterim[5], rseq[5], 0.2, 5)  
  }
})[3]





#
# Transform the PSO outcome into the readable optimal design
# optimDesign <- kStageFreqCrit(
#   nPolarized = optimRes$par[2:nStage],  
#   rProportion = optimRes$par[(nStage + 1):length(optimRes$par)], 
#   nMax = optimRes$par[1], nMin = nMinEachInterim, cliRequirement)
# 
# # The resulting optimal design
# optimDesign$nseq # Sample sizes at each stage (n_1, ..., n_K)
# optimDesign$rseq # Stopping cutoff sizes at each stage (r_1, ..., r_K)
# 
# # Properties of the resulting optimal design
# optimDesign$t1e # Type I error
# optimDesign$t2e # Type II error
# optimDesign$en  # Expected sample size under null hypothesis
# optimDesign$pet_seq # The probabilities of early termination of the trial at each stage
