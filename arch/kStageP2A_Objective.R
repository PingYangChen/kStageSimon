library(Rcpp)
sourceCpp('kstagecli/kStageP2A.cpp')

### ---------------------------------------------------------------------------
### The objective function of the K-stage minimax design
### ---------------------------------------------------------------------------
kStageMinMaxObj <- function(particle, nMin, cliRequirement = NULL) {
  # particle      : vector of the form (nMax, nPolarized, rProportion) with sizes (1, nStage - 1, nStage)
  # nMin          : integer, minimal sample size at each stage
  # cliRequirement: list, requirements of the trial including p0, p1, alpha and beta
  if (is.null(cliRequirement)) {
    p0 <- 0.2
    p1 <- 0.4
    t1eThres <- 0.1
    t2eThres <- 0.1
  } else {
    p0 <- cliRequirement$p0
    p1 <- cliRequirement$p1
    t1eThres <- cliRequirement$alpha
    t2eThres <- cliRequirement$beta
  }
  
  nStage <- as.integer(length(particle)/2)
  nMax <- particle[1]
  nPolarized <- particle[2:nStage]
  rProportion <- particle[(nStage + 1):length(particle)]

  result <- kStageFreqCrit(nPolarized, rProportion, nMax, nMin, cliRequirement)
  # Return the objective function value under null hypothesis 
  if ((result$t1e <= t1eThres) & (result$t2e <= t2eThres)) {
    return(result$en/nMax + nMax)
  } else {
    return(1e8)
  }
}

### ---------------------------------------------------------------------------
### The objective function of the K-stage optimal design
### ---------------------------------------------------------------------------
kStageOptimObj <- function(particle, nMin, cliRequirement = NULL) {
  # particle      : vector of the form (nMax, nPolarized, rProportion) with sizes (1, nStage - 1, nStage)
  # nMin          : integer, minimal sample size at each stage
  # cliRequirement: list, requirements of the trial including p0, p1, alpha and beta
  if (is.null(cliRequirement)) {
    p0 <- 0.2
    p1 <- 0.4
    t1eThres <- 0.1
    t2eThres <- 0.1
  } else {
    p0 <- cliRequirement$p0
    p1 <- cliRequirement$p1
    t1eThres <- cliRequirement$alpha
    t2eThres <- cliRequirement$beta
  }
  
  nStage <- as.integer(length(particle)/2)
  nMax <- particle[1]
  nPolarized <- particle[2:nStage]
  rProportion <- particle[(nStage + 1):length(particle)]
  #
  result <- kStageFreqCrit(nPolarized, rProportion, nMax, nMin, cliRequirement)
  # Return the expected sample size under null hypothesis 
  if ((result$t1e <= t1eThres) & (result$t2e <= t2eThres)) {
    return(result$en)
  } else {
    return(1e8)
  }
}

### ---------------------------------------------------------------------------
### The function computing the expected sample size under null hypothesis for 
###  the input PSO particle 
### ---------------------------------------------------------------------------
kStageFreqCrit <- function(nPolarized, rProportion, nMax, nMin, cliRequirement = NULL) {

  if (is.null(cliRequirement)) {
    p0 <- 0.2
    p1 <- 0.4
    t1eThres <- 0.1
    t2eThres <- 0.1
  } else {
    p0 <- cliRequirement$p0
    p1 <- cliRequirement$p1
    t1eThres <- cliRequirement$alpha
    t2eThres <- cliRequirement$beta
  }
  
  nStage <- length(nPolarized) + 1
  nseq <- get_cohort(nPolarized, nStage, nMax, nMin)
  rseq <- get_cutoff(rProportion, nseq)
  #
  return(kStageP2A_Cpp(p0, p1, nseq, rseq))
}