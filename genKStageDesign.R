library(globpso)
### Import two local R files "util.R" and "kStageP2A_Objective.R"
setwd(here::here())
source("util.R")
source("kStageP2A_Objective.R")

#' Using Particle Swarm Optimization Algorithms for Finding Multi-stage Clinical Trial Design
#'
#' @param designType The type of the clinical trial design. 
#' Input "optimal" for finding the optimal design. Input "minimax" for finding the minimax design.
#' @param nStage The number of stages of the clinical trial design.
#' @param cliRequirement The list of required trial information. 
#' For example, `list(p0 = 0.2, p1 = 0.4, alpha = 0.1, beta = 0.1)`.
#' `p0`: response rate in the null hypothesis
#' `p1`: response rate in the alternative hypothesis
#' `alpha`: upper bound of type I error
#' `beta`: upper bound of type II error
#' @param OmegaL The lower bound of the total sample size.
#' @param OmegaU The upper bound of the total sample size. 
#' @param m1 The minimal sample size at the first stage. The default value is `10`.
#' @param mk The minimal sample size at stages after the first stage. The default value is `1`.
#' @param qk The minimum incremental critical value at each stage. The default value is `0`.
#' @param psoSetting The list of PSO configurations. For example, `list(nSwarm = 32, maxIter = 100)`.
#' `nSwarm`: The swarm size.
#' `maxIter`: The number of iterations.
#' @param seed Set random seed for reproducibility. The default is `NULL`.
#' @param verbose Set `TRUE` for printing the PSO updating progress and the final results.
#' @return An List.
#' \describe{
#' \item{$result$nseq}{ The sample size enrolled at the kth stage (n_1, ..., n_K).}
#' \item{$result$rseq}{ The incremental critical value at the kth stage (r_1, ..., r_K).}
#' \item{$result$t1e}{ Type I error of the resulting design.}
#' \item{$result$t2e}{ Type II error of the resulting design.}
#' \item{$result$en}{ Expected sample size under null hypothesis.}
#' \item{$result$pet_seq}{ The probabilities of early termination of the trial at each stage.}
#' \item{$elapse}{ Computing time.}
#' }
genKStageDesign <- function(designType = "optimal", nStage = 3, 
                            cliRequirement = list(p0 = 0.2, p1 = 0.4, alpha = 0.1, beta = 0.1),
                            OmegaL = 30, OmegaU = 70, m1 = 10, mk = 1, qk = 0,
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
  
  ### Performing PSO Algorithm to find the clinical trial designs
  if (designType == "optimal") {
    
    ### Run PSO for finding optimal design
    res <- globpso(objFunc = kStageOptimObj, PSO_INFO = algSetting, 
                   lower = lower, upper = upper, seed = seed, verbose = verbose,
                   nMin = nMinEachInterim, rMin = qk, cliRequirement = cliRequirement)  

  } else if (designType == "minimax") {
    
    ### Run PSO for finding minimax design
    res <- globpso(objFunc = kStageMinMaxObj, PSO_INFO = algSetting, 
                   lower = lower, upper = upper, seed = seed, verbose = verbose,
                   nMin = nMinEachInterim, rMin = qk, cliRequirement = cliRequirement)
  
  }
  
  ### Transform the PSO outcome into the readable design format
  resDesign <- kStageFreqCrit(
    nPolarized = res$par[2:nStage],  
    rProportion = res$par[(nStage + 1):length(res$par)], 
    nMax = res$par[1], nMin = nMinEachInterim, rMin = qk, cliRequirement = cliRequirement)
  

  if (verbose) {
    cat(
      sprintf("---- %s Design Problem Setting -----------", ifelse(designType == 'optimal', "Optimal", "Minimax")), '\n',
      sprintf("Trial Requirements:      (p0, p1) = (%.3f, %.3f)", cliRequirement$p0, cliRequirement$p1), '\n',
      sprintf("Error Constraints : (alpha, beta) = (%.3f, %.3f)", cliRequirement$alpha, cliRequirement$beta), '\n', '\n'
    )
    cat(
      sprintf("---- Resulting %s Design -----------------", ifelse(designType == 'optimal', "Optimal", "Minimax")), '\n',
      "rk/nk:        ",
      paste0(sprintf("%d/%d", resDesign$rseq, resDesign$nseq), collapse = ", "), '\n',
      sprintf("Type I Error : t1e = %.4f", resDesign$t1e), '\n',
      sprintf("Type II Error: t2e = %.4f", resDesign$t2e), '\n',
      sprintf("Expected Sample Size: E(N|H0) = %.4f", resDesign$en), '\n',
      sprintf("Maximum Sample Size :       n = %d", resDesign$nseq[nStage]), '\n',
      "Probabilites of Early Termination \n",
      paste0(sprintf("- Stage %d: %.4f", 1:(nStage - 1), resDesign$pet_seq), collapse = "\n "), '\n', '\n'
    )
    cat(
      "---- PSO Performance --------------------------", "\n",
      sprintf("Swarm Size: %d", algSetting$nSwarm), "\n",
      sprintf("Iterations: %d", algSetting$maxIter), "\n",
      sprintf("Elapse Time: %.2f seconds", res$cputime), '\n', '\n'
    )
  }
  
  cat(
    "Output values", "\n", "-------------", "\n",
    "$result$nseq   :", "The sample size enrolled at the kth stage (n_1, ..., n_K).", "\n",
    "$result$rseq   :", "The incremental critical value at the kth stage (r_1, ..., r_K).", "\n",
    "$result$t1e    :", "Type I error of the resulting design.", "\n",
    "$result$t2e    :", "Type II error of the resulting design.", "\n",
    "$result$en     :",  "Expected sample size under null hypothesis", "\n",
    "$result$pet_seq:", "The probabilities of early termination of the trial at each stage", "\n",
    "$elapse        :", "Computing time.", "\n"
  )
    
  return(list(
    "result" = resDesign,
    "elapse" = res$cputime # computing time
  ))
  
}