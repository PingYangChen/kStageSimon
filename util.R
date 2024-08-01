
### ---------------------------------------------------------------------------
### The function transforms the inputting PSO particle's part of 
###  proportions vector to the readable stopping cutoff sizes at each
###  stage of the trial 
### ---------------------------------------------------------------------------
get_cutoff <- function(rProportion, nseq) {
  n_stage <- length(nseq)
  nEachInterim <- c(nseq[1], diff(nseq))
  rEachInterim <- rProportion*nEachInterim
  cutoffs <- round(cumsum(rEachInterim))
  return(cutoffs)
}

### ---------------------------------------------------------------------------
### The function transforms the inputting PSO particle's part of 
###  polarized sample sizes to the readable sample sizes at each
###  stage of the trial 
### ---------------------------------------------------------------------------
get_cohort <- function(w_polarized, n_stage, max_n = 50, min_n = 0) {
  # w_polarized = nPolarized; n_stage = nStage; max_n = nMax; min_n = nMin
  wcumsin <- wcos <- numeric(n_stage)
  wcumsin[1] <- 1; wcumsin[2:n_stage] <- cumprod(sin(w_polarized))
  wcos[1:(n_stage-1)] <- cos(w_polarized); wcos[n_stage] <- 1
  wt <-	(wcumsin*wcos)^2
  if (length(min_n) == 1) {
    assigned_sample <- max_n - min_n*n_stage
    fixed_sample <- rep(min_n, n_stage)
  } else {
    assigned_sample <- max_n - sum(min_n)
    fixed_sample <- min_n
  }
  if (assigned_sample < 0) { stop("min_n is too large") }
  nobs.seq = round(cumsum(wt*assigned_sample + fixed_sample))
  return(nobs.seq)
}

### ---------------------------------------------------------------------------
### The function transforms the readable sample sizes at each stage 
###  of the trial to the inputting PSO particle's part of polarized
###  sample sizes
### ---------------------------------------------------------------------------
cohort2polar <- function(nobs.seq, min_n = 0) {
  n <- length(nobs.seq)
  xx <- nobs.seq
  n_interims <- c(nobs.seq[1], diff(xx)) - min_n
  x <- sqrt(n_interims/max(xx))
  ang <- numeric(n-1)
  for (i in 1:(n-1)) {
    ang[i] <- acos(x[i]/sqrt(sum((x[i:n]^2))))
    if (i < (n-1)) {
      ang[i] <- pi - ang[i]
    }
  }
  return(ang)
}
