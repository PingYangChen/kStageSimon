

cumbin_r <- function(x, n, r, nK, rK, p) {
  #print(x)
  dprob <- 1
  for (i in 1:length(x)) {
    dprob <- dprob*dbinom(x[i], n[i], p)
  }
  out <- dprob*pbinom(rK - sum(x), nK, p)
  return(out)
}

pet_recursive_r <- function(n, r, nK, rK, p, K = 2, k = 1, x = NULL) {
  if (K < 2) {
    stop('avaliable for K > 1, if K = 1, use pbinom(r, n, p))')
  } else {
    if (k < K) {
      if (k == 1) {
        bddpar <- 0
        x <- c(r[1] + 1, rep(0, K - 2))
      } else {
        bddpar <- sum(x[1:(k-1)])
      } 
      #print(x)
      lowbdd <- max(c(0, r[k] + 1 - bddpar))
      uppbdd <- min(c(n[k], rK - bddpar))
      #cat(sprintf("%d, %d - %d\n", k, lowbdd, uppbdd))
      #print(c(k, lowbdd, uppbdd))
      out <- 0
      for (v in lowbdd:uppbdd) {
        #cat(sprintf("%d, %d\n", k, v))
        xin <- x
        xin[k] <- v
        out <- out + pet_recursive_r(n, r, nK, rK, p, k = k+1, K = K, x = xin)
      }
      return(out)
    } else {
      return(cumbin_r(x, n, r, nK, rK, p))  
    }
  }
}

kStageP2A_r <- function(p0, p1, nseq, rseq) {
  nstage <- length(nseq)
  if (nstage < 2) {
    stop('Number of stages should be at least 2')
  }
  rej_seq_0 <- numeric(nstage)
  rej_seq_1 <- numeric(nstage)
  for (i in 1:nstage) {
    if (i == 1) {
      rej_seq_0[i] <- pbinom(rseq[1], nseq[1], p0)
      rej_seq_1[i] <- pbinom(rseq[1], nseq[1], p1)
    } else {
      loc <- 1:(i-1)
      rej_seq_0[i] <- pet_recursive_r(n = nseq[loc], r = rseq[loc], nK = nseq[i], rK = rseq[i], p = p0, K = i, k = 1)
      rej_seq_1[i] <- pet_recursive_r(n = nseq[loc], r = rseq[loc], nK = nseq[i], rK = rseq[i], p = p1, K = i, k = 1)  
    }
  }
  pet_seq <- cumsum(rej_seq_0[1:(nstage - 1)])
  rej_drug_0 <- sum(rej_seq_0)
  rej_drug_1 <- sum(rej_seq_1)
  t1e <- 1 - rej_drug_0
  t2e <- rej_drug_1
  en <- sum(c(1, (1 - pet_seq)) * nseq)
  return(list(
    rej_seq_0 = rej_seq_0, rej_seq_1 = rej_seq_1,
    t1e = t1e, t2e = t2e, en = en, pet_seq = pet_seq
  ))
}