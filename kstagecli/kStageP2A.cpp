
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

void vecPrintf(const arma::vec &v)
{
  for (arma::uword i = 0; i < v.n_elem; i++) Rprintf("%4.4f\t", v(i));
  Rprintf("\n\n");
}

// [[Rcpp::export]]
double cumbin(const arma::ivec x, const arma::ivec n, const arma::ivec r, 
                 const int nK, const int rK, const double p) {
  arma::uword ns = x.n_elem; 
  //vecPrintf(x);
  arma::vec log_dprob_vec(ns, fill::zeros);
  for (arma::uword i = 0; i < ns; i++) {
    log_dprob_vec(i) = R::dbinom(x(i), n(i), p, 1);
  }
  double log_dprob = arma::accu(log_dprob_vec);
  double out = std::exp(log_dprob + R::pbinom(rK - arma::accu(x), nK, p, 1, 1));
  return out;
}

// [[Rcpp::export]]
double pet_recursive(int k, arma::ivec x, arma::ivec n, arma::ivec r, 
                     int nK, int rK, double p, int K) {
  /*
  x must be assigned as a zero vector of length K-1 when first calling pet_recursive
  */
  int bddpar = 0;
  int kcur = k;
  if (K < 2) {
    Rcpp::stop("avaliable for K > 1, if K = 1, use pbinom(r, n, p))");
  } else { 
    if (kcur < K) {
      if (kcur == 1) {
        x(0) = r(0) + 1;
      } else {
        bddpar = arma::accu(x.head(kcur));
      }
      int lowbdd = 0;
      int uppbdd = n(kcur-1);
      if ((r(kcur-1) + 1 - bddpar) > lowbdd) {
        lowbdd = r(kcur-1) + 1 - bddpar;
      }
      if ((rK - bddpar) < uppbdd) {
        uppbdd = rK - bddpar;
      } 
      //Rprintf("%d, %d - %d\n", kcur, lowbdd, uppbdd);
      int v;
      int kplus1 = kcur + 1;
      arma::vec out_new(uppbdd - lowbdd + 1, fill::zeros);
      for (v = lowbdd; v < (uppbdd + 1); v++) {
        //Rprintf("%d, %d\n", kcur, v);
        arma::ivec xin = x;
        xin(kcur-1) = v;
        out_new[v-lowbdd] = pet_recursive(k = kplus1, x = xin, n = n, r = r, nK = nK, rK = rK, p = p, K = K);
      }
      double out = arma::accu(out_new);
      return out;
    } else {
      double out = cumbin(x, n, r, nK, rK, p);
      return out;
    } 
  }
}

// [[Rcpp::export]]
List kStageP2A_Cpp(double p0, double p1, arma::ivec nseq, arma::ivec rseq)
{
  arma::uword nstage = nseq.n_elem;
  arma::ivec nEachInterim(nstage, fill::zeros);
  for (arma::uword i = 0; i < nstage; i++) {
    if (i == 0) {
      nEachInterim(i) = nseq(i);
    } else {
      nEachInterim(i) = nseq(i) - nseq(i-1);
    }
  }
  if (nstage < 2) {
    stop("Number of stages should be at least 2");
  }
  arma::vec rej_seq_0(nstage);
  arma::vec rej_seq_1(nstage);
  for (arma::uword i = 0; i < nstage; i++) {
    if (i == 0) {
      rej_seq_0(i) = std::exp(R::pbinom(rseq(0), nEachInterim(0), p0, 1, 1));
      rej_seq_1(i) = std::exp(R::pbinom(rseq(0), nEachInterim(0), p1, 1, 1));
    } else {
      arma::ivec xx(i, fill::zeros);
      rej_seq_0(i) = pet_recursive(1, xx, nEachInterim.head(i), rseq.head(i), nEachInterim(i), rseq(i), p0, i+1);
      rej_seq_1(i) = pet_recursive(1, xx, nEachInterim.head(i), rseq.head(i), nEachInterim(i), rseq(i), p1, i+1);
    }
  }
  
  arma::vec pet_seq = arma::cumsum(rej_seq_0.head(nstage-1));
  double rej_drug_0 = arma::accu(rej_seq_0);
  double rej_drug_1 = arma::accu(rej_seq_1);
  double t1e = 1 - rej_drug_0;
  double t2e = rej_drug_1;
  
  arma::vec prob_tmp(nstage, fill::zeros);
  for (arma::uword i = 0; i < nstage; i++) {
    if (i == 0) {
      prob_tmp(i) = 1; 
    } else { 
      prob_tmp(i) = 1 - pet_seq(i-1); 
    } 
  }
  double en = arma::accu(prob_tmp % nEachInterim);
  
  List result = List::create(
    Named("nseq") = nseq,
    Named("rseq") = rseq,
    Named("t1e") = t1e, 
    Named("t2e") = t2e, 
    Named("en") = en, 
    Named("pet_seq") = pet_seq
  );
  return result;
}



