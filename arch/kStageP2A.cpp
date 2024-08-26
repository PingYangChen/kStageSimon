#include <Rcpp.h>
#include <stdlib.h>
using namespace Rcpp;
using namespace std;

void vecPrintf(const Rcpp::NumericVector &v)
{
  for (int i = 0; i < v.length(); i++) Rprintf("%4.4f\t", v(i));
  Rprintf("\n\n");
}

// [[Rcpp::export]]
double cumbin(Rcpp::NumericVector x, Rcpp::NumericVector n, Rcpp::NumericVector r, 
              int nK, int rK, double p) {
  //double dprob = 1;
  double log_dprob = 0;
  int ns = x.length(); 
  //vecPrintf(x);
  for (int i = 0; i < ns; i++) {
    //dprob = dprob*R::dbinom(x[i], n[i], p, 0);
    log_dprob += R::dbinom(x[i], n[i], p, 1);
  }
  //double out = dprob*R::pbinom(rK - Rcpp::sum(x), nK, p, 1, 0);
  double out = std::exp(log_dprob + R::pbinom(rK - Rcpp::sum(x), nK, p, 1, 1));
  return(out);
}

// [[Rcpp::export]]
double pet_recursive(int k, Rcpp::NumericVector x, Rcpp::NumericVector n, Rcpp::NumericVector r, 
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
        x[0] = r[0] + 1;
      } else {
        for (int i = 0; i < (kcur - 1); i++) {
          bddpar += x[i];
        }
      }
      int lowbdd = 0;
      int uppbdd = n[kcur-1];
      if ((r[kcur-1] + 1 - bddpar) > lowbdd) {
        lowbdd = r[kcur-1] + 1 - bddpar;
      }
      if ((rK - bddpar) < uppbdd) {
        uppbdd = rK - bddpar;
      }
      //Rprintf("%d, %d - %d\n", kcur, lowbdd, uppbdd);
      double out = 0.0;
      int v;
      int kplus1 = kcur + 1;
      for (v = lowbdd; v < (uppbdd + 1); v++) {
        //Rprintf("%d, %d\n", kcur, v);
        Rcpp::NumericVector xin = Rcpp::clone(x);
        xin[kcur-1] = v;
        double out_new = pet_recursive(k = kplus1, x = xin, n = n, r = r, nK = nK, rK = rK, p = p, K = K);
        out += out_new;
      }
      return out;
    } else {
      double out = cumbin(x, n, r, nK, rK, p);
      return out;
    }
  }
}

// [[Rcpp::export]]
List kStageP2A_Cpp(double p0, double p1, Rcpp::NumericVector nseq, Rcpp::NumericVector rseq)
{
  int nstage = nseq.length();
  Rcpp::NumericVector nEachInterim(nstage);
  for (int i = 0; i < nstage; i++) {
    if (i == 0) {
      nEachInterim[i] = nseq[i];
    } else {
      nEachInterim[i] = nseq[i] - nseq[i-1];
    }
  }
  if (nstage < 2) {
    stop("Number of stages should be at least 2");
  }
  Rcpp::NumericVector rej_seq_0(nstage);
  Rcpp::NumericVector rej_seq_1(nstage);
  for (int i = 0; i < nstage; i++) {
    if (i == 0) {
      rej_seq_0[i] = std::exp(R::pbinom(rseq[0], nEachInterim[0], p0, 1, 1));
      rej_seq_1[i] = std::exp(R::pbinom(rseq[0], nEachInterim[0], p1, 1, 1));
    } else {
      Rcpp::IntegerVector loc = Rcpp::seq(0, i-1);
      Rcpp::NumericVector xx(i);
      rej_seq_0[i] = pet_recursive(1, xx, nEachInterim[loc], rseq[loc], nEachInterim[i], rseq[i], p0, i+1);
      rej_seq_1[i] = pet_recursive(1, xx, nEachInterim[loc], rseq[loc], nEachInterim[i], rseq[i], p1, i+1);
    }
  }
  
  //Rcpp:IntegerVector pidx = ;
  Rcpp::NumericVector pet_seq = Rcpp::cumsum(rej_seq_0[Rcpp::seq(0, nstage - 2)]);
  double rej_drug_0 = Rcpp::sum(rej_seq_0);
  double rej_drug_1 = Rcpp::sum(rej_seq_1);
  double t1e = 1 - rej_drug_0;
  double t2e = rej_drug_1;
  
  Rcpp::NumericVector prob_tmp(nstage);
  for (int i = 0; i < nstage; i++) {
    if (i == 0) {
      prob_tmp[i] = 1; 
    } else {
      prob_tmp[i] = 1 - pet_seq[i-1]; 
    }
  }
  double en = Rcpp::sum(prob_tmp*nEachInterim);
  
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



