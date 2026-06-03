#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]

Rcpp::List suofitglmm(arma::fvec &tauVec, float score, float AI, float tol){


float Dtau = score/AI;

arma::fvec tau0 = tauVec;
tauVec(1) = tau0(1) + Dtau;

for(int i=0; i<tauVec.n_elem; ++i) {
if (tauVec(i) < tol){
tauVec(i) = 0;
}
}

float step = 1.0;
while (tauVec(1) < 0.0){

step = step*0.5;
tauVec(1) = tau0(1) + step * Dtau;

}

for(int i=0; i<tauVec.n_elem; ++i) {
if (tauVec(i) < tol){
tauVec(i) = 0;
}
}

return  Rcpp::List::create(Rcpp::Named("tau") = tauVec);
}
