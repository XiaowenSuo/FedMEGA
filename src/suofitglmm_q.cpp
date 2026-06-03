#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List suofitglmm_q(arma::fvec& tauVec, arma::fvec& score, arma::fmat& AI, float tol){

arma::uvec zeroVec = (tauVec < tol);
arma::fvec Dtau = arma::solve(AI, score);
arma::fvec tau0 = tauVec;
tauVec = tau0 + Dtau;

tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros(); 
float step = 1.0;

while (tauVec(0) < 0.0 || tauVec(1)  < 0.0){      	
step = step*0.5;
tauVec = tau0 + step * Dtau;     	
tauVec.elem( find(zeroVec % (tauVec < tol)) ).zeros();     	
}

tauVec.elem( find(tauVec < tol) ).zeros();
return Rcpp::List::create(Rcpp::Named("tau") = tauVec);
}


