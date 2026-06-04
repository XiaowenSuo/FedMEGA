#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::fmat suoscore(arma::fmat& genoN, arma::fmat& genoM,arma::fmat& genoS,arma::fmat& PGenoN, arma::fmat& PGenoM,arma::fmat& PGenoS){
 arma::fmat res;

 res=genoN.t()*PGenoN+genoM.t()*PGenoM+genoS.t()*PGenoS;
 
 return(res);
}
