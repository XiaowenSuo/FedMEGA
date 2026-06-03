#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::fmat getxsigmaix(const arma::fmat& Xmat, const arma::fmat& sigmaiX) {

arma::fmat xsigmaix = Xmat.t()*sigmaiX;

 return(xsigmaix);
}
