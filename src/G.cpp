#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List G_G(arma::fmat& XXVX_inv, arma::fmat& XV){

arma::fmat XXVX_invXV=XXVX_inv* XV;

return Rcpp::List::create(Rcpp::Named("XXVX_invXV") =XXVX_invXV);
}
