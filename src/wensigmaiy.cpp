#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::fvec wensigmaiY(arma::fmat& invsig, arma::fvec& Yvec){
 arma::fvec sigmaiY;

 //sigmaiX=invsig*Xmat;
 sigmaiY=invsig*Yvec;
 return(sigmaiY);
}
