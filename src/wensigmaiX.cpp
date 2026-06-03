#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::fmat wensigmaiX(arma::fmat& invsig, arma::fmat& Xmat){
 arma::fmat sigmaiX;

 sigmaiX=invsig*Xmat;
 //sigmaiY=invsig*Yvec;
 return(sigmaiX);
}
