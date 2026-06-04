#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]

arma::mat invsig(arma::fvec& wVec, arma::fvec& tauVec){

arma::fvec diagVec;
diagVec = tauVec(1) + tauVec(0)/wVec;

unsigned int N= wVec.size();
for(unsigned int i=0; i< N; i++){
if(diagVec(i) < 1e-4){
diagVec(i) = 1e-4;
}
}

int size = diagVec.size();
arma::mat diagMat = arma::zeros(size, size);

for (int i = 0; i < size; ++i) {
diagMat(i, i) = diagVec(i);
}

arma::mat invsig = arma::inv_sympd(diagMat);

return(invsig);
}
