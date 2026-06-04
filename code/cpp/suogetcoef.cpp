#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]

Rcpp::List suogetcoeff(arma::fvec& wVec, arma::fvec& tauVec, arma::fvec& Yvec, const arma::mat& Xmat){

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
arma::vec sigmaiY =invsig * Yvec;
//return(sigmaiY);
arma::mat sigmaiX = invsig * Xmat;
 
arma::mat xsigmaix = Xmat.t()*sigmaiX;

arma::mat cov;
try {
  cov = arma::inv_sympd(arma::symmatu(xsigmaix));
} catch (const std::exception& e) {
  cov = arma::pinv(arma::symmatu(xsigmaix));
  Rcpp::Rcout << "inv_sympd failed, inverted with pinv" << std::endl;
}


arma::mat A = arma::eye(size, size);
arma::mat sigmaiXt = sigmaiX.t();

arma::vec PY1 = sigmaiY - sigmaiX * (cov * (sigmaiXt * Yvec));
arma::vec APY = A * PY1;
float YPAPY = dot(PY1, APY);

arma::mat sigmaiXcov = sigmaiX*cov;
arma::mat sigmaiXcovT = sigmaiXcov.t();
arma::mat P = invsig - sigmaiX * sigmaiXcovT;
arma::vec PAPY = P*APY;
float AI = dot(APY, PAPY);
float trace=arma::trace(P);
float score=YPAPY-trace;


arma::vec sigmaiXtY = sigmaiXt * Yvec;
arma::vec alpha = cov * sigmaiXtY;
arma::vec eta = Yvec - tauVec(0) * (sigmaiY - sigmaiX * alpha) / wVec;

return Rcpp::List::create(Rcpp::Named("alpha") = alpha, Rcpp::Named("eta") = eta, Rcpp::Named("YPAPY") = YPAPY, Rcpp::Named("APY") = APY, 
Rcpp::Named("PY")=PY1, Rcpp::Named("P")=P,Rcpp::Named("AI")=AI,Rcpp::Named("trace") = trace, Rcpp::Named("score") = score,Rcpp::Named("cov") = cov,Rcpp::Named("invsig") = invsig, 
Rcpp::Named("sigmaiY") = sigmaiY,Rcpp::Named("sigmaiX") = sigmaiX,Rcpp::Named("xsigmaix") = xsigmaix);
}






