#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]

Rcpp::List server(arma::mat& invsig ,const arma::mat& sigmaiX, const arma::mat& xsigmaix,const arma::mat& sigmaiY, arma::fvec& Yvec,arma::fvec& wVec, arma::fvec& tauVec){
arma::mat cov;
try {
  cov = arma::inv_sympd(arma::symmatu(xsigmaix));
} catch (const std::exception& e) {
  cov = arma::pinv(arma::symmatu(xsigmaix));
  Rcpp::Rcout << "inv_sympd failed, inverted with pinv" << std::endl;
}

int size = Yvec.size();
arma::mat A = arma::eye(size, size);
arma::mat sigmaiXt = sigmaiX.t();

arma::vec PY1 = sigmaiY - sigmaiX * (cov * (sigmaiXt * Yvec));
arma::vec APY = A * PY1;
float YPAPY = dot(PY1, APY);

arma::mat sigmaiXcov = sigmaiX*cov;
arma::mat sigmaiXcovT = sigmaiXcov.t();
arma::mat P = invsig - sigmaiX * sigmaiXcovT;
arma::vec PAPY_1 = invsig * APY;
arma::vec PAPY = PAPY_1 - sigmaiX * (cov * (sigmaiXt * PAPY_1));

float AI = dot(APY, PAPY);
float trace=arma::trace(P);
float score=YPAPY-trace;

arma::vec sigmaiXtY = sigmaiXt * Yvec;
arma::vec alpha = cov * sigmaiXtY;
arma::vec eta = Yvec - tauVec(0) * (sigmaiY - sigmaiX * alpha) / wVec;

return Rcpp::List::create(Rcpp::Named("alpha") = alpha, Rcpp::Named("eta") = eta, Rcpp::Named("YPAPY") = YPAPY, Rcpp::Named("APY") = APY,
Rcpp::Named("PY")=PY1, Rcpp::Named("P")=P,Rcpp::Named("AI")=AI,Rcpp::Named("trace") = trace, Rcpp::Named("score") = score,Rcpp::Named("cov") = cov);
}
