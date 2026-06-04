#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::export]]

Rcpp::List server_q(arma::mat& invsig ,const arma::mat& sigmaiX, const arma::mat& xsigmaix,const arma::mat& sigmaiY, arma::fvec& Yvec, arma::fvec& wVec, arma::fvec& tauVec){

arma::mat cov;
try {
  cov = arma::inv_sympd(arma::symmatu(xsigmaix));
} catch (const std::exception& e) {
  cov = arma::pinv(arma::symmatu(xsigmaix));
  std::cout << "inv_sympd failed, inverted with pinv" << std::endl;
}

int size = Yvec.size();
arma::mat A = arma::eye(size, size);
arma::mat sigmaiXt = sigmaiX.t();

arma::vec PY1 = sigmaiY - sigmaiX * (cov * (sigmaiXt * Yvec));
arma::vec APY = A * PY1;
float YPAPY = dot(PY1, APY);
arma::vec A0PY = PY1;
float YPA0PY = dot(PY1, A0PY);

arma::mat sigmaiXcov = sigmaiX*cov;
arma::mat sigmaiXcovT = sigmaiXcov.t();
arma::mat P = invsig - sigmaiX * sigmaiXcovT;
arma::vec PAPY = P*APY;
arma::mat AI(2,2);
arma::vec PA0PY = P * A0PY;
AI(0,0) =  dot(A0PY, PA0PY);
AI(1,1) = dot(APY, PAPY);
AI(0,1) = dot(A0PY, PAPY);
AI(1,0) = AI(0,1);

arma::vec traVec(2);
traVec(1) = arma::trace(P);
traVec(0) = arma::trace(P);

float score0 = YPA0PY - traVec(0);
float score1 = YPAPY - traVec(1);
arma::vec scoreVec(2);
scoreVec(0) = score0;
scoreVec(1) = score1;

arma::vec sigmaiXtY = sigmaiXt * Yvec;
arma::vec alpha = cov * sigmaiXtY;
arma::vec eta = Yvec - tauVec(0) * (sigmaiY - sigmaiX * alpha) / wVec;

return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
Rcpp::Named("eta") = eta,
Rcpp::Named("YPAPY") = YPAPY,
Rcpp::Named("YPA0PY") = YPA0PY,
Rcpp::Named("APY") =APY,
Rcpp::Named("PY")=PY1,
Rcpp::Named("P")=P,
Rcpp::Named("AI")=AI,
Rcpp::Named("trace") = traVec,
Rcpp::Named("score") = scoreVec,
Rcpp::Named("invsig") = invsig,
Rcpp::Named("sigmaiY") = sigmaiY,
Rcpp::Named("sigmaiX") = sigmaiX,
Rcpp::Named("xsigmaix") = xsigmaix);
}
