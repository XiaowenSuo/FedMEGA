#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Struct to parallelize the computation
struct SuoFedWorker : public Worker {
  const arma::mat& genoM2;
  const arma::mat& genoN2;
  const arma::mat& genoS2;
  const arma::mat& XXVX_invXVM;
  const arma::mat& XXVX_invXVN;
  const arma::mat& XXVX_invXVS;
  const arma::vec& V;
  const arma::vec& varRatio_null;
  const int N_m; // Number of rows for M
  const int N_n; // Number of rows for N
  const int N_s; // Number of rows for S
  arma::vec& result;

  // Constructor
  SuoFedWorker(const arma::mat& genoM2, const arma::mat& genoN2, const arma::mat& genoS2,
               const arma::mat& XXVX_invXVM, const arma::mat& XXVX_invXVN, const arma::mat& XXVX_invXVS,
               const arma::vec& V, const arma::vec& varRatio_null,
               const int N_m, const int N_n, const int N_s, arma::vec& result)
    : genoM2(genoM2), genoN2(genoN2), genoS2(genoS2),
      XXVX_invXVM(XXVX_invXVM), XXVX_invXVN(XXVX_invXVN), XXVX_invXVS(XXVX_invXVS),
      V(V), varRatio_null(varRatio_null), N_m(N_m), N_n(N_n), N_s(N_s), result(result) {}

  // Operator to perform the computation in parallel
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      arma::mat GM0 = genoM2.col(j);
      arma::mat GN0 = genoN2.col(j);
      arma::mat GS0 = genoS2.col(j);

      arma::mat ZGM = XXVX_invXVM * GM0;
      arma::mat ZGN = XXVX_invXVN * GN0;
      arma::mat ZGS = XXVX_invXVS * GS0;
      arma::mat ZG = ZGM + ZGN + ZGS;

      // Extract submatrices based on N_m, N_n, N_s
      ZGN = ZG.submat(0, 0, N_n - 1, 0);
      ZGM = ZG.submat(N_n, 0, N_m + N_n - 1, 0);
      ZGS = ZG.submat(N_m + N_n, 0, N_m + N_n + N_s - 1, 0);

      arma::mat GM = GM0 - ZGM;
      arma::mat GN = GN0 - ZGN;
      arma::mat GS = GS0 - ZGS;
      arma::mat G = arma::join_vert(GN, GM, GS);

      arma::mat GWG = V.t() * (G % G);
      arma::vec VAR = varRatio_null % GWG;

      // Store the result
      result(j) = arma::as_scalar(VAR);
    }
  }
};

// [[Rcpp::export]]
arma::vec parallel_fedgg(const int& S, const arma::mat& genoM2, const arma::mat& genoN2, const arma::mat& genoS2,
                         const arma::mat& XXVX_invXVM, const arma::mat& XXVX_invXVN, const arma::mat& XXVX_invXVS,
                         const arma::vec& V, const arma::vec& varRatio_null,
                         const int N_m, const int N_n, const int N_s) {
  // Initialize result vector
  arma::vec result(S, arma::fill::zeros);

  // Create worker
  SuoFedWorker worker(genoM2, genoN2, genoS2, XXVX_invXVM, XXVX_invXVN, XXVX_invXVS,
                      V, varRatio_null, N_m, N_n, N_s, result);

  // Use parallelFor to perform the computation
  parallelFor(0, S, worker);

  return result;
}
