#include <RcppEigen.h>
#include <cmath>   // exp
#include <random>  // normal distribution
#include <vector>  // vector class
#include <pcg/pcg_random.hpp> // pcg prng

#include "armmr_types.h"
#include "pcg.h"  // pcg seeding

using namespace Rcpp;





//' Simulate populations with competition.
//'
//' @param max_t Number of time points to simulate.
//' @param N0 Vector of starting population abundances, one for each species.
//' @param r Vector of growth rates, one for each species.
//' @param alpha Matrix of intra- and inter-specific density dependences.
//' @param sigma Standard deviation of process error.
//'
//'
//' @export
//'
//' @examples
//' sim_pops(10, c(10, 10), c(0.5, 0.5), matrix(rep(1e-3, 4), 2, 2), 0.25)
//'
//[[Rcpp::export]]
MatrixXd sim_pops(const uint32& max_t,
                  const Map<RowVectorXd> N0,
                  const Map<RowVectorXd> r,
                  const Map<MatrixXd> alpha,
                  const double& sigma) {

    uint32 n_spp = N0.size();
    if (r.size() != n_spp || alpha.rows() != n_spp || alpha.cols() != n_spp) {
        stop("N0, r, and alpha lengths must be the same.");
    }
    pcg32 eng = seeded_pcg();

    MatrixXd N(max_t + 1, n_spp);
    for (uint32 i = 0; i < n_spp; i++) N(0,i) = N0(i);

    std::normal_distribution<double> distr;
    if (sigma > 0) distr = std::normal_distribution<double>(0.0, sigma);

    RowVectorXd rnd = RowVectorXd::Zero(n_spp);
    RowVectorXd tmp(n_spp);

    for (uint32 t = 0; t < max_t; t++) {

        if (sigma > 0) {
            for (uint32 i = 0; i < rnd.size(); i++) rnd(i) = distr(eng);
        }

        tmp = r - N.row(t) * alpha + rnd;

        N.row(t+1) = N.row(t).array() * tmp.array().exp();

    }

    return N;

}

