#include <RcppEigen.h>
#include <pcg/pcg_random.hpp> // pcg prng
#include <random>  // normal distribution
#include <vector>  // vector class
#ifdef _OPENMP
#include <omp.h>  // OpenMP
#endif


#include "armmr_types.h"
#include "pcg.h"  // pcg seeding

using namespace Rcpp;




// /*
//  One population simulated using AR1 process.
//  N in this case is already logged
//  */
// //[[Rcpp::export]]
// VectorXd sim_pop_ar(const VectorXd& X, const double& N0,
//                      const double& b0, const double& b1, const double& rho,
//                      const double& sigma) {
//     uint32 n_gen = X.n_elem - 1;
//     VectorXd N(n_gen + 1);
//     N(0) = N0;
//
//     std::normal_distribution<double> rnorm_distr(0.0, sigma);
//     pcg32 engine = seeded_pcg();
//
//     for (uint32 t = 0; t < n_gen; t++) {
//         N(t+1) = b0 + b1 * X(t+1) + rho * ( N(t) - b0 - b1 * X(t) );
//         N(t+1) += rnorm_distr(engine);
//     }
//
//     return N;
// }






/*
 Check that input parameters have the correct dimensionality.
 For more info, see description of `sim_pops_ar` below.
 */
void check_dims(const Map<MatrixXd> X,
                const Map<MatrixXd> N0_mat,
                const Map<MatrixXd> b0_mat,
                const Map<MatrixXd> b1_mat,
                const Map<MatrixXd> rho_mat,
                const std::vector<MatrixXd>& vcv_cube) {

  uint32 n_locs = X.cols();
  if (N0_mat.cols() != n_locs) {
    stop("N0_mat should have the same number of columns as does X");
  }
  if (b0_mat.cols() != n_locs) {
    stop("b0_mat should have the same number of columns as does X");
  }
  if (b1_mat.cols() != n_locs) {
    stop("b1_mat should have the same number of columns as does X");
  }
  if (rho_mat.cols() != n_locs) {
    stop("rho_mat should have the same number of columns as does X");
  }
  if (vcv_cube.size() != n_locs) {
    stop("vcv_cube should have the same number of slices as X has columns");
  }

  uint32 n_spp = N0_mat.rows();
  if (N0_mat.rows() != n_spp) {
    stop("N0_mat should have the same number of rows as does N0");
  }
  if (b0_mat.rows() != n_spp) {
    stop("b0_mat should have the same number of rows as does N0");
  }
  if (b1_mat.rows() != n_spp) {
    stop("b1_mat should have the same number of rows as does N0");
  }
  if (rho_mat.rows() != n_spp) {
    stop("rho_mat should have the same number of rows as does N0");
  }
  if (vcv_cube[0].rows() != n_spp || vcv_cube[0].cols() != n_spp) {
    stop("vcv_cube should have the same number of rows and columns as N0 has rows");
  }
  for (uint32 i = 1; i < vcv_cube.size(); i++) {
    if (vcv_cube[i].rows() != n_spp || vcv_cube[i].cols() != n_spp) {
      stop("All items in vcv_cube should have the same dimensions.");
    }
  }

  return;
}


/*
 Turn variance-covariance matrix into matrix that can be used to generate random
 normal deviates with same covariance structure

 "a vector of independent normal random variables,
 when multiplied by the transpose of the Cholesky deposition of [vcv] will
 have covariance matrix equal to [vcv]."
 */
std::vector<MatrixXd> make_chol_decomp(const std::vector<MatrixXd>& vcv_cube) {

  std::vector<MatrixXd> chol_decomp(vcv_cube.size(), vcv_cube.front());
  for (uint32 i = 0; i < vcv_cube.size(); i++) {
    Eigen::LLT<MatrixXd> tmp(vcv_cube[i]);
    // no need to transpose bc `matrixL` is transpose of `matrixU`
    // and `matrixU` is typical output from R's `chol` fxn
    chol_decomp[i] = tmp.matrixL();
  }

  return chol_decomp;
}


/*
 For testing in R
 */
//[[Rcpp::export]]
std::vector<MatrixXd> make_chol_decomp_cpp(const List& vcv_cube) {

  std::vector<MatrixXd> vcv_cube_(vcv_cube.size());
  for (uint32 i = 0; i < vcv_cube.size(); i++) {
    vcv_cube_[i] = as<MatrixXd>(vcv_cube[i]);
  }

  std::vector<MatrixXd> chol_decomp = make_chol_decomp(vcv_cube_);

  return chol_decomp;
}




//' Multiple populations simulated using AR1 process.
//'
//' Input and output N values are logged.
//' All input matrices other than `X` and `vcv_cube` should have rows associated
//' with a given species and columns associated with a given location.
//' See descriptions for `X` and `vcv_cube`.
//'
//' @param X Matrix of environmental variable.
//'     It should have rows associated with a given time point and
//'     columns associated with a given location.
//' @param N0_mat Matrix of starting population abundances (`log(# individuals)`)
//'     by species and location.
//' @param b0_mat Matrix of \eqn{\beta_0} (the population-abundance intercept) values
//'     by species and location.
//' @param b1_mat Matrix of \eqn{\beta_1} (the effect of \eqn{X} on \eqn{N}) values
//'     by species and location.
//' @param rho_mat Matrix of growth rates by species and location.
//' @param vcv_cube Cube representing variance-covariance matrices for process error
//'     among species, one matrix for each location.
//'     It should have rows and columns associated with a given species,
//'     and slices associated with a given location.
//' @param obs_sigma Vector of standard deviations of observation error for each species.
//' @param n_threads Number of cores to use. Defaults to 1.
//'
//'
//'
//' @return A 3-dimensional array.
//' The output will have rows associated with a given time point,
//' columns associated with a given species, and
//' slices associated with a given location.
//'
//'
//'
//'
//' @export
//'
//' @examples
//' n_spp <- 3
//' max_t <- 10
//' n_locs <- 2
//' X <- matrix(rlnorm(max_t * n_locs), max_t, n_locs)
//' N0 <- matrix(log(10), n_spp, n_locs)
//' b0 <- matrix(log(100), n_spp, n_locs)
//' b1 <- matrix(0.1, n_spp, n_locs)
//' rho <- matrix(0.2, n_spp, n_locs)
//' vcv <- diag(n_spp)
//' vcv[lower.tri(vcv)] <- vcv[upper.tri(vcv)] <- 0.1
//' vcv <- replicate(n_locs, vcv, simplify = FALSE)
//' obs <- rep(0.1, n_spp)
//' sim_pops_ar(X, N0, b0, b1, rho, vcv, obs)
//'
//[[Rcpp::export]]
std::vector<MatrixXd> sim_pops_ar(const Map<MatrixXd> X,
                                  const Map<MatrixXd> N0_mat,
                                  const Map<MatrixXd> b0_mat,
                                  const Map<MatrixXd> b1_mat,
                                  const Map<MatrixXd> rho_mat,
                                  const List& vcv_cube,
                                  const Map<VectorXd> obs_sigma,
                                  const uint32& n_threads = 1) {

  std::vector<MatrixXd> vcv_cube_(vcv_cube.size());
  for (uint32 i = 0; i < vcv_cube.size(); i++) {
    vcv_cube_[i] = as<MatrixXd>(vcv_cube[i]);
  }

  check_dims(X, N0_mat, b0_mat, b1_mat, rho_mat, vcv_cube_);

  // For random number generator
  const std::vector<std::vector<uint64>> seeds = mc_seeds(n_threads);
  // For turning ~N(0,1) to multivariate with given covariance matrix
  const std::vector<MatrixXd> chol_decomp(make_chol_decomp(vcv_cube_));

  const uint32 n_time(X.rows());
  const uint32 n_locs(X.cols());
  const uint32 n_spp(N0_mat.rows());

  std::vector<MatrixXd> N(n_locs, MatrixXd(n_time, n_spp));
  for (uint32 i = 0; i < n_locs; i++) N[i].row(0) = N0_mat.col(i).transpose();

#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(n_threads) if(n_threads > 1)
{
#endif

  std::vector<uint64> active_seeds;

  // Write the active seed per core or just write one of the seeds.
#ifdef _OPENMP
  uint32 active_thread = omp_get_thread_num();
  active_seeds = seeds[active_thread];
#else
  active_seeds = seeds[0];
#endif

  pcg32 engine = seeded_pcg(active_seeds);
  // Random normal distribution:
  std::normal_distribution<double> rnorm_distr(0.0, 1.0);

  // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (uint32 loc = 0; loc < n_locs; loc++) {
    const MatrixXd& cd_loc(chol_decomp[loc]);
    const VectorXd& b0s(b0_mat.col(loc));
    const VectorXd& b1s(b1_mat.col(loc));
    const VectorXd& rhos(rho_mat.col(loc));
    const VectorXd& Xs(X.col(loc));
    MatrixXd& Ns(N[loc]);
    RowVectorXd proc_rnd(n_spp);
    RowVectorXd obs_rnd(n_spp);
    for (uint32 t = 0; t < n_time - 1; t++) {
      for (uint32 sp = 0; sp < n_spp; sp++) {
        Ns(t+1,sp) = b0s(sp) + b1s(sp) * Xs(t+1) + rhos(sp) * (
          Ns(t,sp) - b0s(sp) - b1s(sp) * Xs(t));
        proc_rnd(sp) = rnorm_distr(engine);
        obs_rnd(sp) = rnorm_distr(engine) * obs_sigma(sp);
      }
      // Generate variance and covariance:
      proc_rnd = cd_loc * proc_rnd.transpose();
      Ns.row(t+1) += proc_rnd;
      Ns.row(t+1) += obs_rnd;
    }
  }

#ifdef _OPENMP
}
#endif


return N;
}



// //' Melt a cube into a single data frame.
// //'
// //' @param C Three-dimensional array that you want to melt into a two-dimensional
// //'     data frame.
// //'
// //' @return A melted data frame.
// //'
// //' @noRd
// //'
// //[[Rcpp::export]]
// DataFrame melt_cube(const arma::cube& C) {
//
//     uint32 C_rows = C.rows();
//     MatrixXd M(C.rows() * C.n_slices, C.cols() + 1);
//     for (uint32 i = 0; i < C.n_slices; i++) {
//         M(arma::span(i * C_rows, (i+1) * C_rows - 1), arma::span(0)).fill(i+1);
//         M(arma::span(i * C_rows, (i+1) * C_rows - 1),
//           arma::span(1, C.cols())) = C.slice(i);
//     }
//
//     return M;
// }
//
//
// //' Generate parameter values for simulations.
// //'
// //' Generate parameter values to simulate multi-location, multi-species time series data.
// //'
// //'
// //' @param n_time Number of time steps.
// //' @param n_loc Number of locations.
// //' @param n_spp Number of species.
// //' @param mean_b0 Mean for the b0 parameter relating X to N.
// //' @param mean_b1 Mean for the b1 parameter relating X to N.
// //' @param mean_rho Mean for the rho parameter relating X to N.
// //'     This parameter is on the inverse logit scale.
// //' @param sigma_b0 Standard deviation for the b0 parameter relating X to N.
// //' @param sigma_b1 Standard deviation for the b1 parameter relating X to N.
// //' @param sigma_rho Standard deviation for the rho parameter relating X to N.
// //'     This parameter is on the inverse logit scale.
// //' @param sigma_eps Standard deviation for the epsilon parameter.
// //' @param sigma_obs Standard deviation for observation error. Defaults to 0.
// //' @param corr_method Method for determining correlations between species.
// //'     Options are "none", "phylo", or "random". Defaults to "none".
// //'
// //'
// //' @export
// //'
// //' @examples
// //' generate_pars(10, 2, 3,
// //'               mean_b0 = log(100),
// //'               mean_b1 = 0.1,
// //'               mean_rho = 0.25,
// //'               sigma_b0 = 0.1,
// //'               sigma_b1 = 0.1,
// //'               sigma_rho = 0.1,
// //'               sigma_eps = 0.1)
// //'
// //'
// //[[Rcpp::export]]
// List generate_pars(const uint32& n_time,
//                    const uint32& n_loc,
//                    const uint32& n_spp,
//                    const double& mean_b0,
//                    const double& mean_b1,
//                    const double& mean_rho,
//                    const double& sigma_b0,
//                    const double& sigma_b1,
//                    const double& sigma_rho,
//                    const double& sigma_eps,
//                    const double& sigma_obs = 0,
//                    const std::string& corr_method = "none") {
//
//     std::normal_distribution<double> rnorm_distr(0.0, 1.0);
//     pcg32 eng = seeded_pcg();
//
//     // Set up matrices
//     // Environmental variables
//     MatrixXd X_(n_time, n_loc, arma::fill::zeros);
//     VectorXd rnd(n_time);
//     for (double& r : rnd) r = rnorm_distr(eng);
//     X_.each_col() += rnd;
//     // Mean abundance in average environment
//     MatrixXd b0_mat_(n_spp, n_loc, arma::fill::zeros);
//     rnd.set_size(n_spp);
//     for (double& r : rnd) r = rnorm_distr(eng) * sigma_b0 + mean_b0;
//     b0_mat_.each_col() += rnd;
//     // Response to environmental variables
//     MatrixXd b1_mat_(n_spp, n_loc, arma::fill::zeros);
//     for (double& r : rnd) r = rnorm_distr(eng) * sigma_b1 + mean_b1;
//     b1_mat_.each_col() += rnd;
//     // Autoregressive parameter
//     MatrixXd rho_mat_(n_spp, n_loc, arma::fill::zeros);
//     for (double& r : rnd) r = rnorm_distr(eng) * sigma_rho + mean_rho;
//     rnd = arma::exp(rnd) / (1 + arma::exp(rnd));
//     rho_mat_.each_col() += rnd;
//     // Observation error
//     VectorXd obs_sigma_(n_spp);
//     obs_sigma_.fill(sigma_obs);
//
//     // Initial population size (set to b0)
//     MatrixXd N0_mat_(b0_mat_);
//
//     MatrixXd vcv_(n_spp, n_spp, arma::fill::eye);
//     arma::cube vcv_cube_(n_spp, n_spp, n_loc, arma::fill::zeros);
//     if (corr_method == "phylo") {
//         Environment ape = Environment::namespace_env("ape");
//         Function rcoal = ape["rcoal"];
//         Function vcv_phylo = ape["vcv.phylo"];
//         SEXP phylo = rcoal(n_spp);
//         SEXP vcv_SEXP_ = vcv_phylo(phylo);
//         vcv_ = as<MatrixXd>(vcv_SEXP_);
//         vcv_ /= vcv_.diag()(0);
//     } else if (corr_method == "random") {
//         uint32 n_corrs = ((n_spp - 1) / 2) * (1 + (n_spp - 1));
//         VectorXd rnd_phy(n_corrs);
//         for (double& r : rnd_phy) r = runif_ab(eng, -1, 1);
//         for (uint32 i = 0, rnd_i = 0; i < n_spp; i++) {
//             for (uint32 j = i+1; j < n_spp; j++, rnd_i++) {
//                 vcv_(i,j) = rnd_phy(rnd_i);
//                 vcv_(j,i) = rnd_phy(rnd_i);
//             }
//         }
//     } else if (corr_method == "none") {
//         ;
//     } else {
//         stop("corr_method must be none, phylo, or random.");
//     }
//     // Going from correlations to covariances:
//     vcv_ *= (sigma_eps * sigma_eps);
//     vcv_cube_.each_slice() += vcv_;
//
//     List out = List::create(_["X"] = X_,
//                             _["N0_mat"] = N0_mat_,
//                             _["b0_mat"] = b0_mat_,
//                             _["b1_mat"] = b1_mat_,
//                             _["rho_mat"] = rho_mat_,
//                             _["vcv_cube"] = vcv_cube_,
//                             _["obs_sigma"] = obs_sigma_);
//
//     return out;
// }
