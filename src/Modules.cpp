#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_mod) {


    class_<rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> >("model_armm")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_namespace::model_armm, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_mod) {


    class_<rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> >("model_armm_ss")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_namespace::model_armm_ss, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_b_mod) {


    class_<rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> >("model_armm_ss_b")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_b_namespace::model_armm_ss_b, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_beta_mod) {


    class_<rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> >("model_armm_ss_beta")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_beta_namespace::model_armm_ss_beta, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_catlogit_mod) {


    class_<rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> >("model_armm_ss_catlogit")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_catlogit_namespace::model_armm_ss_catlogit, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_hlnp_mod) {


    class_<rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> >("model_armm_ss_hlnp")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_hlnp_namespace::model_armm_ss_hlnp, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_hnb_mod) {


    class_<rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> >("model_armm_ss_hnb")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_hnb_namespace::model_armm_ss_hnb, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_hp_mod) {


    class_<rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> >("model_armm_ss_hp")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_hp_namespace::model_armm_ss_hp, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_lnb_mod) {


    class_<rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> >("model_armm_ss_lnb")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_lnb_namespace::model_armm_ss_lnb, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_lnp_mod) {


    class_<rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> >("model_armm_ss_lnp")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_lnp_namespace::model_armm_ss_lnp, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_nb_mod) {


    class_<rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> >("model_armm_ss_nb")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_nb_namespace::model_armm_ss_nb, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_p_mod) {


    class_<rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> >("model_armm_ss_p")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_p_namespace::model_armm_ss_p, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_priors_mod) {


    class_<rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> >("model_armm_ss_priors")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_priors_namespace::model_armm_ss_priors, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_reps_mod) {


    class_<rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> >("model_armm_ss_reps")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_reps_namespace::model_armm_ss_reps, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_timehack_mod) {


    class_<rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> >("model_armm_ss_timehack")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_timehack_namespace::model_armm_ss_timehack, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_zilnp_mod) {


    class_<rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> >("model_armm_ss_zilnp")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_zilnp_namespace::model_armm_ss_zilnp, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_zinb_mod) {


    class_<rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> >("model_armm_ss_zinb")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_zinb_namespace::model_armm_ss_zinb, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4armm_ss_zip_mod) {


    class_<rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> >("model_armm_ss_zip")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_armm_ss_zip_namespace::model_armm_ss_zip, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4mm_mod) {


    class_<rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> >("model_mm")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_mm_namespace::model_mm, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
