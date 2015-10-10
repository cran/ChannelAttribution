#ifndef header
#define header 

using namespace std;
using namespace Rcpp;

string to_string(T pNumber);

RcppExport SEXP heuristic_models_cpp(SEXP Dy_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p);
RcppExport SEXP markov_model_cpp(SEXP Dy_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP var_null_p, SEXP nsim_p, SEXP n_boot_p, SEXP n_single_boot_p, SEXP out_trans_matrx_p);

#endif