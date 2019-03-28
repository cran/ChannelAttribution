#ifndef header
#define header 

using namespace std;
using namespace Rcpp;

string to_string(T pNumber);
vector<long int> split_string(const string &s, unsigned long int order);

RcppExport SEXP heuristic_models_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP sep_p);
RcppExport SEXP markov_model_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP var_null_p, SEXP order_p, SEXP nsim_p, SEXP max_step_p, SEXP out_more_p, SEXP sep_p);
RcppExport SEXP choose_order_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_null_p, SEXP max_order_p, SEXP sep_p, SEXP ncore_p, SEXP roc_npt_p)
RcppExport SEXP markov_model_mp_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP var_null_p, SEXP order_p, SEXP nsim_start_p, SEXP max_step_p, SEXP out_more_p, string sep_p, SEXP ncore_p, SEXP nfold_p, SEXP seed_p, SEXP conv_par_p, SEXP rate_step_sim_p, SEXP verbose_p);
transition_matrix_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_null_p, SEXP order_p, SEXP sep_p, SEXP flg_equal_p);

#endif
