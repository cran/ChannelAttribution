#ifndef header
#define header 

using namespace std;
using namespace Rcpp;

string to_string(T pNumber);
bool regex_search_v0(string string0, string exp0);	
bool regex_search_v1(string string0, vector<string> vexp0);
vector<string> split_string(string s, string sep);

RcppExport SEXP heuristic_models_cpp(SEXP Dy_p, SEXP var_channel_p, SEXP var_conv_p, SEXP var_value_p);
RcppExport SEXP markov_model_cpp(SEXP Dy_p, SEXP var_channel_p, SEXP var_conv_p, SEXP var_value_p, SEXP nsim_p, SEXP n_boot_p, SEXP n_single_boot_p, SEXP out_trans_matrx_p);


#endif