//ChannelAttribution: Markov model for online multi-channel attribution
//Copyright (C) 2015 - 2020  Davide Altomare and David Loris <http://www.channelattribution.net>
//
//ChannelAttribution is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//ChannelAttribution is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with ChannelAttribution.  If not, see <http://www.gnu.org/licenses/>.

#ifndef header
#define header 

using namespace std;
using namespace Rcpp;

string to_string(T pNumber);
vector<long int> split_string(const string &s, unsigned long int order);

RcppExport SEXP heuristic_models_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP sep_p);
RcppExport SEXP choose_order_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_null_p, SEXP max_order_p, SEXP sep_p, SEXP ncore_p, SEXP roc_npt_p);
RcppExport SEXP markov_model_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP var_null_p, SEXP order_p, SEXP nsim_start_p, SEXP max_step_p, SEXP out_more_p, string sep_p, SEXP ncore_p, SEXP nfold_p, SEXP seed_p, SEXP conv_par_p, SEXP rate_step_sim_p, SEXP verbose_p);
transition_matrix_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_null_p, SEXP order_p, SEXP sep_p, SEXP flg_equal_p);

#endif
