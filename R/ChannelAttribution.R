
heuristic_models=function(Dy, var_channel, var_conv, var_value){

 res=.Call("heuristic_models_cpp", Dy, var_channel, var_conv, var_value)
 return(as.data.frame(res)) 

}	


markov_model=function(Dy, var_channel, var_conv, var_value, nsim=2000000, n_boot=1000000, n_single_boot=30, out_trans_matrx=0){

 res=.Call("markov_model_cpp", Dy, var_channel, var_conv, var_value, nsim, n_boot, n_single_boot, out_trans_matrx)
 return(as.data.frame(res)) 

}	
 
