
heuristic_models=function(Data, var_path, var_conv, var_value=NULL, sep=">"){

 if(length(sep)>1){stop("Separator must have length 1")}

 if(is.null(var_value)){var_value="0"}

 res=.Call("heuristic_models_cpp", Data, var_path, var_conv, var_value, sep)
 
 return(as.data.frame(res)) 

}	


markov_model=function(Data, var_path, var_conv, var_value=NULL, var_null=NULL, order=1, nsim=NULL, max_step=NULL, out_more=FALSE, sep=">", seed=NULL){
 
 if(length(sep)>1){stop("Separator must have length 1")}

 if(is.null(var_value)){var_value="0"}
 if(is.null(var_null)){var_null="0"}
 if(is.null(nsim)){nsim=0}
 if(is.null(max_step)){max_step=0}
 if(!is.null(seed)){set.seed(seed)}

 res=.Call("markov_model_cpp", Data, var_path, var_conv, var_value, var_null, order, nsim, max_step, out_more, sep)
 
 if(out_more==FALSE){
  return(as.data.frame(res)) 
 }else{
  return(list(result=as.data.frame(res$result),transition_matrix=as.data.frame(res$transition_matrix),removal_effects=as.data.frame(res$removal_effects)))
 }

}	
 
