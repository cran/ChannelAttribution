
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


choose_order=function(Data, var_path, var_conv, var_null, max_order=10, sep=">", ncore=Inf, roc_npt=100, plot=TRUE){
 
 if(length(sep)>1){stop("Separator must have length 1")}

 if(ncore==Inf){
  ncore=max_order
 }else if(ncore<1){
  ncore=1
 }

 res=.Call("choose_order_cpp", Data, var_path, var_conv, var_null, max_order, sep, ncore, roc_npt)
 
 ck=res$auc$order[res$auc$order!=0]
 res$auc$order=res$auc$order[ck]
 res$auc$auc=res$auc$auc[ck]
 res$auc$pauc=res$auc$pauc[ck]
 
 best_order=res$auc$order[res$auc$pauc==max(res$auc$pauc)]
 
 if(best_order==max_order){
  print(paste0("Suggested order not found. Try increasing max_order."))
 }else{
  print(paste0("Suggested order: ", res$auc$order[res$auc$pauc==max(res$auc$pauc)]))
 }
 
 if(plot=="TRUE"){
  plot(res$auc$order,res$auc$pauc,type="l",xlab="order",ylab="penalized auc",main="PENALIZED AUC")
 }
 
 res[['suggested_order']]=best_order
 
 return(res)
 
}

markov_model_mp=function(Data, var_path, var_conv, var_value=NULL, var_null=NULL, order=1, nsim_start=1e5, max_step=NULL, out_more=FALSE, sep=">", ncore=Inf, nfold=10, seed=0, conv_par=0.05, rate_step_sim=1.5, verbose=TRUE){
 
 if(length(sep)>1){stop("Separator must have length 1")}

 if(is.null(var_value)){var_value="0"}
 if(is.null(var_null)){var_null="0"}
 if(is.null(max_step)){max_step=0}
 if(!is.null(seed)){set.seed(seed)}

 res=.Call("markov_model_mp_cpp", Data, var_path, var_conv, var_value, var_null, order, nsim_start, max_step, out_more, sep, ncore, nfold, seed, conv_par, rate_step_sim,verbose)
 
 if(out_more==FALSE){
  return(as.data.frame(res)) 
 }else{
  return(list(result=as.data.frame(res$result),transition_matrix=as.data.frame(res$transition_matrix),removal_effects=as.data.frame(res$removal_effects)))
 }

}
 
 
transition_matrix=function(Data, var_path, var_conv, var_null, order=1, sep=">", flg_equal=TRUE){
 
 if(length(sep)>1){stop("Separator must have length 1")}

 res=.Call("transition_matrix_cpp", Data, var_path, var_conv, var_null, order, sep, flg_equal)
 
 return(list(channels=data.frame(id=1:length(res$channels),channel_name=res$channels),transition_matrix=as.data.frame(res$transition_matrix)))
  
}
 
 
auto_markov_model=function(Data, var_path, var_conv, var_null, var_value=NULL, max_order=10, roc_npt=100, plot=FALSE, nsim_start=1e5, max_step=NULL, out_more=FALSE, sep=">", ncore=Inf, nfold=10, seed=0, conv_par=0.05, rate_step_sim=1.5, verbose=TRUE){
 
 if(length(sep)>1){stop("Separator must have length 1")}

 if(is.null(var_value)){var_value="0"}
 if(is.null(var_null)){var_null="0"}
 if(is.null(max_step)){max_step=0}
 if(!is.null(seed)){set.seed(seed)}
 
 order=choose_order(Data, var_path, var_conv, var_null, max_order=max_order, sep=sep, ncore=ncore, roc_npt=roc_npt, plot=plot)
 order=res[['suggested_order']]
 
 res=markov_model_mp(Data, var_path, var_conv, var_value=var_value, var_null=var_null, order=order, nsim_start=nsim_start, max_step=max_step, out_more=out_more, sep=sep, ncore=ncore, nfold=nfold, seed=seed, conv_par=conv_par, rate_step_sim=rate_step_sim, verbose=verbose)

 if(out_more==FALSE){
  return(as.data.frame(res)) 
 }else{
  return(list(result=as.data.frame(res$result),transition_matrix=as.data.frame(res$transition_matrix),removal_effects=as.data.frame(res$removal_effects)))
 }

}
