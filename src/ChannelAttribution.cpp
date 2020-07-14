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


#define __GXX_EXPERIMENTAL_CXX0X__ 1

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <numeric>
#include <thread>

#include <RcppArmadillo.h>

#define ARMA_USE_CXX11
#define ARMA_64BIT_WORD

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif
 
#ifndef END_RCPP
#define END_RCPP
#endif

using namespace std;
using namespace Rcpp;
//using namespace RcppThread;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

template <typename T>
string to_string(T pNumber)
{
 ostringstream oOStrStream;
 oOStrStream << pNumber;
 return oOStrStream.str();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP heuristic_models_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP sep_p)
{
	
 BEGIN_RCPP

 //inp.a
 
 List Data(Data_p);
 
 CharacterVector var_path_0(var_path_p);
 string var_path = Rcpp::as<string>(var_path_0);
 
 CharacterVector var_conv_0(var_conv_p);
 string var_conv = Rcpp::as<string>(var_conv_0);

 CharacterVector var_value_0(var_value_p); 
 string var_value = Rcpp::as<string>(var_value_0);
 
 CharacterVector sep_0(sep_p); 
 string sep = Rcpp::as<string>(sep_0);
 
 //inp.b
 
 bool flg_var_value;
 flg_var_value=0;
 if(var_value.compare("0")!=0){
  flg_var_value=1;
 } 
	
 CharacterVector vy0 = Data[var_path];
 vector<string> vy = Rcpp::as<vector<string> >(vy0);

 NumericVector vc0 = Data[var_conv];
 vector<unsigned long int> vc = Rcpp::as<vector<unsigned long int> >(vc0);

 vector<double> vv;
 if(flg_var_value==1){
  NumericVector vv0 = Data[var_value];
  vv = Rcpp::as<vector<double> >(vv0);
 }
 
 unsigned long int i,j,k,lvy,ssize;
 bool cfirst;
 unsigned long int start_pos,end_pos;
 unsigned long int nchannels;
 string s,channel,channel_first,channel_last;
  
 lvy=(unsigned long int) vy.size();
 nchannels=0;
 
 map<string,unsigned long int> mp_channels;
 vector<string> vchannels;
	  	
 map<string,double> mp_first_conv;
 map<string,double> mp_first_val;	
 map<string,double> mp_last_conv;
 map<string,double> mp_last_val;
 map<string,double> mp_linear_conv;
 map<string,double> mp_linear_val;
 map<string,double> mp0_linear_conv;
 map<string,double> mp0_linear_val;

 vector<string> vchannels_unique;
 double nchannels_unique;
 string kchannel;
 unsigned long int n_path_length;

 for(i=0;i<lvy;i++){
	 	 
  s=vy[i];
  
  s+=sep[0];
  ssize=(unsigned long int) s.size();
  channel="";
  j=0;
  nchannels_unique=0;
  vchannels_unique.clear();
  
  n_path_length=0;
  mp0_linear_conv.clear();
  mp0_linear_val.clear();
  
  start_pos=0;
  end_pos=0;  
  
  while(j<ssize){  
         
   cfirst=1;
   while(s[j]!=sep[0]){
	if(cfirst==0){   
     if(s[j]!=' '){
	  end_pos=j;	 
	 }
    }else if((cfirst==1) & (s[j]!=' ')){
	 cfirst=0;
	 start_pos=j;
	 end_pos=j;
	}
    ++j;     
   }
   
   if(cfirst==0){
    channel=s.substr(start_pos,(end_pos-start_pos+1));
   
    if(mp_channels.find(channel) == mp_channels.end()){
	 mp_channels[channel]=nchannels;
	 vchannels.push_back(channel);
	 ++nchannels;
	
     mp_first_conv[channel]=0;
	 mp_last_conv[channel]=0;
	 mp_linear_conv[channel]=0;
	 mp0_linear_conv[channel]=0;
	 
	 if(flg_var_value==1){
	  mp_first_val[channel]=0;	
	  mp_last_val[channel]=0;
	  mp_linear_val[channel]=0;
	  mp0_linear_val[channel]=0;
	 }		 
	 
	}
	 	 
    //lista canali unici
    if(nchannels_unique==0){
     vchannels_unique.push_back(channel);
	 ++nchannels_unique;
    }else if(find(vchannels_unique.begin(),vchannels_unique.end(),channel)==vchannels_unique.end()){
	 vchannels_unique.push_back(channel);
	 ++nchannels_unique;
    }

 	mp0_linear_conv[channel]=mp0_linear_conv[channel]+vc[i];
    if(flg_var_value==1){
	 mp0_linear_val[channel]=mp0_linear_val[channel]+vv[i]; 
    }
	++n_path_length;
   
    channel_last=channel;
  
   }//end cfirst
  
   channel="";
   ++j;
    
  }//end while j
   
  channel_first=vchannels_unique[0];
  mp_first_conv[channel_first]=mp_first_conv[channel_first]+vc[i];
 
  mp_last_conv[channel_last]=mp_last_conv[channel_last]+vc[i];
 
  //linear
  for(k=0;k<nchannels_unique;k++){
    kchannel=vchannels_unique[k];
    mp_linear_conv[kchannel]=mp_linear_conv[kchannel]+(mp0_linear_conv[kchannel]/n_path_length);
  }
  
  if(flg_var_value==1){
   mp_first_val[channel_first]=mp_first_val[channel_first]+vv[i];   
   mp_last_val[channel_last]=mp_last_val[channel_last]+vv[i];  
   for(k=0;k<nchannels_unique;k++){
    kchannel=vchannels_unique[k];
    mp_linear_val[kchannel]=mp_linear_val[kchannel]+(mp0_linear_val[kchannel]/n_path_length); 
   }
  }	  
 
  
 }//end for i
 
 vector<double> vfirst_conv(nchannels);
 vector<double> vlast_conv(nchannels);
 vector<double> vlinear_conv(nchannels); 
 
 vector<double> vfirst_val(nchannels);
 vector<double> vlast_val(nchannels);
 vector<double> vlinear_val(nchannels);
 
 for(k=0;k<nchannels;k++){
  kchannel=vchannels[k];	 
  vfirst_conv[k]=mp_first_conv[kchannel];
  vlast_conv[k]=mp_last_conv[kchannel];
  vlinear_conv[k]=mp_linear_conv[kchannel];
  
  if(flg_var_value==1){
   vfirst_val[k]=mp_first_val[kchannel];
   vlast_val[k]=mp_last_val[kchannel];
   vlinear_val[k]=mp_linear_val[kchannel];
  }
  
 }
 
 if(flg_var_value==1){
  return List::create(Named("channel_name")=vchannels, Named("first_touch_conversions") = vfirst_conv, Named("first_touch_value") = vfirst_val, Named("last_touch_conversions") = vlast_conv, Named("last_touch_value") = vlast_val, Named("linear_touch_conversions") = vlinear_conv, Named("linear_touch_value") = vlinear_val);
 }else{
  return List::create(Named("channel_name")=vchannels, Named("first_touch") = vfirst_conv, Named("last_touch") = vlast_conv, Named("linear_touch") = vlinear_conv);
 }
 
 END_RCPP

}//end heuristic_models_cpp


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Classe funzione di ripartizione

class Fx
{
 
 SpMat<unsigned long int> S;
 SpMat<unsigned long int> S0;
 SpMat<unsigned long int> S1;
 vector<unsigned long int> lrS0;
 vector<unsigned long int> lrS; 
 unsigned long int non_zeros,nrows,val0,lval0,i,j,k,s0,lrs0i;
 
 public:
  Fx(unsigned long int nrow0,unsigned long int ncol0): S(nrow0,ncol0), S0(nrow0,ncol0), S1(nrow0,ncol0), lrS0(nrow0,0), lrS(nrow0,0), non_zeros(0), nrows(nrow0) {}
  void init(unsigned long int, unsigned long int);
  void add(unsigned long int, unsigned long int,unsigned long int);
  void cum();
  unsigned long int sim(unsigned long int, double);
  double pconv(unsigned long int, unsigned long int);   
  List tran_matx(vector<string>);      
};  

void Fx::init(unsigned long int nrow1, unsigned long int ncol1)
{
  S.reset();
  S.set_size(nrow1,ncol1);
  
  S0.reset();
  S0.set_size(nrow1,ncol1);
  
  S1.reset();
  S1.set_size(nrow1,ncol1);  
  
  lrS0.clear();
  lrS0.resize(nrow1);
  
  lrS.clear();
  lrS.resize(nrow1);
  
  non_zeros=0;
  nrows=nrow1;
} 


void Fx::add(unsigned long int ichannel_old, unsigned long int ichannel, unsigned long int vxi)
{
    
  val0=S(ichannel_old,ichannel); //riempire f.p. transizione con vxi
  if(val0==0){
   lval0=lrS0[ichannel_old];
   S0(ichannel_old,lval0)=ichannel;
   lrS0[ichannel_old]=lval0+1;
   ++non_zeros;
  }
  S(ichannel_old,ichannel)=val0+vxi; 
  
} 

void Fx::cum()
{

 for(i=0;i<nrows;i++){
  lrs0i=lrS0[i];
  if(lrs0i>0){
   S1(i,0)=S(i,S0(i,0));
   for(j=1;j<lrs0i;j++){
    S1(i,j)=S1(i,j-1)+S(i,S0(i,j));
   }
   lrS[i]=S1(i,lrs0i-1);
  }   
 }    
  
}


unsigned long int Fx::sim(unsigned long int c, double uni) 
{
  
 s0=floor(uni*lrS[c]+1);
 
 for(k=0; k<lrS0[c]; k++){   
  if(S1(c,k)>=s0){return(S0(c,k));}
 }

 return 0;
    
}


List Fx::tran_matx(vector<string> vchannels) 
{

 unsigned long int mij,sm3;
 vector<string> vM1(non_zeros);
 vector<string> vM2(non_zeros);
 vector<double> vM3(non_zeros);
 vector<double> vsm;
 vector<unsigned long int> vk;
 
 k=0;
 for(i=0;i<nrows;i++){
  sm3=0;
  for(j=0;j<lrS0[i];j++){
   mij=S(i,S0(i,j));
   if(mij>0){   
      vM1[k]=vchannels[i];
      vM2[k]=vchannels[S0(i,j)];
      vM3[k]=mij;
      sm3=sm3+mij;
      ++k;
    }
  }
  
  vsm.push_back(sm3);
  vk.push_back(k);
  
 }//end for
 
 unsigned long int w=0;
 for(k=0;k<non_zeros;k++){
  if(k==vk[w]){++w;}
  vM3[k]=vM3[k]/vsm[w]; 
 }
  
 List res=List::create(Named("channel_from")=vM1, Named("channel_to") = vM2, Named("transition_probability") = vM3); 

 return(res);
 
}


double Fx::pconv(unsigned long int ichannel, unsigned long int nchannels)
{

 double res=0;
 for(k=(nchannels-2);k<nchannels;k++){ //considero solo il canale conversion ed il canale null
  res=res+(double) S(ichannel,k);
 }
 
 if(res>0){
  res=(double) S(ichannel,nchannels-2)/res;
 }else{
  res=0;	 
 }
 
 return(res);
  
} 



vector<long int> split_string(const string &s, unsigned long int order) {
    
	char delim=' ';
	vector<long int> result(order,-1);
    stringstream ss (s);
    string item;

	unsigned long int h=0;
    while (getline (ss, item, delim)) {
		result[h]=stoi(item);
		h=h+1;
    }
		
    return result;
}



// void print(auto &input)
// {
	// for (unsigned long int i = 0; i < (unsigned long int) input.size(); i++) {
		// std::cout << input.at(i) << ' ';
	// }
// }


vector<unsigned long int> bounds(unsigned long int parts, unsigned long int mem) {
    vector<unsigned long int>bnd;
    unsigned long int delta = mem / parts;
    unsigned long int reminder = mem % parts;
    unsigned long int N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (unsigned long int i = 0; i < parts; ++i) {
        N2 = N1 + delta;
        if (i == parts - 1)
            N2 += reminder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}


void W_choose_order_1(vector<string> vy, unsigned long int lvy, vector<unsigned long int> vc, vector<unsigned long int> vn, unsigned long int roc_npt, unsigned long int nchannels, vector<unsigned long int> &vorder, vector<double> &vuroc, vector<double> &vuroc_corr, List &L_roc, unsigned long int from_W, unsigned long int to_W)
{


 for(unsigned long int order = (from_W+1); order < (to_W+1); order++){
  
  string s,channel,path;
  unsigned long int nchannels_sim,i,ssize,ichannel_old,j,nc,start_pos,end_pos,start_pos_last,ichannel,npassi,vci,vni,vpi,k,h;  
  vector<string> vchannels_sim;
  map<string,unsigned long int> mp_channels_sim; 
  vector<string> vy2(lvy);  
  bool flg_next_path,flg_next_null;
  vector<double> vprev(lvy);
  double min_prev,max_prev,pauc,nnodes;
  
  vector<double> vth(roc_npt);
  unsigned long int tp,fn,tn,fp;
  double th,tpr,fpr,tpr_old,fpr_old,auc;

  vector<unsigned long int> vlastc(lvy);  
  
  //output
  
  vector<double> vtpr(roc_npt+1);
  vector<double> vfpr(roc_npt+1);
  
  nnodes=exp(lgamma(nchannels-3+order-1+1)-lgamma(order+1)-lgamma(nchannels-3-1+1));
  
  unsigned long int np=0;
  for(i=0;i<lvy; i++){
   np=np+vc[i]+vn[i];
  }
  
  if(nnodes<lvy){
  
   nchannels_sim=0;
    
   mp_channels_sim["(start)"]=0;
   vchannels_sim.push_back("(start)");
   ++nchannels_sim;
        
   for(i=0;i<lvy;i++){
   
    s=vy[i];
    ssize=(unsigned long int) s.size();
   
    channel="";
    ichannel_old=0;
    path="";
    j=0; 
    nc=0;
    start_pos=0;
    end_pos=0;
    start_pos_last=0; 
       
    while(j<ssize){ 
     
     while((j<ssize) & (nc<order)){
      
  	  if(s[j]==' '){
  	   nc=nc+1;
  	   if(nc==1){
  	    start_pos=j;
  	   }
  	  }else{
  	   end_pos=j;	
  	  }
  	
      ++j;
     
     } //while(s[j]!=' ') 
     
     channel=s.substr(start_pos_last,(end_pos-start_pos_last+1));
     
     if(mp_channels_sim.find(channel) == mp_channels_sim.end()){
      mp_channels_sim[channel]=nchannels_sim;
      vchannels_sim.push_back(channel);
      ++nchannels_sim;
     }
     
     ichannel=mp_channels_sim[channel];
     
     vlastc[i]=ichannel;
     
     if(ichannel!=ichannel_old){
      path+=to_string(ichannel);
      path+=" ";
     }   
     
     ichannel_old=ichannel;
   
     if((end_pos+1)==ssize){break;}; 
     j=start_pos+1;
     start_pos_last=j;   
     nc=0;
   
    }//end while(j<ssize)
   
    vy2[i]="0 "+path+"e";
  	
   }//end for i
   
   mp_channels_sim["(conversion)"]=nchannels_sim;
   vchannels_sim.push_back("(conversion)");    
   ++nchannels_sim;
    
   mp_channels_sim["(null)"]=nchannels_sim;
   vchannels_sim.push_back("(null)");     
   ++nchannels_sim;
      
   //Stima matrice transizione
     
   npassi=0;
   Fx S(nchannels_sim,nchannels_sim);  
   
   for(i=0;i<lvy;i++){
      
    flg_next_path=0;
    flg_next_null=0;
                             
    s=vy2[i];
    s+=" ";
    ssize= (unsigned long int) s.size();
   
    channel="";
    ichannel_old=0;
    ichannel=0;
   
    j=0;
    npassi=0;
   
    vci=vc[i];
    vni=vn[i];
         
    vpi=vci+vni;
    
    while((j<ssize) & (flg_next_path==0)){
       
     while(s[j]!=' '){
   
      if(j<ssize){
       channel+=s[j];
      }
      j=j+1;
     }
                
     if(channel[0]!='0'){//se non è il channel start
     
      if(channel[0]=='e'){ //stato finale
      
       ++npassi;
       
       if(vci>0){ //se ci sono conversion
        ichannel=nchannels_sim-2;
	    S.add(ichannel_old,ichannel,vci);
       
        if(vni>0){
         flg_next_null=1;
        }else{
         flg_next_path=1;
        }
               
       }
      
       if(((vni>0) | (flg_next_null==1)) & (flg_next_path==0)){ //se non ci sono conversion
        ichannel=nchannels_sim-1;
        S.add(ichannel_old,ichannel,vni);
      
        flg_next_path=1;
       }
      
      }else{ //stato non finale
       
       if(vpi>0){
        ichannel=atol(channel.c_str());
	    S.add(ichannel_old,ichannel,vpi);
	   }
     
      }
     
      if(flg_next_path==0){
       ++npassi;
      }
     }else{ //stato iniziale
    
      ichannel=0;
    
     }
   
     if(flg_next_path==0){
      ichannel_old=ichannel;
     }
      
     if(flg_next_path==0){
      channel="";
      j=j+1;   
     }
    
    }//end while j<size
       
   }//end for
           
   //fit
    
   for(i=0;i<lvy;i++){	   
    vprev[i]=S.pconv(vlastc[i],nchannels_sim);  
   } 
   
   min_prev=1;
   for(i=0;i<lvy;i++){	   
    min_prev=min(min_prev,vprev[i]);  
   } 
   
   max_prev=0;
   for(i=0;i<lvy;i++){	   
    max_prev=max(max_prev,vprev[i]);  
   } 
     
   for(k=0;k<roc_npt;k++){
    vth[k]=min_prev+(k*(max_prev-min_prev)/(roc_npt-1));	
   }
     
   auc=0; 
   tpr_old=0;
   fpr_old=0;
      
   vtpr[0]=0;
   vfpr[0]=0;  
   h=1;	
	 
   for(k=(roc_npt-1);k>=0 && k<roc_npt;k--){	   
     
    tp=0,fn=0,tn=0,fp=0;   
    th=vth[k];
    
    for(i=0;i<lvy;i++){
     
     if((vprev[i]>=th) & (vc[i]>0)){
  	 tp=tp+vc[i];   
     }else if((vprev[i]<th) & (vc[i]>0)){
      fn=fn+vc[i];
     }
     
     vni=vn[i];
	 
	 if((vprev[i]<th) & (vni>0)){
  	 tn=tn+vni;   
     }else if((vprev[i]>=th) & (vni>0)){
  	 fp=fp+vni;   
     }
       
    }
   
    tpr=(double)tp/(double)(tp+fn);
    fpr=(double)fp/(double)(fp+tn);
    
    auc=auc+((fpr-fpr_old)*tpr_old)+(((fpr-fpr_old)*(tpr-tpr_old))/2);
       
    vtpr[h]=tpr;
    vfpr[h]=fpr;
    
    tpr_old=tpr;
    fpr_old=fpr;
    
    h=h+1;
    
   }//end for k
   
   vtpr[roc_npt]=1;
   vfpr[roc_npt]=1; 
   auc=auc+((1-fpr_old)*tpr_old)+(((1-fpr_old)*(1-tpr_old))/2);
       
   pauc=(double)(1-((1-auc)*((np-1)/(np-nnodes-1))));
   if((pauc<0) | (pauc>1)){
    pauc=0;	  
   }
   
   vuroc[order-1]=auc;
   vuroc_corr[order-1]=pauc;
   vorder[order-1]=order;
   
   L_roc[to_string(order)]=Rcpp::List::create(Rcpp::Named("fpr")=vfpr,Rcpp::Named("tpr")=vtpr);
  
  }//end if(nnodes<lvy)
 
 }
  
}


RcppExport SEXP choose_order_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_null_p, SEXP max_order_p, SEXP sep_p, SEXP ncore_p, SEXP roc_npt_p)
{
 
 BEGIN_RCPP

 //inp.a
  
 List Data(Data_p);
 
 CharacterVector var_path_0(var_path_p);
 string var_path = Rcpp::as<string>(var_path_0);
 
 CharacterVector var_conv_0(var_conv_p);
 string var_conv = Rcpp::as<string>(var_conv_0);
  
 CharacterVector var_null_0(var_null_p); 
 string var_null = Rcpp::as<string>(var_null_0);
 
 NumericVector max_order_0(max_order_p); 
 unsigned long int max_order = Rcpp::as<unsigned long int>(max_order_0);
  
 CharacterVector sep_0(sep_p); 
 string sep = Rcpp::as<string>(sep_0);

 NumericVector ncore_0(ncore_p); 
 unsigned long int ncore = Rcpp::as<unsigned long int>(ncore_0);
 
 
 NumericVector roc_npt_0(roc_npt_p); 
 unsigned long int roc_npt = Rcpp::as<unsigned long int>(roc_npt_0);
 
 //inp.b 
   
 CharacterVector vy0 = Data[var_path];
 vector<string> vy = Rcpp::as<vector<string> >(vy0);
  
 NumericVector vc0 = Data[var_conv];
 vector<unsigned long int> vc = Rcpp::as<vector<unsigned long int> >(vc0);
  
 vector<unsigned long int> vn;
 NumericVector vn0 = Data[var_null];
 vn = Rcpp::as<vector<unsigned long int> >(vn0);
   
 unsigned long int i,j,lvy,ssize;
 unsigned long int nchannels,npassi;
 bool cfirst;
 unsigned long int start_pos,end_pos;
 string s,channel,path;
 map<string,unsigned long int> mp_channels,mp_channels_sim;
 map<unsigned long int,unsigned long int> mp_npassi;
 vector<unsigned long int> vnpassi;
    
 lvy=(unsigned long int) vy.size();
   
 //////////////////////
 //CODIFICA DA ONE STEP 
 //////////////////////
      
 string channel_test="";
 string channel_old;
 
 //unsigned long int order;
 
 nchannels=0;
   
 mp_channels["(start)"]=0;
 vector<string> vchannels;
 vchannels.push_back("(start)");     
 ++nchannels;
 
 //ricodifica nomi canale in interi
   
 for(i=0;i<lvy;i++){
	        
  s=vy[i];
  s+=sep[0];
  ssize=(unsigned long int) s.size();
  channel="";
  path="";
  j=0;
  npassi=0;
  start_pos=0;
  end_pos=0;
   
  //medium.touch
  
  while(j<ssize){  
         
   cfirst=1;
   while(s[j]!=sep[0]){
    if(cfirst==0){   
     if(s[j]!=' '){
      end_pos=j;     
     }
    }else if((cfirst==1) & (s[j]!=' ')){
     cfirst=0;
     start_pos=j;
     end_pos=j;
    }
    ++j;     
   }
   
   if(cfirst==0){
    channel=s.substr(start_pos,(end_pos-start_pos+1));
   
    if(mp_channels.find(channel) == mp_channels.end()){
     mp_channels[channel]=nchannels;
     vchannels.push_back(channel);
     ++nchannels;
    }
			 
    if(npassi==0){
	 path="";
    }else{
     path+=" ";
    }
     
    path+=to_string(mp_channels[channel]);
    ++npassi;      
           
   }//if end_pos
   
   channel="";
   ++j;
   
  }//end while channel
    
  vy[i]=path;
  ++npassi;
  
 }//end for
  
 mp_channels["(conversion)"]=nchannels;
 vchannels.push_back("(conversion)");    
 ++nchannels;
  
 mp_channels["(null)"]=nchannels;
 vchannels.push_back("(null)");     
 ++nchannels;
 
 //riconduco a order=1
 
 List L_roc, L_uroc;
 vector<double> vuroc(max_order);
 vector<double> vuroc_corr(max_order);
 vector<unsigned long int> vorder(max_order);
 
 vector<unsigned long int> limits = bounds(ncore, max_order);
  
 if(ncore==1){
  
  W_choose_order_1(vy, lvy, vc, vn, roc_npt, nchannels, ref(vorder), ref(vuroc), ref(vuroc_corr), ref(L_roc), limits[0], limits[1]);	  
 
 }else{
  
  vector<thread> threads(ncore);

  //Launch ncore threads:	
  for(unsigned long int td=0; td<ncore; td++){
   threads[td]=thread(W_choose_order_1, vy, lvy, vc, vn, roc_npt, nchannels, ref(vorder), ref(vuroc), ref(vuroc_corr), ref(L_roc), limits[td], limits[td+1]);
  }

  //Join the threads with the main thread
  for(auto &t : threads){
   t.join();
  }  
  
 }
 
 L_uroc=Rcpp::List::create(Rcpp::Named("order")=vorder,Rcpp::Named("auc")=vuroc,Rcpp::Named("pauc")=vuroc_corr);
  
 return(List::create(Named("roc")=L_roc, Named("auc")=L_uroc));
 
 END_RCPP
 
} 	

void W_markov_model_1(unsigned long int seed, vector<double> v_vui, unsigned long int lvy, vector<unsigned long int> vc, unsigned long int nch0, vector<double> vv, unsigned long int nfold, unsigned long int nsim_start,  map< unsigned long int, vector<long int> > mp_channels_sim_inv, unsigned long int max_npassi, unsigned long int nchannels_sim, unsigned long int order, bool flg_var_value, unsigned long int nchannels, Fx S, Fx fV, vector<unsigned long int> &nconv, vector<double> &ssval, vector< vector<double> > &T, vector< vector<double> > &TV, vector< vector<double> > &V, vector< vector<double> > &VV, vector<double> &v_inc_path, unsigned long int from_W, unsigned long int to_W)
{
 
 long int id0;
 unsigned long int i0,c,npassi0,k0,c_last=0,n_inc_path;
 vector<bool> C(nchannels,0);
 double sval0=0,sn,sm;
 bool flg_exit;
 
 for(unsigned long int run = from_W; run < to_W; run++){
	 
   mt19937 generator(seed+run);  
   uniform_real_distribution<double> distribution(0,1);
   
   n_inc_path=0;   

   for(i0=0; i0<(unsigned long int) (nsim_start/nfold); i0++){
	          	        
     c=0;
     npassi0=0;
 	    
     for(k0=0; k0<nchannels; k0++){ //svuoto il vettore del flag canali visitati
      C[k0]=0;
     }
    
     C[c]=1; //assegno 1 al channel start
     
 	 flg_exit=0;
 	
     while((npassi0<=max_npassi) & (flg_exit==0)){ //interrompo quando raggiungo il massimo numero di passi
      
      c=S.sim(c,distribution(generator));
     
      if(c==nchannels_sim-2){ //se ho raggiunto lo stato conversion interrompo
       flg_exit=1;
      }else if(c==nchannels_sim-1){ //se ho raggiunto lo stato null interrompo
       flg_exit=1;
      }
      
      if(flg_exit==0){
        if(order==1){
          C[c]=1; //flaggo con 1 il canale visitato   
        }else{       
         for(k0=0; k0<order; k0++){
		   id0=mp_channels_sim_inv[c][k0];
           if(id0>=0){
            C[id0]=1;
           }else{
            break;     
           }
          }
        }
         
        c_last=c; //salvo il canale visitato
        ++npassi0;
      }
    	
     }//end while npassi0 
 	
     
     if(c==nchannels_sim-2){ //solo se ho raggiunto la conversion assegno +1 ai canali interessati (se ho raggiunto il max numero di passi è come se fossi andato a null)
      
     nconv[run]=nconv[run]+1;//incremento le conversion
      
     //genero per il canale c_last un valore di conversion "sval0"
     if(flg_var_value==1){
      sval0=v_vui[fV.sim(c_last,distribution(generator))];
     }   
         
     ssval[run]=ssval[run]+sval0;
       
     for (k0=0; k0<nchannels; k0++){
       if(C[k0]==1){
         T[run][k0]=T[run][k0]+1;
        if(flg_var_value==1){
          V[run][k0]=V[run][k0]+sval0;
        }
       }
      }
    
     }else if(c!=nchannels_sim-1){ //se non ho raggiunto neanche lo stato NULL
	   
       n_inc_path=n_inc_path+1;	   
		 	 
	 }//end if conv
           
   }//end for i
   
   v_inc_path[run]=(double) n_inc_path/ (double) i0;
       
   T[run][0]=0; //pongo channel start = 0
   T[run][nchannels-2]=0; //pongo channel conversion = 0 
   T[run][nchannels-1]=0; //pongo channel null = 0 
   
   sn=0; 
   for(k0=0;k0<lvy; k0++){
    sn=sn+vc[k0];
   }
    
   sm=0;
   for(k0=0;k0<nchannels-1; k0++){
    sm=sm+T[run][k0];
   }
   
   for (k0=1; k0<(nch0+1); k0++){
    if(sm>0){
     TV[run][k0-1]=(T[run][k0]/sm)*sn;
    }
   }
     
   if(flg_var_value==1){
  
    V[run][0]=0; //pongo channel start = 0
    V[run][nchannels-2]=0; //pongo channel conversion = 0 
    V[run][nchannels-1]=0; //pongo channel null = 0 
      
    sn=0;
    for(k0=0;k0<lvy; k0++){
     sn=sn+vv[k0];
    }
    
    sm=0;
    for(k0=0;k0<nchannels-1; k0++){
     sm=sm+V[run][k0];
    }
      
    for(k0=1; k0<(nch0+1); k0++){
     if(sm>0){
      VV[run][k0-1]=(V[run][k0]/sm)*sn;
     }
    }
     
   }
		
 }//end for run

}

string f_print_perc(double num){ 
 
 string res;
 if(num>=1){
  res=to_string((double)(floor(num*10000)/100)).substr(0,6);    
 }else if(num>=0.1){ 
  res=to_string((double)(floor(num*10000)/100)).substr(0,5); 
 }else{
  res=to_string((double)(floor(num*10000)/100)).substr(0,4);    	   
 } 
 return(res);
}


RcppExport SEXP markov_model_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_value_p, SEXP var_null_p, SEXP order_p, SEXP nsim_start_p, SEXP max_step_p, SEXP out_more_p, SEXP sep_p, SEXP ncore_p, SEXP nfold_p, SEXP seed_p, SEXP conv_par_p, SEXP rate_step_sim_p, SEXP verbose_p)
{ 
	 	
 BEGIN_RCPP

 //inp.a
  
 List Data(Data_p);
 
 CharacterVector var_path_0(var_path_p);
 string var_path = Rcpp::as<string>(var_path_0);
 
 CharacterVector var_conv_0(var_conv_p);
 string var_conv = Rcpp::as<string>(var_conv_0);
 
 CharacterVector var_value_0(var_value_p); 
 string var_value = Rcpp::as<string>(var_value_0);
 
 CharacterVector var_null_0(var_null_p); 
 string var_null = Rcpp::as<string>(var_null_0);
 
 NumericVector order_0(order_p); 
 unsigned long int order = Rcpp::as<unsigned long int>(order_0);
 
 NumericVector nsim_start_0(nsim_start_p); 
 unsigned long int nsim_start = Rcpp::as<unsigned long int>(nsim_start_0);
 
 NumericVector max_step_0(max_step_p); 
 unsigned long int max_step = Rcpp::as<unsigned long int>(max_step_0); 
 
 NumericVector out_more_0(out_more_p); 
 unsigned long int out_more = Rcpp::as<unsigned long int>(out_more_0);
 
 CharacterVector sep_0(sep_p); 
 string sep = Rcpp::as<string>(sep_0);
 
 NumericVector ncore_0(ncore_p); 
 unsigned long int ncore = Rcpp::as<unsigned long int>(ncore_0);

 NumericVector nfold_0(nfold_p); 
 unsigned long int nfold = Rcpp::as<unsigned long int>(nfold_0); 
 
 NumericVector seed_0(seed_p); 
 unsigned long int seed = Rcpp::as<unsigned long int>(seed_0);  

 NumericVector conv_par_0(conv_par_p); 
 double conv_par = Rcpp::as<double>(conv_par_0);  
 
 NumericVector rate_step_sim_0(rate_step_sim_p); 
 double rate_step_sim = Rcpp::as<double>(rate_step_sim_0);   
 
 NumericVector verbose_0(verbose_p); 
 unsigned long int verbose = Rcpp::as<unsigned long int>(verbose_0);

 //inp.b 
  
 bool flg_var_value;
 flg_var_value=0;
 if(var_value.compare("0")!=0){
  flg_var_value=1;
 }
 
 bool flg_var_null;
 flg_var_null=0;
 if(var_null.compare("0")!=0){
  flg_var_null=1;
 }
 
 CharacterVector vy0 = Data[var_path];
 vector<string> vy = Rcpp::as<vector<string> >(vy0);
  
 NumericVector vc0 = Data[var_conv];
 vector<unsigned long int> vc = Rcpp::as<vector<unsigned long int> >(vc0);
 
 vector<double> vv;
 if(flg_var_value==1){
  NumericVector vv0 = Data[var_value];
  vv = Rcpp::as<vector<double> >(vv0);
 }
 
 vector<unsigned long int> vn;
 if(flg_var_null==1){
  NumericVector vn0 = Data[var_null];
  vn = Rcpp::as<vector<unsigned long int> >(vn0);
 } 
  
 unsigned long int i,j,k,lvy,ssize;
 unsigned long int nchannels,nchannels_sim,npassi;
 bool cfirst;
 unsigned long int start_pos,end_pos;
 string s,channel,path;
 map<string,unsigned long int> mp_channels,mp_channels_sim;
 map< unsigned long int, vector<long int> > mp_channels_sim_inv;
 map<unsigned long int,unsigned long int> mp_npassi;
 vector<unsigned long int> vnpassi;
    
 lvy=(unsigned long int) vy.size();
   
 //////////////////////
 //CODIFICA DA ONE STEP 
 //////////////////////
    
 //mappa dei conversion value
 unsigned long int l_vui=0;
 map<double,unsigned long int> mp_vui;
 vector<double> v_vui;
 double vui;

 vector<string> rchannels;
 string channel_j;
  
 nchannels=0;
 nchannels_sim=0;
 
 vector<string> vy2(lvy);
 
 mp_channels["(start)"]=0;
 vector<string> vchannels;
 vchannels.push_back("(start)");     
 ++nchannels;

 vector<string> vchannels_sim;

 //definizione mappa conversion value
 if(flg_var_value==1){
  for(i=0;i<lvy;i++){
   if(vc[i]>0){
    vui=vv[i]/vc[i];
    if(mp_vui.find(vui)==mp_vui.end()){
     mp_vui[vui]=l_vui;
     v_vui.push_back(vui);
     ++l_vui;    
    }
   }
  }
 }
 
 //ricodifica nomi canale in interi
   
 for(i=0;i<lvy;i++){
    
  s=vy[i];
  s+=sep[0];
  ssize=(unsigned long int) s.size();
  channel="";
  path="";
  j=0;
  npassi=0;
  start_pos=0;
  end_pos=0;
     
  //medium.touch
  
  while(j<ssize){  
         
   cfirst=1;
   while(s[j]!=sep[0]){ 
    if(cfirst==0){   
     if(s[j]!=' '){
      end_pos=j;     
     }
    }else if((cfirst==1) & (s[j]!=' ')){
     cfirst=0;
     start_pos=j;
     end_pos=j;
    }
    ++j;     
   }
   
   if(cfirst==0){
    channel=s.substr(start_pos,(end_pos-start_pos+1));
   
    if(mp_channels.find(channel) == mp_channels.end()){
     mp_channels[channel]=nchannels;
     vchannels.push_back(channel);
     ++nchannels;
    }
			 
    if(npassi==0){
	 path="";
    }else{
     path+=" ";
    }
     
    path+=to_string(mp_channels[channel]);
    ++npassi;      
           
   }//if end_pos
   
   channel="";
   ++j;
   
  }//end while channel
    
  vy[i]=path;
  ++npassi;
 
 }//end for

 mp_channels["(conversion)"]=nchannels;
 vchannels.push_back("(conversion)");    
 ++nchannels;
  
 mp_channels["(null)"]=nchannels;
 vchannels.push_back("(null)");     
 ++nchannels;
 
 //riconduco a order=1
    	
 nchannels_sim=0;
   
 mp_channels_sim["(start)"]=0;
 vchannels_sim.push_back("(start)");
 ++nchannels_sim;
   
 unsigned long int ichannel,ichannel_old,vpi,vci,vni; 
 bool flg_next_path,flg_next_null;
 unsigned long int start_pos_last,nc;
 
 vector<long int> vtmp(order);
   
 for(i=0;i<lvy;i++){
  
   s=vy[i];
   ssize=(unsigned long int) s.size();
  
   channel="";
   ichannel_old=0;
   path="";
   j=0; 
   nc=0;
   start_pos=0;
   end_pos=0;
   start_pos_last=0; 
      
   while(j<ssize){ 
    
    while((j<ssize) & (nc<order)){
     
 	if(s[j]==' '){
 	 nc=nc+1;
 	 if(nc==1){
 	  start_pos=j;
 	 }
 	}else{
 	 end_pos=j;	
 	}
 	
     ++j;
    
    } //while(s[j]!=' ') 
    
    channel=s.substr(start_pos_last,(end_pos-start_pos_last+1));
    	
    if(mp_channels_sim.find(channel) == mp_channels_sim.end()){
     mp_channels_sim[channel]=nchannels_sim;
	 
	 vtmp=split_string(channel,order);
	 mp_channels_sim_inv[nchannels_sim]=vtmp;
    
	 vchannels_sim.push_back(channel);
     ++nchannels_sim;
    }
    
    ichannel=mp_channels_sim[channel];
        
    if(ichannel!=ichannel_old){
     path+=to_string(ichannel);
     path+=" ";
    }   
    
    ichannel_old=ichannel;

    if((end_pos+1)==ssize){break;}; 
    j=start_pos+1;
    start_pos_last=j;   
    nc=0;
  
   }//end while(j<ssize)
  
   vy2[i]="0 "+path+"e";
 	
 }//end for i
  
 mp_channels_sim["(conversion)"]=nchannels_sim;
 vchannels_sim.push_back("(conversion)");    
 ++nchannels_sim;
   
 mp_channels_sim["(null)"]=nchannels_sim;
 vchannels_sim.push_back("(null)");     
 ++nchannels_sim;
    
 /////////////////////////////////////////////////////
 //CREAZIONE DELLE MATRICI FUNZIONALI ALLE SIMULAZIONI
 ////////////////////////////////////////////////////
  
 string channel_old;
  
 Fx S(nchannels_sim,nchannels_sim);
  
 Fx fV(nchannels_sim,l_vui);
 
 unsigned long int max_npassi; 
 max_npassi=0; 
 
 for(i=0;i<lvy;i++){
     
  flg_next_path=0;
  flg_next_null=0;
                            
  s=vy2[i];
  s+=" ";
  ssize= (unsigned long int) s.size();
  
  channel="";
  ichannel_old=0;
  ichannel=0;
 
  j=0;
  npassi=0;
  
  vci=vc[i];
  if(flg_var_null==1){
   vni=vn[i];
  }else{
   vni=0;
  }      
  vpi=vci+vni;
   
  while((j<ssize) & (flg_next_path==0)){
      
   while(s[j]!=' '){
  
    if(j<ssize){
     channel+=s[j];
    }
    j=j+1;
   }
               
   if(channel[0]!='0'){//se non è il channel start
   
    if(channel[0]=='e'){ //stato finale
    
     ++npassi;
     
     if(vci>0){ //se ci sono conversion
      ichannel=nchannels_sim-2;
      S.add(ichannel_old,ichannel,vci);
    
      if(flg_var_value==1){
       vui=vv[i]/vci;
	   fV.add(ichannel_old,mp_vui[vui],vci);
      }

      if(vni>0){
       flg_next_null=1;  
      }else{
       flg_next_path=1;
      }
             
     }
    
     if(((vni>0) | (flg_next_null==1)) & (flg_next_path==0)){ //se non ci sono conversion
      ichannel=nchannels_sim-1;
      S.add(ichannel_old,ichannel,vni);
      flg_next_path=1;
     }
    
    }else{ //stato non finale
     
     if(vpi>0){
      ichannel=atol(channel.c_str());
      S.add(ichannel_old,ichannel,vpi);
     }
   
    }
   
    if(flg_next_path==0){
     ++npassi;
    }
   }else{ //stato iniziale
   
    ichannel=0;
   
   }
  
   if(flg_next_path==0){
    ichannel_old=ichannel;
   }
     
   if(flg_next_path==0){
    channel="";
    j=j+1;   
   }
   
  }//end while j<size
   
  max_npassi=max(max_npassi,npassi);
   
 }//end for 
 
 //out matrice di transizione
 
 List res_mtx; 
 
 if(out_more==1){
   res_mtx=S.tran_matx(vchannels_sim);
 }
   
 //f.r. transizione
 S.cum(); 
    
 //f.r. conversion value
 if(flg_var_value==1){
  fV.cum();
 }
    
 //SIMULAZIONI
     
 if(max_step>0){
  max_npassi=max_step;
 }
   
 if(nsim_start==0){
  nsim_start=1e5;
 }
       
 unsigned long int run;
 unsigned long int nch0; 
 double sn=0;
  
 vector< vector<double> > T(nfold,vector<double>(nchannels,0));
 vector< vector<double> > V(nfold,vector<double>(nchannels,0));
 
 nch0=nchannels-3;

 vector< vector<double> > TV(nfold,vector<double>(nch0,0));
 vector<double> rTV(nch0,0);
 
 vector< vector<double> > VV(nfold,vector<double>(nch0,0));
 vector<double> rVV(nch0,0);
 
 vector<unsigned long int> nconv(nfold,0);
 vector<double> ssval(nfold,0);
 
 vector<double> TV_fin(nch0);
 vector<double> VV_fin(nch0);
 vector<double> vtmp1(nfold);
 vector<double> v_res_conv(nfold);
 
 vector<double> v_inc_path(nfold);
    
 double max_res_conv=numeric_limits<double>::infinity();
 double min_res_conv;
 double res_conv;
 
 unsigned long int id_mz,run_min_res_conv=0;
  
 while(max_res_conv>conv_par){
	
  min_res_conv=numeric_limits<double>::infinity(); 
	
  vector<unsigned long int> limits = bounds(ncore, nfold);

  if(ncore==1){
	  
	W_markov_model_1(seed, v_vui, lvy, vc, nch0, vv, nfold, nsim_start, mp_channels_sim_inv, max_npassi, nchannels_sim, order, flg_var_value, nchannels, S, fV, ref(nconv), ref(ssval), ref(T), ref(TV), ref(V), ref(VV), ref(v_inc_path), limits[0], limits[1]);
	  
  }else{
  
   vector<thread> threads(ncore);
      
   //Launch ncore threads:	  
   for(unsigned long int td=0; td<ncore; td++){
    threads[td]=thread(W_markov_model_1, seed, v_vui, lvy, vc, nch0, vv, nfold, nsim_start, mp_channels_sim_inv, max_npassi, nchannels_sim, order, flg_var_value, nchannels, S, fV, ref(nconv), ref(ssval), ref(T), ref(TV), ref(V), ref(VV), ref(v_inc_path), limits[td], limits[td+1]);
   }
      
   //Join the threads with the main thread
   for(auto &t : threads){
  	t.join();
   } 
	    
  }
  
  sn=0; 
  for(k=0;k<lvy; k++){
   sn=sn+vc[k];
  }
   
  for(k=0; k<nch0; k++){ 
   for(run=0; run<nfold; run++){
    vtmp1[run]=TV[run][k];
   }
   sort(vtmp1.begin(), vtmp1.end(), std::greater<double>());
   id_mz=(unsigned long int) (nfold/2);
   if(nfold % 2 == 0){
	TV_fin[k]=(vtmp1[id_mz-1]+vtmp1[id_mz])/2;   
   }else{
	TV_fin[k]=vtmp1[id_mz];   
   }
  }
   
  if(flg_var_value==1){ 

   for(k=0; k<nch0; k++){ 
    for(run=0; run<nfold; run++){
     vtmp1[run]=VV[run][k];
    }
    sort(vtmp1.begin(), vtmp1.end(), std::greater<double>());
    id_mz=(unsigned long int) (nfold/2);
    if(nfold % 2 == 0){
	 VV_fin[k]=(vtmp1[id_mz-1]+vtmp1[id_mz])/2;   
    }else{
	 VV_fin[k]=vtmp1[id_mz];   
    }
   }  
  
  } 
     
  max_res_conv=0;
  for(run=0; run<nfold; run++){
	  
   res_conv=0;
   for(k=0; k<nch0; k++){
    res_conv=res_conv+abs(TV[run][k]-TV_fin[k]);
   }
   res_conv=res_conv/sn;
   if(res_conv>max_res_conv){
    max_res_conv=res_conv;
   }
   if(res_conv<min_res_conv){
    min_res_conv=res_conv;
	run_min_res_conv=run;
   }
   
  }
      
  if(verbose==1){
   if(max_res_conv>conv_par){
    Rcout <<  "Number of simulations: "+ to_string(nsim_start) + " - Reaching convergence (wait...): " + f_print_perc(max_res_conv) + "% > " + f_print_perc(conv_par) + "%" << endl;
   }else{
	Rcout << endl;
    Rcout << "Number of simulations: "+ to_string(nsim_start) + " - Convergence reached: " + f_print_perc(max_res_conv) + "% < " + f_print_perc(conv_par) + "%" << endl;  
   }
  }
  
  nsim_start=nsim_start*rate_step_sim;

 }//end while(res_conv>conv_par)
     
 vector<string> vchannels0(nch0);
 for(k=1; k<(nch0+1); k++){
  vchannels0[k-1]=vchannels[k];
 }
  
 double succ_path=0;
 for(k=0; k<nfold; k++){
  succ_path=succ_path+(1-v_inc_path[k])/nfold;
 }	  
 if(verbose==1){ 
  Rcout << endl;
  Rcout << "Percentage of simulated paths that successfully end before maximum number of steps (" + to_string(max_npassi) + ") is reached: " << f_print_perc(succ_path) + "%" << endl;
  Rcout << endl;
 }
 
 if(flg_var_value==1){ 
 
  if(out_more==0){
  
   return List::create(Named("channel_name")=vchannels0, Named("total_conversion") = TV[run_min_res_conv], Named("total_conversion_value") = VV[run_min_res_conv] );
  
  }else{
	
   //removal effects conversion	
   for(k=1; k<nch0; k++){
    rTV[k-1]=T[run_min_res_conv][k]/nconv[run_min_res_conv];
   } 

   //removal effects conversion value	
   for(k=1; k<(nch0+1); k++){
    rVV[k-1]=V[run_min_res_conv][k]/ssval[run_min_res_conv];
   } 
      
   // List res1=List::create(Named("channel_name")=vchannels0, Named("total_conversions") = TV[run_min_res_conv], Named("total_conversion_value") = VV[run_min_res_conv], Named("total_conversions_run") = TV, Named("total_conversion_value_run") = VV );
   List res1=List::create(Named("channel_name")=vchannels0, Named("total_conversions") = TV[run_min_res_conv], Named("total_conversion_value") = VV[run_min_res_conv]);
   List res3=List::create(Named("channel_name")=vchannels0, Named("removal_effects_conversion") = rTV, Named("removal_effects_conversion_value") = rVV);
   return List::create(Named("result") = res1, Named("transition_matrix")=res_mtx, Named("removal_effects") = res3);
     
  } 
 
 }else{
     
  if(out_more==0){
  
   return List::create(Named("channel_name")=vchannels0, Named("total_conversions") = TV[run_min_res_conv]);
  
  }else{
	  
   //removal effects conversion	
   for(k=1; k<(nch0+1); k++){
    rTV[k-1]=T[run_min_res_conv][k]/nconv[run_min_res_conv];
   } 
   
   // List res1=List::create(Named("channel_name")=vchannels0, Named("total_conversions") = TV[run_min_res_conv], Named("total_conversions_run") = TV);
   List res1=List::create(Named("channel_name")=vchannels0, Named("total_conversions") = TV[run_min_res_conv]);
   List res3=List::create(Named("channel_name")=vchannels0, Named("removal_effects") = rTV);
   return List::create(Named("result") = res1, Named("transition_matrix")=res_mtx, Named("removal_effects") = res3);
  
  } 
    
 }
  
 END_RCPP

}	


RcppExport SEXP transition_matrix_cpp(SEXP Data_p, SEXP var_path_p, SEXP var_conv_p, SEXP var_null_p, SEXP order_p, SEXP sep_p, SEXP flg_equal_p)
{
 
 BEGIN_RCPP

 //inp.a
  
 List Data(Data_p);
 
 CharacterVector var_path_0(var_path_p);
 string var_path = Rcpp::as<string>(var_path_0);
 
 CharacterVector var_conv_0(var_conv_p);
 string var_conv = Rcpp::as<string>(var_conv_0);
  
 CharacterVector var_null_0(var_null_p); 
 string var_null = Rcpp::as<string>(var_null_0);
 
 NumericVector order_0(order_p); 
 unsigned long int order = Rcpp::as<unsigned long int>(order_0);
  
 CharacterVector sep_0(sep_p); 
 string sep = Rcpp::as<string>(sep_0);

 NumericVector flg_equal_0(flg_equal_p); 
 unsigned long int flg_equal = Rcpp::as<unsigned long int>(flg_equal_0);
 
 //inp.b 
   
 bool flg_var_null;
 flg_var_null=0;
 if(var_null.compare("0")!=0){
  flg_var_null=1;
 }
 
 CharacterVector vy0 = Data[var_path];
 vector<string> vy = Rcpp::as<vector<string> >(vy0);
  
 NumericVector vc0 = Data[var_conv];
 vector<unsigned long int> vc = Rcpp::as<vector<unsigned long int> >(vc0);
  
 vector<unsigned long int> vn;
 if(flg_var_null==1){
  NumericVector vn0 = Data[var_null];
  vn = Rcpp::as<vector<unsigned long int> >(vn0);
 } 
  
 unsigned long int i,j,lvy,ssize;
 unsigned long int nchannels,nchannels_sim,npassi;
 bool cfirst;
 unsigned long int start_pos,end_pos;
 string s,channel,path;
 map<string,unsigned long int> mp_channels,mp_channels_sim;
 map< unsigned long int, vector<long int> > mp_channels_sim_inv;
 map<unsigned long int,unsigned long int> mp_npassi;
 vector<unsigned long int> vnpassi;
    
 lvy=(unsigned long int) vy.size();
   
 //////////////////////
 //CODIFICA DA ONE STEP 
 //////////////////////
    
 //mappa dei conversion value
 map<double,unsigned long int> mp_vui;
 vector<double> v_vui;

 vector<string> rchannels;
 string channel_j;
  
 nchannels=0;
 nchannels_sim=0;
 
 vector<string> vy2(lvy);
 
 mp_channels["(start)"]=0;
 vector<string> vchannels;
 vchannels.push_back("(start)");     
 ++nchannels;

 vector<string> vchannels_sim;
 
 //ricodifica nomi canale in interi
  
 for(i=0;i<lvy;i++){
       
  s=vy[i];
  s+=sep[0];
  ssize=(unsigned long int) s.size();
  channel="";
  path="";
  j=0;
  npassi=0;
  start_pos=0;
  end_pos=0;
   
  //medium.touch
  
  while(j<ssize){  
         
   cfirst=1;
   while(s[j]!=sep[0]){ 
    if(cfirst==0){   
     if(s[j]!=' '){
      end_pos=j;     
     }
    }else if((cfirst==1) & (s[j]!=' ')){
     cfirst=0;
     start_pos=j;
     end_pos=j;
    }
    ++j;     
   }
   
   if(cfirst==0){
    channel=s.substr(start_pos,(end_pos-start_pos+1));
   
    if(mp_channels.find(channel) == mp_channels.end()){
     mp_channels[channel]=nchannels;
     vchannels.push_back(channel);
     ++nchannels;
    }
			 
    if(npassi==0){
	 path="";
    }else{
     path+=" ";
    }
     
    path+=to_string(mp_channels[channel]);
    ++npassi;      
           
   }//if end_pos
   
   channel="";
   ++j;
   
  }//end while channel
    
  vy[i]=path;
  ++npassi;
 
 }//end for

 mp_channels["(conversion)"]=nchannels;
 vchannels.push_back("(conversion)");    
 ++nchannels;
  
 mp_channels["(null)"]=nchannels;
 vchannels.push_back("(null)");     
 ++nchannels;
 
 //riconduco a order=1
    
 nchannels_sim=0;
   
 mp_channels_sim["(start)"]=0;
 vchannels_sim.push_back("(start)");
 ++nchannels_sim;
   
 unsigned long int ichannel,ichannel_old,vpi,vci,vni; 
 bool flg_next_path,flg_next_null;
 unsigned long int start_pos_last,nc;
 
 vector<long int> vtmp(order);

 for(i=0;i<lvy;i++){ 
   
   s=vy[i];
   
   ssize=(unsigned long int) s.size();
  
   channel="";
   ichannel_old=0;
   path="";
   j=0; 
   nc=0;
   start_pos=0;
   end_pos=0;
   start_pos_last=0; 
      
   while(j<ssize){ 
    
    while((j<ssize) & (nc<order)){
     
 	if(s[j]==' '){
 	 nc=nc+1;
 	 if(nc==1){
 	  start_pos=j;
 	 }
 	}else{
 	 end_pos=j;	
 	}
 	
     ++j;
    
    } //while(s[j]!=' ') 
    
    channel=s.substr(start_pos_last,(end_pos-start_pos_last+1));
    	
    if(mp_channels_sim.find(channel) == mp_channels_sim.end()){
     mp_channels_sim[channel]=nchannels_sim;
	 
	 vtmp=split_string(channel,order);
	 mp_channels_sim_inv[nchannels_sim]=vtmp;
    
	 vchannels_sim.push_back(channel);
     ++nchannels_sim;
    }
    
    ichannel=mp_channels_sim[channel];
        
	if(flg_equal==0){	
     if(ichannel!=ichannel_old){
      path+=to_string(ichannel);
      path+=" ";
	  ichannel_old=ichannel;
     }   
    }else{
     path+=to_string(ichannel);
     path+=" ";	
	}
    
    if((end_pos+1)==ssize){break;}; 
    j=start_pos+1;
    start_pos_last=j;   
    nc=0;
  
   }//end while(j<ssize)
  
   vy2[i]="0 "+path+"e";
   
 }//end for i
  
 mp_channels_sim["(conversion)"]=nchannels_sim;
 vchannels_sim.push_back("(conversion)");    
 ++nchannels_sim;
   
 mp_channels_sim["(null)"]=nchannels_sim;
 vchannels_sim.push_back("(null)");     
 ++nchannels_sim;
  
 // /////////////////////////////////////////////////////
 // //CREAZIONE DELLE MATRICI FUNZIONALI ALLE SIMULAZIONI
 // ////////////////////////////////////////////////////
 
 string channel_old;
  
 Fx S(nchannels_sim,nchannels_sim);
     
 for(i=0;i<lvy;i++){
     
  flg_next_path=0;
  flg_next_null=0;
                            
  s=vy2[i];
  s+=" ";
  ssize= (unsigned long int) s.size();
  
  channel="";
  ichannel_old=0;
  ichannel=0;
 
  j=0;
  npassi=0;
  
  vci=vc[i];
  if(flg_var_null==1){
   vni=vn[i];
  }else{
   vni=0;
  }      
  vpi=vci+vni;
   
  while((j<ssize) & (flg_next_path==0)){
      
   while(s[j]!=' '){
  
    if(j<ssize){
     channel+=s[j];
    }
    j=j+1;
   }
               
   if(channel[0]!='0'){//se non è il channel start
    
     if(channel[0]=='e'){ //stato finale
     
      ++npassi;
      
      if(vci>0){ //se ci sono conversion
       ichannel=nchannels_sim-2;
	   S.add(ichannel_old,ichannel,vci);
    
       if(vni>0){
        flg_next_null=1;  
       }else{
        flg_next_path=1;
       }
              
      }
     
      if(((vni>0) | (flg_next_null==1)) & (flg_next_path==0)){ //se non ci sono conversion
       ichannel=nchannels_sim-1;
       S.add(ichannel_old,ichannel,vni);
       flg_next_path=1;
      }
     
     }else{ //stato non finale
      
      if(vpi>0){
       ichannel=atol(channel.c_str());
	   S.add(ichannel_old,ichannel,vpi);
	  }
    
     }
    
     if(flg_next_path==0){
      ++npassi;
     }
    }else{ //stato iniziale
   
     ichannel=0;
   
    }
  
    if(flg_next_path==0){
     //channel_old=channel;
     ichannel_old=ichannel;
    }
     
   if(flg_next_path==0){
    channel="";
    j=j+1;   
   }
   
  }//end while j<size
      
 }//end for 
 
 unsigned long int k,nch0;
 nch0=nchannels-3;
 
 vector<string> vchannels0(nch0);
 for(k=1; k<(nch0+1); k++){
  vchannels0[k-1]=vchannels[k];
 }
 
 
 List res=List::create(Named("transition_matrix")=S.tran_matx(vchannels_sim), Named("channels") = vchannels0);
  
 return(res);
 
 END_RCPP
 
} 	
