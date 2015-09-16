#define __GXX_EXPERIMENTAL_CXX0X__ 1

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <numeric>

//#include <Rcpp.h>
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
using namespace arma;

template <typename T>
string to_string(T pNumber)
{
 ostringstream oOStrStream;
 oOStrStream << pNumber;
 return oOStrStream.str();
}


bool regex_search_v0(string string0, string exp0){	
	
 long int ls0=(long int) string0.size();
 long int le0=(long int) exp0.size();
 long int i,j,k,ck;
 
 i=0;
 while(i<ls0){
  j=0;
  ck=0;
  
  for(k=0;k<le0;k++){
   if(string0[i+k]==exp0[j+k]){
    ++ck;
   }else{
    break;	  
   }
   if(ck==le0){return(1);}
  }
  
  ++i;
 }
 
 return(0); 

}

	
bool regex_search_v1(string string0, vector<string> vexp0){	
	
 long int ls0=(long int) string0.size();
 long int le0;
 long int i,j,k,z,ck;
 string exp0;
 
 long int lvexp0=(long int) vexp0.size();
 
 for(z=0;z<lvexp0;z++){
 
  exp0=vexp0[z];
  le0=(long int) exp0.size();
  
  i=0;
  while(i<ls0){
   j=0;
   ck=0;
   
   for(k=0;k<le0;k++){
    if(string0[i+k]==exp0[j+k]){
     ++ck;
    }else{
     break;	  
    }
    if(ck==le0){return(1);}
   }
   
   ++i;
  }
 
 }// end for z
 
 return(0); 

}	


vector<string> split_string(string s, string sep){

 long int j,ls;
 
 ls=(long int)s.size();
 string r;
 vector<string> v;
 
 j=0;
 while(j<ls){
  
  while((s[j]!=sep[0]) & (j<ls)){
   r+=s[j];
   ++j;   
  }
	
  v.push_back(r);	
  r="";
  ++j;  
 }
 
 return(v);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP heuristic_models_cpp(SEXP Dy_p, SEXP var_channel_p, SEXP var_conv_p, SEXP var_value_p)
{
	
 BEGIN_RCPP

 //inp.a
 
 List Dy(Dy_p);
 
 CharacterVector var_channel_0(var_channel_p);
 CharacterVector var_conv_0(var_conv_p);
 CharacterVector var_value_0(var_value_p); 
  
 string var_channel = Rcpp::as<string>(var_channel_0);
 string var_conv = Rcpp::as<string>(var_conv_0);
 string var_value = Rcpp::as<string>(var_value_0);
  
 //inp.b
	
 CharacterVector vy0 = Dy[var_channel];
 vector<string> vy = Rcpp::as<vector<string> >(vy0);

 NumericVector vc0 = Dy[var_conv];
 vector<long long int> vc = Rcpp::as<vector<long long int> >(vc0);

 NumericVector vv0 = Dy[var_value];
 vector<double> vv = Rcpp::as<vector<double> >(vv0);
 
 long long int i,j,k,lvy,ssize;
 long long int nchannels;
 string s,channel,channel_first,channel_last;
  
 lvy=(long long int) vy.size();
 nchannels=0;
 
 map<string,long long int> mp_channels;
 vector<string> vchannels;
	  	
 map<string,double> mp_first_conv;
 map<string,double> mp_first_val;	
 map<string,double> mp_last_conv;
 map<string,double> mp_last_val;
 map<string,double> mp_linear_conv;
 map<string,double> mp_linear_val;
 
 vector<string> vchannels_linear;
 double nchannels_linear;
 string kchannel;
  
 for(i=0;i<lvy;i++){
	 	 
  s=vy[i];
  
  s+=" >";
  ssize=(long long int) s.size();
  channel="";
  j=0;
  nchannels_linear=0;
  vchannels_linear.clear();
   
  while(j<ssize){  
  
    if((j>0) & (ssize>1)){
   
     if((s[j-1]=='>') & (s[j]==' ')){
	  j=j+1;
	 }
	 if((s[j]==' ') & (s[j+1]=='>')){
	  j=j+2;
	  break;
	 }
	
	}

    while(s[j]!='>'){
	 channel+=s[j];
     ++j;	
    }

  
    if(mp_channels.find(channel) == mp_channels.end()){
	 mp_channels[channel]=nchannels;
	 vchannels.push_back(channel);
	 ++nchannels;
	
     mp_first_conv[channel]=0;
	 mp_first_val[channel]=0;	
	 mp_last_conv[channel]=0;
	 mp_last_val[channel]=0;
	 mp_linear_conv[channel]=0;
	 mp_linear_val[channel]=0;
    }
	 	 
    //lista canali unici
    if(nchannels_linear==0){
     vchannels_linear.push_back(channel);
	 ++nchannels_linear;
    }else if(find(vchannels_linear.begin(),vchannels_linear.end(),channel)==vchannels_linear.end()){
	 vchannels_linear.push_back(channel);
	 ++nchannels_linear;
    }
	
    channel_last=channel;
  
    channel="";
    ++j;
    
  }//end while j
   
  channel_first=vchannels_linear[0];
  mp_first_conv[channel_first]=mp_first_conv[channel_first]+vc[i];
  mp_first_val[channel_first]=mp_first_val[channel_first]+vv[i];   
 
  mp_last_conv[channel_last]=mp_last_conv[channel_last]+vc[i];
  mp_last_val[channel_last]=mp_last_val[channel_last]+vv[i];  
 
  //linear
  for(k=0;k<nchannels_linear;k++){
    kchannel=vchannels_linear[k];
    mp_linear_conv[kchannel]=mp_linear_conv[kchannel]+(vc[i]/nchannels_linear);
    mp_linear_val[kchannel]=mp_linear_val[kchannel]+(vv[i]/nchannels_linear); 
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
  vfirst_val[k]=mp_first_val[kchannel];
  vlast_conv[k]=mp_last_conv[kchannel];
  vlast_val[k]=mp_last_val[kchannel];
  vlinear_conv[k]=mp_linear_conv[kchannel];
  vlinear_val[k]=mp_linear_val[kchannel];
 }
 
 return DataFrame::create(Named("channel_name")=vchannels, Named("first_conversion") = vfirst_conv, Named("first_value") = vfirst_val, Named("last_conversion") = vlast_conv, Named("last_value") = vlast_val, Named("linear_conversion") = vlinear_conv, Named("linear_value") = vlinear_val);
 
 END_RCPP

}//end heuristic_models_cpp


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


RcppExport SEXP markov_model_cpp(SEXP Dy_p, SEXP var_channel_p, SEXP var_conv_p, SEXP var_value_p, SEXP nsim_p, SEXP n_boot_p, SEXP n_single_boot_p)
{
	
 	
 BEGIN_RCPP

 //inp.a
 
 List Dy(Dy_p);
 
 CharacterVector var_channel_0(var_channel_p);
 
 CharacterVector var_conv_0(var_conv_p);
 CharacterVector var_value_0(var_value_p); 
  
 NumericVector nsim_0(nsim_p); 
 NumericVector n_boot_0(n_boot_p); 
 NumericVector n_single_boot_0(n_single_boot_p); 
 
 string var_channel = Rcpp::as<string>(var_channel_0);
   
 string var_conv = Rcpp::as<string>(var_conv_0);
 string var_value = Rcpp::as<string>(var_value_0);
	 
 long long int nsim = Rcpp::as<long long int>(nsim_0);
 long long int n_boot = Rcpp::as<long long int>(n_boot_0);
 long long int n_single_boot = Rcpp::as<long long int>(n_single_boot_0);

 //inp.b 
 

 CharacterVector vy0 = Dy[var_channel];
 vector<string> vy = Rcpp::as<vector<string> >(vy0);
 
 NumericVector vc0 = Dy[var_conv];
 vector<long long int> vc = Rcpp::as<vector<long long int> >(vc0);
 
 NumericVector vv0 = Dy[var_value];
 vector<double> vv = Rcpp::as<vector<double> >(vv0);
  
 long long int i,j,k,lvy,ssize;
 long long int nchannels,npassi;
 string s,channel,path;
 map<string,long long int> mp_channels;
 map<long long int,long long int> mp_npassi;
 vector<long long int> vnpassi;
  
 //cout << "Processed 1/4" << endl;
 
 lvy=(long long int) vy.size();
   
 //////////////////////
 //CODIFICA DA ONE STEP 
 //////////////////////
 
 nchannels=0;
 
 vector<string> vy2(lvy);
 
 mp_channels["(start)"]=0;
 vector<string> vchannels;
 vchannels.push_back("(start)");	 
 ++nchannels;
 
 map<long long int,long long int> W; 
 long long int lvnpassi;
 
 lvnpassi=0;

   	
 for(i=0;i<lvy;i++){
	 
  s=vy[i];
  s+=" >";
  ssize=(long long int) s.size();
  channel="";
  path="";
  j=0;
  npassi=0;
  
  //medium.touch
  
  while(j<ssize){  
   	   
   if((j>0) & (ssize>1)){
   
     if((s[j-1]=='>') & (s[j]==' ')){
	  j=j+1;
	 }
	 if((s[j]==' ') & (s[j+1]=='>')){
	  j=j+2;
	  break;
	 }
	
   }

   while(s[j]!='>'){
	 channel+=s[j];
     ++j;	
   }

	
   if(mp_channels.find(channel) == mp_channels.end()){
    mp_channels[channel]=nchannels;
    vchannels.push_back(channel);
    ++nchannels;
   }
	
   if(npassi==0){
    path="0 ";	
   }else{
    path+=" ";
   }
    
   path+=to_string(mp_channels[channel]);
  
   channel="";
   ++j;
   ++npassi;
    
  }//end while j
    
  vy2[i]=path+" c"; //aggiungo lo stato finale di conversion 
  
  if(W.find(npassi)==W.end()){
   W[npassi]=0;  
   vnpassi.push_back(npassi);
   ++lvnpassi;
  }
  W[npassi]=W[npassi]+vc[i];
 
 
 }//end for
    
 ++nchannels; //aggiungo canale conversion	 
 //inserisco lo stato/channel conversion
 mp_channels["(conversion)"]=(nchannels-1);
 vchannels.push_back("(conversion)");	 

 //cout << "Processed 2/4" << endl;
  
 /////////////////////////////////////////////////////
 //CREAZIONE DELLE MATRICI FUNZIONALI ALLE SIMULAZIONI
 ////////////////////////////////////////////////////
   
 SpMat<long long int> M(nchannels,nchannels);
 SpMat<long long int> S0(nchannels,nchannels);
 vector<long long int> lrS0(nchannels,0);

 long long int ichannel_old;
 long long int ichannel,non_zeros;
 long long int val0,lval0;
  
 non_zeros=0;
 
 //conv val
 map<long long int,vector<double> > mp_conv_val;
 map<long long int,long long int> mp_l_conv_val;
    
 for(i=0;i<lvy;i++){

  s=vy2[i];
  s+=" ";
  ssize= (long long int) s.size();
  channel="";
  ichannel_old=-1;
  ichannel=-1;
  j=0;
  
  while(j<ssize){
	  
   while(s[j]!=' '){
  
    if(j<ssize){
     channel+=s[j];
    }
    j=j+1;
  
   }

   if(channel[0]!='c'){   
    ichannel=atol(channel.c_str());   
   }else{
    ichannel=nchannels-1;
   }
	
   if(ichannel_old!=-1){
	
    val0=M(ichannel_old,ichannel);
     
	 if(val0==0){
	  lval0=lrS0[ichannel_old];
      S0(ichannel_old,lval0)=ichannel;
	  lrS0[ichannel_old]=lval0+1;
      ++non_zeros; 	 
     }
	 
	 M(ichannel_old,ichannel)=val0+vc[i];	 
   
   }
      
   ichannel_old=ichannel;
   channel="";
   
   j=j+1;   
   
  }
  
  
  //conv val
  if(mp_l_conv_val.find(ichannel)==mp_l_conv_val.end()){
   mp_conv_val[ichannel]=vector<double>();
   mp_l_conv_val[ichannel]=0;   
  }
  mp_conv_val[ichannel].push_back(vv[i]/vc[i]);
  mp_l_conv_val[ichannel]=mp_l_conv_val[ichannel]+1;
  
 }//end for 
  

 //distribuzione della matrice di transizione --> cumulato M --> S
  
 SpMat<long long int> S(nchannels,nchannels); 
 vector<long long int> lrS(nchannels,0);
 long long int lrs0;
 
 for(i=0;i<nchannels;i++){
  
  lrs0=lrS0[i];
  
  for(j=0;j<lrs0;j++){
    
   if(j==0){
    S(i,j)=M(i,S0(i,j));
   }else{
	S(i,j)=S(i,j-1)+M(i,S0(i,j)); 
   }
  
  }   
  
  if(j>0){
   lrS[i]=S(i,j-1);   
  }
   
 }	 
 
 //distribuzione del numero di passi --> cumulato W --> sp
 
 vector<long long int>sp(lvnpassi);
 long long int lsp;
 
 for(i=0;i<lvnpassi;i++){
  if(i==0){
   sp[i]=W[i];   
  }else{
   sp[i]=sp[i-1]+W[i];	
  }
	
 }
 
 lsp=sp[lvnpassi-1];
 
 //distribuzione numeri uniformi
 
 double vu,nuf;
 nuf=1000000;
 NumericVector vunif=runif(nuf);
 
 //distribuzione del max numero di passi
  
 map<long long int,long long int> Wm; 
 long long int lvnpassim=0;
 vector<long long int> vnpassim;
 long long int s0,smax;
 vu=0;
 
 for (i = 0; i < n_boot; i++){
  smax=0;
  for (j = 0; j < n_single_boot; j++){
   //genero il numero di passi 
   if(vu>=nuf){vunif=runif(nuf);vu=0;}
   s0=floor(vunif[vu]*lsp+1);
   ++vu;
   for (k = 0; k < lvnpassi; k++){
    if(sp[k]>=s0){s0=vnpassi[k];break;}
   }
   smax=max(smax,s0); 
  }
  if(Wm.find(smax)==Wm.end()){
	Wm[smax]=0;
	vnpassim.push_back(smax);
	++lvnpassim;
  }
  Wm[smax]=Wm[smax]+1;  
 }
 
 vector<long long int>spm(lvnpassim);
 long long int lspm;
 
 for(i=0;i<lvnpassim;i++){
  if(i==0){
   spm[i]=Wm[i];   
  }else{
   spm[i]=spm[i-1]+Wm[i];	
  }
	
 }

 lspm=spm[lvnpassim-1];
  
 //cout << "Processed 3/4" << endl;
 
 //SIMULAZIONI
  
 long long int s1,c,nconv;
 vector<bool> C(nchannels);
 vector<double> T(nchannels);
 vector<double> V(nchannels);
   
 nconv=0;
 
 for (i = 0; i < nsim; i++){
	 	   
  c=0;
  npassi=0;
  
  //svuoto C
  for (k = 0; k < nchannels; k++){
   C[k]=0;
  }

  C[c]=1; //assegno 1 al channel start
   
  //genero il max numero di passi 
  if(vu>=nuf){vunif=runif(nuf);vu=0;}
  s0=floor(vunif[vu]*lspm+1);
  ++vu;
  for (k = 0; k < lvnpassim; k++){
   if(spm[k]>=s0){s0=vnpassim[k];break;}
  }
  
  while(npassi<=s0){ //interrompo quando raggiungo il massimo numero di passi
   
   //estraggo l'indice del channel
   if(vu>=nuf){vunif=runif(nuf);vu=0;}
   s1=floor(vunif[vu]*lrS[c]+1);
   ++vu;
   for (k = 0; k < lrS0[c]; k++){
    if(S(c,k)>=s1){c=S0(c,k);break;}
   }
  
   C[c]=1;
 
   if(c==(nchannels-1)){ //se ho raggiunto lo stato conversion interrompo
    goto go_1;	
   }
     
   ++npassi;
 
  }//end while npassi 
 
  go_1:
 
  if(c==(nchannels-1)){ //solo se ho raggiunto la conversion assegno +1 ai canali interessati
   for (k = 0; k < nchannels; k++){
    if(C[k]==1){
	 T[k]=T[k]+1;
	 //conv val
	 if(vu>=nuf){vunif=runif(nuf);vu=0;}
     s0=floor(vunif[vu]*mp_l_conv_val[c]);
	 ++vu;
	 
	 V[k]=V[k]+mp_conv_val[c][s0];
    }
   }
   ++nconv;
  }
 	
 }//end for i

 double sn=accumulate(vc.begin(), vc.end(), 0.0);
   
 T[0]=0; //pongo channel start = 0
 T[nchannels-1]=0; //pongo channel conversion = 0 
  
 double sm=accumulate(T.begin(), T.end(), 0.0);
 vector<double> TV(nchannels-2,0);
 
 for (k = 1; k < nchannels-1; k++){
  if(sm>0){
   TV[k-1]=(T[k]/sm)*sn;
  }
 }
  
 V[0]=0; //pongo channel start = 0
 V[nchannels-1]=0; //pongo channel conversion = 0 
   
 sn=accumulate(vv.begin(), vv.end(), 0.0);
 
 sm=accumulate(V.begin(), V.end(), 0.0);
 vector<double> VV(nchannels-2,0);
 
 for (k = 1; k < nchannels-1; k++){
  if(sm>0){
   VV[k-1]=(V[k]/sm)*sn;
  }
 }
    
 vector<string> vchannels0(nchannels-2);
 for (k = 1; k < nchannels-1; k++){
  vchannels0[k-1]=vchannels[k];
 }
 

 //cout << "Processed 4/4" << endl; 
  
 return DataFrame::create(Named("channel_name")=vchannels0, Named("total_conversion") = TV, Named("total_conversion_value") = VV );
 
 END_RCPP

}		