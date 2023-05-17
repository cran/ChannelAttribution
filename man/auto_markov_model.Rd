\name{auto_markov_model}


\alias{auto_markov_model}


\title{Automatic Markov Model.}

\description{Estimate a Markov model from customer journey data after automatically choosing a suitable order. It requires paths that do not lead to conversion as input.}


\usage{
auto_markov_model(Data, var_path, var_conv, var_null, var_value=NULL, 
             max_order=10, roc_npt=100, plot=FALSE, nsim_start=1e5, 
             max_step=NULL, out_more=FALSE, sep=">", 
             ncore=1, nfold=10, seed=0, conv_par=0.05, rate_step_sim=1.5, 
             verbose=TRUE, flg_adv=TRUE)
}

\arguments{
  \item{Data}{data.frame containing customer journeys data.}
  \item{var_path}{column name containing paths.}
  \item{var_conv}{column name containing total conversions.}
  \item{var_null}{column name containing total paths that do not lead to conversions.}
  \item{var_value}{column name containing total conversion value.}  
  \item{max_order}{maximum Markov Model order considered.}
  \item{roc_npt}{number of points used for approximating roc and auc.}
  \item{plot}{if TRUE, a plot with penalized auc with respect to order will be displayed.}
  \item{nsim_start}{minimum number of simulations used in computation.}
  \item{max_step}{maximum number of steps for a single simulated path. if NULL, it is the maximum number of steps found into Data.}
  \item{out_more}{if TRUE, transition probabilities between channels and removal effects will be shown.}
  \item{sep}{separator between the channels.}
  \item{ncore}{number of threads used in computation.}
  \item{nfold}{how many repetitions are used to verify if convergence is reached at each iteration.}
  \item{seed}{random seed. Giving this parameter the same value over different runs guarantees that results will not vary.}
  \item{conv_par}{convergence parameter for the algorithm. The estimation process ends when the percentage of variation of the results over different repetitions is less than convergence parameter.}
  \item{rate_step_sim}{number of simulations used at each iteration is equal to the number of simulations used at previous iteration multiplied by rate_step_sim.}
  \item{verbose}{if TRUE, additional information about process convergence will be shown.}
  \item{flg_adv}{if TRUE, ChannelAttribution Pro banner is printed.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel.}

\author{
  Davide Altomare (\email{info@channelattribution.io}).
}


\examples{

\dontrun{

library(ChannelAttribution)

data(PathData) 

auto_markov_model(Data, "path", "total_conversions", "total_null")

}

}

\keyword{markov model}
\keyword{markov graph}
