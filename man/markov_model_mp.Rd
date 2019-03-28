\name{markov_model_mp}


\alias{markov_model_mp}


\title{Multiprocessing Markov model function.}

\description{This function is a multiprocessing version of markov_model function.}



\usage{
markov_model_mp(Data, var_path, var_conv, var_value=NULL, var_null=NULL, 
             order=1, nsim_start=1e5, max_step=NULL, out_more=FALSE, sep=">", 
             ncore=Inf, nfold=10, seed=0, conv_par=0.05, rate_step_sim=1.5, 
             verbose=TRUE)
}


\arguments{
  \item{Data}{data.frame containing customer journeys data.}
  \item{var_path}{column name containing paths.}
  \item{var_conv}{column name containing total conversions.}
  \item{var_value}{column name containing total conversion value.}
  \item{var_null}{column name containing total paths that do not lead to conversions.}
  \item{order}{Markov Model order.}
  \item{nsim_start}{minimum number of simulations used in computation.}
  \item{max_step}{maximum number of steps for a single simulated path.}
  \item{out_more}{if TRUE returns the transition probabilities between channels and removal effects.}
  \item{sep}{separator between the channels.}
  \item{ncore}{number of threads used in computation. Default is number of CPUs available.}
  \item{nfold}{how many repetitions are used to verify if convergence is reached at each iteration.}
  \item{seed}{random seed. Giving this parameter the same value over different runs guarantees that results will not vary.}
  \item{conv_par}{convergence parameter for the algorithm. The estimation process ends when the percentage of variation of the results over different repetions is less than convergence parameter.}
  \item{rate_step_sim}{number of simulations used at each iteration is equal to the number of simulations used at previous iteration multiplied by rate_step_sim.}
  \item{verbose}{show additional information about process convergence.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel.}

\author{
  Davide Altomare (\email{davide.altomare@gmail.com}).
}


\examples{

\dontrun{

library(ChannelAttribution)

data(PathData) 

markov_model_mp(Data,"path","total_conversions", var_null="total_null")

}

}

\keyword{markov model}
\keyword{markov graph}
