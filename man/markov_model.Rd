\name{markov_model}


\alias{markov_model}


\title{Markov model for the online attribution problem.}

\description{Estimate a k-order Markov model from customer journey data. Differently from markov_model, this function iterates estimation until convergence is reached and enables multiprocessing.}



\usage{
markov_model(Data, var_path, var_conv, var_value=NULL, var_null=NULL, 
             order=1, nsim_start=1e5, max_step=NULL, out_more=FALSE, sep=">", 
             ncore=1, nfold=10, seed=0, conv_par=0.05, rate_step_sim=1.5, 
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
  \item{max_step}{maximum number of steps for a single simulated path. if NULL, it is the maximum number of steps found into Data.}
  \item{out_more}{if TRUE, transition probabilities between channels and removal effects will be returned.}
  \item{sep}{separator between the channels.}
  \item{ncore}{number of threads used in computation.}
  \item{nfold}{how many repetitions are used to verify if convergence has been reached at each iteration.}
  \item{seed}{random seed. Giving this parameter the same value over different runs guarantees that results will not vary.}
  \item{conv_par}{convergence parameter for the algorithm. The estimation process ends when the percentage of variation of the results over different repetions is less than convergence parameter.}
  \item{rate_step_sim}{number of simulations used at each iteration is equal to the number of simulations used at previous iteration multiplied by rate_step_sim.}
  \item{verbose}{if TRUE, additional information about process convergence will be shown.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel.}

\author{
  Davide Altomare (\email{info@channelattribution.net}).
}


\examples{

\dontrun{

library(ChannelAttribution)

data(PathData) 

#Estimate a Makov model using total conversions
markov_model(Data, var_path="path", "total_conversions")

#Estimate a Makov model using total conversions and revenues
markov_model(Data, "path", "total_conversions",
var_value="total_conversion_value")

#Estimate a Makov model using total conversions, revenues and paths that do not lead to conversions
markov_model(Data, "path", "total_conversions",
var_value="total_conversion_value", var_null="total_null")

#Estimate a Makov model returning transition matrix and removal effects
markov_model(Data, "path", "total_conversions",
var_value="total_conversion_value", var_null="total_null", out_more=TRUE)

#Estimate a Markov model using 4 threads
markov_model(Data, "path", "total_conversions",
var_value="total_conversion_value", ncore=4)

}

}

\keyword{markov model}
\keyword{markov graph}
