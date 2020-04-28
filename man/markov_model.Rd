\name{markov_model}


\alias{markov_model}


\title{Markov model for the online attribution problem.}

\description{Estimate a k-order Markov model from customer journey data.
}


\usage{
markov_model(Data, var_path, var_conv, var_value=NULL, var_null=NULL, 
             order=1, nsim=NULL, max_step=NULL, out_more=FALSE, sep=">", 
             seed=NULL,verbose=TRUE)
}


\arguments{
  \item{Data}{data.frame containing paths and conversions.}
  \item{var_path}{column name containing paths.}
  \item{var_conv}{column name containing total conversions.}
  \item{var_value}{column name containing total conversion value.}
  \item{var_null}{column name containing total paths that do not lead to conversions.}
  \item{order}{Markov Model order.}
  \item{nsim}{total simulations from transition matrix.}
  \item{max_step}{maximum number of steps for a single simulated path. if NULL, it is the maximum number of steps found into Data.}
  \item{out_more}{if TRUE, transition probabilities between channels and removal effects will be shown.}
  \item{sep}{separator between the channels.}
  \item{seed}{random seed. Giving to this parameter the same value over different runs guarantee that results will not vary.}
  \item{verbose}{if FALSE, warning message will be disabled.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel.}

\author{
  Davide Altomare (\email{davide.altomare@gmail.com}).
}


\examples{

\dontrun{

library(ChannelAttribution)

data(PathData) 

markov_model(Data, "path", "total_conversions")
markov_model(Data, "path", "total_conversions", var_value="total_conversion_value")
markov_model(Data,"path","total_conversions", var_value="total_conversion_value", 
             var_null="total_null")
markov_model(Data, "path", "total_conversions", var_value="total_conversion_value", 
             var_null="total_null", out_more=TRUE)

}

}

\keyword{markov model}
\keyword{markov graph}
