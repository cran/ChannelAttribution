\name{markov_model}


\alias{markov_model}


\title{Markov model for the online attribution problem}
\description{Estimate a k-order Markov model from customer journey data.
}


\usage{
markov_model(Data, var_path, var_conv, var_value=NULL, var_null=NULL, 
             order=1, nsim=NULL, max_step=NULL, out_more=FALSE, sep=">", seed=NULL)
}


\arguments{
  \item{Data}{data.frame containing paths and conversions.}
  \item{var_path}{name of the the column containing paths.}
  \item{var_conv}{name of the column containing total conversions.}
  \item{var_value}{name of the column containing total conversion value.}
  \item{var_null}{name of the column containing total paths that do not lead to conversions.}
  \item{order}{order of Markov Model.}
  \item{nsim}{total simulations from transition matrix.}
  \item{max_step}{maximum number of steps for a single simulated path.}
  \item{out_more}{if TRUE returns the transition probabilities between channels and removal effects.}
  \item{sep}{separator between the channels.}
  \item{seed}{random seed. giving to this parameter the same value over different runs guarantee that results will not vary.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel.}

\author{
  Davide Altomare (\email{davide.altomare@gmail.com}).
}


\examples{

data(PathData) 

#uncomment the following lines to run the examples

#markov_model(Data, "path", "total_conversions")
#markov_model(Data, "path", "total_conversions", 
#             var_value="total_conversion_value")
#markov_model(Data,"path","total_conversions",
#             var_value="total_conversion_value", var_null="total_null")
#markov_model(Data, "path", "total_conversions", 
#             var_value="total_conversion_value", var_null="total_null", out_more=TRUE)



}

\keyword{markov model}
\keyword{markov graph}
