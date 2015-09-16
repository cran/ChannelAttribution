\name{markov_model}


\alias{markov_model}


\title{Markov Model for the online attribution problem}
\description{Estimate a first-order Markov model from customer journey data.
}


\usage{
markov_model(Dy, var_channel, var_conv, var_value, 
             nsim = 2e+06, n_boot = 1e+06, n_single_boot = 30, out_trans_matrx = 0)
}


\arguments{
  \item{Dy}{data.frame containing paths and conversions.}
  \item{var_channel}{name of the the column containing paths.}
  \item{var_conv}{name of the column containing total conversions.}
  \item{var_value}{name of the column containing total conversion value.}
  \item{nsim}{number of total simulations from transition matrix.}
  \item{n_boot}{number of bootstrap samples from the empirical distribution of the maxium number of steps between channels.}
  \item{n_single_boot}{dimension of the bootstrap sample.}
  \item{out_trans_matrx}{if set to 1 return the transition matrix between channels.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel.}

\author{
  Davide Altomare (\email{davide.altomare@gmail.com}).
}


\examples{

data(PathData) 

markov_model(Dy, "channel_path", "total_conversions", "total_conversion_value")

}

\keyword{online attribution}
\keyword{attribution problem}
\keyword{channel attribution}
\keyword{multi channel funnel}
\keyword{markov graph}
