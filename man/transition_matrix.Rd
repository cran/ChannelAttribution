\name{transition_matrix}


\alias{transition_matrix}


\title{Transition matrix.}

\description{Estimate a k-order transition matrix from customer journey data.}



\usage{
transition_matrix(Data, var_path, var_conv, var_null, order=1, sep=">", 
                  flg_equal=TRUE)
}


\arguments{
  \item{Data}{data.frame containing customer journeys data.}
  \item{var_path}{column name containing paths.}
  \item{var_conv}{column name containing total conversions.}
  \item{var_null}{column name containing paths that do not lead to conversions.}
  \item{order}{Markov Model order.}
  \item{sep}{separator between the channels.}
  \item{flg_equal}{if TRUE, transitions from a channel to itself will be considered.}
}


\value{An object of \code{class} \code{List} containing a dataframe with channel names and a dataframe with the estimated transition matrix.}

\author{
  Davide Altomare (\email{info@channelattribution.net}).
}


\examples{

\dontrun{

library(ChannelAttribution)

data(PathData) 

transition_matrix(Data, var_path="path", var_conv="total_conversions",
                  var_null="total_null", order=1, sep=">", flg_equal=TRUE)

transition_matrix(Data, var_path="path", var_conv="total_conversions",
                  var_null="total_null", order=3, sep=">", flg_equal=TRUE)

}

}

\keyword{markov model}
\keyword{markov graph}
