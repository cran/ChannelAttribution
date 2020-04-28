\name{heuristic_models}


\alias{heuristic_models}


\title{Heuristic models for the online attribution problem.}

\description{Estimate theree heuristic models (first-touch, last-touch and linear) from customer journey data.
}


\usage{
heuristic_models(Data, var_path, var_conv, var_value=NULL, sep=">")
}


\arguments{
  \item{Data}{data.frame containing paths and conversions.}
  \item{var_path}{column name containing paths.}
  \item{var_conv}{column name containing total conversions.}
  \item{var_value}{column name containing total conversion value.}
  \item{sep}{separator between the channels.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel for each model.}

\author{
  Davide Altomare (\email{davide.altomare@gmail.com}).
}


\examples{

\dontrun{

library(ChannelAttribution)

data(PathData) 

heuristic_models(Data,"path","total_conversions")
heuristic_models(Data,"path","total_conversions",var_value="total_conversion_value")

}

}

\keyword{first touch}
\keyword{last touch}
\keyword{linear touch}



