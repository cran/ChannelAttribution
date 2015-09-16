\name{heuristic_models}


\alias{heuristic_models}


\title{Heuristic models for the online attribution problem}
\description{Estimate theree heuristic models (first-touch, last-touch, linear) from customer journey data.
}


\usage{
heuristic_models(Dy, var_channel, var_conv, var_value)
}


\arguments{
  \item{Dy}{data.frame containing paths and conversions.}
  \item{var_channel}{name of the the column containing paths.}
  \item{var_conv}{name of the column containing total conversions.}
  \item{var_value}{name of the column containing total conversion value.}
}


\value{An object of \code{class} \code{data.frame} with the estimated number of conversions and the estimated conversion value attributed to each channel for each model.}

\author{
  Davide Altomare (\email{davide.altomare@gmail.com}).
}


\examples{

data(PathData) 

heuristic_models(Dy, "channel_path", "total_conversions", "total_conversion_value")

}


\keyword{online attribution}
\keyword{channel attribution}
\keyword{multi channel funnel}
\keyword{attribution problem}
\keyword{first touch}
\keyword{last touch}
\keyword{linear touch}



