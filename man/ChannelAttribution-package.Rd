\name{ChannelAttribution-package}
\alias{ChannelAttribution-package}
\alias{ChannelAttribution}
\docType{package}

\title{
Markov Model for Online Multi-Channel Attribution
}

\description{
Advertisers use a variety of online marketing channels to reach consumers and they want to know the degree each channel contributes to their marketing success. This is called online multi-channel attribution problem. In many cases, advertisers approach this problem through some simple heuristics methods that do not take into account any customer interactions and often tend to underestimate the importance of small channels in marketing contribution. 
This package provides a function that approaches the attribution problem in a probabilistic way. It uses a k-order Markov representation to identify structural correlations in the customer journey data. This would allow advertisers to give a more reliable assessment of the marketing contribution of each channel. 
The approach basically follows the one presented in Eva Anderl, Ingo Becker, Florian v. Wangenheim, Jan H. Schumann (2014). Differently for them, we solved the estimation process using stochastic simulations. In this way it is also possible to take into account conversion values and their variability in the computation of the channel importance. 
The package also contains a function that estimates three heuristic models (first-touch, last-touch and linear-touch approach) for the same problem.
}

\details{
\tabular{ll}{
Package: \tab ChannelAttribution\cr
Type: \tab Package\cr
Version: \tab 2.0.4 \cr
Date: \tab 2021-03-29 \cr
License: \tab GPL (>= 2)\cr
}
 Package contains functions for channel attribution in web marketing. 
}

\author{
Davide Altomare, David Loris

Maintainer Davide Altomare <info@channelattribution.net>
}

\references{

\href{https://www.channelattribution.net}{ChannelAttribution Official Website}

\doi{10.2139/ssrn.2343077}{Eva Anderl, Ingo Becker, Florian v. Wangenheim, Jan H. Schumann. Mapping the Customer Journey (2014)}

}

\keyword{marketing attribution}
\keyword{web marketing}
\keyword{multi channel marketing}
\keyword{web statistics}
\keyword{channel marketing}
\keyword{online attribution}
\keyword{customer journey}
\keyword{multi channel funnel}
\keyword{channel attribution}

