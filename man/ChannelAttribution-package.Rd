\name{ChannelAttribution-package}
\alias{ChannelAttribution-package}
\alias{ChannelAttribution}
\docType{package}

\title{
Markov Model for the Online Multi-Channel Attribution Problem
}

\description{
Advertisers use a variety of online marketing channels to reach consumers and they want to know the degree each channel contributes to their marketing success. It's called the online multi-channel attribution problem. In many cases, advertisers approach this problem through some simple heuristics methods that do not take into account any customer interactions and often tend to underestimate the importance of small channels in marketing contribution. 
This package provides a function that approaches the attribution problem in a probabilistic way. It uses a k-order Markov representation to identifying structural correlations in the customer journey data. This would allow advertisers to give a more reliable assessment of the marketing contribution of each channel. 
The approach basically follows the one presented in Eva Anderl, Ingo Becker, Florian v. Wangenheim, Jan H. Schumann (2014). Differently for them, we solved the estimation process using stochastic simulations. In this way it is also possible to take into account conversion values and their variability in the computation of the channel importance. 
The package also contains a function that estimates three heuristic models (first-touch, last-touch and linear-touch approach) for the same problem.
For customizations or interest in other statistical methodologies for web data analysis please contact <davide.altomare@gmail.com>.
}

\details{
\tabular{ll}{
Package: \tab ChannelAttribution\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-10-09\cr
License: \tab GPL (>= 2)\cr
}
 Package contains two functions: markov_model which estimates a k-order Markov model and heuristic_model which estimates three heuristic models (first-touch, last-touch and linear-touch) from customer journey data. 
}

\author{

Davide Altomare and David Loris

Maintainer Davide Altomare <davide.altomare@gmail.com>

}

\references{

Davide Altomare, David Loris (2015). \href{http://www.slideshare.net/adavide1982/markov-model-for-the-multichannel-attribution-problem}{Markov Model for the Online Multi-Channel Attribution Problem}. 

Eva Anderl, Ingo Becker, Florian v. Wangenheim, Jan H. Schumann. Mapping the Customer Journey (2014). A Graph-Based Framework for Online Attribution Modeling.

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

