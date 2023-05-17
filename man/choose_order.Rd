\name{choose_order}


\alias{choose_order}


\title{Choose order for Markov model.}

\description{Find the minimum Markov Model order that gives a good representation of customers' behaviour for data considered. It requires paths that do not lead to conversion as input. Minimum order is found maximizing a penalized area under ROC curve. 
}


\usage{

choose_order(Data, var_path, var_conv, var_null, max_order=10, sep=">", 
             ncore=1, roc_npt=100, plot=TRUE, flg_adv=TRUE)

}


\arguments{
  \item{Data}{data.frame containing customer journeys.}
  \item{var_path}{column name of Data containing paths.}
  \item{var_conv}{column name of Data containing total conversions.}
  \item{var_null}{column name of Data containing total paths that do not lead to conversion.}
  \item{max_order}{maximum Markov Model order considered.}
  \item{sep}{separator between channels.}
  \item{ncore}{number of threads used in computation.}
  \item{roc_npt}{number of points used for approximating roc and auc.}
  \item{plot}{if TRUE, a plot with penalized auc with respect to order will be displayed.}
  \item{flg_adv}{if TRUE, ChannelAttribution Pro banner is printed.}
}


\value{An object of \code{class} \code{List} with the estimated roc, auc and penalized auc.}

\author{
  Davide Altomare (\email{info@channelattribution.io}).
}


\examples{

\dontrun{

library(ChannelAttribution)

data(PathData) 

res=choose_order(Data, var_path="path", var_conv="total_conversions",
                 var_null="total_null")

#plot auc and penalized auc	   
	   
plot(res$auc$order,res$auc$auc,type="l",xlab="order",ylab="pauc",main="AUC")
lines(res$auc$order,res$auc$pauc,col="red")
legend("right", legend=c("auc","penalized auc"), 
       col=c("black","red"),lty=1)

}
}

\keyword{choose markov model order}
\keyword{choose markov graph order}
