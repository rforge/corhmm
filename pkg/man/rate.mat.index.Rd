\name{rate.mat.index}
\alias{rate.mat.index}
\title{Rate matrix index}
\description{Prints the rate matrix index to allow the user to manipulate the underlying model}
\usage{
rate.mat.index(rate.cat, hrm=TRUE, ntraits=NULL, charnum=NULL, model=c("ER", "SYM", "ARD"), par.drop=NULL, par.eq=NULL)
}       
\arguments{
\item{rate.cat}{specifies the number of rate categories in the HRM.}
\item{hrm}{a logical indicating whether the underlying model is the hidden rates model (HRM). The default is \code{FALSE}.} 
\item{ntraits}{specifies the number of traits in the data file if the underlying model is not the HRM.}
\item{charnum}{specifies the number of characters in the data file used in rayDISC.}
\item{model}{if the model is not HRM, specifies the underlying model.}
\item{par.drop}{a vector of transitions to be dropped from the model. Use \code{rate.mat.index} to see what correspond to which transition.}
\item{par.eq}{a vector of transitions pairs to be set equal. Use \code{rate.mat.index} to see what correspond to which transition.}
}
\details{
Outputs the full index of the rate parameters that are optimized. The intention is that a user might want to see how the matrix is designed prior to an analysis and perhaps drops a few parameters beforehand due to some hypothesis that he or she might have.
}
\value{
}
\examples{
}
\references{
}
\author{Jeremy M. Beaulieu}
\keyword{models}
