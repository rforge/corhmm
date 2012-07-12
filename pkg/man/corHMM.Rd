\name{corHMM}
\alias{corHMM}
\title{Hidden Rates Model}
\description{Estimates hidden rates underlying the evolution of a binary character}
\usage{
corHMM(phy, data, rate.cat, node.states=c("joint", "marginal","scaled"), optim.method=c("subplex"), 
p=NULL, par.drop=NULL, par.eq=NULL, root.p=NULL, ip=NULL, nstarts=10, n.cores=NULL, lb=0.00001, ub=100)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{rate.cat}{specifies the number of rate categories in the HRM.} 
\item{node.states}{method used to calculate ancestral states at internal nodes (see Details).}
\item{optim.method}{method used to perform optimization. The default is \code{subplex}.} 
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.}
\item{par.drop}{a vector of transitions to be dropped from the model. Use \code{rate.mat.index} to see what correspond to which transition.}
\item{par.eq}{a vector of transitions pairs to be set equal. Use \code{rate.mat.index} to see what correspond to which transition.}
\item{root.p}{a vector used to fix the probabilities at the root.}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{nstarts}{the number of random restarts to be performed. The default is \code{nstarts=10}.}
\item{n.cores}{the number of processor cores to spread out the random restarts.}
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0.00001}.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=100}.}
}
\details{
The function takes a tree and a trait file and estimates the transitions between a binary character using the hidden rate model. The first column of the trait file must contain the species labels to match to the tree, with the second corresponding to the binary trait of interest. Can test models that assume 1, 2, 3, 4, or 5 rate categories underlying the observed data. The user can fix the root state probabilities by supplying a vector to the \code{root.p}, otherwise, the program assumes the marginal probability for the root.. 
}
\value{
\code{corHMM} returns an object of class \code{corHMM}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$rate.cat}{The number of rate categories specified.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the transition rates.}
\item{$solution.se}{a matrix containing the approximate standard errors of the transition rates. The standard error is calculated as the diagonal of the inverse of the Hessian matrix.}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint.}
\item{$opts}{Internal settings of the likelihood search}
\item{$data}{User-supplied dataset.}
\item{$phy}{User-supplied tree.}
\item{$states}{The likeliest states at each internal node.}
\item{$tip.states}{The likeliest states at each tip.}
\item{$iterations}{The number of iterations used by the optimization routine.}
\item{$eigval}{The eigenvalues from the decomposition of the Hessian of the likelihood function. If any \code{eigval<0} then one or more parameters were not optimized during the likelihood search}
\item{$eigvect}{The eigenvectors from the decomposition of the Hessian of the likelihood function is returned}
}
\examples{
}
\references{
}
\author{Jeremy M. Beaulieu}
\keyword{models}
