\name{corDISC}
\alias{corDISC}
\title{Correlated evolution binary traits}
\description{Fits a model of correlated evolution between two or three binary traits}
\usage{
corDISC(phy,data, ntraits=2, rate.mat=NULL, model=c("ER","SYM","ARD"), 
node.states=c("joint", "marginal", "scaled"), p=NULL, root.p=NULL, ip=NULL, 
lb=0, ub=100, diagn=TRUE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{ntraits}{specifies the number of traits to included in the analysis.}
\item{rate.mat}{a user-supplied rate matrix index of parameters to be optimized.}
\item{model}{specifies the underlying model.} 
\item{node.states}{method used to calculate ancestral states at internal nodes (see Details).}
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.}
\item{root.p}{a vector used to fix the probabilities at the root.}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0}.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=100}.}
\item{diagn}{logical indicating whether diagnostic tests should be performed. The default is \code{TRUE}.}
}
\details{
The function takes a tree and a trait file and estimates the transitions between a two or three binary characters. The first column of the trait file must contain the species labels to match to the tree, with the second, and so on, corresponding to the binary traits of interest. The user can fix the root state probabilities by supplying a vector to the \code{root.p}. Also, the input phylogeny need not be bifurcating as the algorithm implemented can easily handle multifucations. Note that rayDISC can be used to evaluate the same models are in corDISC. The difference is that with rayDISC the rate matrix would have to be manipulated in order to remove parameters associated with dual transitions. 
}
\value{
\code{corDISC} returns an object of class \code{corDISC}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$ntraits}{The number of traits specified.}
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
# Not run
# Load tree and data
# data(primates)
# Obtain the fit for two binary characters
# pp<-corDISC(tree,trait,ntraits=2,model="ARD",node.states="marginal")
# pp
}
\references{
Beaulieu J.M., and M.J. Donoghue 2013. Fruit evolution and diversification in campanulid angiosperms. Evolution In revision.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
