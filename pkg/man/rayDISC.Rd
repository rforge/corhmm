\name{rayDISC}
\alias{rayDISC}
\title{Evolution of categorical traits}
\description{Fits a model of evolution for categorical traits, allowing for multi-state characters, polymorphisms, missing data, and incompletely resolved trees}
\usage{
rayDISC(phy,data, ntraits=1, charnum=1, model=c("ER","SYM","ARD"), node.states=c("joint", "marginal", "scaled"), 
p=NULL, par.drop=NULL, par.eq=NULL, root.p=NULL, ip=NULL, lb=0, ub=100)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{ntraits}{specifies the number of traits to included in the analysis.}
\item{charnum}{specified the character to analyze.}
\item{model}{specifies the underlying model.} 
\item{node.states}{method used to calculate ancestral states at internal nodes.}
\item{p}{a vector of transition rates. Allows the user to calculate the likelihood given a specified set of parameter values to specified as fixed and calculate the likelihood.}
\item{par.drop}{a vector of transitions to be dropped from the model. Use \code{rate.mat.index} to see what correspond to which transition.}
\item{par.eq}{a vector of transitions pairs to be set equal. Use \code{rate.mat.index} to see what correspond to which transition.}
\item{root.p}{a vector used to fix the probabilities at the root.}
\item{ip}{initial values used for the likelihood search. Can be a single value or a vector of unique values for each parameter. The default is \code{ip=1}.}
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0}.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=100}.}
}
\details{
The function takes a tree and a trait file and estimates the transitions among character states. The first column of the trait file must contain the species labels to match to the tree, with the second, third, fourth, and so on, corresponding to the traits of interest. Use the \code{charnum} variable to select the trait for analysis. The user can fix the root state probabilities by supplying a vector to the \code{root.p}, otherwise, the program assumes the marginal probability for the root. 
}
\value{
\code{rayDISC} returns an object of class \code{rayDISC}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$ntraits}{The number of traits specified.}
\item{$solution}{a matrix containing the maximum likelihood estimates of the transition rates.}
\item{$solution.se}{a matrix containing the approximate standard errors of the transition rates. The standard error is calculated as the diagonal of the inverse of the Hessian matrix.}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint.}
\item{$opts}{Internal settings of the likelihood search.}
\item{$data}{User-supplied dataset.}
\item{$phy}{User-supplied tree.}
\item{$states}{The likeliest states at each internal node.}
\item{$tip.states}{The likeliest states at each tip.}
\item{$iterations}{The number of iterations used by the optimization routine.}
\item{$eigval}{The eigenvalues from the decomposition of the Hessian of the likelihood function. If any \code{eigval<0} then one or more parameters were not optimized during the likelihood search.}
\item{$eigvect}{The eigenvectors from the decomposition of the Hessian of the likelihood function is returned.}
\item{$bound.hit}{A logical for diagnosing if rate parameters were constrained by \code{lb} or \code{ub} values during optimization.}
\item{$message.tree}{A list of taxa which were listed in the data matrix, but were not present in the passed \code{phylo} object. These taxa will be excluded from the analysis. \code{message.tree} is null if all taxa in data are included in tree.}
\item{$message.data}{A list of taxa which were present in the passed \code{phylo} object, but lacked data in the passed data matrix. These taxa will be coded as missing data (all states equally likely). \code{message.data} is null if all taxa in tree have entries in data matrix.}
}
\author{Jeffrey C. Oliver}
\keyword{models}
