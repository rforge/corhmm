\name{ancRECON}
\alias{ancRECON}
\title{Ancestral state reconstruction}
\description{Infers ancestral states based on a set of model parameters}
\usage{
ancRECON(phy,data, p, method=c("joint", "marginal", "scaled"), hrm=TRUE, 
rate.cat, ntraits=NULL, charnum=NULL, rate.mat=NULL,
model=c("ER", "SYM", "ARD"), root.p=NULL)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{data}{a data matrix containing species information (see Details).}
\item{p}{a vector of transition rates to be used to estimate ancestral states.}
\item{method}{method used to calculate ancestral states at internal nodes (see Details).}
\item{hrm}{a logical indicating whether the underlying model is the hidden rates model (HRM). The default is \code{FALSE}.} 
\item{rate.cat}{specifies the number of rate categories in the HRM.} 
\item{ntraits}{specifies the number of traits in the data file if the underlying model is not the HRM.}
\item{charnum}{specifies the number of characters in the data file used in rayDISC.}
\item{rate.mat}{a user-supplied rate matrix index of parameters to be optimized.}
\item{model}{if the model is not HRM, specifies the underlying model.}
\item{root.p}{a vector used to fix the probabilities at the root.}
}
\details{
This function computes the marginal, joint, or scaled likelihoods of internal nodes for a given set of transition rates. The tree does not have to be bifurcating in order for the analyses to be carried out. The algorithm based on Pupko et al. (2000, 2002) is used to calculate the \code{joint} estimates of ancestral states. The \code{marginal} method was originally based on a description of an algorithm by Yang (2006). The basic idea of Yang (2006) is that the tree is rerooted on each internal node, with the marginal likelihood being the probabilities of observing the tips states given that the the focal node is the root. However, this takes a ton of time as the number of nodes increases. But more importantly, this does not work easily when the model contains asymmetric rates. Here we use the same dynamic programming algorithm as Mesquite, which calculates the marginal probability at a node using an additional up and down pass of the tree and is time linear with the number of species. If \code{scaled}, the function uses the same algorithm from ace(), though, here the trees do not need to be bifurcating. Note that the \code{scaled} method of ace() is simply the conditional likelihoods of observing everything at or above the focal node. 
}
\value{
For the \code{joint}, a vector of likeliest states at internal nodes and tips. For either \code{marginal} or \code{scaled}, a matrix of the probabilities of each state for each internal node are returned. 
}
\examples{
## Load tree and trait
data(primetes)
## Obtain the marginal reconstruction for a set of parameters:
param<-c(0.05254515,100.00000000,0.01451868,0.01687257,0.006741277,0.000000000,0.01938886,51.68536462)
states<-ancRECON(tree,trait,p=param,method="marginal",hrm=FALSE,ntraits=2,model="ARD")
## Put likeliest states on the tree:
pr<-apply(states$lik.anc.states,1,which.max)
tree$node.label <- pr
}
\author{Jeremy M. Beaulieu and Jeffrey C. Oliver}
\keyword{reconstructions}
