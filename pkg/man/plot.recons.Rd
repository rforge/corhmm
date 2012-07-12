\name{plot.recons}
\alias{plot.recons}
\title{Plot ancestral state reconstructions}
\description{Plots maximum likelihood ancestral state estimates on tree}
\usage{
plot.recons(phy, likelihoods, piecolors=NULL, cex=0.5, file=NULL, height=11, width=8.5, showtiplabels=TRUE,title=NULL)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{likelihoods}{likelihoods for ancestral states (see Details).}
\item{piecolors}{a vector of colors for states.}
\item{cex}{specifies the size of the symbols to plot on tree.}
\item{file}{filename to which a pdf is saved.}
\item{height}{height of plot.}
\item{width}{width of plot.}
\item{showtiplabels}{a logical indicating whether to draw tip labels to tree.  The default is \code{TRUE}.}
\item{title}{an optional title for the plot.}
}
\details{
Plots ancestral state estimates on provided tree.  The \code{likelihoods} can be the \code{states} of an object of class \code{rayDISC} or class \code{corDISC}, or the \code{lik.anc} of an object of class \code{ace} (from the \code{ape} package).
}
\value{
}
\examples{
}
\references{
}
\author{Jeffrey C. Oliver}
\keyword{plot}
