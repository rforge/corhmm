#EVOLUTION OF DISCRETE TRAITS, ALLOWING POLYMORPHIC AND MISSING STATES

#written by Jeremy M. Beaulieu & Jeffrey C. Oliver

#Takes a tree and a trait file and estimates the rate of transitions among states,
#The first column of the trait file must contain the species labels 
#to match to the tree, with any additional columns corresponding to traits of interest. 
#Can test three models: ER = equal rates;SYM = forw/back equal, ARD = all 
#rates unequal. Can return three different types of most likely ancestral states:
#joint, marginal, or scaled.

require(ape)
require(nloptr)
require(numDeriv)
require(expm)
require(corpcor)
require(phangorn)
#require(multicore)
source("ancRECON.one.trait.R")

#corDISC<-function(phy,data, ntraits=2, model=c("ER","SYM","ARD"), method=c("joint", "marginal", "scaled"), nstarts=10, n.cores=NULL, p=NULL, par.drop=NULL, par.eq=NULL, root.p=NULL, ip=NULL){
rayDISC<-function(phy,data, ntraits=1,charnum=1, model=c("ER","SYM","ARD"), method=c("joint", "marginal", "scaled"), nstarts=10, n.cores=NULL, p=NULL, par.drop=NULL, par.eq=NULL, root.p=NULL, ip=NULL){

	#Creates the data structure and orders the rows to match the tree
	phy$edge.length[phy$edge.length==0]=1e-5

	# Checks to make sure phy & data have same taxa.  Fixes conflicts (see match.tree.data function).
	matching <- match.tree.data(phy,data) 
	data <- matching$data
	phy <- matching$phy

	# Won't perform reconstructions on invariant characters
	if(nlevels(as.factor(data[,charnum+1])) <= 1){
		obj <- NULL
		obj$loglik <- NULL
		obj$diagnostic <- paste("Character ",charnum," is invariant.  Reconstructions stopped.",sep="")
		return(obj)
	} else {
		# Still need to make sure second level isn't just an ambiguity
		lvls <- as.factor(data[,charnum+1])
		if(nlevels(as.factor(data[,charnum+1])) == 2 && length(which(lvls == "?"))){
			obj <- NULL
			obj$loglik <- NULL
			obj$diagnostic <- paste("Character ",charnum," is invariant.  Reconstructions stopped.",sep="")
			return(obj)
		}
	}

#	data <- data.frame(data[,charnum+1],data[,charnum+1],row.names=data[,1]) # added character twice, because at least two columns are necessary
#	data <- data[phy$tip.label,] # this might have already been done by match.tree.data	

	workingData <- data.frame(data[,charnum+1],data[,charnum+1],row.names=data[,1]) # added character twice, because at least two columns are necessary
	workingData <- workingData[phy$tip.label,] # this might have already been done by match.tree.data	

	counts <- table(workingData[,1])
	levels <- levels(workingData[,1])
	cols <- as.factor(workingData[,1])
	cat("State distribution in data:\n")
	cat("States:",levels,"\n",sep="\t")
	cat("Counts:",counts,"\n",sep="\t")
	#Have to collect this here. When you reorder, the branching time function is not correct:
	tl<-max(branching.times(phy))
	#Some initial values for use later - will clean up
	k <- 1 # Only one trait allowed
	factored <- factorData(workingData) # just factoring to figure out how many levels (i.e. number of states) in data.
	nl <- ncol(factored)
	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	
	model=model
	par.drop=par.drop
	par.eq=par.eq
	root.p=root.p	
	nstarts=nstarts
	ip=ip

	model.set.final<-rate.cat.set.oneT(phy=phy,data=workingData,model=model,par.drop=par.drop,par.eq=par.eq)
	lower = rep(0.0, model.set.final$np)
	upper = rep(100, model.set.final$np)

	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)

	if(!is.null(p)){
		cat("Calculating likelihood from a set of fixed parameters", "\n")
		out<-NULL
		out$solution<-p
		out$objective<-dev.raydisc(out$solution,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
		loglik <- -out$objective
		est.pars<-out$solution
	}
	if(is.null(ip)){
		cat("Initializing...", "\n")
		#Sets parameter settings for random restarts by taking the parsimony score and dividing
		#by the total length of the tree
		model.set.init<-rate.cat.set.oneT(phy=phy,data=workingData,model="ER",par.drop=par.drop,par.eq=par.eq)
		opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)
		dat<-as.matrix(workingData)
		dat<-phyDat(dat,type="USER", levels=levels(as.factor(workingData[,1])))
		par.score<-parsimony(phy, dat, method="fitch")
		tl <- sum(phy$edge.length)
		mean = par.score/tl
		ip<-rexp(1, 1/mean)
		lower.init = rep(0, model.set.init$np)
		upper.init = rep(100, model.set.init$np)
		init = nloptr(x0=rep(ip, length.out = model.set.init$np), eval_f=dev.raydisc, lb=lower.init, ub=upper.init, opts=opts, phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,root.p=root.p)
		cat("Finished. Begin thorough search...", "\n")
		lower = rep(0, model.set.final$np)
		upper = rep(100, model.set.final$np)
		out <- nloptr(x0=rep(init$solution, length.out = model.set.final$np), eval_f=dev.raydisc, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
		loglik <- -out$objective
		est.pars<-out$solution
	}
	#If a user-specified starting value(s) is supplied: is this redundant with if(!is.null(p)) conditional above?
	else{
		cat("Begin subplex optimization routine -- Starting value(s):", ip, "\n")
		opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)
		out = nloptr(x0=rep(ip, length.out = model.set.final$np), eval_f=dev.raydisc, lb=lower, ub=upper, opts=opts)
		loglik <- -out$objective
		est.pars<-out$solution
	}
	
	#Starts the summarization process:
	cat("Finished. Inferring ancestral states using", method, "reconstruction.","\n")
	
	lik.anc <- ancRECON.one.trait(phy, data, est.pars, method=method, model=model, charnum=charnum,par.drop=par.drop, par.eq=par.eq, root.p=root.p)
	if(method == "marginal" || method == "scaled"){
		pr<-apply(lik.anc$lik.anc.states,1,which.max)
		phy$node.label <- pr
		tip.states <- NULL
	}
	if(method == "joint"){
		phy$node.label <- lik.anc$lik.anc.states
		tip.states <- lik.anc$lik.tip.states
	}

	cat("Finished. Performing diagnostic tests.", "\n")
	
	#Approximates the Hessian using the numDeriv function
	h <- hessian(func=dev.raydisc, x=est.pars, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
	solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
	solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix], dim(model.set.final$index.matrix))
	
	rownames(solution) <- rownames(solution.se) <- c("0","1")
	colnames(solution) <- colnames(solution.se) <- c("0","1")
	if(is.character(method)){
		if (method == "marginal"){
			colnames(lik.anc$lik.anc.states) <- c("0","1")
		}
	}

	hess.eig <- eigen(h,symmetric=TRUE)
	eigval<-signif(hess.eig$values,2)
	eigvect<-round(hess.eig$vectors, 2)
	obj = list(loglik = loglik, AIC = -2*loglik+2*model.set.final$np,AICc = -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1))),ntraits=1, solution=solution, solution.se=solution.se, index.mat=model.set.final$index.matrix, opts=opts, data=data, phy=phy, states=lik.anc$lik.anc.states, tip.states=tip.states, iterations=out$iterations, eigval=eigval, eigvect=eigvect) 
	class(obj)<-"raydisc"
	return(obj)
}


#Print function
print.raydisc<-function(x,...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")
	
	param.est<- x$solution
	cat("Rates\n")
	print(param.est)
	cat("\n")
	
	if(any(x$eigval<0)){
		index.matrix <- x$index.mat
		#If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
		if (any(x$eigval<0)) {
			cat("The objective function may be at a saddle point", "\n")
		}
	}
	else{
		cat("Arrived at a reliable solution","\n")
	}
}

dev.raydisc<-function(p,phy,liks,Q,rate,root.p){

	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	TIPS <- 1:nb.tip
	comp <- numeric(nb.tip + nb.node)
	phy <- reorder(phy, "pruningwise")
	#Obtain an object of all the unique ancestors
	anc <- unique(phy$edge[,1])

	if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)

	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	for (i  in seq(from = 1, length.out = nb.node)) {
		#the ancestral node at row i is called focal
		focal <- anc[i]
		#Get descendant information of focal
		desRows<-which(phy$edge[,1]==focal)
		desNodes<-phy$edge[desRows,2]
		v <- 1

		for (desIndex in sequence(length(desRows))){
			v<-v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
		}
		comp[focal] <- sum(v)
		liks[focal, ] <- v/comp[focal]
	}
	root <- nb.tip + 1L	
	if (is.na(sum(log(comp[-TIPS])))){return(1000000)}
	else{
		#If root.p!=0 then will fix root probabilities according to FitzJohn et al 2009 Eq. 10.
		if (is.null(root.p)){
			-sum(log(comp[-TIPS]))
		}
		else{				
			-sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,]))
		}
	}	
}

rate.cat.set.oneT<-function(phy,data,model,par.drop,par.eq){
	
	k <- 1
	factored <- factorData(data)

#	nl=2 #We need to accommodate traits >2 states
	nl <- ncol(factored)

	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode

	#Using one of three pre-defined models instead of user-defined model
	if (is.character(model)) {
		#rate is a matrix of rate categories (not actual rates)
		rate <- matrix(NA, nl, nl) #An nl x nl matrix, filled with NA values
		#np is the number of parameters in the rate matrix
		#Equal rates model, one parameter, all rate categories enumerated as '1'
		if (model == "ER") np <- rate[] <- 1
		#All rates different model, number of different parameters just nl x (nl-1)
		#rates are enumerated up to np (e.g. for a two-state character there are 2 rate categories)
		if (model == "ARD") {
			np <- nl*(nl - 1)
			rate[col(rate) != row(rate)] <- 1:np
		}
		#Symmetrical model, like ARD, with half the values
		if (model == "SYM") {
			np <- nl * (nl - 1)/2
			sel <- col(rate) < row(rate)
			rate[sel] <- 1:np
			#Use transpose of the rate category matrix to finish enumerating:
			rate <- t(rate)
			rate[sel] <- 1:np
		}
	} else { #Using user-defined rate matrix
		if (ncol(model) != nrow(model))
			stop("the matrix given as `model' is not square")
		if (ncol(model) != nl)
			stop("the matrix `model' must have as many rows as the number of categories in `x'")
		#Model matrix is OK, and categories should already be enumerated
		rate <- model
		np <- max(rate,na.rm=TRUE)
	}
	index.matrix <- rate
	tmp <- cbind(1:nl, 1:nl)
	rate[tmp] <- 0
	rate[rate == 0] <- np + 1 # Got to make sure to assign np + 1 to empty rate categories

	stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip
	for(column in 1:nl){
		stateTable <- cbind(stateTable,factored[,column])
	}
	colnames(stateTable) <- colnames(factored)

	ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
	liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
	rownames(liks) <- NULL

	Q <- matrix(0, nl^k, nl^k)

	obj$np<-np
	obj$rate<-rate
	obj$index.matrix<-index.matrix
	obj$liks<-liks
	obj$Q<-Q
	
	return(obj)

}

#########################
#    match.tree.data    #
#########################
# Compares a tree and data to make sure they include the same taxa
# Taxa which are in the tree, but not the data matrix, are added to the matrix and coded as missing data.
# Any taxa in the data matrix which are not in the tree are removed from the matrix
# The function returns an object with three parts:
#	$phy: the tree
#	$data: the matrix, omitting taxa not in tree and taxa that were present in the tree but not in the matrix
#	$message.data: a brief message explaining modifications (if any) to the data
#	$message.tree: a brief message explaining modificatoins (if any) to the tree
match.tree.data <- function(phy, data){
	matchobj <- NULL
	matchobj$phy <- phy
	matchobj$data <- data
	matchobj$message.data <- NULL
	matchobj$message.tree <- NULL
	# First look at data matrix to see if each taxon in matrix is also in tree
	missing.fromtree <- NULL
	for(datarow in 1:length(data[,1])){
		if(is.na(match(data[datarow,1],phy$tip.label))){
			missing.fromtree <- c(missing.fromtree,datarow)
		}
	}
	if(length(missing.fromtree) > 0){ # At least one taxa is listed in the matrix, but is not in the tree
		# Make message so user knows taxa have been removed
		matchobj$message.data <- "The following taxa in the data matrix were not in the tree and were excluded from analysis: "
		first <- TRUE
		for(toRemove in 1:length(missing.fromtree)){
			if(first){
				matchobj$message.data <- paste(matchobj$message.data,as.character(data[missing.fromtree[toRemove],1]),sep="")
				first <- FALSE
			} else { #not the first one, so add leading comma
				matchobj$message.data <- paste(matchobj$message.data,", ",as.character(data[missing.fromtree[toRemove],1]),sep="")
			}
		}
		matchobj$data <- data[-missing.fromtree,] # omits those data rows which have no match in the tree
		for(datacol in 2:length(matchobj$data[1,])){
			matchobj$data[,datacol] <- factor(matchobj$data[,datacol]) # have to use factor to remove any factors not present in the final dataset
		}
	}

	missing.taxa <- NULL
	for(tip in 1:length(phy$tip.label)){
		if(is.na(match(phy$tip.label[tip],matchobj$data[,1]))){
			if(is.null(matchobj$message.tree)){ # The first missing taxon
				missing.taxa <- as.character(phy$tip.label[tip])
				matchobj$message.tree <- "The following taxa were in the tree but did not have corresponding data in the data matrix.  They are coded as missing data for subsequent analyses: "
			} else { # not the first missing taxon, add with leading comma
				missing.taxa <- paste(missing.taxa,", ",as.character(phy$tip.label[tip]),sep="")
			}
			# missing taxa will be coded as having missing data "?"
			addtaxon <- as.character(phy$tip.label[tip])
			numcols <- length(matchobj$data[1,])
			newrow <- matrix(as.character("\x3F"),1,numcols) # absurd, but it works
			newrow[1,1] <- addtaxon
			newrowdf <- data.frame(newrow)
			colnames(newrowdf) <- colnames(matchobj$data)
			matchobj$data <- rbind(matchobj$data,newrowdf)
		}
	}
	rownames(matchobj$data) <- matchobj$data[,1] # Use first column (taxon names) as row names
	matchobj$data <- matchobj$data[matchobj$phy$tip.label,] # Sort by order in tree
	rownames(matchobj$data) <- NULL # remove row names after sorting
	if(!is.null(missing.taxa)){
		matchobj$message.tree <- paste(matchobj$message.tree,missing.taxa,sep="")
	}
	return(matchobj)
}

##############
#  findAmps  #
##############
# A function to find positions of ampersands for separating different states.  
# Will allow character state to be greater than one character long.
findAmps <- function(string){
	if(!is.character(string)) return(NULL)
	locs <- NULL # Will hold location values
	for(charnum in 1:nchar(as.character(string))){
		if(substr(string,charnum,charnum) == "&"){
			locs <- c(locs,charnum)
		}
	}
	return(locs)
}

##############
# factorData #
##############
# Function to make factored matrix as levels are discovered.
factorData <- function(data,whichchar=1){
	charcol <- whichchar+1
#	charcol <- whichchar
	factored <- NULL # will become the matrix.  Starts with no data.
	lvls <- NULL
	numrows <- length(data[,charcol])
	missing <- NULL

	for(row in 1:numrows){
		currlvl <- NULL
		levelstring <- as.character(data[row,charcol])
		ampLocs <- findAmps(levelstring)
		if(length(ampLocs) == 0){ #No ampersands, character is monomorphic
			currlvl <- levelstring
			if(currlvl == "?" || currlvl == "-"){ # Check for missing data
				missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
			}
			else { # Not missing data
				if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
					if(length(factored) == 0){ # Matrix is empty, need to create it
						factored <- matrix(0,numrows,1)
						colnames(factored) <- currlvl
#						rownames(factored) <- data[,1] # data object should already have taxa names assinged to rownames
						rownames(factored) <- rownames(data)
					} else { # matrix already exists, but need to add a column for the new level
						zerocolumn <- rep(0,numrows)
						factored <- cbind(factored, zerocolumn)
						colnames(factored)[length(factored[1,])] <- currlvl
					}
					lvls <- c(lvls,currlvl) # add that level to the list
				} # already found this level in another state.  Set the value to one
					whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
					factored[row,whichlvl] <- 1
			}
		} else { #At least one ampersand found, polymorphic character
			start <- 1
			numlvls <- length(ampLocs)+1
			for(part in 1:numlvls){
				# Pull out level from levelstring
				if(part <= length(ampLocs)){ # Haven't reached the last state
					currlvl <- substr(levelstring,start,(ampLocs[part]-1)) # pull out value between start and the location-1 of the next ampersand
				} else { # Final state in list
					currlvl <- substr(levelstring,start,nchar(levelstring)) # pull out value between start and the last character of the string
				}
				if(currlvl == "?" || currlvl == "-"){ # Missing data, but polymorphic?
					missing <- c(missing,row) # add to list of taxa with missing values, will fill in entire row later
				}
				else { # Not missing data
					if(length(which(lvls == currlvl)) == 0){# encountered a level not seen yet
						if(length(factored) == 0){ # Matrix is empty, need to create it
							factored <- matrix(0,numrows,1)
							colnames(factored) <- currlvl
#							rownames(factored) <- data[,1] # data object should already have taxa names assinged to rownames
							rownames(factored) <- rownames(data)
						} else { # matrix already exists, but need to add a column for the new level
							zerocolumn <- rep(0,numrows)
							factored <- cbind(factored, zerocolumn)
							colnames(factored)[length(factored[1,])] <- currlvl
						}
						lvls <- c(lvls,currlvl) # add that level to the list
					} # already found this level in another state.  Set the value to one
						whichlvl <- which(lvls == currlvl) # this index number should correspond to the column number of the state
						factored[row,whichlvl] <- 1
					start <- ampLocs[part] + 1
				}
			}
		}
	}
	# Need to deal with any rows with missing data; fill in NA for all columns for that row
	for(missingrows in 1:length(missing)){
		for(column in 1:length(factored[1,])){
			factored[missing[missingrows],column] <- 1 # All states equally likely
		}
	}
	factored <- factored[,order(colnames(factored))]
	return(factored)
}

###################
#   plot.recons   #
###################
# For plotting ancestral state reconstructions with pie-charts for relative likelihoods
# 'likelihoods' arguement is the $states object from corDISC or rayDISC or the $lik.anc
# object from ancestral state reconstructions, from ape (ace$lik.anc), corEvol
# (corEvol$lik.anc), polyStates (polyStates$lik.anc)
plot.recons <- function(phy, likelihoods, piecolors=NULL, cex=0.5, file=NULL, height=11, width=8.5, showtiplabels=TRUE,title=NULL){
	if(is.null(piecolors)){
		piecolors=c("white","black","red","yellow","forestgreen","blue","coral","aquamarine","darkorchid","gold")
	}
	if(!is.null(file)){
		pdf(file, height=height, width=width,useDingbats=FALSE)
	}
	plot(phy, cex=cex, show.tip.label=showtiplabels)

	if(!is.null(title)){
		title(main=title)
	}
	nodelabels(pie=likelihoods,piecol=piecolors, cex=cex)
	states <- colnames(likelihoods)
	legend(x="topleft", states, cex=0.8, pt.bg=piecolors,col="black",pch=21);

	if(!is.null(file)){
		dev.off()
	}
}

