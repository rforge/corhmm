#COVARION MODEL OF ANCESTRAL STATE RECONSTRUCTION OF A SINGLE BINARY TRAIT

#written by Jeremy M. Beaulieu, though the basic algorithm is a modified version of ace() 
#found in the 'ape' package written and maintained by Emmanual Paradis.

#Takes a tree and a trait file and estimates the rate of dependent transitions between
#a binary character. The first column of the trait file must contain the species labels 
#to match to the tree, with the second corresponding to the binary trait of interest. 
#Can test models that assume 1, 2, 3, 4, or 5 rate categories underlying the observed data. 
#Setting the Probs=TRUE will output the relative probabilities of the ancestral state at 
#each node, and the user can fix the root state probabilities by supplying a vector to the
#root.p=c() option -- NOTE: At the moment fixing the root state can causing issues with 
#Hessian calculation and the program fails. In all other cases, the program assumes the marginal 
#probability for the root. The user can also drop particular parameters from the model by 
#specifying the parameter number in a vector supplied to par.drop=c(), also the user can fix
#particular parameters to have the same rate. However, at the moment both of these features are 
#not documented and may be difficult to use (it is a work in progress). 

#A plethora of dependencies:
require(ape)
require(nloptr)
require(numDeriv)
require(expm)
require(corpcor)
require(phangorn)
require(multicore)
source("recon.joint.R")
source("recon.marginal.R")

corHMM<-function(phy, data, rate.cat, nstarts=10, n.cores=NULL, node.states=c("joint", "marginal"), p=NULL, par.drop=NULL, par.eq=NULL, root.p=NULL, ip=NULL){
	
	#Creates the data structure and orders the rows to match the tree. 
	phy$edge.length[phy$edge.length<=1e-5]=1e-5
	data <- data.frame(data[,2], data[,2],row.names=data[,1])
	data <- data[phy$tip.label,]
	#Have to collect this here. When you reorder, the branching time function is not correct:
	tl<-max(branching.times(phy))
	
	#Some initial values for use later
	k=2
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	
	par.drop=par.drop
	par.eq=par.eq
	root.p=root.p
	
	#Builds the rate matrix based on the specified rate.cat. Not exactly the best way
	#to go about this, but it is the best I can do for now -- it works, so what me worry?
	if (rate.cat == 1){
		rate <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- 2
		
		index<-matrix(TRUE,k*rate.cat,k*rate.cat)
		diag(index) <- FALSE
		rate[index] <- 1:np
		index.matrix <- rate
		
		diag(rate)<-0
		rate[rate == 0] <- np + 1
	}
	if (rate.cat == 2){
		rate <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- 8
		tmp <- cbind(1:(k*rate.cat), (k*rate.cat):1)
		tmp2 <- cbind(1:(k*rate.cat), 1:(k*rate.cat))
		
		index <- matrix(TRUE,k*rate.cat,k*rate.cat)
		diag(index) <- FALSE
		index[tmp] <- FALSE
		index[tmp2] <- FALSE
		rate[index] <- 1:np
		#If par.drop is not null will adjust the rate matrix
		if(!is.null(par.drop)==TRUE){
			for(i in 1:length(par.drop)){
				tmp3 <- which(rate==par.drop[i], arr.ind=T)
				index[tmp3] <- FALSE
				rate[tmp3] <- 0
			}
			np <- np-length(par.drop)
			rate[index] <- 1:np
		}
		#If par.eq is not null then pairs of parameters are set equal to each other.
		if(!is.null(par.eq)==TRUE){
			for (i  in seq(from = 1, by = 2, length.out = length(par.eq)/2)) {
				j<-i+1
				tmp3 <- which(rate==par.eq[j], arr.ind=T)
				index[tmp3] <- FALSE
				rate[tmp3] <- 0
				np <- np-1
				rate[index] <- 1:np
				rate[tmp3] <- par.eq[i]
			}
		}
		index.matrix <- rate
		index.matrix[index.matrix == 0] = NA
		
		rate[tmp] <- 0
		rate[tmp2] <- 0
		rate[rate == 0] <- np + 1
	}
	if (rate.cat == 3){
		rate <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- 14
		tmp <-	c(4,5,6,3,5,6,2,6,1,5,1,2,4,1,2,3)
		tmp2 <- c(1,1,1,2,2,2,3,3,4,4,5,5,5,6,6,6)
		tmp3 <- cbind(tmp,tmp2)
		
		index <- matrix(TRUE,(k*rate.cat),(k*rate.cat))
		diag(index) <- FALSE
		index[tmp3] <- FALSE			
		rate[index] <- 1:np
		#If par.drop is not null will adjust the rate matrix
		if(!is.null(par.drop)==TRUE){
			for(i in 1:length(par.drop)){
				tmp4 <- which(rate==par.drop[i], arr.ind=T)
				index[tmp4] <- FALSE
				rate[tmp4] <- 0
			}
			np <- np-length(par.drop)
			rate[index] <- 1:np
		}
		#If par.eq is not null then pairs of parameters are set equal to each other.
		if(!is.null(par.eq)==TRUE){
			for (i  in seq(from = 1, by = 2, length.out = length(par.eq)/2)) {
				j <- i+1
				tmp4 <- which(rate==par.eq[j], arr.ind=T)
				index[tmp4] <- FALSE
				rate[tmp4] <- 0
				np <- np-1
				rate[index] <- 1:np
				rate[tmp4] <- par.eq[i]
				par.eq<-par.eq-1
			}
		}		
		index.matrix <- rate
		index.matrix[index.matrix == 0] = NA
		
		diag(rate) <- 0
		rate[tmp3] <- 0
		rate[rate == 0] <- np + 1
	}
	if (rate.cat == 4){
		rate <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- 20
		tmp <- c(4,5,6,7,8,3,5,6,7,8,2,6,7,8,1,5,7,8,1,2,4,8,1,2,3,7,1,2,3,4,6,1,2,3,4,5)
		tmp2 <-c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,8,8,8,8,8)
		tmp3 <- cbind(tmp,tmp2)
		
		index <- matrix(TRUE,(k*rate.cat),(k*rate.cat))
		diag(index) <- FALSE
		index[tmp3] <- FALSE			
		rate[index] <- 1:np
		#If par.drop is not null will adjust the rate matrix
		if(!is.null(par.drop)==TRUE){
			for(i in 1:length(par.drop)){
				tmp4 <- which(rate==par.drop[i], arr.ind=T)
				index[tmp4] <- FALSE
				rate[tmp4] <- 0
			}
			np <- np-length(par.drop)
			rate[index] <- 1:np
		}
		#If par.eq is not null then pairs of parameters are set equal to each other.
		if(!is.null(par.eq)==TRUE){
			for (i  in seq(from = 1, by = 2, length.out = length(par.eq)/2)) {
				j<-i+1
				tmp4 <- which(rate==par.eq[j], arr.ind=T)
				index[tmp4] <- FALSE
				rate[tmp4] <- 0
				np <- np-1
				rate[index] <- 1:np
				rate[tmp4] <- par.eq[i]
			}
		}		
		index.matrix <- rate
		index.matrix[index.matrix == 0] = NA
		
		diag(rate) <- 0
		rate[tmp3]<-0
		rate[rate == 0] <- np + 1
	}
	if (rate.cat == 5){
		rate <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- 26
		tmp <- c(4,5,6,7,8,9,10,3,5,6,7,8,9,10,2,6,7,8,9,10,1,5,7,8,9,10,1,2,4,8,9,10,1,2,3,7,9,10,1,2,3,4,6,10,1,2,3,4,5,8,9,1,2,3,4,5,6,8,1,2,3,4,5,6,7)
		tmp2 <-c(1,1,1,1,1,1, 1,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,8,9,9,9,9,9,9,9,10,10,10,10,10,10,10)		
		tmp3 <- cbind(tmp,tmp2)
		
		index <- matrix(TRUE,(k*rate.cat),(k*rate.cat))
		diag(index) <- FALSE
		index[tmp3] <- FALSE			
		rate[index] <- 1:np
		#If par.drop is not null will adjust the rate matrix
		if(!is.null(par.drop)==TRUE){
			for(i in 1:length(par.drop)){
				tmp4 <- which(rate==par.drop[i], arr.ind=T)
				index[tmp4] <- FALSE
				rate[tmp4] <- 0
			}
			np <- np-length(par.drop)
			rate[index] <- 1:np
		}
		#If par.eq is not null then pairs of parameters are set equal to each other.
		if(!is.null(par.eq)==TRUE){
			for (i  in seq(from = 1, by = 2, length.out = length(par.eq)/2)) {
				j<-i+1
				tmp4 <- which(rate==par.eq[j], arr.ind=T)
				index[tmp4] <- FALSE
				rate[tmp4] <- 0
				np <- np-1
				rate[index] <- 1:np
				rate[tmp4] <- par.eq[i]
			}
		}		
		index.matrix <- rate
		index.matrix[index.matrix == 0] = NA
		
		diag(rate) <- 0
		rate[tmp3] <- 0
		rate[rate == 0] <- np + 1
	}
	#Makes a matrix of tip states and empty cells corresponding 
	#to ancestral nodes during the optimization process.	
	x <- data[,1]
	TIPS <- 1:nb.tip
	
	for(i in 1:nb.tip){
		if(is.na(x[i])){x[i]=2}
	}
	if (rate.cat == 1){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,1]=1}
			if(x[i]==1){liks[i,2]=1}
			if(x[i]==2){liks[i,1:2]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 2){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3)]=1}
			if(x[i]==1){liks[i,c(2,4)]=1}
			if(x[i]==2){liks[i,1:4]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 3){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5)]=1}
			if(x[i]==1){liks[i,c(2,4,6)]=1}
			if(x[i]==2){liks[i,1:6]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 4){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8)]=1}
			if(x[i]==2){liks[i,1:8]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 5){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7,9)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8,10)]=1}
			if(x[i]==2){liks[i,1:10]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	phy <- reorder(phy, "pruningwise")
	comp <- numeric(nb.tip + nb.node)
	#Obtain an object of all the unique ancestors
	anc <- unique(phy$edge[,1])
	#Generalized ace() function that allows analysis to be carried out when there are polytomies:
	dev <- function(p) {
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
			#Loops through all descendants of focal (how we deal with polytomies):
			for (desIndex in sequence(length(desRows))){
				v<-v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
			}
			#Sum the likelihoods:
			comp[focal] <- sum(v)
			#Divide each likelihood by the sum to obtain probabilities:
			liks[focal, ] <- v/comp[focal]
		}
		#Specifies the root:
		root <- nb.tip + 1L
		#If any of the logs have NAs restart search:
		if (is.na(sum(log(comp[-TIPS])))){return(1000000)}
		else{
			if (is.null(root.p)){
				-sum(log(comp[-TIPS]))
			}
			#root.p!=NULL, will fix root probabilities according to FitzJohn et al 2009 Eq. 10.
			else{				
				#Interesting development -- must have a non-zero rate in order to properly fix the root!
				-sum(log(comp[-TIPS])) + log(sum(root.p * liks[root,]))
			}
		}	
	}
	
	#Sets the bounds on the parameter search:
	lower = rep(0.00001, np)
	upper = rep(100, np)
	
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25, "xtol_rel"=.Machine$double.eps^0.25)
	
	if(!is.null(p)){
		cat("Calculating likelihood from a set of fixed parameters", "\n")
		out<-NULL
		out$solution<-p
		out$objective<-dev(out$solution)
	}
	else{	   
		#If a user-specified starting value(s) is not supplied this begins loop through a set of randomly chosen starting values:
		if(is.null(ip)){
			
			cat("Begin thorough optimization search -- performing", nstarts, "random restarts", "\n")
			
			#If the analysis is to be run a single processor:
			if(is.null(n.cores)){
				#Sets parameter settings for random restarts by taking the parsimony score and dividing
				#by the total length of the tree
				dat<-as.matrix(data)
				dat<-phyDat(dat,type="USER", levels=c("0","1"))
				par.score<-parsimony(phy, dat, method="fitch")/2
				mean = par.score/tl
				starts<-rexp(np, mean)
				ip = starts
				out = nloptr(x0=rep(starts, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)			
				#Initializes a logfile, tmp, of the likelihood for different starting values. A quasi check-point in case computer gets disrupted during an analysis
				tmp = matrix(,1,ncol=(1+np))
				tmp[,1] = out$objective
				tmp[,2:(np+1)] = starts
				for(i in 2:nstarts){
					starts<-rexp(np, mean)
					out.alt = nloptr(x0=rep(starts, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)
					tmp[,1] = out.alt$objective
					tmp[,2:(np+1)] = starts
					if(out.alt$objective < out$objective){
						out = out.alt
						ip = starts
					}
					else{
						out = out
						ip = ip
					}
				}
			}
			#If the analysis is to be run on multiple processors:
			else{
				#Sets parameter settings for random restarts by taking the parsimony score and dividing
				#by the total length of the tree
				dat<-as.matrix(data)
				dat<-phyDat(dat,type="USER", levels=c("0","1"))
				par.score<-parsimony(phy, dat, method="fitch")
				mean = par.score/tl
				if(mean<0.1){
					mean=0.1
				}
				random.restart<-function(nstarts){
					tmp = matrix(,1,ncol=(1+np))
					starts<-rexp(np, mean)
					out = nloptr(x0=rep(starts, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)
					tmp[,1] = out$objective
					tmp[,2:(np+1)] = out$solution
					tmp
				}
				restart.set<-mclapply(1:nstarts,random.restart, mc.cores=n.cores)
				#Finds the best fit within the restart.set list
				best.fit<-which.min(unlist(lapply(1:nstarts,function(i) lapply(restart.set[[i]][,1],min))))
				#Generates an object to store results from restart algorithm:
				out<-NULL
				out$objective=unlist(restart.set[[best.fit]][,1])
				out$solution=unlist(restart.set[[best.fit]][,2:(np+1)])
			}
		}
		#If a user-specified starting value(s) is supplied:
		else{
			cat("Begin subplex optimization routine -- Starting value(s):", ip, "\n")
			ip=ip
			out = nloptr(x0=rep(ip, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)
		}
	}
	
	#Starts the summarization process:
	cat("Finished. Performing diagnostic tests.", "\n")
	
	obj$loglik <- -out$objective	
	obj$AIC <- -2*obj$loglik+2*np
	obj$AICc <- -2*obj$loglik+(2*np*(nb.tip/(nb.tip-np-1)))
	#Approximates the Hessian using the numDeriv function
	h <- hessian(x=out$solution, func=dev)
	#Initiates the summary process
	if (rate.cat == 1){
		obj$Param.est<- matrix(out$solution[index.matrix], dim(index.matrix))
		#Calculates the standard error of the parameter by taking the sqrt of diagonal of the inverse of the Hessian: 
		obj$Param.SE <- matrix(sqrt(diag(pseudoinverse(h)))[index.matrix], dim(index.matrix))
		rownames(obj$Param.est) <- rownames(obj$Param.SE) <- c("(0)","(1)")
		colnames(obj$Param.est) <- colnames(obj$Param.SE) <- c("(0)","(1)")			
		#Initiates user-specified reconstruction method:
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				lik.anc <- recon.marginal(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				colnames(lik.anc$lik.anc.states) <- c("P(0)","P(1)")
				write.table(lik.anc$lik.anc.states, file="Anc.EstimatesHMM1cat.xls", quote=FALSE, sep="\t")
				pr<-apply(lik.anc$lik.anc.states,1,which.max)
				phy$node.label <- pr
				write.tree(phy, file="AncReconStatesHMM1cat.tre", append=TRUE)
				phy$node.label <- 1:nb.node
				write.tree(phy, file="AncReconKey.tre")
			}
			if (node.states == "joint"){
				lik.anc <- recon.joint(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL,par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				phy$node.label <- lik.anc$lik.anc.states
				write.tree(phy, file="AncReconStatesHMM1cat.tre", append=TRUE)
			}
		}
	}
	if (rate.cat == 2){
		obj$Param.est<- matrix(out$solution[index.matrix], dim(index.matrix))
		obj$Param.SE <- matrix(sqrt(diag(pseudoinverse(h)))[index.matrix], dim(index.matrix))
		rownames(obj$Param.est) <- rownames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
		colnames(obj$Param.est) <- colnames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){		
				lik.anc <- recon.marginal(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
				write.table(lik.anc$lik.anc.states, file="Anc.EstimatesHMM2cat.xls", quote=FALSE, sep="\t")
				pr<-apply(lik.anc$lik.anc.states,1,which.max)
				phy$node.label<-pr
				write.tree(phy, file="AncReconStatesHMM2cat.tre", append=TRUE)
				phy$node.label<-1:nb.node
				write.tree(phy, file="AncReconKey.tre")
			}
			if (node.states == "joint"){
				lik.anc <- recon.joint(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				write.table(cbind(row.names(data), lik.anc$lik.tip.states), file="Tipstates.HMM2cat.xls", quote=FALSE, row.names=F, sep="\t")
				phy$node.label <- lik.anc$lik.anc.states
				write.tree(phy, file="AncReconStatesHMM2cat.tre", append=TRUE)
			}			
		}
	}
	if (rate.cat == 3){
		obj$Param.est<- matrix(out$solution[index.matrix], dim(index.matrix))
		obj$Param.SE <- matrix(sqrt(diag(pseudoinverse(h)))[index.matrix], dim(index.matrix))
		rownames(obj$Param.est) <- rownames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
		colnames(obj$Param.est) <- colnames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){		
				lik.anc <- recon.marginal(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
				write.table(lik.anc$lik.anc.states, file="Anc.EstimatesHMM3cat.xls", quote=FALSE, sep="\t")
				pr<-apply(lik.anc$lik.anc.states,1,which.max)
				phy$node.label <- pr
				write.tree(phy, file="AncReconStatesHMM3cat.tre", append=TRUE)
				phy$node.label <- 1:nb.node
				write.tree(phy, file="AncReconKey.tre")
			}
			if (node.states == "joint"){
				lik.anc <- recon.joint(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				write.table(cbind(row.names(data),lik.anc$lik.tip.states), file="Tipstates.HMM3cat.xls", row.names=F, quote=FALSE, sep="\t")
				phy$node.label <- lik.anc$lik.anc.states
				write.tree(phy, file="AncReconStatesHMM3cat.tre", append=TRUE)
			}			
		}
	}
	if (rate.cat == 4){
		obj$Param.est <- matrix(out$solution[index.matrix], dim(index.matrix))
		obj$Param.SE <- matrix(sqrt(diag(pseudoinverse(h)))[index.matrix], dim(index.matrix))
		rownames(obj$Param.est) <- rownames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
		colnames(obj$Param.est) <- colnames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){	
				lik.anc <- recon.marginal(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
				write.table(lik.anc$lik.anc.states, file="Anc.EstimatesHMM4cat.xls", quote=FALSE, sep="\t")
				pr<-apply(lik.anc$lik.anc.states,1,which.max)
				phy$node.label <- pr
				write.tree(phy, file="AncReconStatesHMM4cat.tre", append=TRUE)
				phy$node.label <- 1:nb.node
				write.tree(phy, file="AncReconKey.tre")
			}
			if (node.states == "joint"){
				lik.anc <- recon.joint(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				write.table(cbind(row.names(data), lik.anc$lik.tip.states), file="Tipstates.HMM4cat.xls", row.names=F, quote=FALSE, sep="\t")
				phy$node.label <- lik.anc$lik.anc.states
				write.tree(phy, file="AncReconStatesHMM4cat.tre", append=TRUE)
			}			
		}
	}
	if (rate.cat == 5){
		obj$Param.est<- matrix(out$solution[index.matrix], dim(index.matrix))
		obj$Param.SE <- matrix(sqrt(diag(pseudoinverse(h)))[index.matrix], dim(index.matrix))
		rownames(obj$Param.est) <- rownames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
		colnames(obj$Param.est) <- colnames(obj$Param.SE) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){	
				lik.anc <- recon.marginal(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
				write.table(lik.anc$lik.anc.states, file="Anc.EstimatesHMM5cat.xls", quote=FALSE, sep="\t")
				pr<-apply(lik.anc$lik.anc.states,1,which.max)
				phy$node.label <- pr
				write.tree(phy, file="AncReconStatesHMM5cat.tre", append=TRUE)
				phy$node.label <- 1:nb.node
				write.tree(phy, file="AncReconKey.tre")
			}
			if (node.states == "joint"){
				lik.anc <- recon.joint(phy, data, out$solution, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				write.table(cbind(row.names(data), lik.anc$lik.tip.states), file="Tipstates.HMM5cat.xls", row.names=F, quote=FALSE, sep="\t")
				phy$node.label <- lik.anc$lik.anc.states
				write.tree(phy, file="AncReconStatesHMM5cat.tre", append=TRUE)
			}			
		}
	}
	
	obj$iterations <- out$iterations
	#The eigendecomposition of the Hessian matrix to assess whether or not the function has found the minimum
	hess.eig <- eigen(h,symmetric=TRUE)
	obj$eigval <- signif(hess.eig$values,2)	
	#If all the eigenvalues are positive then the function has found the maximum likelihood solution -- eigen ratio might be a good thing to add: 
	#any parameter with a ratio exceeding 5000 being considered to be not very identifiable?
	if(any(obj$eigval<0)){
		obj$diagnostic <- 'The objective function may be at a saddle point -- check eigenvectors'
		obj$eigvect <- round(hess.eig$vectors, 2)
		obj$index.matrix <- index.matrix
	}
	else{obj$diagnostic<-'Arrived at a reliable solution'}
	
	obj
}

