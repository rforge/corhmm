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
require(rgenoud)
source("recon.joint.R")
source("recon.marginal.R")

corHMM<-function(phy, data, rate.cat, node.states=c("joint", "marginal"), method=c("genoud", "subplex"), p=NULL, par.drop=NULL, par.eq=NULL, root.p=NULL, ip=NULL, nstarts=10, n.cores=NULL){
	
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
	rate.cat=rate.cat
	par.drop=par.drop
	par.eq=par.eq
	root.p=root.p
	nstarts=nstarts
	ip=ip
	
	model.set.final<-rate.cat.set(phy=phy,data=data,rate.cat=rate.cat,par.drop=par.drop,par.eq=par.eq)
	lower = rep(0.00001, model.set.final$np)
	upper = rep(1000, model.set.final$np)
	
	if(method=="genoud"){
		if(!is.null(p)){
			cat("Calculating likelihood from a set of fixed parameters", "\n")
			out<-NULL
			est.pars<-p
			out$objective<-dev.corhmm(est.pars,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=model.set.final$root.p)
		}
		else{
			cat("Initializing...", "\n")
			model.set.init<-rate.cat.set(phy=phy,data=data,rate.cat=1,par.drop=NULL,par.eq=c(1,2))
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)
			dat<-as.matrix(data)
			dat<-phyDat(dat,type="USER", levels=c("0","1"))
			par.score<-parsimony(phy, dat, method="fitch")
			mean = par.score/tl
			if(mean<0.1){
				mean=0.1
			}
			ip<-rexp(1, mean)
			lower = rep(0.00001, model.set.init$np)
			upper = rep(1000, model.set.init$np)
			init = nloptr(x0=rep(ip, length.out = model.set.init$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,root.p=root.p)
			
			cat("Finished. Begin thorough search...", "\n")
			Domains<-cbind(lower,upper)
			starting.values=rep(init$solution,length.out = model.set.final$np)
			out<-genoud(fn=dev.corhmm, starting.values=starting.values, nvars=model.set.final$np, print.level=0, boundary.enforcement=2, Domains=Domains, wait.generations=20, max.generations=100, pop.size=1000, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			loglik <- -out$value
			est.pars<-out$par
		}
	}

	if(method=="subplex"){
		opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)
		if(!is.null(p)){
			cat("Calculating likelihood from a set of fixed parameters", "\n")
			out<-NULL
			out$solution<-p
			out$objective<-dev.corhmm(out$solution,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=model.set.final$root.p)
			loglik <- -out$objective
			est.pars<-out$solution
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
					starts<-rexp(model.set.final$np, mean)
					ip = starts
					out = nloptr(x0=rep(ip, length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)			
					#Initializes a logfile, tmp, of the likelihood for different starting values. A quasi check-point in case computer gets disrupted during an analysis
					tmp = matrix(,1,ncol=(1+model.set.final$np))
					tmp[,1] = out$objective
					tmp[,2:(model.set.final$np+1)] = starts
					for(i in 2:nstarts){
						starts<-rexp(model.set.final$np, mean)
						out.alt = nloptr(x0=rep(starts, length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
						tmp[,1] = out.alt$objective
						tmp[,2:(model.set.final$np+1)] = starts
						if(out.alt$objective < out$objective){
							out = out.alt
							ip = starts
						}
						else{
							out = out
							ip = ip
						}
					}
					loglik <- -out$objective
					est.pars<-out$solution
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
						tmp = matrix(,1,ncol=(1+model.set.final$np))
						starts<-rexp(model.set.final$np, mean)
						out = nloptr(x0=rep(starts, length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
						tmp[,1] = out$objective
						tmp[,2:(model.set.final$np+1)] = out$solution
						tmp
					}
					restart.set<-mclapply(1:nstarts,random.restart, mc.cores=n.cores)
					#Finds the best fit within the restart.set list
					best.fit<-which.min(unlist(lapply(1:nstarts,function(i) lapply(restart.set[[i]][,1],min))))
					#Generates an object to store results from restart algorithm:
					out<-NULL
					out$objective=unlist(restart.set[[best.fit]][,1])
					out$solution=unlist(restart.set[[best.fit]][,2:(model.set.final$np+1)])
					loglik <- -out$objective
					est.pars<-out$solution
				}
			}
			#If a user-specified starting value(s) is supplied:
			else{
				cat("Begin subplex optimization routine -- Starting value(s):", ip, "\n")
				ip=ip
				out = nloptr(x0=rep(ip, length.out = model.set.final$np), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
				loglik <- -out$objective
				est.pars<-out$solution
			}
		}
	}

	#Starts the summarization process:
	cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
	
	if (node.states == "marginal"){
		lik.anc <- recon.marginal(phy, data, est.pars, hrm=TRUE, rate.cat, ntraits=NULL, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
		pr<-apply(lik.anc$lik.anc.states,1,which.max)
		phy$node.label <- pr
		phy$node.label <- 1:nb.node
		tip.states <- NULL
	}
	if (node.states == "joint"){
		lik.anc <- recon.joint(phy, data, est.pars, hrm=TRUE, rate.cat, ntraits=NULL,par.drop=par.drop, par.eq=par.eq, root.p=root.p)
		phy$node.label <- lik.anc$lik.anc.states
		tip.states <- lik.anc$lik.tip.states
	}
	
	cat("Finished. Performing diagnostic tests.", "\n")
	
	#Approximates the Hessian using the numDeriv function
	h <- hessian(func=dev.corhmm, x=est.pars, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
	solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
	solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix], dim(model.set.final$index.matrix))

	if (rate.cat == 1){
		rownames(solution) <- rownames(solution.se) <- c("(0)","(1)")
		colnames(solution) <- colnames(solution.se) <- c("(0)","(1)")			
		#Initiates user-specified reconstruction method:
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				colnames(lik.anc$lik.anc.states) <- c("P(0)","P(1)")
			}
		}
	}
	if (rate.cat == 2){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){		
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
			}
		}
	}
	if (rate.cat == 3){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){		
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
			}
		}
	}
	if (rate.cat == 4){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){	
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
			}
		}
	}
	if (rate.cat == 5){
		rownames(solution) <- rownames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
		colnames(solution) <- colnames(solution.se) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){	
				colnames(lik.anc$lik.anc.states) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
			}
		}
	}
	hess.eig <- eigen(h,symmetric=TRUE)
	eigval<-signif(hess.eig$values,2)
	eigvect<-round(hess.eig$vectors, 2)
	obj = list(loglik = loglik, AIC = -2*loglik+2*model.set.final$np,AICc = -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1))),rate.cat=rate.cat,solution=solution, solution.se=solution.se, index.mat=model.set.final$index.matrix, opts=opts, data=data, phy=phy, states=lik.anc$lik.anc.states, tip.states=tip.states, iterations=out$iterations, eigval=eigval, eigvect=eigvect) 
	class(obj)<-"corhmm"
	return(obj)
}

#Print function
print.corhmm<-function(x,...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,x$rate.cat,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","Rate.cat","ntax")
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

#Generalized ace() function that allows analysis to be carried out when there are polytomies:
dev.corhmm <- function(p,phy,liks,Q,rate,root.p) {
	
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

rate.cat.set<-function(phy,data,rate.cat,par.drop,par.eq){
	
	k=2
	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	obj$rate.cat<-rate.cat
	#Builds the rate matrix based on the specified rate.cat. Not exactly the best way
	#to go about this, but it is the best I can do for now -- it works, so what me worry?
	if (rate.cat == 1){
		rate <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- 2
		index<-matrix(TRUE,k*rate.cat,k*rate.cat)
		diag(index) <- FALSE
		rate[index] <- 1:np
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
	obj$np<-np
	obj$rate<-rate
	obj$index.matrix<-index.matrix
	obj$liks<-liks
	obj$Q<-Q
	
	obj
}


