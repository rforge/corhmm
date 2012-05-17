#CORRELATED EVOLUTION OF TWO or THREE BINARY TRAITS

#written by Jeremy M. Beaulieu

#Takes a tree and a trait file and estimates the rate of dependent transitions between
#two or three binary characters. The first column of the trait file must contain the species labels 
#to match to the tree, with the second and third columns corresponding to x and y, 
#respectively. Can test three models: ER = equal rates;SYM = forw/back equal, ARD = all 
#rates unequal. Setting the Probs=TRUE will output the relative probabilities of the 
#ancestral state at each node.

require(ape)
require(nloptr)
require(numDeriv)
require(expm)
require(corpcor)
require(phangorn)
require(multicore)
source("recon.joint.R")
source("recon.marginal.R")

corDISC<-function(phy,data, ntraits=2, model=c("ER","SYM","ARD"), node.states=c("joint", "marginal"), nstarts=10, n.cores=NULL, p=NULL, par.drop=NULL, par.eq=NULL, root.p=NULL, ip=NULL){
	
	#Creates the data structure and orders the rows to match the tree
	phy$edge.length[phy$edge.length==0]=1e-5
	
	if(ntraits==2){
		data<-data.frame(data[,2], data[,3], row.names=data[,1])
		data<-data[phy$tip.label,]
	}
	if(ntraits==3){
		data<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
		data<-data[phy$tip.label,]
	}
	#Have to collect this here. When you reorder, the branching time function is not correct:
	tl<-max(branching.times(phy))
	#Some initial values for use later - will clean up
	k=ntraits
	nl=2
	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	
	par.drop=par.drop
	par.eq=par.eq
	root.p=root.p	
	
	if(ntraits==2){
		
		if (is.character(model)) {
			
			rate <- matrix(NA, nl^k, nl^k)
			
			if (model == "ER"){
				np <- 1
				tmp <- cbind(1:(nl^k), (nl^k):1)
				tmp2 <- cbind(1:(nl^k), 1:(nl^k))
				
				index<-matrix(TRUE,nl^k,nl^k)
				diag(index)<-FALSE
				index[tmp]<-FALSE
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
			if (model == "SYM") {
				np <- 4
				tmp <- cbind(1:(nl^k), (nl^k):1)
				tmp2 <- cbind(1:(nl^k), 1:(nl^k))
				
				index<-matrix(TRUE,nl^k,nl^k)
				diag(index)<-FALSE
				index[tmp]<-FALSE
				
				rate[index][c(1,2,4,6)] <- rate[index][c(3,5,7,8)] <- 1:np
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
			
			if (model == "ARD") {
				np <- 8
				tmp <- cbind(1:(nl^k), (nl^k):1)
				tmp2 <- cbind(1:(nl^k), 1:(nl^k))
				
				index<-matrix(TRUE,nl^k,nl^k)
				diag(index)<-FALSE
				index[tmp]<-FALSE			
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
		} else {
			if (ncol(model) != nrow(model))
			stop("the matrix given as 'model' is not square")
			if (ncol(model) != nl)
			stop("the matrix 'model' must have as many rows
				 as the number of categories in `x'")
			rate <- model
			np <- max(rate)
		}
		
		x<-data[,1]
		y<-data[,2]
		
		liks <- matrix(0, nb.tip + nb.node, nl^k)
		TIPS <- 1:nb.tip
		
		for(i in 1:nb.tip){
			if(is.na(x[i])){x[i]=2 & y[i]=2}
		}
		for(i in 1:nb.tip){
			if(x[i]==0 & y[i]==0){liks[i,1]=1}
			if(x[i]==0 & y[i]==1){liks[i,2]=1}
			if(x[i]==1 & y[i]==0){liks[i,3]=1}
			if(x[i]==1 & y[i]==1){liks[i,4]=1}
			if(x[i]==2 & y[i]==2){liks[i,1:4]=1}
		}
	}
	if(ntraits==3){
		if (is.character(model)) {
			rate <- matrix(NA, nl^k, nl^k)
			
			if (model == "ER"){
				np <- 1
				tmp <- cbind(1:(nl^k), (nl^k):1)
				tmp2 <- cbind(1:(nl^k), 1:(nl^k))
				col1 <- c(5:7,3:4,8,2,4,8,2,3,8,1,6,7,1,5,7,1,5,6,2,3,4,5,5,5,2,3,8,8,8,6,7)
				col2 <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3),rep(8,3),2,3,8,5,5,5,6,7,8,8)
				tmp3 <- cbind(col1, col2)	
				
				index<-matrix(TRUE,nl^k,nl^k)
				index[tmp2]<-FALSE
				index[tmp]<-FALSE
				index[tmp3]<-FALSE
				
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
				
				rate[tmp] <- 0
				rate[tmp2] <- 0
				rate[tmp3] <- 0
				rate[rate == 0] <- np + 1
			}
			
			if (model == "SYM") {
				np <- 12
				tmp <- cbind(1:(nl^k), (nl^k):1)
				tmp2 <- cbind(1:(nl^k), 1:(nl^k))
				col1 <- c(5:7,3:4,8,2,4,8,2,3,8,1,6,7,1,5,7,1,5,6,2,3,4)
				col2 <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3),rep(8,3))
				tmp3 <- cbind(col1, col2)	
				
				index<-matrix(TRUE,nl^k,nl^k)
				index[tmp2]<-FALSE
				index[tmp]<-FALSE
				index[tmp3]<-FALSE
				rate[index][c(1,2,3,5,6,8,9,11,12,15,18,21)] <- rate[index][c(4,7,10,13,16,14,19,17,20,22,23,24)] <- 1:np
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
				
				rate[tmp] <- 0
				rate[tmp2] <- 0
				rate[tmp3] <- 0
				rate[rate == 0] <- np + 1	
			}
			
			if (model == "ARD") {
				np <- 24
				tmp <- cbind(1:(nl^k), (nl^k):1)
				tmp2 <- cbind(1:(nl^k), 1:(nl^k))
				col1 <- c(5:7,3:4,8,2,4,8,2,3,8,1,6,7,1,5,7,1,5,6,2,3,4)
				col2 <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3),rep(8,3))
				tmp3 <- cbind(col1, col2)	
				
				index<-matrix(TRUE,nl^k,nl^k)
				index[tmp2]<-FALSE
				index[tmp]<-FALSE
				index[tmp3]<-FALSE			
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
				
				rate[tmp] <- 0
				rate[tmp2] <- 0
				rate[tmp3] <- 0
				rate[rate == 0] <- np + 1
			}
		} 
		else {
			if (ncol(model) != nrow(model))
			stop("the matrix given as 'model' is not square")
			if (ncol(model) != nl)
			stop("the matrix `model' must have as many rows
				 as the number of categories in `x'")
			rate <- model
			np <- max(rate)
		}
		x<-data[,1]
		y<-data[,2]
		z<-data[,3]
		
		liks <- matrix(0, nb.tip + nb.node, nl^k)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(is.na(x[i])){x[i]=2 & y[i]=2 & z[i]=2}
		}
		for(i in 1:nb.tip){
			if(x[i]==0 & y[i]==0 & z[i]==0){liks[i,1]=1}
			if(x[i]==1 & y[i]==0 & z[i]==0){liks[i,2]=1}
			if(x[i]==0 & y[i]==1 & z[i]==0){liks[i,3]=1}
			if(x[i]==0 & y[i]==0 & z[i]==1){liks[i,4]=1}
			if(x[i]==1 & y[i]==1 & z[i]==0){liks[i,5]=1}
			if(x[i]==1 & y[i]==0 & z[i]==1){liks[i,6]=1}
			if(x[i]==0 & y[i]==1 & z[i]==1){liks[i,7]=1}
			if(x[i]==1 & y[i]==1 & z[i]==1){liks[i,8]=1}
			if(x[i]==2 & y[i]==2 & z[i]==2){liks[i,1:8]=1}
		}
	}

	phy <- reorder(phy, "pruningwise")
	Q <- matrix(0, nl^k, nl^k)
	# from Rich FitzJohn - attenuates underflow problems.
	comp <- numeric(nb.tip + nb.node) #Storage...
	anc <- unique(phy$edge[,1])
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
	
	lower = rep(0, np)
	upper = rep(100, np)
	
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25, "xtol_rel"=.Machine$double.eps^0.25)

	if(!is.null(p)){
		cat("Calculating likelihood from a set of fixed parameters", "\n")
		out<-NULL
		out$solution<-p
		out$objective<-dev(out$solution)
	}
	else{	   
		if(is.null(ip)){
			cat("Begin thorough optimization search -- performing", nstarts, "random restarts", "\n")
			#If the analysis is to be run a single processor:
			if(is.null(n.cores)){
				#Sets parameter settings for random restarts by taking the parsimony score and dividing
				#by the total length of the tree
				dat<-as.matrix(data)
				dat<-phyDat(dat,type="USER", levels=c("0","1"))
				par.score<-parsimony(phy, dat, method="fitch")
				mean = par.score/tl
				if(mean<0.1){
					mean=0.1
				}			
				starts<-rexp(np, mean)
				ip = starts
				out = nloptr(x0=rep(starts, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)			
				tmp = matrix(,1,ncol=(1+np))
				tmp[,1] = out$objective
				tmp[,2:(np+1)] = out$solution
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
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25, "xtol_rel"=.Machine$double.eps^0.25)
			out = nloptr(x0=rep(ip, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts)
		}
	}
	
	cat("Finished. Performing diagnostic tests.", "\n")
	obj$loglik <- -out$objective
	obj$AIC <- -2*obj$loglik+2*np
	obj$AICc <- -2*obj$loglik+(2*np*(nb.tip/(nb.tip-np-1)))
	
	h <- hessian(x=out$solution, func=dev)
	
	#Initiates user-specified reconstruction method:
	if(ntraits==2){
		obj$Param.est<- matrix(out$solution[index.matrix], dim(index.matrix))
		#Calculates the standard error of the parameter by taking the sqrt of diagonal of the inverse of the Hessian: 
		obj$Param.SE <- matrix(sqrt(diag(pseudoinverse(h)))[index.matrix], dim(index.matrix))
		rownames(obj$Param.est) <- rownames(obj$Param.SE) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
		colnames(obj$Param.est) <- colnames(obj$Param.SE) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				lik.anc <- recon.marginal(phy, data, out$solution, hrm=FALSE, rate.cat=NULL, ntraits=2, model=model, par.drop, par.eq, root.p=root.p)
				colnames(lik.anc$lik.anc.states) <-  c("P(0,0)","P(0,1)","P(1,0)","P(1,1)")
				write.table(lik.anc$lik.anc.states, file="Anc.EstimatesDISCRETE.xls", quote=FALSE, sep="\t")
				pr<-apply(lik.anc$lik.anc.states,1,which.max)
				phy$node.label <- pr
				write.tree(phy, file="AncReconStatesDISCRETE.tre", append=TRUE)
				phy$node.label <- 1:nb.node
				write.tree(phy, file="AncReconKey.tre")
			}
			if (node.states == "joint"){
				lik.anc <- recon.joint(phy, data, out$solution, hrm=FALSE, rate.cat=NULL, ntraits=2, model=model, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				write.table(cbind(row.names(data), lik.anc$lik.tip.states), file="Tipstates.DISCRETE.xls", row.names=F, quote=FALSE, sep="\t")
				phy$node.label <- lik.anc$lik.anc.states
				write.tree(phy, file="AncReconStatesDISCRETE.tre", append=TRUE)
			}
		}
	}
	if(ntraits==3){
		obj$Param.est<- matrix(out$solution[index.matrix], dim(index.matrix))
		#Calculates the standard error of the parameter by taking the sqrt of diagonal of the inverse of the Hessian: 
		obj$Param.SE <- matrix(sqrt(diag(pseudoinverse(h)))[index.matrix], dim(index.matrix))
		rownames(obj$Param.est) <- rownames(obj$Param.SE) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
		colnames(obj$Param.est) <- colnames(obj$Param.SE) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
		hess.eig <- eigen(h,symmetric=TRUE)
		obj$eigval <- signif(hess.eig$values,2)
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				lik.anc <- recon.marginal(phy, data, out$solution, hrm=FALSE, rate.cat=NULL, ntraits=3, model=model, par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				colnames(lik.anc$lik.anc.states) <-  c("P(0,0,0)","P(1,0,0)","P(0,1,0)","P(0,0,1)","P(1,1,0)","P(1,0,1)","P(0,1,1)","P(1,1,1)")
				write.table(lik.anc$lik.anc.states, file="Anc.EstimatesDISCRETE.xls", quote=FALSE, sep="\t")
				pr<-apply(lik.anc$lik.anc.states,1,which.max)
				phy$node.label <- pr
				write.tree(phy, file="AncReconStatesDISCRETE.tre", append=TRUE)
				phy$node.label <- 1:nb.node
				write.tree(phy, file="AncReconKey.tre")
			}
			if (node.states == "joint"){
				lik.anc <- recon.joint(phy, data, out$solution, hrm=FALSE, rate.cat=NULL, ntraits=3, model=model,par.drop=par.drop, par.eq=par.eq, root.p=root.p)
				write.table(cbind(row.names(data), lik.anc$lik.tip.states), file="Tipstates.DISCRETE.xls", row.names=F, quote=FALSE, sep="\t")
				phy$node.label <- lik.anc$lik.anc.states
				write.tree(phy, file="AncReconStatesDISCRETE.tre", append=TRUE)
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
	}
	else{obj$diagnostic<-'Arrived at a reliable solution'}
	
	obj
}

