#CORRELATED EVOLUTION OF TWO or THREE BINARY TRAITS

#written by Jeremy M. Beaulieu

corDISC<-function(phy, data, ntraits=2, rate.mat=NULL, model=c("ER","SYM","ARD"), node.states=c("joint", "marginal", "scaled"), p=NULL, root.p=NULL, ip=NULL, lb=0, ub=100){
	
	#Creates the data structure and orders the rows to match the tree
	phy$edge.length[phy$edge.length==0]=1e-5
	
	if(ntraits==2){
		data.sort<-data.frame(data[,2], data[,3], row.names=data[,1])
	}
	if(ntraits==3){
		data.sort<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
	}
	data.sort<-data.sort[phy$tip.label,]
	#Some initial values for use later - will clean up
	k=ntraits
	nl=2
	
	# Check to make sure values are reasonable (i.e. non-negative)
	if(ub < 0){
		ub <- 100
	}
	if(lb < 0){
		lb <- 0
	}
	if(ub < lb){ # This user really needs help
		ub <- 100
		lb <- 0
	}

	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	
	ntraits=ntraits
	model=model
	root.p=root.p	
	ip=ip
	
	model.set.final<-rate.mat.set(phy,data.sort,ntraits,model=model)
	if(!is.null(rate.mat)){
		rate <- rate.mat
		rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
		model.set.final$rate <- rate
		model.set.final$index.matrix <- rate.mat
	}
	
	lower = rep(lb, model.set.final$np)
	upper = rep(ub, model.set.final$np)
	
	opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)
	
	if(!is.null(p)){
		cat("Calculating likelihood from a set of fixed parameters", "\n")
		out<-NULL
		out$solution<-p
		out$objective<-dev.cordisc(out$solution,phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
		loglik <- -out$objective
		est.pars<-out$solution
	}
	else{
		if(is.null(ip)){
			cat("Initializing...", "\n")
			#If the analysis is to be run a single processor:
			#Sets parameter settings for random restarts by taking the parsimony score and dividing
			#by the total length of the tree
			model.set.init<-rate.mat.set(phy,data.sort,ntraits,model="ER")
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)
			dat<-as.matrix(data.sort)
			dat<-phyDat(dat,type="USER", levels=c("0","1"))
			par.score<-parsimony(phy, dat, method="fitch")
			tl <- sum(phy$edge.length)
			mean = par.score/tl
			ip<-rexp(1, 1/mean)
			lower.init = rep(lb, model.set.init$np)
			upper.init = rep(ub, model.set.init$np)
			init = nloptr(x0=rep(ip, length.out = model.set.init$np), eval_f=dev.cordisc, lb=lower.init, ub=upper.init, opts=opts, phy=phy,liks=model.set.init$liks,Q=model.set.init$Q,rate=model.set.init$rate,root.p=root.p)
			cat("Finished. Begin thorough search...", "\n")
			lower = rep(lb, model.set.final$np)
			upper = rep(ub, model.set.final$np)	
			out = nloptr(x0=rep(init$solution, length.out = model.set.final$np), eval_f=dev.cordisc, lb=lower, ub=upper, opts=opts, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
			loglik <- -out$objective
			est.pars<-out$solution
		}
		#If a user-specified starting value(s) is supplied:
		else{
			cat("Begin subplex optimization routine -- Starting value(s):", ip, "\n")
			opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.25)
			out = nloptr(x0=rep(ip, length.out = model.set.final$np), eval_f=dev.cordisc, lb=lower, ub=upper, opts=opts)
			loglik <- -out$objective
			est.pars<-out$solution
		}
	}
	#Starts the summarization process:
	cat("Finished. Inferring ancestral states using", node.states, "reconstruction.","\n")
	
	lik.anc <- ancRECON(phy, data, est.pars, hrm=FALSE, rate.cat=NULL, ntraits=ntraits, method=node.states, model=model, root.p=root.p)
	if(node.states == "marginal" || node.states == "scaled"){
		pr<-apply(lik.anc$lik.anc.states,1,which.max)
		phy$node.label <- pr
		tip.states <- NULL
	}
	if(node.states == "joint"){
		phy$node.label <- lik.anc$lik.anc.states
		tip.states <- lik.anc$lik.tip.states
	}
	
	cat("Finished. Performing diagnostic tests.", "\n")
	
	#Approximates the Hessian using the numDeriv function
	h <- hessian(func=dev.cordisc, x=est.pars, phy=phy,liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p)
	solution <- matrix(est.pars[model.set.final$index.matrix], dim(model.set.final$index.matrix))
	solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[model.set.final$index.matrix], dim(model.set.final$index.matrix))
	
	if(ntraits==2){
		rownames(solution) <- rownames(solution.se) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
		colnames(solution) <- colnames(solution.se) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				colnames(lik.anc$lik.anc.states) <-  c("P(0,0)","P(0,1)","P(1,0)","P(1,1)")
			}
		}
	}
	if(ntraits==3){
		rownames(solution) <- rownames(solution.se) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
		colnames(solution) <- colnames(solution.se) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
		if (is.character(node.states)) {
			if (node.states == "marginal"){
				colnames(lik.anc$lik.anc.states) <-  c("P(0,0,0)","P(1,0,0)","P(0,1,0)","P(0,0,1)","P(1,1,0)","P(1,0,1)","P(0,1,1)","P(1,1,1)")
			}
		}
	}
	hess.eig <- eigen(h,symmetric=TRUE)
	eigval<-signif(hess.eig$values,2)
	eigvect<-round(hess.eig$vectors, 2)
	obj = list(loglik = loglik, AIC = -2*loglik+2*model.set.final$np,AICc = -2*loglik+(2*model.set.final$np*(nb.tip/(nb.tip-model.set.final$np-1))),ntraits=ntraits, solution=solution, solution.se=solution.se, index.mat=model.set.final$index.matrix, opts=opts, data=data.sort, phy=phy, states=lik.anc$lik.anc.states, tip.states=tip.states, iterations=out$iterations, eigval=eigval, eigvect=eigvect) 
	class(obj)<-"cordisc"
	return(obj)
}


#Print function
print.cordisc<-function(x,...){
	
	ntips=Ntip(x$phy)
	output<-data.frame(x$loglik,x$AIC,x$AICc,x$ntraits,ntips, row.names="")
	names(output)<-c("-lnL","AIC","AICc","N.traits","ntax")
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


dev.cordisc<-function(p,phy,liks,Q,rate,root.p){
	
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

rate.mat.set<-function(phy,data.sort,ntraits,model){

	k=ntraits
	nl=2
	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode
	
	if(ntraits==2){
		rate<-rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
		index.matrix<-rate
		rate[is.na(rate)]<-max(rate,na.rm=T)+1

		x<-data.sort[,1]
		y<-data.sort[,2]
		
		liks <- matrix(0, nb.tip + nb.node, nl^k)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(is.na(x[i])){
				x[i]=2
				y[i]=2
			}
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
		rate.mat<-rate.mat.maker(hrm=FALSE,ntraits=ntraits,model=model)
		x<-data.sort[,1]
		y<-data.sort[,2]
		z<-data.sort[,3]
		
		liks <- matrix(0, nb.tip + nb.node, nl^k)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(is.na(x[i])){
				x[i]=2
				y[i]=2
				z[i]=2
			}
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
	Q <- matrix(0, nl^k, nl^k)
	
	obj$np<-max(rate)-1
	obj$rate<-rate
	obj$index.matrix<-index.matrix
	obj$liks<-liks
	obj$Q<-Q
	
	obj
}


