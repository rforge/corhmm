#JOINT RECONSTRUCTION OF ANCESTRAL STATES

#written by Jeremy M. Beaulieu

#Algorithm is based on Pupko et al (2000). Trees do not need to be bifurcating. Also, the code is written
#so that it can be used as a separate function from corHMM. All that is required is a tree, trait, and a vector
#of estimated parameter values and the user is provided the joint ancestral reconstruction.

recon.joint <- function(phy, data, p, rate.cat, par.drop=NULL, par.eq=NULL, root.p=NULL){
	
	#Note: Does not like zero branches at the tips. Here I extend these branches by just a bit:
	phy$edge.length[phy$edge.length<=1e-5]=1e-5
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
	
	if (rate.cat == 1){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		TIPS <- 1:nb.tip
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,1]=1}
			if(x[i]==1){liks[i,2]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 2){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3)]=1}
			if(x[i]==1){liks[i,c(2,4)]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 3){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5)]=1}
			if(x[i]==1){liks[i,c(2,4,6)]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 4){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8)]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 5){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,c(1,3,5,7,9)]=1}
			if(x[i]==1){liks[i,c(2,4,6,8,10)]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	phy <- reorder(phy, "pruningwise")
	#Number of columns should be equal to the number of states.
	comp<-matrix(0,nb.tip + nb.node,ncol(liks))
	lik.states<-numeric(nb.tip + nb.node)
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)

	for (i  in seq(from = 1, length.out = nb.node)) {
		#The ancestral node at row i is called focal:
		focal <- anc[i]
		#Get descendant information of focal:
		desRows<-which(phy$edge[,1]==focal)
		#Get node information for each descendant:
		desNodes<-phy$edge[desRows,2]
		#Initiates a loop to check if any nodes are tips:
		for (desIndex in sequence(length(desRows))){
			#If a tip calculate C_y(i) for the tips and stores in liks matrix:
			if(any(desNodes[desIndex]==phy$edge[,1])==FALSE){
				liks[desNodes[desIndex],] <- expm(Q * phy$edge.length[i], method=c("Ward77")) %*% liks[desNodes[desIndex],]
				#Divide by the sum of the liks to deal with underflow issues:
				liks[desNodes[desIndex],] <- liks[desNodes[desIndex],]/sum(liks[desNodes[desIndex],])
				#Collects the likeliest state at the tips:
				comp[desNodes[desIndex],] = which.max(liks[desNodes[desIndex],])
			}
		}
		#Collects t_z, or the branch subtending focal:
		tz<-phy$edge.length[which(phy$edge[,2] == focal)]	
		if(length(tz)==0){
			#The focal node is the root, calculate P_k:
			if(is.null(root.p)){
				root.p=1
				for (desIndex in sequence(length(desRows))){
					#This is the basic marginal calculation:
					root.p <- root.p * expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
				}
				#Divide by the sum of the liks to deal with underflow issues:
				liks[focal, ] <- root.p/sum(root.p)
			}
			else{
				liks[focal, ] <- root.p
			}
		}
		else{
			#Calculates P_ij(t_z):
			Pij <- expm(Q * tz, method=c("Ward77"))
			#Calculates L_z(i):
			v<-c(rep(1, k*rate.cat))
			for (desIndex in sequence(length(desRows))){
				v = v * liks[desNodes[desIndex],]
			}
			#Finishes L_z(i):
			L <- t(Pij) * v
			#Collects which is the highest likelihood and which state it corresponds to:
			liks[focal,] <- apply(L, 2, max)
			comp[focal,] <- apply(L, 2, which.max)
			#Divide by the sum of the liks to deal with underflow issues:
			liks[focal,] <- liks[focal,]/sum(liks[focal,])
		}
	}
	#If the state at the root is not specified root will just be the joint estimate:
	if (is.null(root.p)){
		root <- nb.tip + 1L
		lik.states[root] <- which.max(liks[root,])
		N <- dim(phy$edge)[1]
		for(i in N:1){
			des <- phy$edge[i,2]
			tmp <- which.max(liks[des,])
			lik.states[des] <- comp[des,tmp]
		}
	}
	#If the state at the root is specified then will start up-pass based on the maximum of the user-defined root probabilities:
	else{
		root <- nb.tip + 1L
		lik.states[root] <- which.max(root.p)
		N <- dim(phy$edge)[1]
		for(i in N:1){
			des <- phy$edge[i,2]
			tmp <- which.max(liks[des,])
			lik.states[des] <- comp[des,tmp]
		}
	}
	#Outputs likeliest tip states
	obj$lik.tip.states <- lik.states[TIPS]
	#Outputs likeliest node states
	obj$lik.anc.states <- lik.states[-TIPS]
	
	obj
}


