#MARGINAL RECONSTRUCTION OF ANCESTRAL STATE

#written by Jeremy M. Beaulieu

#Algorithm is based on ace(), though trees do not need to be bifurcating. Also, the code is written
#so that it can be used as a separate function from corHMM. All that is required is a tree, trait, and a vector
#of estimated parameter values and the user is provided the marginal ancestral reconstruction.

recon.marginal <- function(phy, data, p, rate.cat, par.drop=NULL, par.eq=NULL, root.p=NULL){
	
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
		tmp <-	c(3,5,6,4,6,1,4,5,2,3,6,1,3,1,2,4)
		tmp2 <- c(1,1,1,2,2,3,3,3,4,4,4,5,5,6,6,6)
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
		tmp <- c(3,4,6,7,8,4,5,7,8,1,5,6,8,1,2,5,6,7,2,3,4,7,8,1,3,4,8,1,2,4,5,1,2,3,5,6)
		tmp2 <-c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8)
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
		tmp <- c(3,4,5,7,8,9,10,4,5,6,8,9,10,1,5,6,7,9,10,1,2,6,7,8,10,1,2,3,6,7,8,9,2,3,4,5,8,9,10,1,3,4,5,9,10,1,2,4,5,6,10,1,2,3,5,6,7,1,2,3,4,6,7,8)
		tmp2 <-c(1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,10,10)		
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
			if(x[i]==0){liks[i,1:2]=1}
			if(x[i]==1){liks[i,3:4]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 3){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,1:3]=1}
			if(x[i]==1){liks[i,4:6]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 4){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,1:4]=1}
			if(x[i]==1){liks[i,5:8]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	if (rate.cat == 5){
		liks <- matrix(0, nb.tip + nb.node, k*rate.cat)
		for(i in 1:nb.tip){
			if(x[i]==0){liks[i,1:5]=1}
			if(x[i]==1){liks[i,6:10]=1}
		}
		Q <- matrix(0, k*rate.cat, k*rate.cat)
	}
	phy <- reorder(phy, "pruningwise")
	#Number of columns should be equal to the number of states.
	comp <- numeric(nb.tip + nb.node)
	lik.states<-numeric(nb.tip + nb.node)
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	#These probabilities are just for example and taken from Table 1 in Pupko et al 2000.
	#Remember to add in exponentiation of Q to get probabilities for real datasets:

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
	obj$liks.anc.states<-liks[-TIPS, ]
	
	obj
}


