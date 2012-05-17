#MARGINAL RECONSTRUCTION OF ANCESTRAL STATE

#written by Jeremy M. Beaulieu

#Algorithm is based on ace(), though trees do not need to be bifurcating. Also, the code is written
#so that it can be used as a separate function from corHMM. All that is required is a tree, trait, and a vector
#of estimated parameter values and the user is provided the marginal ancestral reconstruction.

recon.marginal <- function(phy, data, p, hrm=TRUE, rate.cat, ntraits=NULL, model=c("ER", "SYM", "ARD"), par.drop=NULL, par.eq=NULL, root.p=NULL){
	
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
	if(hrm==TRUE){
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
	}
	if(hrm==FALSE){
		if(ntraits==2){
			k=2
			nl=2
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
			Q <- matrix(0, nl^k, nl^k)
		}
		if(ntraits==3){
			k=3
			nl=2
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
			} else {
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
			Q <- matrix(0, nl^k, nl^k)
		}
	}

	phy <- reorder(phy, "pruningwise")
	#Number of columns should be equal to the number of states.
	comp <- numeric(nb.tip + nb.node)
	lik.states<-numeric(nb.tip + nb.node)
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	
	#The same algorithm as in the main function. See comments in corHMM.R for details:
	for (i  in seq(from = 1, length.out = nb.node)) {
		#the ancestral node at row i is called focal
		focal <- anc[i]
		#Get descendant information of focal
		desRows<-which(phy$edge[,1]==focal)
		desNodes<-phy$edge[desRows,2]
		v <- 1
		for (desIndex in sequence(length(desRows))){
			v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
		}
		comp[focal] <- sum(v)
		liks[focal, ] <- v/comp[focal]
	}
	obj$lik.anc.states <- liks[-TIPS, ]
	
	obj
}


