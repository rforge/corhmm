#RATE MATRIX INDEX

rate.mat.index<-function(rate.cat, hrm=TRUE, ntraits=NULL, nstates=NULL, model=c("ER", "SYM", "ARD"), par.drop=NULL, par.eq=NULL){
	
	k=2

	par.drop=par.drop
	par.eq=par.eq
	
	#Builds the rate matrix based on the specified rate.cat. Not exactly the best way
	#to go about this, but it is the best I can do for now -- it works, so what me worry?
	if(hrm==TRUE){
		if (rate.cat == 1){
			rate <- matrix(NA, k*rate.cat, k*rate.cat)
			np <- 2
			
			index<-matrix(TRUE,k*rate.cat,k*rate.cat)
			diag(index) <- FALSE
			rate[index] <- 1:np
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
			index.matrix <- rate # TODO: Should this be index.matrix[index.matrix == 0] = NA

			rownames(index.matrix) <- c("(0,R1)","(1,R1)")
			colnames(index.matrix) <- c("(0,R1)","(1,R1)")		
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
			rownames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")
			colnames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)")		
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
			rownames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")
			colnames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)")		
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
			rownames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")
			colnames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)")		
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
			rownames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")
			colnames(index.matrix) <- c("(0,R1)","(1,R1)","(0,R2)","(1,R2)","(0,R3)","(1,R3)","(0,R4)","(1,R4)","(0,R5)","(1,R5)")		
		}
	}
	if(hrm==FALSE){
		#Imported from Jeffs rayDISC -- will clean up later, but for now, it works fine:
		if(ntraits==1){
			k <- 1
			nl <- nstates
			obj <- NULL
			
			#Using one of three pre-defined models instead of user-defined model
			if (is.character(model)) {
				#rate is a matrix of rate categories (not actual rates)
				rate <- matrix(NA, nl, nl) #An nl x nl matrix, filled with NA values
				tmp2 <- cbind(1:(nl^k), 1:(nl^k)) # For setting diagonals
				index<-matrix(TRUE,nl^k,nl^k)
				diag(index) <- FALSE
				#Equal rates model, one parameter, all rate categories enumerated as '1'
				if (model == "ER") {
					np <- 1 #np is the number of parameters in the rate matrix
					rate[index] <- 1:np
					# TODO: par.drop doesn't work
					#If par.drop is not null will adjust the rate matrix
#					if(!is.null(par.drop)==TRUE){
#						for(i in 1:length(par.drop)){
#							tmp3 <- which(rate==par.drop[i], arr.ind=T)
#							index[tmp3] <- FALSE
#							rate[tmp3] <- 0
#						}
#						np <- np-length(par.drop)
#						rate[index] <- 1:np
#					}
				}
				#All rates different model, number of different parameters just nl x (nl-1)
				#rates are enumerated up to np (e.g. for a two-state character there are 2 rate categories)
				if (model == "ARD") {
					np <- nl*(nl - 1)
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
				}
				#Symmetrical model, like ARD, with half the values
				if (model == "SYM") {
					np <- nl * (nl - 1)/2
					sel <- col(rate) < row(rate)
					rate[sel] <- 1:np
					#Use transpose of the rate category matrix to finish enumerating:
					rate <- t(rate)
					rate[sel] <- 1:np
					#If par.drop is not null will adjust the rate matrix
					if(!is.null(par.drop)==TRUE){
						for(i in 1:length(par.drop)){
							tmp3 <- which(rate==par.drop[i], arr.ind=T)
							index[tmp3] <- FALSE
							rate[tmp3] <- 0
							decrement <- which(rate > par.drop[i],arr.ind=TRUE) # rate categories above the one to be dropped
							rate[decrement] <- rate[decrement] - 1
						}
						np <- np-length(par.drop)
#						rate[index] <- 1:np # TODO: this renumbering doesn't work with symmetric rate matrix
					}
					#If par.eq is not null then pairs of parameters are set equal to each other.
					if(!is.null(par.eq)==TRUE){
						for (i  in seq(from = 1, by = 2, length.out = length(par.eq)/2)) {
							j<-i+1
							tmp3 <- which(rate==par.eq[j], arr.ind=T)
							if(length(tmp3) > 0){
								index[tmp3] <- FALSE
								rate[tmp3] <- 0
								np <- np-1
#								rate[index] <- 1:np  # TODO: this renumbering doesn't work with symmetric rate matrix
								rate[tmp3] <- par.eq[i]
								decrement <- which(rate > par.drop[i],arr.ind=TRUE) # rate categories above the one to be dropped
								rate[decrement] <- rate[decrement] - 1
							}
						}
					}
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
			index.matrix[index.matrix == 0] = NA

			rownames(index.matrix) <- c(0:(nl-1))
			colnames(index.matrix) <- c(0:(nl-1))
		}
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
				}
			} 
			rownames(index.matrix) <- c("(0,0)","(0,1)","(1,0)","(1,1)")
			colnames(index.matrix) <- c("(0,0)","(0,1)","(1,0)","(1,1)")		
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
					
				}
				rownames(index.matrix) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")
				colnames(index.matrix) <- c("(0,0,0)","(1,0,0)","(0,1,0)","(0,0,1)","(1,1,0)","(1,0,1)","(0,1,1)","(1,1,1)")	
			} 
		}
		
	}
	index.matrix
}
	
	
