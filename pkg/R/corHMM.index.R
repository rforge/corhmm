#corHMM INDEX MATRIX

#Outputs the full index of the rate parameters that are optimized. The intention is that
#a user might want to see how the matrix is designed prior to an analysis and perhaps drops
#a few parameters beforehand due to some hypothesis that he or she might have.

corHMM.index<-function(rate.cat, par.drop=NULL, par.eq=NULL){
	k=2
	if (rate.cat == 1){
		rate <- matrix(NA, k*rate.cat, k*rate.cat)
		np <- 2
		
		index<-matrix(TRUE,k*rate.cat,k*rate.cat)
		diag(index) <- FALSE
		rate[index] <- 1:np
		index.matrix <- rate
		rownames(index.matrix) <- c("(0)","(1)")
		colnames(index.matrix) <- c("(0)","(1)")
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
	index.matrix
}

