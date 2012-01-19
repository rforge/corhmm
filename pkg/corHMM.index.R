#corHMM INDEX MATRIX

#Outputs the full index of the rate parameters that are optimized. The intention is that
#a user might want to see how the matrix is designed prior to an analysis and perhaps drops
#may be a few parameters a priori due to some hypothesis that he or she might have.

corHMM.index<-function(rate.cat){

	k=2
	nl=2

	if (rate.cat == 1){
		rate <- matrix(NA, nl, nl)
		np <- 2
		
		index<-matrix(TRUE,nl,nl)
		diag(index) <- FALSE
		rate[index] <- 1:np
		index.matrix <- rate
		rownames(index.matrix) <- c("(0)","(1)")
		colnames(index.matrix) <- c("(0)","(1)")
	}
	if (rate.cat == 2){
		rate <- matrix(NA, nl^k, nl^k)
		np <- 8
		tmp <- cbind(1:(nl^k), (nl^k):1)
		tmp2 <- cbind(1:(nl^k), 1:(nl^k))
		
		index <- matrix(TRUE,nl^k,nl^k)
		diag(index) <- FALSE
		index[tmp] <- FALSE
		index[tmp2] <- FALSE			
		rate[index] <- 1:np
		index.matrix <- rate
		rownames(index.matrix) <- c("(0,R1)","(0,R2)","(1,R1)","(1,R2)")
		colnames(index.matrix) <- c("(0,R1)","(0,R2)","(1,R1)","(1,R2)")
	}
	if (rate.cat == 3){
		rate <- matrix(NA, nl*k+2, nl*k+2)
		np <- 14
		tmp <-	c(3,5,6,4,6,1,4,5,2,3,6,1,3,1,2,4)
		tmp2 <- c(1,1,1,2,2,3,3,3,4,4,4,5,5,6,6,6)
		tmp3 <- cbind(tmp,tmp2)
		
		index <- matrix(TRUE,(nl*k+2),(nl*k+2))
		diag(index) <- FALSE
		index[tmp3] <- FALSE			
		rate[index] <- 1:np
		index.matrix <- rate
		rownames(index.matrix) <- c("(0,R1)","(0,R2)","(0,R3)","(1,R1)","(1,R2)","(1,R3)")
		colnames(index.matrix) <- c("(0,R1)","(0,R2)","(0,R3)","(1,R1)","(1,R2)","(1,R3)")
	}
	if (rate.cat == 4){
		rate <- matrix(NA, nl*k+4, nl*k+4)
		np <- 20
		tmp <- c(3,4,6,7,8,4,5,7,8,1,5,6,8,1,2,5,6,7,2,3,4,7,8,1,3,4,8,1,2,4,5,1,2,3,5,6)
		tmp2 <-c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8)
		tmp3 <- cbind(tmp,tmp2)
		
		index <- matrix(TRUE,(nl*k+4),(nl*k+4))
		diag(index) <- FALSE
		index[tmp3] <- FALSE			
		rate[index] <- 1:np
		index.matrix <- rate
		rownames(index.matrix) <- c("(0,R1)","(0,R2)","(0,R3)","(0,R4)","(1,R1)","(1,R2)","(1,R3)","(1,R4)")
		colnames(index.matrix) <- c("(0,R1)","(0,R2)","(0,R3)","(0,R4)","(1,R1)","(1,R2)","(1,R3)","(1,R4)")
	}
	if (rate.cat == 5){
		rate <- matrix(NA, nl*k+6, nl*k+6)
		np <- 26
		tmp <- c(3,4,5,7,8,9,10,4,5,6,8,9,10,1,5,6,7,9,10,1,2,6,7,8,10,1,2,3,6,7,8,9,2,3,4,5,8,9,10,1,3,4,5,9,10,1,2,4,5,6,10,1,2,3,5,6,7,1,2,3,4,6,7,8)
		tmp2 <-c(1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,10,10)
		tmp3 <- cbind(tmp,tmp2)
		
		index <- matrix(TRUE,(nl*k+6),(nl*k+6))
		diag(index) <- FALSE
		index[tmp3] <- FALSE			
		rate[index] <- 1:np
		index.matrix <- rate
		rownames(index.matrix) <- c("(0,R1)","(0,R2)","(0,R3)","(0,R4)","(0,R5)","(1,R1)","(1,R2)","(1,R3)","(1,R4)","(1,R5)")
		colnames(index.matrix) <- c("(0,R1)","(0,R2)","(0,R3)","(0,R4)","(0,R5)","(1,R1)","(1,R2)","(1,R3)","(1,R4)","(1,R5)")
	}
	index.matrix
}

