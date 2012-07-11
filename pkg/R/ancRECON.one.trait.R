#MARGINAL RECONSTRUCTION OF ANCESTRAL STATES FOR GIVEN MATRIX

#written by Jeremy M. Beaulieu & Jeffrey C. Oliver

#Algorithm was based on a description by Yang (2006), and trees do not need to be bifurcating. The basic idea of Yang (2006) 
#is that the tree is rerooted on each internal node, with the marginal likelihood being the probabilities of observing the tips 
#states given that the the focal node is the root. However, this takes a ton of time as the number of nodes increases. But more 
#importantly, this does not work when the model contains asymmetric rates. Here we use a dynamic programming algorithm that 
#calculates the marginal probability at a node during an additional up and down passes of the tree. This approach was 
#developed with significant influence of the Mesquite modules zMargLikeCateg and AsymmModel.  The results are identical to those 
#from Mesquite and diversitree. Also, the code is written so that it can be used as a separate function from corHMM. All that is 
#required is a tree, trait, and a vector of estimated parameter values and the user is provided the marginal ancestral reconstruction.


ancRECON.one.trait <- function(phy, data, p, method=c("joint", "marginal", "scaled"), model=c("ER", "SYM", "ARD"), charnum=1,par.drop=NULL, par.eq=NULL, root.p=NULL){

	#Note: Does not like zero branches at the tips. Here I extend these branches by just a bit:
	phy$edge.length[phy$edge.length<=1e-5]=1e-5

	data <- data.frame(data[,charnum+1],data[,charnum+1],row.names=data[,1]) # added character twice, because at least two columns are necessary
	data <- data[phy$tip.label,] # this might have already been done by match.tree.data	

	#Some initial values for use later

	obj <- NULL
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	
	par.drop=par.drop
	par.eq=par.eq
	root.p=root.p
	
	k <- 1
	factored <- factorData(data)

	nl <- ncol(factored)
	obj <- NULL
	nb.tip<-length(phy$tip.label)
	nb.node <- phy$Nnode

	#Using one of three pre-defined models instead of user-defined model
	if (is.character(model)) {
		#rate is a matrix of rate categories (not actual rates)
		rate <- matrix(NA, nl, nl) #An nl x nl matrix, filled with NA values
		#np is the number of parameters in the rate matrix
		#Equal rates model, one parameter, all rate categories enumerated as '1'
		if (model == "ER") np <- rate[] <- 1
		#All rates different model, number of different parameters just nl x (nl-1)
		#rates are enumerated up to np (e.g. for a two-state character there are 2 rate categories)
		if (model == "ARD") {
			np <- nl*(nl - 1)
			rate[col(rate) != row(rate)] <- 1:np
		}
		#Symmetrical model, like ARD, with half the values
		if (model == "SYM") {
			np <- nl * (nl - 1)/2
			sel <- col(rate) < row(rate)
			rate[sel] <- 1:np
			#Use transpose of the rate category matrix to finish enumerating:
			rate <- t(rate)
			rate[sel] <- 1:np
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
	tmp <- cbind(1:nl, 1:nl)
	rate[tmp] <- 0
	rate[rate == 0] <- np + 1 # Got to make sure to assign np + 1 to empty rate categories

	stateTable <- NULL # will hold 0s and 1s for likelihoods of each state at tip
	for(column in 1:nl){
		stateTable <- cbind(stateTable,factored[,column])
	}
	colnames(stateTable) <- colnames(factored)

	ancestral <- matrix(0,nb.node,nl) # all likelihoods at ancestral nodes will be 0
	liks <- rbind(stateTable,ancestral) # combine tip likelihoods & ancestral likelihoods
	rownames(liks) <- NULL

	Q <- matrix(0, nl^k, nl^k)
	Q[] <- c(p, 0)[rate]
	diag(Q) <- -rowSums(Q)
	phy <- reorder(phy, "pruningwise")
	TIPS <- 1:nb.tip
	anc <- unique(phy$edge[,1])

	if(method=="joint"){
		lik.states<-numeric(nb.tip + nb.node)
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
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
					liks[desNodes[desIndex],] <- expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
					#Divide by the sum of the liks to deal with underflow issues:
					liks[desNodes[desIndex],] <- liks[desNodes[desIndex],]/sum(liks[desNodes[desIndex],])
					#Collects the likeliest state at the tips:
					comp[desNodes[desIndex],] <- which.max(liks[desNodes[desIndex],])
				}
			}
			#Collects t_z, or the branch subtending focal:
			tz<-phy$edge.length[which(phy$edge[,2] == focal)]	
			if(length(tz)==0){
				#The focal node is the root, calculate P_k:
				if(is.null(root.p)){
					root.state=1
					for (desIndex in sequence(length(desRows))){
						#This is the basic marginal calculation:
						root.state <- root.state * expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks[desNodes[desIndex],]
					}
					#Divide by the sum of the liks to deal with underflow issues:
					liks[focal, ] <- root.state/sum(root.state)
				}
				else{
					liks[focal, ] <- root.p
				}
			}
			else{
				#Calculates P_ij(t_z):
				Pij <- expm(Q * tz, method=c("Ward77"))
				v<-c(rep(1, nl^k))
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
		root <- nb.tip + 1L
		lik.states[root] <- which.max(liks[root,])
		N <- dim(phy$edge)[1]
		for(i in N:1){
			des <- phy$edge[i,2]
			tmp <- which.max(liks[des,])
			lik.states[des] <- comp[des,tmp]
		}
		#Outputs likeliest tip states
		obj$lik.tip.states <- lik.states[TIPS]
		#Outputs likeliest node states
		obj$lik.anc.states <- lik.states[-TIPS]
	}

	if(method=="marginal"){
		#A temporary likelihood matrix so that the original does not get written over:
		liks.down<-liks
		# root equilibrium frequencies
		equil.root <- NULL
		for(i in 1:ncol(Q)){
			posrows <- which(Q[,i] >= 0)
			rowsum <- sum(Q[posrows,i])
			poscols <- which(Q[i,] >= 0)
			colsum <- sum(Q[i,poscols])
			equil.root <- c(equil.root,rowsum/(rowsum+colsum))
		}
		#A transpose of Q for assessing probability of j to i, rather than i to j:
		tranQ <- matrix(0, nl^k, nl^k)
		tranQ<-t(Q)

		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		#The first down-pass: The same algorithm as in the main function to calculate the conditional likelihood at each node:
		for (i  in seq(from = 1, length.out = nb.node)) {
			#the ancestral node at row i is called focal
			focal <- anc[i]
			#Get descendant information of focal
			desRows<-which(phy$edge[,1]==focal)
			desNodes<-phy$edge[desRows,2]
			v <- 1
			for (desIndex in sequence(length(desRows))){
				v <- v*expm(Q * phy$edge.length[desRows[desIndex]], method=c("Ward77")) %*% liks.down[desNodes[desIndex],]
			}
			comp[focal] <- sum(v)
			liks.down[focal, ] <- v/comp[focal]
		}

		root <- nb.tip + 1L
		#Enter the root defined root probabilities if they are supplied by the user:
		if(!is.null(root.p)){
			root <- nb.tip + 1L	
			liks.down[root, ]<-root.p
		}
		#The up-pass 
		liks.up<-liks
		states<-apply(liks,1,which.max)
		N <- dim(phy$edge)[1]
		comp <- numeric(nb.tip + nb.node)
		for(i in length(anc):1){
			focal <- anc[i]
			if(!focal==root){
				#Gets mother and sister information of focal:
				focalRow<-which(phy$edge[,2]==focal)
				motherRow<-which(phy$edge[,1]==phy$edge[focalRow,1])
				motherNode<-phy$edge[focalRow,1]
				desNodes<-phy$edge[motherRow,2]
				sisterNodes<-desNodes[(which(!desNodes==focal))]
				sisterRows<-which(phy$edge[,2]%in%sisterNodes==TRUE)
				#If the mother is not the root then you are calculating the probability of the being in either state.
				#But note we are assessing the reverse transition, j to i, rather than i to j, so we transpose Q to carry out this calculation:
				if(motherNode!=root){
					v <- expm(tranQ * phy$edge.length[which(phy$edge[,2]==motherNode)], method=c("Ward77")) %*% liks.up[motherNode,]
				}
				#If the mother is the root then just use the marginal. This can also be the prior, which I think is the equilibrium frequency. 
				#But for now we are just going to use the marginal at the root -- it is unclear what Mesquite does.
				else{
					v <- equil.root
				}
				#Now calculate the probability that each sister is in either state. Sister can be more than 1 when the node is a polytomy. 
				#This is essentially calculating the product of the mothers probability and the sisters probability:
				for (sisterIndex in sequence(length(sisterRows))){
					v <- v*expm(Q * phy$edge.length[sisterRows[sisterIndex]], method=c("Ward77")) %*% liks.down[sisterNodes[sisterIndex],]
				}
				comp[focal] <- sum(v)
				liks.up[focal,] <- v/comp[focal]
			}
		}
		#The final pass
		liks.final<-liks
		comp <- numeric(nb.tip + nb.node)
		for (i  in seq(from = 1, length.out = nb.node-1)) { # In this final pass, root is never encountered.  But its OK, because root likelihoods are set after loop.
			#the ancestral node at row i is called focal
			focal <- anc[i]
			focalRows<-which(phy$edge[,2]==focal)
			#Now you are asessing the change along the branch subtending the focal by multiplying the probability of 
			#everything at and above focal by the probability of the mother and all the sisters given time t:
			v <- liks.down[focal,]*expm(tranQ * phy$edge.length[focalRows], method=c("Ward77")) %*% liks.up[focal,]

			comp[focal] <- sum(v)
			liks.final[focal, ] <- v/comp[focal]
		}
		#Just add in the marginal at the root calculated on the original downpass or if supplied by the user:
		liks.final[root,] <- liks.down[root,] * equil.root

		root.final <- liks.down[root,] * equil.root
		comproot <- sum(root.final)
		liks.final[root,] <- root.final/comproot

		#Reports just the probabilities at internal nodes:
		obj$lik.anc.states <- liks.final[-TIPS, ]
	}
	if(method=="scaled"){
		comp<-matrix(0,nb.tip + nb.node,ncol(liks))
		#The same algorithm as in the main function. See comments in either corHMM.R or corDISC.R for details:
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
		if(!is.null(root.p)){
			root <- nb.tip + 1L	
			liks[root, ]<-root.p
		}
		obj$lik.anc.states <- liks[-TIPS, ]
	}	
	obj
}
