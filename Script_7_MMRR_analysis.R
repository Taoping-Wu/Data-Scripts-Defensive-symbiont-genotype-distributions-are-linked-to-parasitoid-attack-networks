############################################
### create aphid genetic distance matrix ###
############################################

library(ape)

setwd("C:/Users/wtp25/Desktop/En taro Feng Huan/PhD in QMUL/Research/Distribution of Hamiltonella defensa is shaped by parasitoids attack/Data and script/Data&script")

## read aphid COI fasta file
dna_sequences <- read.dna("Data_4_Aphid sequences.fas", format = "fasta")
## create pairwise genetic distance matrix
aphid_genetic <- dist.gene(dna_sequences, method = "pairwise", pairwise.deletion = FALSE,
                           variance = FALSE)
genetic_distances_matrix <- as.matrix(aphid_genetic)

#the biggest value in this genetic distance file is 1.06, to normalize this, we divided all value by 1.06 manually. 
write.csv(genetic_distances_matrix, file = "genetic distance.csv")







library(ape)
library(vegan)

############
### MMRR ###
############

# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

MMRR<-function(Y,X,nperm=999){
	#compute regression coefficients and test statistics
	nrowsY<-nrow(Y)
	y<-unfold(Y)
	if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
        Xmats<-sapply(X,unfold)
        fit<-lm(y~Xmats)
	coeffs<-fit$coefficients
	summ<-summary(fit)
	r.squared<-summ$r.squared
	tstat<-summ$coefficients[,"t value"]
	Fstat<-summ$fstatistic[1]
	tprob<-rep(1,length(tstat))
	Fprob<-1

	#perform permutations
	for(i in 1:nperm){
		rand<-sample(1:nrowsY)
		Yperm<-Y[rand,rand]
		yperm<-unfold(Yperm)
		fit<-lm(yperm~Xmats)
		summ<-summary(fit)
                Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
                tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
	}

	#return values
	tp<-tprob/(nperm+1)
	Fp<-Fprob/(nperm+1)
	names(r.squared)<-"r.squared"
	names(coeffs)<-c("Intercept",names(X))
	names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
	names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
	names(Fstat)<-"F-statistic"
	names(Fp)<-"F p-value"
	return(list(r.squared=r.squared,
		coefficients=coeffs,
		tstatistic=tstat,
		tpvalue=tp,
		Fstatistic=Fstat,
		Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
	x<-vector()
	for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
	return(x)
}

######################
### Sample loading ###
######################

## Set the workpath
setwd("C:/Users/wtp25/Desktop/En taro Feng Huan/PhD in QMUL/Research/Distribution of Hamiltonella defensa is shaped by parasitoids attack/Data and script/Data&script")

##################
### parasitoid ###
##################

### Parasitoid & Aphid distance as predictors, Hamiltonella as dependent matrix

# read parasitoid matrix, convert into Bray-Curtis dissimilarity matrix
Para <- read.csv("Para_Aphid_matrix_22species.csv", header=T, row.names = 1)
Para.dists <- vegdist(Para, method = "bray") # Bray–Curtis dissimilarity
Para.dists <- as.matrix(Para.dists)

# read aphid genetic distance matrix
Aphid <- read.csv("Aphid_phylogenetic_relatedness_22species.csv", header=T, row.names = 1)
Aphid.dists <- as.matrix(Aphid)

#read Hamiltonella matrix, convert into Bray-Curtis dissimilarity list (as the Y dependent distance matrix is a list)
Ham <- read.csv("Hamiltonella_Aphid_matrix_22species.csv", header=T, row.names = 1)
Ham.dists <- as.matrix(vegdist(Ham, method = "bray")) # Bray–Curtis dissimilarity

# Combine predictor matrices into a list
predictors <- list(X1 = Para.dists, X2 = Aphid.dists)
# Run MMRR
result <- MMRR(Ham.dists, predictors, nperm = 9999)

# View the results
print(result$coefficients)  # Regression coefficients
print(result$tpvalue)      # P-values for each predictor


#############
### Plant ###
#############

### Plant & Aphid distance as predictors, Hamiltonella as dependent matrix

Plant <- read.csv("Plant_Aphid_matrix_31species.csv", header=T, row.names = 1)
Plant.dists <- as.matrix(vegdist(Plant, method = "bray")) # Bray–Curtis dissimilarity

Aphid <- read.csv("Aphid_phylogenetic_relatedness_31species.csv", header=T, row.names = 1)
Aphid.dists <- as.matrix(Aphid)

Ham <- read.csv("Hamiltonella_Aphid_matrix_31species.csv", header=T, row.names = 1)
Ham.dists <- as.matrix(vegdist(Ham, method = "bray")) # Bray–Curtis dissimilarity


# Combine predictor matrices into a list
predictors <- list(X1 = Plant.dists, X2 = Aphid.dists)
# Run MMRR
result <- MMRR(Ham.dists, predictors, nperm = 9999)

# View the results
print(result$coefficients)  # Regression coefficients
print(result$tpvalue)      # P-values for each predictor


###############
### both 22 ###
###############

Para <- read.csv("Para_Aphid_matrix_22species.csv", header=T, row.names = 1)
Para.dists <- vegdist(Para, method = "bray") # Bray–Curtis dissimilarity
Para.dists <- as.matrix(Para.dists)

Plant <- read.csv("Plant_Aphid_matrix_22species.csv", header=T, row.names = 1)
Plant.dists <- vegdist(Plant, method = "bray") # Bray–Curtis dissimilarity
Plant.dists <- as.matrix(Plant.dists)

Aphid <- read.csv("Aphid_phylogenetic_relatedness_22species.csv", header=T, row.names = 1)
Aphid.dists <- as.matrix(Aphid)

Ham <- read.csv("Hamiltonella_Aphid_matrix_22species.csv", header=T, row.names = 1)
Ham.dists <- as.matrix(vegdist(Ham, method = "bray")) # Bray–Curtis dissimilarity

# Combine predictor matrices into a list
predictors <- list(X1 = Para.dists, X2 = Plant.dists, X3 = Aphid.dists)
# Run MMRR
result <- MMRR(Ham.dists, predictors, nperm = 9999)

# View the results
print(result$coefficients)  # Regression coefficients
print(result$tpvalue)      # P-values for each predictor


