############
### MMRR ###
############

rm(list = ls())

# load package

library(ape)
library(vegan)

#the MMRR analysis code is generated from Wang 2013. https://doi.org/10.1111/evo.12134

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

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##################
### parasitoid ###
##################

### Parasitoid & Aphid distance as predictors, Hamiltonella as dependent matrix

### 22 species, full model
# read parasitoid matrix, convert into Bray-Curtis dissimilarity matrix
Para <- read.csv("../Data_in_process/Para_Aphid_matrix_22species.csv", header=T, row.names = 1)
Para.dists_22 <- vegdist(Para, method = "bray") # Bray–Curtis dissimilarity
Para.dists_22 <- as.matrix(Para.dists_22)
Para.dists_22 <- 1-Para.dists_22         #as we are using 1-Bray–Curtis, similarity
# read aphid genetic distance matrix
Aphid <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_22species.csv", header=T, row.names = 1)
Aphid.dists_22 <- as.matrix(Aphid)
Aphid.dists_22 <- 1-Aphid.dists_22      #as we are using 1-Bray–Curtis, similarity
#read Hamiltonella matrix, convert into Bray-Curtis dissimilarity list (as the Y dependent distance matrix is a list)
Ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv", header=T, row.names = 1)
Ham.dists_22 <- as.matrix(vegdist(Ham, method = "bray")) # Bray–Curtis dissimilarity
Ham.dists_22 <- 1-Ham.dists_22                              #as we are using 1-Bray–Curtis, similarity
#read Plant matrix, convert into Bray-Curtis dissimilarity list (as the Y dependent distance matrix is a list)
Plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_22species.csv", header=T, row.names = 1)
Plant.dists_22 <- as.matrix(vegdist(Plant, method = "bray")) # Bray–Curtis dissimilarity
Plant.dists_22 <- 1-Plant.dists_22                              #as we are using 1-Bray–Curtis, similarity

# Combine predictor matrices into a list
predictors <- list(X1 = Para.dists_22, X2 = Aphid.dists_22)
# Run MMRR， Fig 3D
result_22_par_Aph <- MMRR(Ham.dists_22, predictors, nperm = 9999)

# single predictor models
# parasitoid only, same as Mantel's test, Fig 3A
predictors <- list(X1 = Para.dists_22)
result_22_par <- MMRR(Ham.dists_22, predictors, nperm = 9999)

# Aphid relatedness only, same as Mantel's test, , Fig 3B
predictors <- list(X1 = Aphid.dists_22)
result_22_aph <- MMRR(Ham.dists_22, predictors, nperm = 9999)
# predictors interaction model, Fig 3C (predictor is still Aphid.dist_22)
result_22_predictors_interaction <- MMRR(Para.dists_22, predictors, nperm = 9999)

# Combine predictor matrices into a list
predictors <- list(X1 =Para.dists_22, X2 =Plant.dists_22, X3 = Aphid.dists_22)
# Run MMRR
result_22_all <- MMRR(Ham.dists_22, predictors, nperm = 9999)


### 16 species, reduced model

# read parasitoid matrix, convert into Bray-Curtis dissimilarity matrix
Para <- read.csv("../Data_in_process/Para_Aphid_matrix_16species.csv", header=T, row.names = 1)
Para.dists_16 <- vegdist(Para, method = "bray") # Bray–Curtis dissimilarity
Para.dists_16 <- as.matrix(Para.dists_16)
Para.dists_16 <- 1-Para.dists_16         #as we are using 1-Bray–Curtis, similarity
# read aphid genetic distance matrix
Aphid <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv", header=T, row.names = 1)
Aphid.dists_16 <- as.matrix(Aphid)
Aphid.dists_16 <- 1-Aphid.dists_16      #as we are using 1-Bray–Curtis, similarity
#read Hamiltonella matrix, convert into Bray-Curtis dissimilarity list (as the Y dependent distance matrix is a list)
Ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1)
Ham.dists_16 <- as.matrix(vegdist(Ham, method = "bray")) # Bray–Curtis dissimilarity
Ham.dists_16 <- 1-Ham.dists_16                              #as we are using 1-Bray–Curtis, similarity

# Combine predictor matrices into a list
predictors <- list(X1 = Para.dists_16, X2 = Aphid.dists_16)
# Run MMRR, Fig S5C
result_16_par_Aph <- MMRR(Ham.dists_16, predictors, nperm = 9999)

# single predictor models
# parasitoid only, same as Mantel's test, Fig S5A
predictors <- list(X1 = Para.dists_16)
result_16_par <- MMRR(Ham.dists_16, predictors, nperm = 9999)

# Aphid relatedness only, same as Mantel's test, Fig S5D
predictors <- list(X1 = Aphid.dists_16)
result_16_aph <- MMRR(Ham.dists_16, predictors, nperm = 9999)

# predictors interaction model, Fig S5C
result_16_par_predictors_interaction <- MMRR(Para.dists_16, predictors, nperm = 9999)



#############
### Plant ###
#############

### Plant & Aphid distance as predictors, Hamiltonella as dependent matrix
# 31 species, reduced model
Plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_31species.csv", header=T, row.names = 1)
Plant.dists_31 <- as.matrix(vegdist(Plant, method = "bray")) # Bray–Curtis dissimilarity
Plant.dists_31  <- 1-Plant.dists_31                               #as we are using 1-Bray–Curtis, similarity
Aphid <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_31species.csv", header=T, row.names = 1)
Aphid.dists_31 <- as.matrix(Aphid)
Aphid.dists_31 <- 1-Aphid.dists_31                              #as we are using 1-Bray–Curtis, similarity
Ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv", header=T, row.names = 1)
Ham.dists_31 <- as.matrix(vegdist(Ham, method = "bray")) # Bray–Curtis dissimilarity
Ham.dists_31 <- 1-Ham.dists_31                              #as we are using 1-Bray–Curtis, similarity

# Combine predictor matrices into a list
predictors <- list(X1 = Plant.dists_31, X2 = Aphid.dists_31)
# Run MMRR, Fig 3H
result_31_pla_aph <- MMRR(Ham.dists_31, predictors, nperm = 9999)

# single predictor models
# Host plant only, same as Mantel's test, Fig 3E
predictors <- list(X1 = Plant.dists_31)
result_31_pla <- MMRR(Ham.dists_31, predictors, nperm = 9999)

# Aphid relatedness only, same as Mantel's test, Fig 3F
predictors <- list(X1 = Aphid.dists_31)
result_31_aph <- MMRR(Ham.dists_31, predictors, nperm = 9999)

# predictors interaction model, Fig 3G
result_31_pla_predictors_interaction <- MMRR(Plant.dists_31, predictors, nperm = 9999)





### 16 species, reduced model
Plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_16species.csv", header=T, row.names = 1)
Plant.dists_16 <- as.matrix(vegdist(Plant, method = "bray")) # Bray–Curtis dissimilarity
Plant.dists_16  <- 1-Plant.dists_16                               #as we are using 1-Bray–Curtis, similarity
Aphid <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv", header=T, row.names = 1)
Aphid.dists_16 <- as.matrix(Aphid)
Aphid.dists_16 <- 1-Aphid.dists_16                              #as we are using 1-Bray–Curtis, similarity
Ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1)
Ham.dists_16 <- as.matrix(vegdist(Ham, method = "bray")) # Bray–Curtis dissimilarity
Ham.dists_16 <- 1-Ham.dists_16                              #as we are using 1-Bray–Curtis, similarity

# Combine predictor matrices into a list
predictors <- list(X1 = Plant.dists_16, X2 = Aphid.dists_16)
# Run MMRR
result_16_pla_aph <- MMRR(Ham.dists_16, predictors, nperm = 9999)


# single predictor models
# Host plant only, same as Mantel's test
predictors <- list(X1 = Plant.dists_16)
result_16_pla <- MMRR(Ham.dists_16, predictors, nperm = 9999)

# Aphid relatedness only, same as Mantel's test
predictors <- list(X1 = Aphid.dists_16)
result_16_aph <- MMRR(Ham.dists_16, predictors, nperm = 9999)

# predictors interaction model
result_16_pla_predictors_interaction <- MMRR(Plant.dists_16, predictors, nperm = 9999)

# Combine predictor matrices into a list
predictors <- list(X1 =Para.dists_16, X2 =Plant.dists_16, X3 = Aphid.dists_16)
# Run MMRR
result_16_all <- MMRR(Ham.dists_16, predictors, nperm = 9999)

