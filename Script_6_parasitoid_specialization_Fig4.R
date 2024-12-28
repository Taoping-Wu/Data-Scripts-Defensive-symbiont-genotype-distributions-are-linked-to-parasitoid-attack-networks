##########################
# Bipartite Plot: Fig. 4 #
# and H2 analysis        #
##########################

rm(list = ls())


library("bipartite")

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load aphid-Hamiltonella and aphid-parasitoid interaction matrices
ham.aph<-read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv", row.names = 1, check.names = FALSE)
par.aph<-read.csv("../Data_in_process/Para_Aphid_matrix_22species.csv", row.names = 1, check.names = FALSE)

# Check if row names of the two matrices align (e.g., aphid species match in both datasets)
cbind(rownames(par.aph),rownames(ham.aph))

#Example visualization of networks (only for checking the data is correct)
#plotweb(ham.aph)
#plotweb(par.aph)

# Initialize a new matrix to hold Hamiltonella-parasitoid interactions
ham.par<-matrix(nrow=ncol(ham.aph), ncol=ncol(par.aph), data=0)  # Creates a zero matrix
rownames(ham.par) <- colnames(ham.aph)   # Assign row names as Hamiltonella strains
colnames(ham.par) <- colnames(par.aph)   # Assign column names as parasitoid species


#for each aphid species...

#assign Hamiltonella from each strain to parasitoid species at random, proportional to the
#frequency of each parasitoid on that aphid

#i<-1
#j<-1
#k<-1

# Simulate Hamiltonella-parasitoid interactions based on observed data
for (i in 1:nrow(ham.aph)){ #going through the different aphid species
  for (j in 1:ncol(ham.aph)){ #going through each of the Hamiltonella strains
    freq <- ham.aph[i,j] #number of aphids found with that particular strain of Hamiltonella
    if (freq==0) next #skip if none
    par.assignments <- sample(colnames(par.aph), size = freq, replace = T, prob=par.aph[i,]/sum(par.aph[i,])) #assign those "individuals" of Hamiltonella to parasitoid species, proportional to the frequency of those parasitoids on the focal aphid species
    for (k in 1:length(par.assignments)){
      ham.par[j,match(par.assignments[k], colnames(ham.par))] <- ham.par[j,match(par.assignments[k], colnames(ham.par))] + 1 #add those "individuals" into the final generated network
    }
  }
}

# Plotting Hamiltonella-aphid, parasitoid-aphid, and Hamiltonella-parasitoid networks
par(mfrow=c(1,1),mar = c(1.1, 1.1, 1.1, 1.1), font = 3)
plotweb(ham.aph, text.rot=90, y.lim=c(0,2)) # Plot Hamiltonella-aphid network
mtext("(a)", side=1, line=-20, at=-0.15, cex=2)
plotweb(par.aph, text.rot=90, y.lim=c(0,2)) # Plot parasitoid-aphid network
mtext("(b)", side=1, line=-20, at=-0.15, cex=2)
plotweb(ham.par, text.rot=90, y.lim=c(0,2)) # Plot Hamiltonella-parasitoid network
mtext("(c)", side=1, line=-20, at=-0.15, cex=2)


# Calculate total interactions in each matrix
sum(ham.aph) # Total Hamiltonella-aphid interactions
sum(par.aph) # Total parasitoid-aphid interactions
sum(ham.par) # Total Hamiltonella-parasitoid interactions


# Null model analysis for network specialization (H2') for each network, ps: the SES will fluctuate a bit in random processes

#Hamiltonella - Aphid specialization
nulls.ham.aph <- nullmodel(ham.aph, N=1000, method="r2dtable")
outputs.ham.aph<-vector(length=1000)
for (i in 1:1000){
  outputs.ham.aph[i]<-networklevel(nulls.ham.aph[[i]], index="H2")
}
SES.ham.aph <- (networklevel(ham.aph, index="H2")-mean(outputs.ham.aph))/sd(outputs.ham.aph)

#Parasitoid - Aphid specialization
nulls.par.aph <- nullmodel(par.aph, N=1000, method="r2dtable")
outputs.par.aph<-vector(length=1000)
for (i in 1:1000){
  outputs.par.aph[i]<-networklevel(nulls.par.aph[[i]], index="H2")
}
SES.par.aph <- (networklevel(par.aph, index="H2")-mean(outputs.par.aph))/sd(outputs.par.aph)

#Hamiltonella - Parasitoid specialization
nulls.ham.par <- nullmodel(ham.par, N=1000, method="r2dtable")
outputs.ham.par<-vector(length=1000)
for (i in 1:1000){
  outputs.ham.par[i]<-networklevel(nulls.ham.par[[i]], index="H2")
}
SES.ham.par <- (networklevel(ham.par, index="H2")-mean(outputs.ham.par))/sd(outputs.ham.par)


#Plotting the result
par(mfrow=c(3,1),mar = c(4.1, 4.1, 1.1, 2.1))
hist(outputs.ham.aph, xlim=c(0,1), ylim=c(0,350), breaks = c((1:200)/200), xlab=NA, main=NA) #, main="Hamiltonella-aphid network", xlab="Network specialisation (H2')", SES = 153.3
abline(v=networklevel(ham.aph, index="H2"), lty=2)
mtext("(a)", side=1, line=-9, at=-0.15)
hist(outputs.par.aph, xlim=c(0,1), ylim=c(0,350), breaks = c((1:200)/200), xlab=NA, main=NA) #"Parasitoid-aphid network", SES = 34
abline(v=networklevel(par.aph, index="H2"), lty=2)
mtext("(b)", side=1, line=-9, at=-0.15)
hist(outputs.ham.par, xlim=c(0,1), ylim=c(0,350), breaks = c((1:200)/200), xlab="Network specialisation (H2')", main=NA)#"Hamiltonella-parasitoid network", SES = 48.3
abline(v=networklevel(ham.par, index="H2"), lty=2)
mtext("(c)", side=1, line=-9, at=-0.15)

