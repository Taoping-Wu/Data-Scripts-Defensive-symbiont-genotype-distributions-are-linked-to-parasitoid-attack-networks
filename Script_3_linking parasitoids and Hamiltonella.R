rm(list = ls())

#install.packages("bipartite")
library("bipartite")

ham.aph<-read.csv("Hamiltonella_Aphid_matrix_22species.csv", row.names = 1, check.names = FALSE)
par.aph<-read.csv("Para_Aphid_matrix_22species.csv", row.names = 1, check.names = FALSE)

cbind(rownames(par.aph),rownames(ham.aph))

#plotweb(ham.aph)
#plotweb(par.aph)

ham.par<-matrix(nrow=ncol(ham.aph), ncol=ncol(par.aph), data=0)
rownames(ham.par) <- colnames(ham.aph)
colnames(ham.par) <- colnames(par.aph)

#for each aphid species...

#assign Hamiltonella from each strain to parasitoid species at random, proportional to the
#frequency of each parasitoid on that aphid

#i<-1
#j<-1
#k<-1

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

par(mfrow=c(1,1),mar = c(1.1, 1.1, 1.1, 1.1), font = 3)
plotweb(ham.aph, text.rot=90, y.lim=c(0,2))
mtext("(a)", side=1, line=-20, at=-0.15, cex=2)
plotweb(par.aph, text.rot=90, y.lim=c(0,2))
mtext("(b)", side=1, line=-20, at=-0.15, cex=2)
plotweb(ham.par, text.rot=90, y.lim=c(0,2))
mtext("(c)", side=1, line=-20, at=-0.15, cex=2)
sum(ham.aph)
sum(par.aph)
sum(ham.par)

nulls.ham.aph <- nullmodel(ham.aph, N=1000, method="r2dtable")
outputs.ham.aph<-vector(length=1000)
for (i in 1:1000){
  outputs.ham.aph[i]<-networklevel(nulls.ham.aph[[i]], index="H2")
}
SES.ham.aph <- (networklevel(ham.aph, index="H2")-mean(outputs.ham.aph))/sd(outputs.ham.aph)

nulls.par.aph <- nullmodel(par.aph, N=1000, method="r2dtable")
outputs.par.aph<-vector(length=1000)
for (i in 1:1000){
  outputs.par.aph[i]<-networklevel(nulls.par.aph[[i]], index="H2")
}
SES.par.aph <- (networklevel(par.aph, index="H2")-mean(outputs.par.aph))/sd(outputs.par.aph)

nulls.ham.par <- nullmodel(ham.par, N=1000, method="r2dtable")
outputs.ham.par<-vector(length=1000)
for (i in 1:1000){
  outputs.ham.par[i]<-networklevel(nulls.ham.par[[i]], index="H2")
}
SES.ham.par <- (networklevel(ham.par, index="H2")-mean(outputs.ham.par))/sd(outputs.ham.par)

par(mfrow=c(3,1),mar = c(4.1, 4.1, 1.1, 2.1))
hist(outputs.ham.aph, xlim=c(0,1), ylim=c(0,350), breaks = c((1:200)/200), xlab=NA, main=NA) #, main="Hamiltonella-aphid network", xlab="Network specialisation (H2')"
abline(v=networklevel(ham.aph, index="H2"), lty=2)
mtext("(a)", side=1, line=-9, at=-0.15)
hist(outputs.par.aph, xlim=c(0,1), ylim=c(0,350), breaks = c((1:200)/200), xlab=NA, main=NA) #"Parasitoid-aphid network"
abline(v=networklevel(par.aph, index="H2"), lty=2)
mtext("(b)", side=1, line=-9, at=-0.15)
hist(outputs.ham.par, xlim=c(0,1), ylim=c(0,350), breaks = c((1:200)/200), xlab="Network specialisation (H2')", main=NA)#"Hamiltonella-parasitoid network"
abline(v=networklevel(ham.par, index="H2"), lty=2)
mtext("(c)", side=1, line=-9, at=-0.15)