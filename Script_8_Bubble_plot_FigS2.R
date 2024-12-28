################################
# Fig. S2                      #
#                              #
# Bubble plot and co-phylogeny #
################################

rm(list = ls())

# load packages
library(ape)
library(picante)
library(geiger)
library(ade4)

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load phylogeny tree file 
tree <- read.tree("../Rawdata/Data_8_Aphid_species_phylogeny.txt")
tree2 <- read.tree("../Rawdata/Data_9_Hamiltonella_phylogeny.txt")

#ultrametric (for plotting characters suggest not using ultrametric as tree becomes squished)
chronopl(tree, lambda=0.1) -> treeUltra
#check if any tree contains the 0 distance, 
plot(treeUltra, cex=0.9)
tree$tip.label

#ultrametric (for plotting characters suggest not using ultrametric as tree becomes squished)
chronopl(tree2, lambda=0.1) -> treeUltra2
#check if any tree contains the 0 distance, 
plot(treeUltra2, cex=0.9)
tree$tip.label

# load Hamiltonella community across aphid species matrix
phenodata <- read.csv(file = "../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv")

# Logarithmic transformation: Added 1 to all values beforehand to prevent issues with zeros 
# (log(0) is undefined, so 0 remains 0 after transformation).
phenodata[,-1] <- phenodata[,-1]+1
phenodata[,-1] <- log2(phenodata[,-1])

#convert data in dataframe - breaks a matrix so that each column is a distinct object
colnames(phenodata)[1] <- "aphid"
phenodata$aphid <- gsub(" ", "_", phenodata$aphid)
#make data frame for tip labels
aphid<-tree$tip.label

#make data frame with species names in same order as they are in tree
tipdata<-data.frame(aphid, aphid=phenodata$aphid[match(aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X345=phenodata$X345[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X4089=phenodata$X4089[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X93=phenodata$X93[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, N_P1_54=phenodata$N_P1_54[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X93=phenodata$X93[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X201=phenodata$X201[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, N_2_65=phenodata$N_2_65[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X2567=phenodata$X2567[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X3467=phenodata$X3467[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X3562=phenodata$X3562[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X2611=phenodata$X2611[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, N_1_88=phenodata$N_1_88[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X3521=phenodata$X3521[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X5088=phenodata$X5088[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X233=phenodata$X233[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X2578=phenodata$X2578[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X2687=phenodata$X2687[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X3004=phenodata$X3004[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X2184=phenodata$X2184[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X4981=phenodata$X4981[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X1205=phenodata$X1205[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X990=phenodata$X990[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X5032=phenodata$X5032[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X2825=phenodata$X2825[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X3448=phenodata$X3448[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X1014=phenodata$X1014[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X51=phenodata$X51[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X1473=phenodata$X1473[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X573=phenodata$X573[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, N_4_75=phenodata$N_4_75[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, N_4_15=phenodata$N_4_15[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X1485=phenodata$X1485[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X1516=phenodata$X1516[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, N_P4_6=phenodata$N_P4_6[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X231=phenodata$X231[match(tipdata$aphid,phenodata$aphid)])
tipdata<-data.frame(tipdata, X2910=phenodata$X2910[match(tipdata$aphid,phenodata$aphid)])
tipdata

#add names to the rows in our data frame, then make dataframe the default.
rownames(phenodata) <- phenodata[,1]
attach(phenodata)

############### the following two dataset is for the label and for tree only
#######################################
### for Label, the plot is aligned  ###
#######################################
#Match row names of dataframe with tip.label of phylogeny 
phenodata <-phenodata[tree$tip.label,]
#Check overlap of phenotypic traits with tree
name.check(treeUltra, phenodata)
#Save two lists from name.check
name.check(treeUltra, phenodata) -> phenodataOverlap
phenodataOverlap
#Use names in $tree.not.data in drop tip function
drop.tip(treeUltra, phenodataOverlap$tree_not_data) -> ComparativeTree
ComparativeTree$tip.label
#Match row names of dataframe with tip.label of phylogeny 
phenodata <-phenodata[ComparativeTree$tip.label,]

#plot tree  W:1373, H:705
plot(ComparativeTree, use.edge.length = TRUE, node.pos = 1, show.tip.label = TRUE, show.node.label = TRUE,cex=1.3,label.offset =3.8, no.margin=TRUE, edge.color="black")

###add the dash lines 
# Determine the rightmost x coordinate of the tree to start the dashed lines
x_start <- max(ComparativeTree$edge.length) * 2.1 # Adjust as needed

# Set the x coordinate for the labels further to the right
x_end <- x_start * 5.65 # Adjust this value based on your plot dimensions

# The y coordinates are the positions of the tips
y_positions <- 1:Ntip(ComparativeTree)

# Add labels and dashed lines for each tip
for(i in seq_along(y_positions)) {
  y_coord <- y_positions[i]
  
  # Add the label to the right
  #text(x = x_end, y = y_coord, labels = ComparativeTree$tip.label[i], pos = 4, cex = 2)
  
  # Draw a faint dashed line connecting the tip to its label
  segments(x0 = x_start, y0 = y_coord, x1 = x_end, y1 = y_coord, lty = 2, col = "#e0e0e0")
}




##

# Loop over intervals of 0.1 within the distance range
for(x in seq(0.843, 4.442, by = 0.1)) {
  abline(v = x, lty = 2, col = "#e0e0e0") # Add vertical dashed line at each interval
}
#plot tip characters - standard tree
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X4089*1, lwd=0.5, adj=0.6)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X93*1, lwd=0.5, adj=0.7)
tiplabels(pch = 21, bg = "#80d8ff", col="black", cex =phenodata$N_P1_54*1, lwd=0.5, adj=0.8)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X345*1, lwd=0.5, adj=0.9)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X201*1, lwd=0.5, adj=1)

tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X2567*1, lwd=0.5, adj=1.1)
tiplabels(pch = 21, bg = "#80d8ff", col="black", cex =phenodata$N_2_65*1, lwd=0.5, adj=1.2)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X3467*1, lwd=0.5, adj=1.3)
tiplabels(pch = 21, bg = "#80d8ff", col="black", cex =phenodata$N_1_88*1, lwd=0.5, adj=1.4)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X3521*1, lwd=0.5, adj=1.5)

tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X5088*1, lwd=0.5, adj=1.6)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X2611*1, lwd=0.5, adj=1.7)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X3562*1, lwd=0.5, adj=1.8)
tiplabels(pch = 21, bg = "#80d8ff", col="black", cex =phenodata$N_2_83*1, lwd=0.5, adj=1.9)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X1205*1, lwd=0.5, adj=2)

tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X4981*1, lwd=0.5, adj=2.1)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X990*1, lwd=0.5, adj=2.2)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X2825*1, lwd=0.5, adj=2.3)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X5032*1, lwd=0.5, adj=2.4)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X3448*1, lwd=0.5, adj=2.5)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X1014*1, lwd=0.5, adj=2.6)

tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X1473*1, lwd=0.5, adj=2.7)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X51*1, lwd=0.5, adj=2.8)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X573*1, lwd=0.5, adj=2.9)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X231*1, lwd=0.5, adj=3.00)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X2910*1, lwd=0.5, adj=3.1)

tiplabels(pch = 21, bg = "#80d8ff", col="black", cex =phenodata$N_4_15*1, lwd=0.5, adj=3.2)
tiplabels(pch = 21, bg = "#80d8ff", col="black", cex =phenodata$N_4_75*1, lwd=0.5, adj=3.3)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X1516*1, lwd=0.5, adj=3.40)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X1485*1, lwd=0.5, adj=3.5)
tiplabels(pch = 21, bg = "#80d8ff", col="black", cex =phenodata$N_P4_6*1, lwd=0.5, adj=3.60)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X233*1, lwd=0.5, adj=3.7)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X2578*1, lwd=0.5, adj=3.80)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X3004*1, lwd=0.5, adj=3.9)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X2687*1, lwd=0.5, adj=4.00)
tiplabels(pch = 21, bg = "dodgerblue3", col="black", cex =phenodata$X2184*1, lwd=0.5, adj=4.1)





#########################
### phylogeny drawing ###
#########################

# We are using the evolutionary distance phylogeny, so we draw it sepreately, 
# then combine the bubble plot (above) and the phylogeny tree (below) together.

##################
### Aphid tree ###
##################

tree <- read.tree("../Rawdata/Data_8_Aphid_species_phylogeny.txt")

phenodata <- read.csv(file = "../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv")
#convert data in dataframe - breaks a matrix so that each column is a distinct object
colnames(phenodata)[1] <- "aphid"
phenodata$aphid <- gsub(" ", "_", phenodata$aphid)
#make data frame for tip labels
aphid<-tree$tip.label

rownames(phenodata) <- phenodata[,1]
attach(phenodata)

#Match row names of dataframe with tip.label of phylogeny 
phenodata <-phenodata[tree$tip.label,]
#Check overlap of phenotypic traits with tree
name.check(tree, phenodata)
#Save two lists from name.check
name.check(tree, phenodata) -> phenodataOverlap
phenodataOverlap
#Use names in $tree.not.data in drop tip function
drop.tip(tree, phenodataOverlap$tree_not_data) -> ComparativeTree
ComparativeTree$tip.label
#Match row names of dataframe with tip.label of phylogeny 
phenodata <-phenodata[ComparativeTree$tip.label,]

#plot Aphid tree  W:1373, H:705
plot(ComparativeTree, use.edge.length = TRUE, node.pos = 1, show.tip.label = TRUE, show.node.label = TRUE,cex=1.3,label.offset =3.8, no.margin=TRUE, edge.color="black")




######################
## Hamiltonella tree##
######################

tree2 <- read.tree("../Rawdata/Data_9_Hamiltonella_phylogeny.txt")

# load data and convert the Hamiltonella name into row name
phenodata <- read.csv(file = "../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv", row.names = 1)
phenodata <- as.data.frame(t(phenodata))
rownames(phenodata) <- gsub("^X", "", rownames(phenodata))
phenodata <- cbind(Ham = rownames(phenodata), phenodata)

#convert data in dataframe - breaks a matrix so that each column is a distinct object
#make data frame for tip labels
Ham<-tree2$tip.label


#Match row names of dataframe with tip.label of phylogeny 
phenodata <-phenodata[tree2$tip.label,]
#Check overlap of phenotypic traits with tree
name.check(tree2, phenodata)
#Save two lists from name.check
name.check(tree2, phenodata) -> phenodataOverlap
phenodataOverlap
#Use names in $tree.not.data in drop tip function
drop.tip(tree2, phenodataOverlap$tree_not_data) -> ComparativeTree
ComparativeTree$tip.label
#Match row names of dataframe with tip.label of phylogeny 
phenodata <-phenodata[ComparativeTree$tip.label,]

#plot Hamiltonella tree  W:1373, H:705
plot(ComparativeTree, use.edge.length = TRUE, node.pos = 1, show.tip.label = TRUE, show.node.label = TRUE,cex=1.3,label.offset =0.1, no.margin=TRUE, edge.color="black")

