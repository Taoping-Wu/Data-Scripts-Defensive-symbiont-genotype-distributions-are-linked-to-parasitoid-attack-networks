###################################################################
# Plot: Fig. 2&3                                                  #
#                                                                 #
# MMRR analysis figure and aphid species linkage figure           #
#                                                                 #
# Figure Size: W = 720, H =720                                    #
###################################################################

############################################
### create aphid genetic distance matrix ###
############################################

############################
### MMRR analysis Figure ###
############################
rm(list=ls())

library(writexl)
library(vegan)
library(igraph)
library(dplyr)
library(ggplot2)

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




###################
### MMRR Figure ###
###################

# To use a single continuous value in all analysis of this paper, we use similarity (1- Bray-Curtis dissimilarity) in all analysis


########################
### Fig 3A & Fig S5A ###
########################

### aphid distance-Hamiltonella community figure

### Function - plot the similarity figure For Fig 3A & Fig S5A
plot_similarity_Fig3A <- function(aph.para, aph.symb) {

  # Data transformation: Calculate Bray-Curtis dissimilarity
  aph.y.dists <- vegdist(aph.para, method = "bray") # Parasitoid dissimilarity
  aph.x.dists <- vegdist(aph.symb, method = "bray") # Symbiont dissimilarity
  
  # Convert dissimilarity to similarity
  aph.y.similarity <- 1 - aph.y.dists
  aph.x.similarity <- 1 - aph.x.dists
  
  # Combine similarities into a data frame
  data <- data.frame(
    Parasitoid_Similarity = aph.y.similarity,
    Symbiont_Similarity = aph.x.similarity
  )
  
  # Generate the plot using ggplot2
  p <- ggplot(data, aes(x = Parasitoid_Similarity, y = Symbiont_Similarity)) +
    geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1, color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "blue", level = 0.95) +  # Regression line with CI
    labs(
      x = "Parasitoid community similarity",
      y = "Symbiont composition similarity"
    ) +
    ylim(0, 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    )
  
  # Print and return the plot
  print(p)
  return(p)
}

# Fig 3A, load data - Parasitoid community-Hamiltonella all years data,22 species. 
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_22species.csv", header=T, row.names = 1)     #Aphid matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3A(aph.para, aph.symb)
# Fig S5A, load data - Parasitoid community-Hamiltonella all years data,31 species.
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_16species.csv", header=T, row.names = 1)     #Aphid matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3A(aph.para, aph.symb)


##########################
### Fig 3B,F & Fig S5D ###
##########################

### aphid distance-Hamiltonella community figure

### Function - plot the similarity figure For Fig 3B,F & Fig S5D
plot_similarity_Fig3BF <- function(aph, aph.symb) {
  
  # data transformation
  aph.y.dists <- as.dist(aph)                       # Aphid relatedness distance
  aph.x.dists <- vegdist(aph.symb, method = "bray") # Bray–Curtis dissimilarity

  # convert dissimilarity to similarity
  aph.y.similarity <- 1 - aph.y.dists
  aph.x.similarity <- 1 - aph.x.dists

  # For figure plotting, combine two matrix into a dataframe
  data <- data.frame(aph.x.similarity, aph.y.similarity)

  # plotting in ggplot
  p <- ggplot(data, aes(y = aph.x.similarity, x = aph.y.similarity)) +
    geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    labs(x = "Aphid genetic distance", 
         y = "Symbiont composition similarity") +
    ylim(0,1)+ # use it to modify on the top of plot
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 18),
          axis.text.y=element_text(size = 18),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank())
  print(p)
return(p)
}


# Fig 3B, load data - aphid distance-Hamiltonella all years data,22 species. 
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_22species.csv", header=T, row.names = 1)     #Aphid matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3BF(aph, aph.symb)
# Fig 3F, load data - aphid distance-Hamiltonella all years data,31 species.
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_31species.csv", header=T, row.names = 1)     #Aphid matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3BF(aph, aph.symb)
# Fig S5D, load data - aphid distance-Hamiltonella reduced data,16 species. 
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv", header=T, row.names = 1)     #Aphid matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3BF(aph, aph.symb)


########################
### Fig 3C & Fig S5B ###
########################

#aphid distance-Parasitoid community plot

### Function - plot the similarity figure For Fig 3C & Fig S5B
plot_similarity_Fig3C <- function(aph, aph.para) {
  # data transform
  aph.para.dists <- vegdist(aph.para, method = "bray") # Bray–Curtis dissimilarity 
  aph.dists  <- as.dist(aph) #
  
  # convert dissimilarity to similarity
  aph.para.similarity <- 1 - aph.para.dists
  aph.dists.similarity <- 1 - aph.dists

  # For figure plotting, combine two matrix into a dataframe
  data <- data.frame(aph.para.similarity, aph.dists.similarity)
    
  # plotting in ggplot
  p <- ggplot(data, aes(y = aph.para.similarity, x = aph.dists.similarity)) +
    geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    labs(y = "Parasitoid composition similarity", x = "Aphid genetic distance") +
    ylim(0,0.65)+ # use it to modify on the top of plot
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 18),
          axis.text.y=element_text(size = 18),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank())
  print(p)
return(p)
}
# Fig 3C, load data - aphid distance-Parasitoid data, 22 species
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_22species.csv", header=T, row.names = 1)           #Para matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_22species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3C(aph, aph.para)
# Fig S5B, load data - aphid distance-Parasitoid data, 16species, reduced data
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_16species.csv", header=T, row.names = 1)           #Para matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3C(aph, aph.para)

########################
### Fig 3D & Fig S5C ###
########################

# MMRR plot

### Fig 3D, load data - MMRR matrix data, 22 species, full model
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_22species.csv", header=T, row.names = 1)           #plant matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_22species.csv", header=T, row.names = 1)     #symb matrix
aph.ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv", header=T, row.names = 1)     #symb matrix

# data transform
aph.para.dists <- vegdist(aph.para, method = "bray") # Bray–Curtis dissimilarity 
aph.symb.dists <- vegdist(aph.ham, method = "bray")  #
aph.aph.dists  <- as.dist(aph) #

# convert dissimilarity to similarity, data based on Supplementary table S6, MMRR test result
aph.high.similarity <- 1 - (0.367*aph.para.dists + 0.325*aph.aph.dists)
aph.symb.similarity <- 1 - aph.symb.dists

# For figure plotting, combine two matrix into a dataframe
data <- data.frame(aph.high.similarity, aph.symb.similarity)

# plotting in ggplot
p <- ggplot(data, aes(x = aph.high.similarity, y = aph.symb.similarity)) +
  geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
  labs(y = "Hamiltonella composition similarity", x = "0.37 (parasitoid)+ 0.33 (Phylogenetic)") +
  #ylim(0,0.6)+ # use it to modify on the top of plot
  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18),
        axis.text.y=element_text(size = 18),
        panel.grid.major = element_blank(),  # remove major grid lines
        panel.grid.minor = element_blank())
print(p)

#### Fig S5C, load data - MMRR matrix data, 16 species, reduced model
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_16species.csv", header=T, row.names = 1)           #plant matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv", header=T, row.names = 1)     #symb matrix
aph.ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1)     #symb matrix

# data transform
aph.para.dists <- vegdist(aph.para, method = "bray") # Bray–Curtis dissimilarity 
aph.symb.dists <- vegdist(aph.ham, method = "bray")  #
aph.aph.dists  <- as.dist(aph) #

# convert dissimilarity to similarity, data based on Supplementary table S6, MMRR test result
aph.high.similarity <- 1 - (0.385*aph.para.dists + 0.469*aph.aph.dists)
aph.symb.similarity <- 1 - aph.symb.dists

# For figure plotting, combine two matrix into a dataframe
data <- data.frame(aph.high.similarity, aph.symb.similarity)

# plotting in ggplot
p <- ggplot(data, aes(x = aph.high.similarity, y = aph.symb.similarity)) +
  geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
  labs(y = "Hamiltonella composition similarity", x = "0.38 (Parasitoid)+ 0.47 (Phylogenetic)") +
  #ylim(0,0.6)+ # use it to modify on the top of plot
  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18),
        axis.text.y=element_text(size = 18),
        panel.grid.major = element_blank(),  # remove major grid lines
        panel.grid.minor = element_blank())
print(p)



########################
### Fig 3E & Fig S5E ###
########################

# aphid distance-Parasitoid plot

### Function - plot the similarity figure For Fig 3E & Fig S5E
plot_similarity_Fig3E <- function(aph.plant, aph.symb) 
{
  #  Data transform
  aph.high.dists <- vegdist(aph.plant, method = "bray") # Bray–Curtis dissimilarity 
  aph.symb.dists <- vegdist(aph.symb, method = "bray") # Bray–Curtis dissimilarity 
  
  # convert dissimilarity to similarity
  aph.high.similarity <- 1 - aph.high.dists
  aph.symb.similarity <- 1 - aph.symb.dists
  
  # For figure plotting, combine two matrix into a dataframe
  data <- data.frame(aph.high.similarity, aph.symb.similarity)
  # plotting in ggplot
  p <- ggplot(data, aes(y = aph.high.similarity, x = aph.symb.similarity)) +
    geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    labs(x = "Plant composition similarity", y = "Hamiltonella composition similarity") +
    #ylim(0,0.75)+ # use it to modify on the top of plot
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 18),
          axis.text.y=element_text(size = 18),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank())
  print(p)
return(p)
}

# Fig 3E, load data - aphid distance-Parasitoid data, 31 species, full model
aph.plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_31species.csv", header=T, row.names = 1)           #plant matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3E(aph.plant, aph.symb)
# Fig S5E, load data - aphid distance-Parasitoid data, 16 species, reduced model
aph.plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_16species.csv", header=T, row.names = 1)           #plant matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3E(aph.plant, aph.symb)


#########################
### Fig. 3G & Fig S5F ###
#########################


# aphid distance-Plant interaction

### Function - plot the similarity figure For Fig. 3G & Fig S5F
plot_similarity_Fig3G <- function(aph.plant, aph) 
{
  #data transform
  aph.plant.dists <- vegdist(aph.plant, method = "bray") # Bray–Curtis dissimilarity 
  aph.dists <- as.dist(aph) # Bray–Curtis dissimilarity 
  
  # convert dissimilarity to similarity
  aph.plant.similarity <- 1 - aph.plant.dists
  aph.similarity <- 1 - aph.dists

  # For figure plotting, combine two matrix into a dataframe
  data <- data.frame(aph.plant.similarity, aph.similarity)

  # plotting in ggplot
  p <- ggplot(data, aes(y = aph.plant.similarity, x = aph.similarity)) +
    geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    labs(y = "Plant composition similarity", x = "Aphid genetic distance") +
    ylim(0,0.6)+
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 18),
          axis.text.y=element_text(size = 18),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank())
  print(p)
  return(p)
}


# Fig 3G, load data - aphid distance-Plant data, 31 species, full model
aph.plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_31species.csv", header=T, row.names = 1)           #plant matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_31species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3G(aph.plant, aph)
# Fig S5F, load data - aphid distance-Plant data, 16 species, reduced model 
aph.plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_16species.csv", header=T, row.names = 1)           #plant matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv", header=T, row.names = 1)     #symb matrix
p <- plot_similarity_Fig3G(aph.plant, aph)

#########################
### Fig. 3H & Fig S5G ###
#########################

# MMRR mix model, aphid-plant-Hamiltonella

### Fig. 3H, load data - aphid distance-Parasitoid data, 31 species, full model 
aph.plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_31species.csv", header=T, row.names = 1)           #plant matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_31species.csv", header=T, row.names = 1)     #symb matrix
aph.ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv", header=T, row.names = 1)     #symb matrix

# data transform
aph.plant.dists <- vegdist(aph.plant, method = "bray") # Bray–Curtis dissimilarity 
aph.symb.dists <- vegdist(aph.ham, method = "bray")  #
aph.aph.dists  <- as.dist(aph) #

# convert dissimilarity to similarity
aph.high.similarity <- 1 - (0.034*aph.plant.dists + 0.197*aph.aph.dists)
aph.symb.similarity <- 1 - aph.symb.dists

# For figure plotting, combine two matrix into a dataframe
data <- data.frame(aph.high.similarity, aph.symb.similarity)

# plotting in ggplot
p <- ggplot(data, aes(x = aph.high.similarity, y = aph.symb.similarity)) +
  geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
  labs(y = "Hamiltonella composition similarity", x = "0.03 (Plant)+ 0.19 (Phylogenetic)") +
  #ylim(0,0.6)+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18),
        axis.text.y=element_text(size = 18),
        panel.grid.major = element_blank(),  # remove major grid lines
        panel.grid.minor = element_blank())
print(p)

### Fig S5G, load data - aphid distance-Parasitoid data, 16 species, reduced model 
aph.plant <- read.csv("../Data_in_process/Plant_Aphid_matrix_16species.csv", header=T, row.names = 1)           #plant matrix
aph <- read.csv("../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv", header=T, row.names = 1)     #symb matrix
aph.ham <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1)     #symb matrix

# data transform
aph.plant.dists <- vegdist(aph.plant, method = "bray") # Bray–Curtis dissimilarity 
aph.symb.dists <- vegdist(aph.ham, method = "bray")  #
aph.aph.dists  <- as.dist(aph) #

# convert dissimilarity to similarity
aph.high.similarity <- 1 - (0.04*aph.plant.dists + 0.61*aph.aph.dists)
aph.symb.similarity <- 1 - aph.symb.dists

# For figure plotting, combine two matrix into a dataframe
data <- data.frame(aph.high.similarity, aph.symb.similarity)

# plotting in ggplot
p <- ggplot(data, aes(x = aph.high.similarity, y = aph.symb.similarity)) +
  geom_point(alpha = 0.5, shape = 1, size = 4, stroke = 1,color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
  labs(y = "Hamiltonella composition similarity", x = "0.04 (Plant)+ 0.61 (Phylogenetic)") +
  #ylim(0,0.6)+
  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 18),
        axis.text.y=element_text(size = 18),
        panel.grid.major = element_blank(),  # remove major grid lines
        panel.grid.minor = element_blank())
print(p)









