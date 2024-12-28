rm(list =ls())
######################
### Fig 2 & Fig S4 ###
######################

### Figure Size: W = 720, H =419  ###

library(writexl)
library(vegan)
library(igraph)
library(dplyr)
library(ggplot2)

create_aphid_mantel <- function(aph.high, aph.symb, Axis) {
  
  #create matrices for pairwise distances between aphid species in terms of 
  #highsitoid species/ hostplant disimilarity (in code I call it "high"), and symbiont dissimilarity
  aph.high.dists <- vegdist(aph.high, method = "bray") # Bray–Curtis dissimilarity 
  aph.symb.dists <- vegdist(aph.symb, method = "bray") # Bray–Curtis dissimilarity
  
  # convert dissimilarity to similarity
  aph.high.similarity <- 1 - aph.high.dists
  aph.symb.similarity <- 1 - aph.symb.dists
  
  # run Mantel test - this first calculates a standard correlation coefficient between the two matrices, disregarding the fact that they represent pairwise comparisons, so this cannot be intepreted in the usual way
  # the test then repeatedly shuffles one of the matrices at random, recalculating the correlation coefficient each time. This generates a null distribution, which the observed value is compared with to generate a P-value
  mantel.result <- mantel(aph.high.similarity, aph.symb.similarity, method = "pearson", permutations = 9999, na.rm = TRUE)
  ##plotting in base R
  par(mfrow=c(1,2))
  return(mantel.result) 
}




### function that make the aphid-ham high/plant network
create_aphid_network <- function(hostinfo_path, aph.para, aph.symb) {
  # create matrices for pairwise distances between aphid species in terms of highsitoid species disimilarity, and symbiont dissimilarity
  aph.para.dists <- vegdist(aph.para, method = "bray")
  aph.symb.dists <- vegdist(aph.symb, method = "bray")
  
  # similarity
  aph.para.similarity <- 1- aph.para.dists
  aph.symb.similarity <- 1- aph.symb.dists
  
  para_matrix <- as.matrix(aph.para.similarity)
  symb_matrix <- as.matrix(aph.symb.similarity)
  # Sum the two matrices
  sum_matrix <- para_matrix + symb_matrix
  
  
  # Initialize an empty data frame to store the results
  result_para <- data.frame(Species_A = character(),
                            Species_B = character(),
                            Sum_Value = numeric(),
                            para_Value = numeric(),
                            symb_Value = numeric())
  
  # Loop through the matrix to populate the data frame
  for(i in 1:(nrow(sum_matrix))) {
    for(j in i:(ncol(sum_matrix))) {  # Only consider upper triangular part
      if (i != j) {  # Skip diagonal
        new_row <- data.frame(Species_A = rownames(sum_matrix)[i],
                              Species_B = colnames(sum_matrix)[j],
                              Sum_Value = sum_matrix[i, j],
                              para_Value = para_matrix[i, j],
                              symb_Value = symb_matrix[i, j])
        result_para <- rbind(result_para, new_row)
      }
    }
  }
  
  
  
  # pick positive correlated samples as 
  links_para <- result_para %>% 
    filter(para_Value != "0") %>%
    filter(symb_Value != "0") %>% 
    filter(symb_Value != "NaN")
  # 
  
  #read the host (grass/herb/tree) info
  hostinfo <- read.csv(file=hostinfo_path)
  
  
  #create vertices info based on the aphid name and merge host plant info
  Aphidlist <- as.data.frame(unique(c(result_para$Species_A, result_para$Species_B)))
  colnames(Aphidlist) <-"Aphid" 
  vertices <- merge(Aphidlist, hostinfo, by = "Aphid", all.x = T)
  
  # create the igraph matrix
  g <- graph_from_data_frame(d=links_para, vertices =vertices, directed=F)
  
  #define color for Herb, Grass, Tree aphdis.
  colrs <- c("#c0fa4d", "#e7c736", "#eeb0af")
  V(g)$color <- colrs[V(g)$Host]
  samplesize <- rowSums(aph.symb)
  V(g)$size <- log10(samplesize)*5+2
  E(g)$weight <- links_para$Sum_Value
  # Plot the graph
  
  plot(g,edge.width = E(g)$weight*10,
       vertex.label.dist = -0.4, 
       vertex.label.cex = 1.2, 
       vertex.label.font = 3,
       vertex.label.color = "black",
       edge.color = "#bdbdbd")
  plot(g,edge.width = E(g)$weight*10,
       vertex.label = NA,
       edge.color = "#bdbdbd")
  return(result_para)
}



######################################
####Parasitoids tripartite analysis###
######################################

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Figure Size: W = 720, H =419  ###

# Fig 2A
#aphid-para all years data. 
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_22species.csv", header=T, row.names = 1)   #para matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv", header=T, row.names = 1)    #symb matrix
create_aphid_mantel(aph.para, aph.symb)      # Mantel's test
result <- create_aphid_network("../Rawdata/Data_10_aphid_host_info.csv", aph.para, aph.symb)  # network figure
write.csv(result, file = "../Data_in_process/para_similarity_22species.csv") # output for Supplementary table S4

# Fig S4A
#aphid-para only for 2021-2022 year data, Mummies only
aph.para <- read.csv("../Data_in_process/Para_Aphid_matrix_16species.csv", header=T, row.names = 1)         #para matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1) #symb matrix
create_aphid_mantel(aph.para, aph.symb)                              # Mantel's test
result <- create_aphid_network("../Rawdata/Data_10_aphid_host_info.csv", aph.para, aph.symb)  # network figure
write.csv(result, file = "../Data_in_process/para_similarity_16species.csv")

#################################
####plants tripartite analysis###
#################################

### Figure Size: W = 720, H =419  ###

# Fig 2B
# load data - aphid-plant all years data. 
aph.para <- read.csv("../Data_in_process/Plant_Aphid_matrix_31species.csv", header=T, row.names = 1)           #plant matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv", header=T, row.names = 1)     #symb matrix
create_aphid_mantel(aph.para, aph.symb)                      # Mantel's test
result <- create_aphid_network("../Rawdata/Data_10_aphid_host_info.csv", aph.para, aph.symb)    # network figure
write.csv(result, file = "../Data_in_process/plant_similarity_allYear.csv") # output for Supplementary table S4

# Fig S4B
# load data - only for 2021-2022 year data, Mummies only
aph.para <- read.csv("../Data_in_process/Plant_Aphid_matrix_16species.csv", header=T, row.names = 1)        #plant matrix
aph.symb <- read.csv("../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", header=T, row.names = 1) #symb matrix
create_aphid_mantel(aph.para, aph.symb)       # Mantel's test
result <- create_aphid_network("../Rawdata/Data_10_aphid_host_info.csv", aph.para, aph.symb)    # network figure
write.csv(result, file = "../Data_in_process/plant_similarity_2Year.csv")

