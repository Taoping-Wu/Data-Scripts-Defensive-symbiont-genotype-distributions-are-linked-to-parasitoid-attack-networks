rm(list=ls())

#####################################################
## Convert the list from Supplementary table S2 #####
#####################################################

##load necessary packages
library(tidyr)
library(dplyr)
library(ape)

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### function: convert the list into an output data.frame which has aphid(rows) and para/ham/plant (columns)
aggregate_and_pivot <- function(df, names_from_column) {
  # Aggregate data
  New_data <- aggregate(list(numdup = rep(1, nrow(df))), df, length)
  
  # Pivot wider using the specified names_from_column
  library(tidyr) # Ensure tidyr is loaded for pivot_wider
  matrix <- pivot_wider(data = New_data, names_from = {{names_from_column}}, values_from = numdup)
  matrix[is.na(matrix)] <- 0 #convert the NA inside the Matrix into 0 for downstream analysis
  # Convert matrix to a data frame
  matrix <- as.data.frame(matrix)
  rownames(matrix) <- matrix[[1]] # rename the rowname using the first column(aphid species name)
  matrix <- matrix[,-1] # Remove the first column after setting it as row names
  
  return(matrix)
}

# Input: the list that only has parasitoid/Hamiltonella/Plant, then aphid species, output the analysis matrix



#all dataset, read the list of Parasitoid-aphid file #####
df <- read.csv("../Rawdata/Table_S2.csv")         ## read the list
df <- df[,c(8,7)]                                   ## wasp taxa first, aphid species second.
df <- filter(df, !is.na(Wasp.taxa))
matrix_para_22 <- aggregate_and_pivot(df, "Wasp.taxa")   ## call the function, 22 species aphid-parasitoid matrix.
write.csv(matrix_para_22, file = "../Data_in_process/Para_Aphid_matrix_22species.csv")

#all dataset, read the list of Plant-aphid file #####
df <- read.csv("../Rawdata/Table_S2.csv")         ## read the list
df <- df[,c(9,7)]                                   ## Plant first, aphid species second.
matrix_plant_31 <- aggregate_and_pivot(df, "Plant.species")   ## call the function, 31 species aphid-plant matrix.
write.csv(matrix_plant_31, file = "../Data_in_process/Plant_Aphid_matrix_31species.csv")


############################################
####read the list of Hamiltonella Aphid ####
############################################


### Hamiltonella (from 31 aphid species for Plant dataset)
df <- read.csv("../Rawdata/Table_S2.csv")        ## read the list
df <- df[,c(10,7)]                                          ## Hamiltonella species first, aphid species second.
matrix_ham_31 <- aggregate_and_pivot(df, "Hamiltonella.strain")   ## call the function, 31 species aphid-Hamiltonella matrix.
write.csv(matrix_ham_31, file = "../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv")

### filter the 22species from the whole dataset of Hamiltonella.
List_para_22 <- row.names(matrix_para_22)                                   # extract the list of 22 parasitoid species
matrix_ham_31 <- matrix_ham_31 %>% mutate(row_id = rownames(matrix_ham_31)) 
matrix_ham_22 <- matrix_ham_31 %>% filter(row_id %in% List_para_22)         # keep the 22 species from aphid-Hamiltonella matrix.
matrix_ham_22 <- matrix_ham_22[,c(1:36)]
write.csv(matrix_ham_22, file = "../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv")

### filter the 22species from the whole dataset of host plant
List_para_22 <- row.names(matrix_para_22)                                   # extract the list of 22 parasitoid species
matrix_plant_31 <- matrix_plant_31 %>% mutate(row_id = rownames(matrix_plant_31)) 
matrix_plant_22 <- matrix_plant_31 %>% filter(row_id %in% List_para_22)         # keep the 22 species from aphid-Hamiltonella matrix.
matrix_plant_22 <- matrix_plant_22[,c(1:64)]
write.csv(matrix_plant_22, file = "../Data_in_process/Plant_Aphid_matrix_22species.csv")

### filter the 16species from the whole dataset of Hamiltonella.
df <- read.csv("../Rawdata/Table_S2.csv")        ## read the list
df <- filter(df, Data.source == "This study")              ## filter the data from 2022
df <- df[,c(10,7)]                                          ## Hamiltonella species first, aphid species second.
matrix_ham_16 <- aggregate_and_pivot(df, "Hamiltonella.strain")   ## call the convert function, 16 species aphid-Hamiltonella matrix.
write.csv(matrix_ham_16, file = "../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv")

### filter the 16 species of Parasitoid from 22 species dataset. 
List_ham_16 <- row.names(matrix_ham_16)
matrix_para_22 <- matrix_para_22 %>% mutate(row_id = rownames(matrix_para_22))
matrix_para_16 <- matrix_para_22 %>% filter(row_id %in% List_ham_16)
matrix_para_16 <- matrix_para_16[,c(1:45)]
write.csv(matrix_para_16, file = "../Data_in_process/Para_Aphid_matrix_16species.csv")

### filter the 16 species of Host plant matrix from 31 species dataset. 
df <- read.csv("../Rawdata/Table_S2.csv")        ## read the list
df <- filter(df, Data.source == "This study")              ## filter the data from 2022
df_16_plant <- df[,c(9,7)]                                          ## Hamiltonella species first, aphid species second.
matrix_plant_16 <- aggregate_and_pivot(df_16_plant, "Plant.species")   ## call the convert function
matrix_plant_16 <- matrix_plant_16 %>% mutate(row_id = rownames(matrix_plant_16))
matrix_plant_16 <- matrix_plant_16 %>% filter(row_id %in% List_ham_16)
matrix_plant_16 <- matrix_plant_16[,c(1:17)] # delete the row_id column
write.csv(matrix_plant_16, file = "../Data_in_process/Plant_Aphid_matrix_16species.csv")


#####################################
### Aphid relatedness calculation ###
#####################################

# read aphid COI fasta file
dna_sequences <- read.dna("../Rawdata/Data_4_Aphid sequences.fas", format = "fasta")
# create pairwise genetic distance matrix
aphid_genetic <- dist.gene(dna_sequences, method = "pairwise", pairwise.deletion = FALSE,
                           variance = FALSE)
genetic_distances_matrix <- as.matrix(aphid_genetic)

#the biggest distance value in this genetic distance file is 1.06, to normalize the biggest to 1, we divided all value by 106 manually. 
genetic_distances_matrix <- genetic_distances_matrix/106

#This is the original file of Aphid_phylogenetic_relatedness_31species.csv, including the outgroup Adelges_cooleyi
write.csv(genetic_distances_matrix, file = "../Data_in_process/genetic distance.csv")

#delete the outgroup, output the Aphid_phylogenetic_relatedness_31species.csv file for the following analyses
genetic_distances_matrix <- genetic_distances_matrix[-1,-1]
genetic_distances_matrix <- round(genetic_distances_matrix, 9) # Keep only 9 decimal places
write.csv(genetic_distances_matrix, file = "../Data_in_process/Aphid_phylogenetic_relatedness_31species.csv")

# Subset the matrix, keep only the 22 species of Aphid for Parasitoid analysis
List_para_22 <- gsub(" ", "_",List_para_22)
subset_matrix_22 <- genetic_distances_matrix[List_para_22, List_para_22]
write.csv(subset_matrix_22, file = "../Data_in_process/Aphid_phylogenetic_relatedness_22species.csv")

# Subset the matrix, keep only the 22 species of Aphid for Parasitoid analysis
List_ham_16 <- gsub(" ", "_",List_ham_16)
subset_matrix_16 <- genetic_distances_matrix[List_ham_16, List_ham_16]
write.csv(subset_matrix_16, file = "../Data_in_process/Aphid_phylogenetic_relatedness_16species.csv")

