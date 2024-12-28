########################
#       Table 1        #
#                      #
# Biodiversity indices #
########################

### Figure Size: W = 720, H =720 ###

rm(list = ls())

####################################
###Biodiversity indices analysis####
####################################
#load packages
library(tidyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpmisc)
library(lme4)
## Normal analysis

biodiversity_indices_df <- function(diversity) {
  # Apply biodiversity indices calculation to each row
  results <- apply(diversity, 1, function(row) {
    total_individuals <- sum(row)
    richness <- sum(row > 0)
    p <- row / total_individuals
    p <- p[p > 0]
    H <- -sum(p * log(p))
    evenness <- H / log(richness)
    simpson_index <- 1 - sum(p^2)
    
    return(list(
      Richness = richness,
      Total_Individuals = total_individuals,
      Shannon_Index = H,
      Evenness = evenness,
      Simpson_Index = simpson_index
    ))
  })
  
  # Convert the results to a data frame
  results_df <- do.call(rbind, results)
  rownames(results_df) <- rownames(diversity)
  results_df <- as.data.frame(results_df)
  print(results_df)
  
  
  results_df <- data.frame(lapply(results_df, function(x) {
    if(is.list(x)) unlist(x) else x
  }))
  
  return(results_df)
}


###rarefaction analysis####

# Define the main function to process the diversity data
process_biodiversity_data <- function(diversity, N = 7, iterations = 1000) {
  # Inner function for calculating biodiversity indices
  biodiversity_indices_rarefied <- function(row, N, iterations) {
    if(sum(row) < N) {
      return(list(
        Richness = NA,
        Shannon_Index = NA,
        Evenness = NA,
        Simpson_Index = NA
      ))
    }
    
    richness_results <- numeric(iterations)
    shannon_results <- numeric(iterations)
    evenness_results <- numeric(iterations)
    simpson_results <- numeric(iterations)
    
    for (i in 1:iterations) {
      rarefied_row <- rrarefy(row, sample = N)
      
      if(sum(rarefied_row) >= N) {
        total_individuals <- sum(rarefied_row)
        p <- rarefied_row / total_individuals
        p <- p[p > 0]
        
        richness_results[i] <- sum(rarefied_row > 0)
        shannon_results[i] <- -sum(p * log(p))
        evenness_results[i] <- shannon_results[i] / log(richness_results[i])
        simpson_results[i] <- 1 - sum(p^2)
      }
    }
    
    avg_richness <- mean(richness_results, na.rm = TRUE)
    avg_shannon <- mean(shannon_results, na.rm = TRUE)
    avg_evenness <- mean(evenness_results, na.rm = TRUE)
    avg_simpson <- mean(simpson_results, na.rm = TRUE)
    
    return(list(
      Richness = avg_richness,
      Shannon_Index = avg_shannon,
      Evenness = avg_evenness,
      Simpson_Index = avg_simpson
    ))
  }
  
  # Apply the indices function to each row of the diversity matrix/data frame
  results <- apply(diversity, 1, function(row) biodiversity_indices_rarefied(row, N, iterations))
  
  # Convert the list of results to a data frame
  results_df <- do.call(rbind, results)
  rownames(results_df) <- rownames(diversity)
  results_df <- as.data.frame(results_df)
  
  # Ensure all elements are properly formatted
  results_df <- data.frame(lapply(results_df, function(x) if(is.list(x)) unlist(x) else x))
  
  return(results_df)
}

###################
# parasitoid test #
###################

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### biodiversity indices - aphid-para all years data. 

#aphid-Hamiltonella matrix
diversity <- read.csv(file = "../Data_in_process/Hamiltonella_Aphid_matrix_22species.csv", row.names = 1)##load Hamiltonella-aphid matrix
Biodiversity_ham_22 <- biodiversity_indices_df(diversity)                            
colnames(Biodiversity_ham_22) <- paste0(colnames(Biodiversity_ham_22), "_Ham")

# rarefy biodiversity indices - aphid-Hamiltonella all years data.
Biodiversity_ham_22_rarefy <- process_biodiversity_data(diversity)                 ##call rarefy test
colnames(Biodiversity_ham_22_rarefy) <- paste0(colnames(Biodiversity_ham_22_rarefy), "_Ham")



#aphid-parasitoid matrix
diversity <- read.csv(file = "../Data_in_process/Para_Aphid_matrix_22species.csv", row.names = 1)##load Parasitoid-aphid matrix
Biodiversity_para_22 <- biodiversity_indices_df(diversity)                             ##call function
colnames(Biodiversity_para_22) <- paste0(colnames(Biodiversity_para_22), "_Para")


# rarefy biodiversity indices - aphid-parasitoid all years data.
Biodiversity_para_22_rarefy <- process_biodiversity_data(diversity)                  ##call rarefy test
colnames(Biodiversity_para_22_rarefy) <- paste0(colnames(Biodiversity_para_22_rarefy), "_Para")


###aphid-para only for 2021-2022 year data, Mummies only
diversity <- read.csv(file = "../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", row.names = 1) ##load Hamiltonella-aphid matrix
Biodiversity_ham_16 <- biodiversity_indices_df(diversity)                                     ##call function
colnames(Biodiversity_ham_16) <- paste0(colnames(Biodiversity_ham_16), "_Ham")


diversity <- read.csv(file = "../Data_in_process/Para_Aphid_matrix_16species.csv", row.names = 1)        ##load Parasitoid-aphid matrix
Biodiversity_para_16 <- biodiversity_indices_df(diversity)                            ##call function
colnames(Biodiversity_para_16) <- paste0(colnames(Biodiversity_para_16), "_Para")



#Merge the Hamiltonella/ Parasitoid information together

# All year 22 species model
Biodiversity_Aphid_ham_para_22 <- merge(Biodiversity_ham_22,Biodiversity_para_22, by = "row.names")
write.csv(Biodiversity_Aphid_ham_para_22, file = "../Data_in_process/Eco_Aphid_ham_para_22species.csv")

# All year 22 species model, rarefied, NA values removed
Biodiversity_Aphid_ham_para_22_rarefy <- merge(Biodiversity_ham_22_rarefy,Biodiversity_para_22_rarefy, by = "row.names")
# remove NA values
Biodiversity_Aphid_ham_para_22_rarefy <- Biodiversity_Aphid_ham_para_22_rarefy[!is.na(Biodiversity_Aphid_ham_para_22_rarefy$Simpson_Index_Ham), ]
Biodiversity_Aphid_ham_para_22_rarefy <- Biodiversity_Aphid_ham_para_22_rarefy[!is.na(Biodiversity_Aphid_ham_para_22_rarefy$Shannon_Index_Para), ]
write.csv(Biodiversity_Aphid_ham_para_22_rarefy, file = "../Data_in_process/Eco_Aphid_ham_para_22species_rarefy_n7.csv")

# 2 year 16 species model
Biodiversity_Aphid_ham_para_16 <- merge(Biodiversity_ham_16,Biodiversity_para_16, by = "row.names")
write.csv(Biodiversity_Aphid_ham_para_16, file = "../Data_in_process/Eco_Aphid_Ham_Para_16species.csv")



############
#plant test#
############

#aphid-Hamiltonella matrix
diversity <- read.csv(file = "../Data_in_process/Hamiltonella_Aphid_matrix_31species.csv", row.names = 1)##load Hamiltonella-aphid matrix
Biodiversity_ham_31 <- biodiversity_indices_df(diversity)                             ##call function
colnames(Biodiversity_ham_31) <- paste0(colnames(Biodiversity_ham_31), "_Ham")

# rarefy biodiversity indices - aphid-Hamiltonella all years data.
Biodiversity_ham_31_rarefy <- process_biodiversity_data(diversity)                 ##call rarefy test
colnames(Biodiversity_ham_31_rarefy) <- paste0(colnames(Biodiversity_ham_31_rarefy), "_Ham") ##write result

#aphid-plant matrix
diversity <- read.csv(file = "../Data_in_process/Plant_Aphid_matrix_31species.csv", row.names = 1)##load plant-aphid matrix
Biodiversity_plant_31 <- biodiversity_indices_df(diversity)                             ##call function
colnames(Biodiversity_plant_31) <- paste0(colnames(Biodiversity_plant_31), "_Plant") ##write result
# rarefy biodiversity indices - aphid-plant all years data.
Biodiversity_plant_31_rarefy <- process_biodiversity_data(diversity)                  ##call rarefy test
colnames(Biodiversity_plant_31_rarefy) <- paste0(colnames(Biodiversity_plant_31_rarefy), "_Plant") ##write result

###aphid-plant only for 2021-2031 year data, Mummies only
diversity <- read.csv(file = "../Data_in_process/Hamiltonella_Aphid_matrix_16species.csv", row.names = 1) ##load Hamiltonella-aphid matrix
Biodiversity_ham_16 <- biodiversity_indices_df(diversity)                                     ##call function
colnames(Biodiversity_ham_16) <- paste0(colnames(Biodiversity_ham_16), "_Ham")


diversity <- read.csv(file = "../Data_in_process/plant_Aphid_matrix_16species.csv", row.names = 1)        ##load plant-aphid matrix
Biodiversity_plant_16 <- biodiversity_indices_df(diversity)                            ##call function
colnames(Biodiversity_plant_16) <- paste0(colnames(Biodiversity_plant_16), "_Plant") ##write result

#Merge the Hamiltonella/ Plant information together
# All year model 33 species
Biodiversity_Aphid_ham_plant_31 <- merge(Biodiversity_ham_31,Biodiversity_plant_31, by = "row.names")
write.csv(Biodiversity_Aphid_ham_plant_31, file = "../Data_in_process/Eco_Aphid_ham_plant_31species.csv")

# All year model 33 species, rarefied, NA values removed
Biodiversity_Aphid_ham_plant_31_rarefy <- merge(Biodiversity_ham_31_rarefy,Biodiversity_plant_31_rarefy, by = "row.names")
# remove NA values
Biodiversity_Aphid_ham_plant_31_rarefy <- Biodiversity_Aphid_ham_plant_31_rarefy[!is.na(Biodiversity_Aphid_ham_plant_31_rarefy$Simpson_Index_Ham), ]
Biodiversity_Aphid_ham_plant_31_rarefy <- Biodiversity_Aphid_ham_plant_31_rarefy[!is.na(Biodiversity_Aphid_ham_plant_31_rarefy$Shannon_Index_Plant), ]
write.csv(Biodiversity_Aphid_ham_plant_31_rarefy, file = "../Data_in_process/Eco_Aphid_ham_plant_31species_rarefy.csv")

# Two year model 16 species
Biodiversity_Aphid_ham_plant_16 <- merge(Biodiversity_ham_16,Biodiversity_plant_16, by = "row.names")
write.csv(Biodiversity_Aphid_ham_plant_16, file = "../Data_in_process/Eco_Aphid_ham_plant_16species.csv")

##################################
### Biological indices figures ###
##################################

# Set figure theme art parameter
artdata <-  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 14),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14)) 


#####################
###  Parasitoid   ###
#####################

plot_diversity_para <- function(diversity) 
{
  pd <- position_dodge(0.1)
  
  # Richness plot - parasitoid
  p <- ggplot(data = diversity, aes(x = Richness_Para ,y = Richness_Ham))+
    geom_point(alpha = 0.5, shape = 16, size = 4, stroke = 1,color = "#424242") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16),
          axis.text.y=element_text(size = 16),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank()) +
    ylab("Hamiltonella Richness")+
    xlab("Parasitoid Richness")+
    #xlim(-1.2,19)+
    guides(colour = "none", size = "none", shape = "none")
  print(p)
  
  # Shannon index plot - parasitoid
  p <- ggplot(data = diversity, aes(x = Shannon_Index_Para ,y = Shannon_Index_Ham))+
    geom_point(alpha = 0.5, shape = 16, size = 4, stroke = 1,color = "#424242") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16),
          axis.text.y=element_text(size = 16),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank()) +
    ylab("Hamiltonella Shannon's index")+
    xlab("Parasitoid Shannon's index")+
    #xlim(-0.25,2.2)+
    guides(colour = "none", size = "none", shape = "none")
  print(p)
  
  # Simpson index plot - parasitoid
  p <- ggplot(data = diversity, aes(x = Simpson_Index_Para ,y = Simpson_Index_Ham))+
    geom_point(alpha = 0.5, shape = 16, size = 4, stroke = 1,color = "#424242") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16),
          axis.text.y=element_text(size = 16),
          panel.grid.major = element_blank(),  # remove major grid lines
         panel.grid.minor = element_blank()) +
    ylab("Hamiltonella Simpson's index")+
   xlab("Parasitoid Simpson's index")+
   #xlim(-0.1,0.92)+
   guides(colour = "none", size = "none", shape = "none")
  print(p)
  
  # Linear model 
  # Initialize a list to store models and summaries
  model_results <- list()
  
  # Linear model of Hamiltonella/Parasitoid Richness
  model_results$Richness <- lm(diversity$Richness_Ham ~ diversity$Richness_Para)
  print(summary(model_results$Richness))  # Print summary
  
  # Linear model of Hamiltonella/Parasitoid Shannon Index
  model_results$Shannon <- lm(diversity$Shannon_Index_Ham ~ diversity$Shannon_Index_Para)
  print(summary(model_results$Shannon))  # Print summary
  
  # Linear model of Hamiltonella/Parasitoid Simpson Index
  model_results$Simpson <- lm(diversity$Simpson_Index_Ham ~ diversity$Simpson_Index_Para)
  print(summary(model_results$Simpson))  # Print summary
  
}

# load diversity data
# All year Parasitoid 22 species model
diversity <- read.csv(file = "../Data_in_process/Eco_Aphid_ham_para_22species.csv", row.names = 1)
plot_diversity_para(diversity)

# All year Parasitoid 22 species rarefied model
diversity <- read.csv(file = "../Data_in_process/Eco_Aphid_ham_para_22species_rarefy_n7.csv", row.names = 1)
plot_diversity_para(diversity)

# All year Parasitoid 16 species model
diversity <- read.csv(file = "../Data_in_process/Eco_Aphid_Ham_Para_16species.csv")
plot_diversity_para(diversity)







#####################
###   Plant       ###
#####################

plot_diversity_plant <- function(diversity) 
{
  # Richness plot - plant
  p <- ggplot(data = diversity, aes(x = Richness_Plant,y = Richness_Ham))+
    geom_point(alpha = 0.5, shape = 16, size = 4, stroke = 1,color = "#424242") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16),
          axis.text.y=element_text(size = 16),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank()) +
    ylab("Hamiltonella Richness")+
    xlab("Plant Richness")+
    #xlim(0.2,6.1)+
    guides(colour = "none", size = "none", shape = "none")
    
    print(p)
  
  # Shannon Index plot - plant
  p <- ggplot(data = diversity, aes(x = Shannon_Index_Plant ,y = Shannon_Index_Ham))+
    geom_point(alpha = 0.5, shape = 16, size = 4, stroke = 1,color = "#424242") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16),
          axis.text.y=element_text(size = 16),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank()) +
    ylab("Hamiltonella Shannon's index")+
    xlab("Plant Shannon's index")+
    #xlim(-0.2,1.8)+
    guides(colour = "none", size = "none", shape = "none")
    print(p)
    
    # Simpson Index plot - plant
  p <- ggplot(data = diversity, aes(x = Simpson_Index_Plant ,y = Simpson_Index_Ham))+
    geom_point(alpha = 0.5, shape = 16, size = 4, stroke = 1,color = "#424242") +
    geom_smooth(method = "lm", se = TRUE, color = "blue",level = 0.95) +  # add regression line with CI
    theme_bw()+
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 18),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 16),
          axis.text.y=element_text(size = 16),
          panel.grid.major = element_blank(),  # remove major grid lines
          panel.grid.minor = element_blank()) +
    ylab("Hamiltonella Simpson's Index")+
    xlab("Plant Simpson's Index")+
    #xlim(-0.11,0.87)+
    guides(colour = "none", size = "none", shape = "none")
  
    print(p)
  
  # Linear model 
  # Initialize a list to store models and summaries
  model_results <- list()
  # Linear model of Hamiltonella/Parasitoid Richness
  model_results$Richness <- lm(diversity$Richness_Ham ~ diversity$Richness_Plant)
  print(summary(model_results$Richness))  # Print summary

  # Linear model of Hamiltonella/Parasitoid Shannon Index
  model_results$Shannon <- lm(diversity$Shannon_Index_Ham ~ diversity$Shannon_Index_Plant)
  print(summary(model_results$Shannon))  # Print summary

  # Linear model of Hamiltonella/Parasitoid Simpson Index
  model_results$Simpson <- lm(diversity$Simpson_Index_Ham ~ diversity$Simpson_Index_Plant)
  print(summary(model_results$Simpson))  # Print summary
}


# load diversity plant data

# All year Plant 31 species model
diversity <- read.csv(file = "../Data_in_process/Eco_Aphid_ham_plant_31species.csv", row.names = 1)
plot_diversity_plant(diversity)

# All year Plant 31 species model, rarefied
diversity <- read.csv(file = "../Data_in_process/Eco_Aphid_ham_plant_31species_rarefy.csv", row.names = 1)
plot_diversity_plant(diversity)

# Two year Plant 16 species model, rarefied
diversity <- read.csv(file = "../Data_in_process/Eco_Aphid_ham_plant_16species.csv", row.names = 1)
plot_diversity_plant(diversity)
