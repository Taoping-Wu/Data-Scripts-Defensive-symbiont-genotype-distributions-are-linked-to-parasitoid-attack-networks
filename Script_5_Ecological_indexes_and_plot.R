####################################
###Biodiversity indices analysis####
####################################
library(tidyr)
library(dplyr)
library(vegan)
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

library(vegan) # Ensure vegan is loaded for rrarefy

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
      Average_Richness = avg_richness,
      Average_Shannon_Index = avg_shannon,
      Average_Evenness = avg_evenness,
      Average_Simpson_Index = avg_simpson
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

#################
#parasitoid test#
#################

### biodiversity indices - aphid-para all years data. 

#aphid-Hamiltonella matrix
diversity <- read.csv(file = "Hamiltonella_Aphid_matrix_22species.csv", row.names = 1)##load Hamiltonella-aphid matrix
Biodiversity_ham_22 <- biodiversity_indices_df(diversity)                             ##call function
write.csv(Biodiversity_ham_22, file = "Eco_Aphid_ham_22species.csv")                  ##write result

# rarefy biodiversity indices - aphid-Hamiltonella all years data.
Biodiversity_ham_22_rarefy <- process_biodiversity_data(diversity)                 ##call rarefy test
write.csv(Biodiversity_ham_22_rarefy, file = "Eco_Aphid_ham_22species_rarefy.csv") ##write result

#aphid-parasitoid matrix
diversity <- read.csv(file = "Para_Aphid_matrix_22species_all_year.csv", row.names = 1)##load Parasitoid-aphid matrix
Biodiversity_para_22 <- biodiversity_indices_df(diversity)                             ##call function
write.csv(Biodiversity_para_22, file = "Eco_Aphid_para_22species.csv")                 ##write result
# rarefy biodiversity indices - aphid-parasitoid all years data.
Biodiversity_para_22_rarefy <- process_biodiversity_data(diversity)                  ##call rarefy test
write.csv(Biodiversity_para_22_rarefy, file = "Eco_Aphid_para_22species_rarefy.csv") ##write result

###aphid-para only for 2021-2022 year data, Mummies only
diversity <- read.csv(file = "Hamiltonella_Aphid_matrix_16species.csv", row.names = 1) ##load Hamiltonella-aphid matrix
Biodiversity_ham_16 <- biodiversity_indices_df(diversity)                                     ##call function
write.csv(Biodiversity_ham_16, file = "Eco_Aphid_ham_16species.csv")                          ##write result
diversity <- read.csv(file = "Para_Aphid_matrix_16species.csv", row.names = 1)        ##load Parasitoid-aphid matrix
Biodiversity_para_16 <- biodiversity_indices_df(diversity)                            ##call function
write.csv(Biodiversity_para_16, file = "Eco_Aphid_para_16species.csv")                ##write result







############
#plant test#
############

#aphid-Hamiltonella matrix
diversity <- read.csv(file = "Hamiltonella_Aphid_matrix_31species.csv", row.names = 1)##load Hamiltonella-aphid matrix
Biodiversity_ham_31 <- biodiversity_indices_df(diversity)                             ##call function
write.csv(Biodiversity_ham_31, file = "Eco_Aphid_ham_31species.csv")                  ##write result

# rarefy biodiversity indices - aphid-Hamiltonella all years data.
Biodiversity_ham_31_rarefy <- process_biodiversity_data(diversity)                 ##call rarefy test
write.csv(Biodiversity_ham_31_rarefy, file = "Eco_Aphid_ham_31species_rarefy.csv") ##write result

#aphid-plant matrix
diversity <- read.csv(file = "Plant_Aphid_matrix_31species_all_year.csv", row.names = 1)##load plant-aphid matrix
Biodiversity_plant_31 <- biodiversity_indices_df(diversity)                             ##call function
write.csv(Biodiversity_plant_31, file = "Eco_Aphid_plant_31species.csv")                 ##write result
# rarefy biodiversity indices - aphid-plant all years data.
Biodiversity_plant_31_rarefy <- process_biodiversity_data(diversity)                  ##call rarefy test
write.csv(Biodiversity_plant_31_rarefy, file = "Eco_Aphid_plant_31species_rarefy.csv") ##write result

###aphid-plant only for 2021-2031 year data, Mummies only
diversity <- read.csv(file = "Hamiltonella_Aphid_matrix_16species.csv", row.names = 1) ##load Hamiltonella-aphid matrix
Biodiversity_ham_16 <- biodiversity_indices_df(diversity)                                     ##call function
write.csv(Biodiversity_ham_16, file = "Eco_Aphid_ham_16species.csv")                          ##write result
diversity <- read.csv(file = "plant_Aphid_matrix_16species.csv", row.names = 1)        ##load plant-aphid matrix
Biodiversity_plant_16 <- biodiversity_indices_df(diversity)                            ##call function
write.csv(Biodiversity_plant_16, file = "Eco_Aphid_plant_16species.csv")                ##write result

##########################
### Biological indices figures ###
##########################
artdata <-  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 14),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14)) 

library(ggplot2)
library(ggpmisc)
library(tidyr)



### choose one out of three analysis dataset from Table S3
### or, manually combine the dataset derived from above
### e.g. Eco_Aphid_ham_22species.csv + Eco_Aphid_para_22species.csv = Eco_Aphid_ham_para_22species.csv

diversity <- read.csv(file = "Eco_Aphid_ham_para_22species.csv", row.names = 1)
diversity <- read.csv(file = "Eco_Aphid_ham_para_22species_rarefy_n7.csv", row.names = 1)
diversity <- read.csv(file = "Eco_Aphid_Ham_Para_16species.csv")


artdata <-  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 14),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14)) 

pd <- position_dodge(0.1)

ggplot(data = diversity, aes(x = Richness_Para ,y = Richness_Ham))+
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


ggplot(data = diversity, aes(x = Shannon_Index_Para ,y = Shannon_Index_Ham))+
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

ggplot(data = diversity, aes(x = Simpson_Index_Para ,y = Simpson_Index_Ham))+
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

library(lme4)
m1 <- lm(diversity$Richness_Ham ~ diversity$Richness_Para) 
summary(m1)

m1 <- lm(diversity$Shannon_Index_Ham ~ diversity$Shannon_Index_Para) 
summary(m1)

m1 <- lm(diversity$Simpson_Index_Ham ~ diversity$Simpson_Index_Para) 
summary(m1)


#####################
###   Plant       ###
#####################
### choose one out of three analysis dataset from Table S3
### or, manually combine the dataset derived from above
### e.g. Eco_Aphid_ham_31species.csv + Eco_Aphid_plant_31species.csv = Eco_Aphid_ham_para_31species.csv

diversity <- read.csv(file = "Eco_Aphid_ham_plant_16species.csv", row.names = 1)
diversity <- read.csv(file = "Eco_Aphid_ham_plant_31species.csv", row.names = 1)
diversity <- read.csv(file = "Eco_Aphid_ham_plant_31species_rarefy.csv", row.names = 1)

diversity_log <- log(diversity)
diversity_log[diversity_log == -Inf] <- NA

ggplot(data = diversity, aes(x = Richness_Plant,y = Richness_Ham))+
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

ggplot(data = diversity, aes(x = Shannon_Index_Plant ,y = Shannon_Index_Ham))+
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

ggplot(data = diversity, aes(x = Simpson_Index_Plant ,y = Simpson_Index_Ham))+
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

m1 <- lm(diversity$Richness_Ham ~ diversity$Richness_Plant)
summary(m1)
m1 <- lm(diversity$Shannon_Index_Ham ~ diversity$Shannon_Index_Plant) 
summary(m1)
m1 <- lm(diversity$Simpson_Index_Ham ~ diversity$Simpson_Index_Plant) 
summary(m1)
