
###################################################################
# Bar Plot: Fig. 1                                                #
#                                                                 #
# The bar order in Fig. 1 is arranged by proportion,              #
# which is not directly supported in the ggplot2 package in R.    #
# Therefore, we manually reordered it later in Adobe Illustrator. #
###################################################################

# Data_11 and Data_12 used in this script were derived from 
# Para_Aphid_Matrix_22.csv and Plant_Aphid_Matrix_31.csv
# by merging data from OTUs associated with the same parasitoid species
# and plant species belonging.

rm(list=ls())

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# package loading
library(tidyr)
library(dplyr)
library(ggplot2)

# load figure theme art data
artdata <-  theme_bw()+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 14),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14)) 

##########################
### Fig.1A Parasitoids ###
##########################

# load dataset (Species from same Genus were merged in this dataset)
dataset <- read.csv("../Rawdata/Data_11_Parasitoid_species_aphid_22_species_Fig.1.csv", header=T, row.names = 1)


# Normalize the dataset so that each row sums to 1
normalized_dataset <- apply(dataset, 1, function(row) row / sum(row))
normalized_dataset <- t(normalized_dataset) # Transpose the result to match the original data structure
df <- as.data.frame(normalized_dataset)     # Convert the normalized matrix into a data frame

df$Area <- rownames(df)

# Convert the data to a "long" format for easier plotting
df_long <- pivot_longer(df, cols = -Area, names_to = "Species", values_to = "Proportion")

#### cluster methods comparison####

distances <- dist(normalized_dataset) # Compute distance matrix
hc <- hclust(distances) # Perform hierarchical clustering

# Plot the dendrogram to visualize similarity
plot(hc)
area_names_ordered <- rownames(dataset)[hc$order]

# Now, create an ordered factor in df_long using these ordered area names
df_long$Area_ordered <- factor(df_long$Area, levels = area_names_ordered, ordered = TRUE)

# Now, reorder 'df_long' based on this ordered factor
df_long_ordered <- df_long[order(df_long$Area_ordered), ]

# Create colour dataset
# color for parasitoid ##
Manual.Color <- rep("white",length(unique(df_long_ordered$Species)))
names(Manual.Color) <- unique(df_long_ordered$Species)
unique(df_long_ordered$Species)
DistinctColors <- c("#263238" , #Aphelinus sp 8
                    "#455a64",  #Aphelinus varipes
                    
                    "#3e2723", "#e64a19", "#ff5722", "#ff8a65", "#ffccbc", #Aphidius
                    "#e65100", "#f57c00", "#ff9800", "#ffb74d", #Aphidius rhopalosiphi
                    "#ffe0b2", #A. rosae
                    "#ffd600", "#ffff00", "#ffff8d", "#e6ee9c", "#c6ff00",
                    "#76ff03",                       #Euaphidius setiger 18
                    "#1b5e20", "#388e3c", "#4caf50", #Binodoxys19
                    "#81c784", "#c8e6c9", "#80cbc4", #Trioxys 21
                    "#006064", "#0097a7", "#00bcd4", #Ephedrus24
                    "#bbdefb",                       #Falciconus27
                    "#039be5", "#4fc3f7",            #Lysiphlebus28
                    "#1a237e", "#01579b",            #Monoctonus31 32
                    "#4a148c", "#7b1fa2", "#9c27b0", "#ba68c8", "#e1bee7","#880e4f", "#d81b60", #praon
                    #"#3e2723", "#5d4037", "#795548", "#a1887f", 
                    "#d7ccc8"#39
)
## 
Manual.Color["Aphelinus.sp.8"]<-DistinctColors[1]
Manual.Color["Aphelinus.varipes"]<-DistinctColors[2]
Manual.Color["Aphidius.absinthii"]<-DistinctColors[3]
Manual.Color["Aphidius.avenae"]<-DistinctColors[4]
Manual.Color["Aphidius.banksae"]<-DistinctColors[5]
Manual.Color["Aphidius.eadyi"]<-DistinctColors[6]
Manual.Color["Aphidius.ervi"]<-DistinctColors[7]
Manual.Color["Aphidius.funebris"]<-DistinctColors[8]
Manual.Color["Aphidius.matricariae"]<-DistinctColors[9]
Manual.Color["Aphidius.microlophii"]<-DistinctColors[10]
Manual.Color["Aphidius.rhopalosiphi"]<-DistinctColors[11]
Manual.Color["Aphidius.rosae"]<-DistinctColors[12]
Manual.Color["Aphidius.sonchi"]<-DistinctColors[13]
Manual.Color["Aphidius.sp..21"]<-DistinctColors[14]
Manual.Color["Aphidius.sp..23"]<-DistinctColors[15]
Manual.Color["Aphidius.sp..27"]<-DistinctColors[16]
Manual.Color["Aphidius.urticae"]<-DistinctColors[17]
Manual.Color["Euaphidius.setiger"]<-DistinctColors[18]
Manual.Color["Binodoxys.acalephae"]<-DistinctColors[19]
Manual.Color["Binodoxys.angelicae"]<-DistinctColors[20]
Manual.Color["Binodoxys.sp..1"]<-DistinctColors[21]
Manual.Color["Trioxys.sp..176"]<-DistinctColors[22]
Manual.Color["Trioxys.sp..39"]<-DistinctColors[23]
Manual.Color["Trioxys.sp..56"]<-DistinctColors[24]
Manual.Color["Ephedrus.californicus"]<-DistinctColors[25]
Manual.Color["Ephedrus.lacertosus"]<-DistinctColors[26]
Manual.Color["Ephedrus.plagiator"]<-DistinctColors[27]
Manual.Color["Falciconus.pseudoplatani"]<-DistinctColors[28]
Manual.Color["Lysiphlebus.fabarum"]<-DistinctColors[29]
Manual.Color["Lysiphlebus.testaceipes"]<-DistinctColors[30]
Manual.Color["Monoctonus.caricis"]<-DistinctColors[31]
Manual.Color["Monoctonus.leclanti"]<-DistinctColors[32]
Manual.Color["Areopraon.silvestre"]<-DistinctColors[33]
Manual.Color["Praon.barbatum"]<-DistinctColors[34]
Manual.Color["Praon.dorsale"]<-DistinctColors[35]
Manual.Color["Praon.gallicum"]<-DistinctColors[36]
Manual.Color["Praon.sp..37"]<-DistinctColors[37]
Manual.Color["Praon.sp..82"]<-DistinctColors[38]
Manual.Color["Praon.volucre"]<-DistinctColors[39]
Manual.Color["Toxares.deltiger"]<-DistinctColors[40]




# reorder the columns from the bigest proportion to the smallest proportion in each column.
# by reordering the csv sheet in excel, manually
df_long_ordered <- df_long_ordered %>%
  mutate(Species = factor(Species, levels = unique(Species)))

# Plot, when use the "fill = Species" the order just disappear/ when need colur, use this code, when need order, don't use
ggplot(data = df_long_ordered, aes(x = Area_ordered, y = Proportion, fill = Species)) + #, fill = Species
  geom_bar(stat = "identity", position = "stack", colour = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), expand = c(0, 0)) +
  ylab("Abundance (%)") +
  xlab("A")+
  theme_void() + 
  theme(
    axis.text.y = element_text(color = "#2E4053", size = 16, hjust = 1, face = "italic"),
    axis.title.y = element_text(angle = 90, color = "#2E4053", size = 16, face = "bold"),
    axis.ticks.length = unit(0.1, "cm"),
    legend.position = "bottom",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(10, 10, 10, 10),
    legend.text = element_text(size = 10, colour = "#2E4053", face = "italic"),
    legend.title = element_text(face = "bold", margin = ggplot2::margin(5, 5, 5, 5, "pt"), size = 12, colour = "#2E4053"),
    plot.title = element_text(color = "#2E4053", margin = ggplot2::margin(10, 10, 10, 10, "pt"), size = 15, face = "bold")
  ) +
  scale_fill_manual(values = Manual.Color) +
  guides(fill = guide_legend(nrow = 8)) +
  ggtitle("Parasitoid Proportion of Aphid species") +
  coord_flip()



####################
### Fig.1B Plant ###
####################

# load dataset (Species from same Genus were merged in this dataset)
dataset <- read.csv("../Rawdata/Data_12_Plant_genus_aphid_31_species_Fig.1.csv", header=T, row.names = 1)


# Normalize the dataset so that each row sums to 1
normalized_dataset <- apply(dataset, 1, function(row) row / sum(row))
normalized_dataset <- t(normalized_dataset) # Transpose the result to match the original data structure
df <- as.data.frame(normalized_dataset)     # Convert the normalized matrix into a data frame

# Add a column for the area names (original row names)
df$Area <- rownames(df)

# Convert the data to a "long" format for easier plotting
df_long <- pivot_longer(df, cols = -Area, names_to = "Species", values_to = "Proportion")



#### cluster methods comparison####

distances <- dist(normalized_dataset) # Compute distance matrix
hc <- hclust(distances) # Perform hierarchical clustering

# Plot the dendrogram to visualize similarity
plot(hc)
area_names_ordered <- rownames(dataset)[hc$order]

# Now, create an ordered factor in df_long using these ordered area names
df_long$Area_ordered <- factor(df_long$Area, levels = area_names_ordered, ordered = TRUE)


# Now, reorder 'df_long' based on this ordered factor
df_long_ordered <- df_long[order(df_long$Area_ordered), ]

## COLOURS ###
Manual.Color <- rep("white",length(unique(df_long$Species)))
names(Manual.Color) <- unique(df_long$Species)


## color for plant 31 aphid species, 53 plant genus##
DistinctColors <- c("#263238", "#37474f", "#607d8b", #Acer
                    "#424242", "#424242", "#757575", "#9e9e9e", "#e0e0e0", #Anthriscus
                    "#3e2723", "#5d4037", "#bcaaa4", 
                    "#dd2c00", "#ff3d00", "#ff6e40", "#ff9e80", "#bf360c", 
                    "#d84315", "#e64a19", "#f4511e", "#ff5722", "#ff7043", "#ffccbc", #euphorbiae
                    "#ff6f00", "#ffa000", "#ffca28", #Galium
                    "#fff8e1", #geum
                    "#827717", #glyceria
                    "#ffd600", "#ffff00", #Heracleum
                    "#ffff8d", "#eeff41", "#c6ff00","#e6ee9c", "#00796b", #lupinus
                    "#a5d6a7", #Medicago
                    "#006064", "#0097a7", "#00bcd4", "#039be5", "#4fc3f7", "#01579b", "#bbdefb",#rosa
                    "#1a237e",  #Rubus
                    "#7b1fa2", #rumex
                    "#d500f9", "#ba68c8", "#e1bee7","#f8bbd0", "#f48fb1", "#e91e63",#Sonchus
                    "#ffcdd2", "#d50000", "grey"
)
## 

Manual.Color["Abrosia"]<-DistinctColors[1]
Manual.Color["Acanthus"]<-DistinctColors[2]
Manual.Color["Acer"]<-DistinctColors[3]
Manual.Color["Aconitum"]<-DistinctColors[4]
Manual.Color["Alchemilla"]<-DistinctColors[5]
Manual.Color["Amaranthus"]<- DistinctColors[6]
Manual.Color["Ammi"]<- DistinctColors[7]
Manual.Color["Anthriscus"]<- DistinctColors[8]
Manual.Color["Aquilegia"]<-DistinctColors[9]
Manual.Color["Arctium"]<-DistinctColors[10]
Manual.Color["Artemisia"]<-DistinctColors[11]
Manual.Color["Atriplex"]<-DistinctColors[12]
Manual.Color["Avena"]<- DistinctColors[13]
Manual.Color["Bromus"]<- DistinctColors[14]
Manual.Color["Buddleja"]<- DistinctColors[15]
Manual.Color["Centranthus"]<- DistinctColors[16]
Manual.Color["Cirsium"]<- DistinctColors[17]
Manual.Color["Conyza"]<- DistinctColors[18]
Manual.Color["Dactylis"]<-DistinctColors[19]
Manual.Color["Dipsacus"]<-DistinctColors[20]
Manual.Color["Encelia"]<-DistinctColors[21]
Manual.Color["Euphorbia"]<-DistinctColors[22]
Manual.Color["Ferula"]<- DistinctColors[23]
Manual.Color["Festuca"]<- DistinctColors[24]
Manual.Color["Gallium"]<- DistinctColors[25]
Manual.Color["Geum"]<- DistinctColors[26]
Manual.Color["Glyceria"]<- DistinctColors[27]
Manual.Color["Helleborus"]<- DistinctColors[28]
Manual.Color["Heracleum"]<-DistinctColors[29]
Manual.Color["Holcus"]<-DistinctColors[30]
Manual.Color["Hordeum"]<-DistinctColors[31]
Manual.Color["Iris"]<-DistinctColors[32]
Manual.Color["Knautia"]<- DistinctColors[33]
Manual.Color["Lupinus"]<- DistinctColors[34]
Manual.Color["Medicago"]<- DistinctColors[35]
Manual.Color["Onopordum"]<- DistinctColors[36]
Manual.Color["Papaver"]<- DistinctColors[37]
Manual.Color["Pastinaca"]<- DistinctColors[38]
Manual.Color["Phalaris"]<-DistinctColors[39]
Manual.Color["Poa"]<-DistinctColors[40]
Manual.Color["Quercus"]<-DistinctColors[41]
Manual.Color["Rosa"]<-DistinctColors[42]
Manual.Color["Rubus"]<- DistinctColors[43]
Manual.Color["Rumex"]<- DistinctColors[44]
Manual.Color["Sambucus"]<- DistinctColors[45]
Manual.Color["Scirpus"]<- DistinctColors[46]
Manual.Color["Secale"]<- DistinctColors[47]
Manual.Color["Senecio"]<- DistinctColors[48]
Manual.Color["Solanum"]<-DistinctColors[49]
Manual.Color["Sonchus"]<-DistinctColors[50]
Manual.Color["Sphaeralcea"]<-DistinctColors[51]
Manual.Color["Triticum"]<-DistinctColors[52]
Manual.Color["Valerinana"]<- DistinctColors[53]

# Plot, when use the "fill = Species" the order just disappear/ when need colur, use this code, when need order, don't use
ggplot(data = df_long_ordered, aes(x = Area_ordered, y = Proportion, fill = Species)) + #, fill = Species
  geom_bar(stat = "identity", position = "stack", colour = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), expand = c(0, 0)) +
  ylab("Abundance (%)") +
  xlab("A")+
  theme_void() + 
  theme(
    axis.text.y = element_text(color = "#2E4053", size = 16, hjust = 1, face = "italic"),
    axis.title.y = element_text(angle = 90, color = "#2E4053", size = 16, face = "bold"),
    axis.ticks.length = unit(0.1, "cm"),
    legend.position = "bottom",
    legend.justification = c("left", "top"),
    legend.box.just = "left",
    legend.margin = margin(10, 10, 10, 10),
    legend.text = element_text(size = 10, colour = "#2E4053", face = "italic"),
    legend.title = element_text(face = "bold", margin = ggplot2::margin(5, 5, 5, 5, "pt"), size = 12, colour = "#2E4053"),
    plot.title = element_text(color = "#2E4053", margin = ggplot2::margin(10, 10, 10, 10, "pt"), size = 15, face = "bold")
  ) +
  scale_fill_manual(values = Manual.Color) +
  guides(fill = guide_legend(nrow = 8)) +
  ggtitle("Host Plant Proportion of Aphid species") +
  coord_flip()

