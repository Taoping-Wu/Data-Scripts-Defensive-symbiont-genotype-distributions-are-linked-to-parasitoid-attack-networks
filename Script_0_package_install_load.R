###################################
###                             ###
### R packages installation     ###
### renv environment activation ###
###                             ###
###################################

# Install renv package
install.packages("renv")

# activate renv environment
renv::activate()
renv::restore()

# If necessary, set the R package installation pathway to the renv library
#.libPaths(renv::paths$library())

# write all library information to renv.lock file
# renv::snapshot()





######################################
#                                    #
# When The renv package doesn't work #
# Install the following packages     #
#                                    #
######################################


# Overall packages
install.packages("rstudioapi")
install.packages("BiocManager")

# Packages used in Script 1
install.packages("tidyr")
install.packages("dplyr")
install.packages("ape")


# Packages used in Script 2
install.packages("ggplot2")

# Packages used in Script 3
install.packages("vegan")

# Packages used in Script 4
install.packages("igraph")
install.packages("writexl")

# Packages used in Script 6
install.packages("bipartite")

# Packages used in Script 7
install.packages("ggpmisc")
install.packages("lme4")

# Packages used in Script 8
install.packages("picante")
install.packages("geiger")
install.packages("ade4")

# Packages used in Script 9 (Optional)
# When using the BiocManager to install packages, please avoid spaces in paths
# e.g. C:/Desktop/Data_Script
BiocManager::install("dada2")

