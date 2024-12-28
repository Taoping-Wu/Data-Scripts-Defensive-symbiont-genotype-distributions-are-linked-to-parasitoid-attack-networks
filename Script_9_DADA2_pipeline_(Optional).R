#####################################################################################################
## Name: DADA analysis for amplicon DNA                                                            ##
## Info: Quick run through of amplicon DNA samples                                                 ##
## Date are in NCBI bioproject PRJNA1139364                                                        ##
## accession numbers are:                                                                          ##
##      SRX25431825 to SRX25432194 (Data_1_Block 1)                                                ## 
##      SRX25431012 to SRX25431744 (Data_2_Block 2)                                                ##
## Two blocks of data have to run separately                                                       ##
## All of them are paired-end Illumina highthroughput sequencing .fastq file                       ##
#####################################################################################################

rm(list=ls() )

## load necessary library
library(dada2)
library(rentrez)

# IMPORTANT: Two blocks of data have to run separately

## Set the working directory to the Script folder ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################
### Start analysis ###
######################

########
### All fastq file from a SINGLE BLOCK should download into a single filefolder
########


### Optional, test data is included in the R project
Block <- "test"
path <- "../Rawdata/Illumina_test/"

### Block 1: 2021 year data 370 samples, download needed
Block <- "Block1"
path <- "../Rawdata/Illumina_Block1/" # download the Genbank data, extract to fastq files by sraTools, then store them in this folder

### Block 2: 2022 year data 738 samples, download needed
Block <- "Block2"
path <- "../Rawdata/Illumina_Block2/" # download the Genbank data, extract to fastq files by sraTools, then store them in this folder





list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_1.fastq"), `[`, 1)
sample.names <- sapply(strsplit(basename(fnRs), "_2.fastq"), `[`, 1)

# Check sample quality in batches - for the cut off for filtering
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

#plotQualityProfile(fnRs[1:3])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filtering - play with truncLen and maxEE - Error rate
?filterAndTrim # for extra options/arguments
# Removed all filter and trim parameters for the single ended

#trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN)
#out <- filterAndTrim(fnFs, filtFs, maxN=0, maxEE=2, verbose=TRUE)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen= c(250, 250), trimLeft = 20,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# check the quality of the filtered reads - change truncLen and maxEE to find best filtering parameters
plotQualityProfile(filtFs[1:4])
plotQualityProfile(filtRs[1:4])

# learning the error rate - dataset specific
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# quick check of estimated error rates - observed (dots) to follow the expected black line from the model, red line is the normal distrbution. Generally decreasing error with increasing quality (left - right)
plotErrors(errF, nominalQ=TRUE)

# Dereplication 
derep_forward <- derepFastq(filtFs, verbose=TRUE)
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
# now use the core algothrim to calculate the unique sequences with the samples - based on the filtered reads and the error rate calculated before
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# number of ASVs from the unique sequence/haplotypes - change number in brackets to investigate the different samples 

dadaFs[[4]]

#merge the reverse and foward sequence
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[5]])
# generate the ASV table 
seqtabM <- makeSequenceTable(mergers)
filepath <- file.path(path, paste0(Block, "_seqtabM.csv"))
write.csv(seqtabM, file = filepath)
# Inspect distribution of sequence lengths
# maybe want to change the truncLen parameter in the filtering step to get less length variation
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# length distribution of chimeras
table(nchar(getSequences(seqtab.nochim)))

head(seqtab.nochim)

# frequency of non-chimeric sequences - chimera frequency is 1 minus this number 
sum(seqtab.nochim)/sum(seqtab)

# track the number of reads passing each filtering step - good for supplementary tables, comparing parameters 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)


filepath <- file.path(path, paste0(Block, "_seqtab.nochim.csv"))
write.csv(seqtab.nochim, file = filepath)

## Outputting the tables:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
filepath <- file.path(path, paste0(Block, "_ASVs.fasta"))
write.csv(asv_fasta, file = filepath)

asv_seqs2 <- colnames(seqtab)
asv_headers2 <- vector(dim(seqtab)[2], mode="character")

for (i in 1:dim(seqtab)[2]) {
  asv_headers2[i] <- paste(">ASV", i, sep="_")
}

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
filepath <- file.path(path, paste0(Block, "_ASVs_counts.csv"))
write.csv(asv_tab, file = filepath)






