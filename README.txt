This R project is data and scripts for paper - Defensive symbiont genotype distributions are linked to parasitoid attack networks

Script for analysis:

    Run all the scripts in order, to get all analyses and results in this paper. 

    Script_0_package_install_load.R: 
        Installs and loads the necessary R packages required for all subsequent scripts.


    Script_1_Extract_from_Table_S2.R
        Extracts information from:
            Table_S2.csv (located in /Rawdata/)    Data_4_Aphid_sequences.fas
        Generates 10 files for further analyses:
            Aphid relatedness distance:
              Aphid_phylogenetic_relatedness_16species.csv
              Aphid_phylogenetic_relatedness_22species.csv
              Aphid_phylogenetic_relatedness_31species.csv
            Hamiltonella - Aphid Matrices:
              Hamiltonella_Aphid_matrix_16species.csv
              Hamiltonella_Aphid_matrix_22species.csv
              Hamiltonella_Aphid_matrix_31species.csv
            Parasitoid - Aphid Matrices:
              Para_Aphid_matrix_16species.csv
              Para_Aphid_matrix_22species.csv
            Plant - Aphid Matrices:
              Plant_Aphid_matrix_16species.csv
              Plant_Aphid_matrix_31species.csv


    Script_2_BarPlot_Fig1.R: 
        Generates Fig. 1 using: 
            Data_11_Parasitoid_genus_aphid_22_species_Fig.1.csv  Data_12_Plant_genus_aphid_31_species_Fig.1.csv


    Script_3_MMRR_analysis.R
        Uses 10 distance/distribution matrices of Parasitoid, Plant, Hamiltonella, and Aphid relationships to compute matrix correlations with the MMRR (Multiple Matrix Regression with Randomization) model.


    Script_4_Species_linkage_Fig2_FigS4.R
        Uses 7 distribution matrices of Parasitoid, Plant, and Hamiltonella relationships with Aphids (excluding Aphid genetic distance matrices) to generate the species linkage diagrams for Fig. 2 and Fig. S4.


    Script_5_MMRR_Fig3_FigS5.R
        Plots the MMRR correlations using the 10 distance/distribution matrices for Parasitoid, Plant, Hamiltonella, and Aphid relationships (Fig. 3 & Fig. S5). 

    Script_6_parasitoid_specialization_Fig4.R
        Computes specialization levels using the H2 Index for:
            Parasitoid-Aphid
            Aphid-Hamiltonella
            Parasitoid-Hamiltonella relationships 
        Generates Figure 4.

    Script_7_Ecologicial_indices&plots_Table1_FigS3.R
        Uses 7 distribution matrices of Parasitoid, Plant, and Hamiltonella relationships with Aphids (excluding Aphid genetic distance matrices) to:
            Calculate ecological indices (Richness, Shannon Index, Simpson Index).
            Use linear models to analyze the relationships between Parasitoid/Plant-Aphid and Aphid-Hamiltonella communities.   
        Outputs results for Table 1 and Figure S3.


    Script_8_Bubble_plot_FigS2.R
        This script using 
            Hamiltonella_Aphid_matrix_31species.csv 
        and two phylogeny trees 
            Data_8_Aphid_species_phylogeny.txt & 
            Data_9_Hamiltonella_phylogeny.txt
        To create a Cophylogeny tree of Hamiltonella strains and aphid species that we identified in this study and previously known strains.
        (Fig. S2)

    Script_9_DADA2_pipeline_(Optional).R 
        Implements the DADA2 pipeline to generate: 
            Data_1_Deep_sequencing_data_block_1.xlsx & 
            Data_2_Deep_sequencing_data_block_2.xlsx
        Requirements:
            Download raw sequencing data from GenBank (BioProject PRJNA1139364) before running this script.
        A 10-sample test dataset is provided in /Rawdata/Illumina_test/.

        Raw FASTQ files Genbank inofrmation:
            The 370 samples of Block1 are SRX25431825 to SRX25432194
            The 738 samples of Block2 are SRX25431012 to SRX25431744     



Raw data and inofrmation:

    Data_1_Deep_sequencing_data_block_1.xlsx & Data_2_Deep_sequencing_data_block_2.xlsx 
        Sample: The sample name from Illumina sequencing, which can later be found on GenBank (BioProject PRJNA1139364).
        Aphid: The detected aphid species for each sample.
        Parasitoid: The detected parasitoid species.
        The remaining columns contain ASV (Amplicon Sequence Variant) data derived from high-throughput barcoding sequencing.

    Data_3_Parasitoid sequences.fas & 
    Data_4_Aphid sequences.fas & 
    Data_5_Hamiltonella sequences.fas,
        FASTA files for:
            Parasitoid species
            Aphid species
            Hamiltonella strains
        These sequences were used for phylogenetic reconstruction and subsequent analyses.

    Data_6_Parasitoid_Aphid_pooled_table.xlsx 
        Contains the OTU (Operational Taxonomic Unit) pooling results for Parasitoid and Aphid species based on 99% (4 base pair) sequence similarity:
            Parasitoid pooling together group & Aphid pooling together group: Original names from Illumina sequencing.
            Original name: Representative sequences for each group.
            Pooled species name: Final species names used in all analyses.
            Pooled sequences: Final representative sequences used in all analyses.

    Data_7_Parasitoid_species_phylogeny.txt & Data_8_Aphid_species_phylogeny.txt & Data_9_Hamiltonella_phylogeny.txt 
        Phylogenetic trees for:
            Parasitoid species
            Aphid species
            Hamiltonella strains
        These trees were generated using the PhyML tool on the ATGC Montpellier platform, original fasta file was Data_3, 4 & 5.

    Data_10_aphid_host_info.csv 
        Provides aphid host information for generating Figure 2 and Figure S4:
        Aphid: Names of aphid species included in this study.
        Host: Host categories:
            1: Herb aphids
            2: Grass aphids
            3: Tree aphids
        Host_category: Detailed descriptions of host categories.

    Data_11_Parasitoid_genus_aphid_22_species_Fig.1.csv & Data_12_Plant_genus_aphid_31_species_Fig.1.csv 
        Reduced matrices used for Figure 1, created by merging data from:
            OTUs associated with the same parasitoid species.
            Plant species belonging to the same genus.

    7 Distribution Matrices (from Script_1_Extract_from_Table_S2.R)
        Hamiltonella-Aphid Matrices:
           Hamiltonella_Aphid_matrix_16species.csv
            Hamiltonella_Aphid_matrix_22species.csv
            Hamiltonella_Aphid_matrix_31species.csv
        Parasitoid-Aphid Matrices:
            Para_Aphid_matrix_16species.csv
            Para_Aphid_matrix_22species.csv
        Plant-Aphid Matrices:
            Plant_Aphid_matrix_16species.csv
            Plant_Aphid_matrix_31species.csv
        Structure:
            Rows represent aphid species.
            Columns represent Hamiltonella strains, parasitoid species, or host plant species linked to each aphid species.
        Usage: These matrices were used to generate Figures 2, 3, 4, Figures S2, S3, S4, S5, and for MMRR tests, ecological indices, and H2 index calculations.

    3 Aphid phylogenetic relatedness matrices derived from Script_1_Extract_from_Table_S2.R
        Aphid_phylogenetic_relatedness_16species.csv
        Aphid_phylogenetic_relatedness_22species.csv
        Aphid_phylogenetic_relatedness_31species.csv
        These matrices provide phylogenetic distance information, showing pairwise genetic distances between the 31 aphid species included in this study..