This repository contains all the scripts and external data files necessary to reproduce the analysis presented in the manuscript "Transcriptional plasticity of virulence genes provides malaria parasites with greater adaptive capacity for avoiding host immunity”. To perform the analysis, download all files to a local folder. A working version of R is required, including the packages loaded in the beginning of each R script. The analysis was performed using Rstudio and R version 4.3.2, available for MasOS, Windows, and Linux operating systems. All software and data should download and/or install in less than five minutes. 
________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
The repository contains:

•	Drop-Seq_matrix_allgenes - A zip file containing all Drop-Seq gene expression matrices.

•	3D7HighSingleA_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “High Single A”.

•	3D7HighSingleB_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “High Single B”.

•	3D7LowManyA_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “Low Many A”.

•	3D7LowManyB_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “Low Many B”.

•	IT4HighSingle_HIVE_100 - A zip file containing the HIVE gene expression matrix for the IT4 sample “High Single”.

•	IT4LowMany_HIVE_100 - A zip file containing the HIVE gene expression matrix for the IT4 sample “Low Many”.
________________________________________________________________________________________________________________________________________________________________________________
•	Top50cellsRelative – An Excel file with the Drop-Seq Top50 cells in terms of UMI count used for the graphs in Figure 4.

•	HIVETop100Relative – An Excel file with the HIVE Top100 cells in terms of UMI count used for the graphs in Figure 6.
________________________________________________________________________________________________________________________________________________________________________________
•	HiveFigures.R – R code for all the analysis in Figure 5, Supplementary Figure 4 and Supplementary Figure 7.

•	Figure4_6.R – R code for the analysis in Figure 4 and 6.

•	Figure5C.R – R code for the analysis in Figure 5C.

•	Generate_50bp_Windows.R – R code for generating 50bp windows for the analysis in Figure S5.

•	Blast_Heatmap.R – R code for completing the analysis in Figure S5.

________________________________________________________________________________________________________________________________________________________________________________
Supplementary files needed for the code to work:

•	Multicopy3D7.txt – List of all clonally-variant genes in 3D7 strain, used for the analysis in Figure 5F and Supplementary Table 6. The list is based on Cortes and Deitsch, 2017 (PMID: 28320828).

•	MulticopyIT4.txt – List of all clonally-variant genes in IT4 strain, used for the analysis in Supplementary Figure 4D, E and Supplementary Table 6.

•	Bozdech.csv – Timecourse data from Bozdech et al. 2003 (PMID: 12929205), used for the analysis in Figure 5C.

•	GeneByLocusTag.csv – Converts GeneIDs from previous nomenclature to updated nomenclature, used for the analysis in Figure 5C.

•	Collapsed_Sums.csv – Sums of all HIVE UMI for the 3D7 samples, used for the analysis in Figure 5C.

