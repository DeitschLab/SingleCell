This repository contains all the scripts and external data files necessary to reproduce the analysis presented in the manuscript "Transcriptional plasticity of virulence genes provides malaria parasites with greater adaptive capacity for avoiding host immunity”. To perform the analysis, download all files to a local folder. A working version of R is required, including the packages loaded in the beginning of each R script. The analysis was performed using Rstudio and R version 4.3.2, available for MasOS, Windows, and Linux operating systems. All software and data should install or download in less than five minutes.

The repository contains:
•	Drop-Seq_matrix_allgenes - A zip file containing all Drop-Seq gene expression matrices.
•	3D7HighSingleA_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “High Single A”.
•	3D7HighSingleB_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “High Single B”.
•	3D7LowManyA_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “Low Many A”.
•	3D7LowManyB_HIVE_100 - A zip file containing the HIVE gene expression matrix for the 3D7 sample “Low Many B”.
•	IT4HighSingle_HIVE_100 - A zip file containing the HIVE gene expression matrix for the IT4 sample “High Single”.
•	IT4LowMany_HIVE_100 - A zip file containing the HIVE gene expression matrix for the IT4 sample “Low Many”.

•	Top50cellsRelative – An Excel file with the Drop-Seq Top50 cells in terms of UMI count used for the graphs in Figure 4.
•	HIVETop100Relative – An Excel file with the HIVE Top100 cells in terms of UMI count used for the graphs in Figure 6.

•	HiveFigures.R – R code for all the analysis in Figure 5, Supplementary Figure 4 and Supplementary Figure 7.
•	Figure4_6.R – R code for the analysis in Figure 4 and 6.

Supplementary files needed for the codes to work:
•	Multicopy3D7.txt – List of all clonally-variant genes in 3D7 strain, used for the analysis in Figure 5F and Supplementary Table 6. The list is based on Cortes and Deitsch, 2017 (PMID: 28320828).
•	MulticopyIT4.txt – List of all clonally-variant genes in IT4 strain, used for the analysis in Supplementary Figure 4D, E and Supplementary Table 6.


