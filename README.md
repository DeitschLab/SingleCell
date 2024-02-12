This repository contains all the scripts and external data files necessary to reproduce the analysis presented in the manuscript "Transcriptional plasticity of virulence genes provides malaria parasites with greater adaptive capacity for avoiding host immunity”. To perform the analysis, download all files to a local folder. A working version of R is required, including the packages loaded in the beginning of each R script. The analysis was performed using Rstudio and R version 4.3.2.

The repository contains:
•	Drop-Seq_matrix_allgenes - A zip file containing all Drop-Seq gene expression matrices.
•	HighSingleA_matrix_allgenes - A zip file containing the HIVE gene expression matrix for the sample “High Single A”.
•	LowManyA_matrix_allgenes - A zip file containing the HIVE gene expression matrix for the sample “Low Many A”.
•	HighSingleB_matrix_allgenes - A zip file containing the HIVE gene expression matrix for the sample “High Single B”.
•	LowManyB_matrix_allgenes - A zip file containing the HIVE gene expression matrix for the sample “Low Many B”.

•	Top50cells – An Excel file with the Drop-Seq Top50 cells in terms of UMI count used for the graphs in Figure 4.
•	HighSingleA_HIVE_Top100 – A csv file with the HIVE Top 100 cells in terms of UMI count for “High Single A” sample, used for the graphs in Figure 6A,C.
•	LowManyA_HIVE_Top100 – A csv file with the HIVE Top 100 cells in terms of UMI count for “Low Many A” sample, used for the graphs in Figure 6B,D.

•	HiveFigures.R – R code for all the analysis in Figure 5 and Supplementary Figure 2.
•	Figure4_6.R – R code for the analysis in Figure 4 and 6.

Supplementary files needed for the codes to work:
•	multicopy.txt – List of all clonally-variant genes used for the analysis in Figure 5E. The list is based on Cortes and Deitsch, 2017 (PMID: 28320828).

