#This Script will pull all reads mapping to var genes and plot the ratio of exon1 to exon2 reads for each sample

setwd("#########")

library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)

gff_file <- "example/PlasmoDB-68_PfalciparumIT.gff"   # Path to the GFF file
bam_file <- "example/your.sorted.bam"   # Path to the BAM file

# vargenes
df <- read.csv("vargenes_Strain.csv", header = F) 
gene_ids <- df$geneIDs

# gff
gff_data <- import(gff_file, format = "gff")

# filter for exons and keep only vars
gff_exons <- gff_data %>%
  as.data.frame() %>%
  filter(type == "exon", gene_id %in% gene_ids) %>%
  dplyr::select(seqnames, start, end, strand, gene_id, ID)


# convert to GRanges
exon_gr <- makeGRangesFromDataFrame(gff_exons, keep.extra.columns = TRUE)

# count exon 1 and exon 2 of a single gene
count_exon_reads <- function(bam_file, exons_gr) {
  # sort by position
  exons_sorted <- exons_gr[order(start(exons_gr))]
  
  # vector for counts
  counts <- numeric(length(exons_sorted))
  
  # for loop exon and count reads
  for (i in seq_along(exons_sorted)) {
    exon <- exons_sorted[i]
    param <- ScanBamParam(which = exon)
    counts[i] <- countBam(bam_file, param = param)$records
  }
  
  # return the counts
  return(counts[1:2])
}

# store results
results <- data.frame(GeneID = character(), Strand = character(),
                      Exon1_Counts = integer(), Exon2_Counts = integer(),
                      Exon1_Length = integer(), Exon2_Length = integer(),
                      stringsAsFactors = FALSE)

# for loop for each gene
for (gene in gene_ids) {
  # exons for the current gene
  gene_exons <- exon_gr[exon_gr$gene_id == gene]
  
  # if only 1 exon (pseudogenes that were left in the list)
  if (length(gene_exons) < 2) {
    warning(paste("Skipping", gene, "under 2 exons"))
    next
  }
  
  # sort exons by start coordinate
  gene_exons <- gene_exons[order(start(gene_exons))]
  
  # calculate lengths
  exon1_length <- width(gene_exons[1])
  exon2_length <- width(gene_exons[2])
  
  # strand information 
  strand_info <- as.character(strand(gene_exons)[1])
  
  # count reads
  counts <- count_exon_reads(bam_file, gene_exons)
  
  # write results to dataframe
  results <- rbind(results, data.frame(GeneID = gene,
                                       Strand = strand_info,
                                       Exon1_Counts = counts[1],
                                       Exon2_Counts = counts[2],
                                       Exon1_Length = exon1_length,
                                       Exon2_Length = exon2_length))
}

print(results)
#Input your strain and sample information
write.csv(results, "Strain_Sample.csv", row.names = FALSE)

############
#Align_Exon1_Exon2_Strands
#This block aligns the strands. Swaps length and counts data so that Exon1 is always first. 
IT4_7B <- read.csv("IT4_7B.csv", header = TRUE)

#If IT4_7B$Strand = - then change IT4_7B$Strand to + and switch IT4_7B$Exon1_Counts with IT4_7B$Exon2_Counts and switch IT4_7B$Exon1_Length with IT4_7B$Exon2_Length
Strand_neg <- IT4_7B$Strand == "-"

# Swap Exon1_Counts and Exon2_Counts only for rows where Strand == "-"
temp_counts <- IT4_7B$Exon1_Counts[Strand_neg]
IT4_7B$Exon1_Counts[Strand_neg] <- IT4_7B$Exon2_Counts[Strand_neg]
IT4_7B$Exon2_Counts[Strand_neg] <- temp_counts

# Swap Exon1_Length and Exon2_Length only for rows where Strand == "-"
temp_length <- IT4_7B$Exon1_Length[Strand_neg]
IT4_7B$Exon1_Length[Strand_neg] <- IT4_7B$Exon2_Length[Strand_neg]
IT4_7B$Exon2_Length[Strand_neg] <- temp_length

# Finally, change Strand from "-" to "+"
IT4_7B$Strand[Strand_neg] <- "+"

write.csv(IT4_7B, "IT4_7B_strands_aligned.csv", row.names = FALSE)

###########
#Plot rpkb

IT4_7B_SA <- read.csv("IT4_7B_strands_aligned.csv", header = TRUE)
#Repeat for other samples of interest to load the strand aligned counts

IT4_7B_SA$Exon1_Length_Kb <- IT4_7B_SA$Exon1_Length * .001
IT4_7B_SA$Exon2_Length_Kb <- IT4_7B_SA$Exon2_Length * .001

IT4_7B_SA$Exon1_rpkb <- IT4_7B_SA$Exon1_Counts / IT4_7B_SA$Exon1_Length_Kb
IT4_7B_SA$Exon2_rpkb <- IT4_7B_SA$Exon2_Counts / IT4_7B_SA$Exon2_Length_Kb

IT4_7B_SA$ratio <- IT4_7B_SA$Exon1_rpkb / IT4_7B_SA$Exon2_rpkb

IT4_7B_SA$ratio

ggplot(IT4_7B_SA, aes(x = seq_along(ratio), y = ratio)) +
  geom_point(color = "blue", size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
  labs(
    x = "Var Gene",
    y = "Ratio Exon1:Exon2",
    title = "Exon1 v Exon2 Alignments"
  ) +
  theme_minimal()

SA_Dataframes <- list(A03_SA, A11_SA, B09_SA, A04_SA, IT4_7D_SA, IT4_7B_SA)

Totals <- data.frame(
  dataframe = paste0("df", seq_along(SA_Dataframes)), # Name of each dataframe
  Exon1_rpkb = sapply(SA_Dataframes, function(df) sum(df$Exon1_rpkb, na.rm = TRUE)),
  Exon2_rpkb = sapply(SA_Dataframes, function(df) sum(df$Exon2_rpkb, na.rm = TRUE))
)

data_long <- Totals %>%
  tidyr::pivot_longer(cols = c(Exon1_rpkb, Exon2_rpkb), 
                      names_to = "Exon", 
                      values_to = "RPKB")

# Calculate proportions within each sample
data_long <- data_long %>%
  group_by(dataframe) %>%
  mutate(Proportion = RPKB / sum(RPKB)) %>%
  ungroup()

# Plot option 1: Stacked bar chart (proportions)
ggplot(data_long, aes(x = dataframe, y = Proportion, fill = Exon)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Exon1 and Exon2 RPKB per Sample",
       x = "Sample", y = "Proportion") +
  theme_bw() +
  scale_fill_manual(values=c("red", "blue"),
                    labels = c("Exon 1", "Exon 2")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() + 
  scale_x_discrete(limits = rev) 
