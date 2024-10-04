setwd("YourWorkingDirectory")

library(Biostrings)
library(seqinr)

split_fasta_with_sliding_window <- function(input_fasta, output_fasta, window_size=50, step_size=1) {
  # Load input fasta
  seq <- readDNAStringSet(input_fasta)
  seq_id <- names(seq)
  seq <- as.character(seq[[1]])
  
  # seq and name lists
  window_sequences <- list()
  window_names <- list()
  
  # Generate windows list
  for (i in seq(1, nchar(seq) - window_size + 1, by = step_size)) {
    window_seq <- substr(seq, i, i + window_size - 1)
    window_id <- paste(seq_id, "window", i, "to", i + window_size - 1, sep = "_")
    
    window_sequences <- c(window_sequences, window_seq)
    window_names <- c(window_names, window_id)
  }
  
  # fasta file
  write.fasta(sequences = window_sequences, names = window_names, file.out = output_fasta)
}

# load files and run
input_fasta <- "GeneOfInterest.fasta"
output_fasta <- "Output_50bp_Windows.fasta"
split_fasta_with_sliding_window(input_fasta, output_fasta, window_size = 50, step_size = 1)

#Take this output and run on NCBI blast: command line tool or web browser version