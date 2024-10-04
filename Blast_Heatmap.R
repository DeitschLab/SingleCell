setwd("YourWorkingDirectory")

library(tidyverse)
library(dplyr)

#Blastn Cleanup --> Find Unique Hits to Gene Models
Blast_Hits_ALL <- read.csv("GeneOfInterest-Alignment-HitTable.csv", header = FALSE)

Blast_Hits_ALL

attributes <- "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
attributes

Blast_Hits <- Blast_Hits_ALL[substr(BlastQueryPrefix$V2, 1, 2) == "XM", ]
Blast_Hits

Unique_Blast_Hits <- Blast_Hits[BlastQueryPrefix$V2 != "XM_(YOUR_TRANSCRIPT_ID)", ]
Unique_Blast_Hits

######Now Search the results table

create_search_terms <- function(start, end) {
  map_chr(start:end, ~ paste0("BlastQueryPrefix_window_", ., "_to_", . + 49))
}

search_terms <- create_search_terms(1, #NumberofWindows)

search_and_output <- function(Unique_Blast_Hits, search_terms) {
  
  output <- tibble(search = character(), result = numeric())
  
  for (term in search_terms) {
    match_found <- Unique_Blast_Hits %>% filter(V1 == term)
    
    if (nrow(match_found) > 0) {
      output <- output %>% add_row(search = term, result = match_found$V12)
    } else {
      output <- output %>% add_row(search = term, result = 36.5)
    }
  }
  
  return(output)
}

output <- search_and_output(Unique_Blast_Hits, search_terms)

print(output)

ggplot(output, aes(x = search, y = result)) +
  geom_point()  # or 
  #geom_line()

###Only top value of the output

top_output <- output %>%
  group_by(search) %>%
  summarize(result = max(result))

#Reorder by window position
top_output
 
ordered_top_output <- top_output %>%
  mutate(
    WINDOW = as.numeric(str_extract(search, "\\d+$"))
  ) %>%
  arrange(WINDOW) %>%
  select(-search)

print(ordered_top_output)

ggplot(ordered_top_output, aes(x = WINDOW, y = result)) +
  geom_point()  # or 

ggplot(ordered_top_output, aes(x = WINDOW, y = 1, fill = result)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

