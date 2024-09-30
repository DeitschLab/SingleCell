#This code generates the graphs in Figures 4 and 6

setwd("/Users/Fra/Desktop") #Set to the folder where all scripts and files are found

install.packages("tidyverse")
install.packages("readxl")
install.packages("xlsx")
install.packages("RColorBrewer")

library(tidyverse)
library(readxl)
library(xlsx)
library(RColorBrewer)
nb.cols <- 70
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

#Read DropSeq Enrichment matrices (Filtered for the top 50 cells in terms of total UMI count)
LowManyA<-read_excel('Top50cellsRelative.xlsx')
LowManyA_Enriched<-read_excel('Top50cellsRelative.xlsx', sheet=2)
HighSingleA<-read_excel('Top50cellsRelative.xlsx', sheet=3)
HighSingleA_Enriched<-read_excel('Top50cellsRelative.xlsx', sheet=4)

#Read HIVE matrices (Filtered for the top 100 cells in terms of total UMI count)
HighSingleA_HIVE <- read_excel('HiveTop100Relative.xlsx')
LowManyA_HIVE<-read_excel('HiveTop100Relative.xlsx', sheet=3)

#Reorganization of the table (shown for one sample as an example)
LowManyAReshaped<-LowManyA_HIVE %>% 
  pivot_longer(
    cols=!Cells,
    names_to="vars",
    values_to="count"
  )

#to distribute vars by amount of reads
LowManyAOrdered <- LowManyAReshaped %>% 
  group_by(Cells) %>% 
  arrange(desc(count),Cells) %>% 
  mutate(
    vars = factor(vars, levels = unique(vars))
  )

#Graphs showing expression of different var genes as relative number of UMIs
ggplot(LowManyAOrdered, aes(x=Cells, y=count, fill=vars))+
  geom_bar (colour="black", linewidth=0.15,stat="identity", show.legend = FALSE)+
  theme(axis.text.x=element_blank(),panel.background = element_blank(), axis.line=element_line(), plot.title = element_text(hjust = 0.5))+
  ggtitle("Low Many HIVE")+
  scale_fill_manual(values = mycolors) +
  ylim(0,0.08)+ labs(y = "Relative UMI")

