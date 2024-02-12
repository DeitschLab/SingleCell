#This code generates the graphs in Figures 4 and 6

setwd("") #Set to the folder where all scripts and files are found

library(tidyverse)
library(readxl)
library(xlsx)
library(RColorBrewer)
nb.cols <- 70
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

#Read DropSeq Enrichment matrices (Filtered for the top 50 cells in terms of total UMI count)
LowManyA<-read_excel('Top50cells.xlsx')
LowManyA_Enriched<-read_excel('Top50cells.xlsx', sheet=2)
HighSingleA<-read_excel('Top50cells.xlsx', sheet=3)
HighSingleA_Enriched<-read_excel('Top50cells.xlsx', sheet=4)

#Read HIVE matrices (Filtered for the top 100 cells in terms of total UMI count)
HighSingleA_HIVE<-read.csv("HighSingleA_HIVE.csv")
LowManyA_HIVE<-read.csv("LowManyA_HIVE.csv")

#Reorganization of the table (shown for one sample as an example)
HSAHIVEReshaped<-HighSingleA_HIVE %>% 
  pivot_longer(
    cols=!Cells,
    names_to="vars",
    values_to="count"
  )

#to distribute vars by amount of reads
HSAHIVEOrdered <- HSAHIVEReshaped %>% 
  group_by(Cells) %>% 
  arrange(desc(count),Cells) %>% 
  mutate(
    vars = factor(vars, levels = unique(vars))
  )

#Graphs showing expression of different var genes as number of UMIs
#Figure 4A,B,C,D and Figure 6A,B
ggplot(HSAHIVEOrdered, aes(x=Cells, y=count, fill=vars))+
  geom_bar (colour="black", linewidth=0.15,stat="identity", show.legend = FALSE)+
  theme(axis.text.x=element_blank(),panel.background = element_blank(), axis.line=element_line(), plot.title = element_text(hjust = 0.5))+
  ggtitle("High Single A HIVE")+
  scale_fill_manual(values = mycolors) +
  ylim(0,100)+ labs(y = "Number of UMI")


#Graphs showing expression of different var genes as percentage of total var UMIs
#Figure 4E,F and Figure 6C,D
ggplot(HSAHIVEOrdered, aes(x=Cells, y=count, fill=vars))+
  geom_bar (position = "fill",colour="black", linewidth=0.15,stat="identity", show.legend = FALSE)+
  theme(axis.text.x=element_blank(),panel.background = element_blank(), axis.line=element_line(), plot.title = element_text(hjust = 0.5))+
  ggtitle("High Single A HIVE") +
  scale_fill_manual(values = mycolors)+ 
  labs(y = "Percentage of UMI")
