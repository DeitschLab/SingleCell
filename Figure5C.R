setwd("YourWorkingDirectory")

#Load Dependencies 
library(ggplot2)
library(reshape2)

#Import Required Data
Bozdech <- read.table("Bozdech.csv", sep = ",", header = TRUE)
GENE_ID_KEY <- read.table("GeneByLocusTag.csv", sep = ",", header = TRUE)

#Import Collapsed Sum of UMI for 3D7 Hive Samples
SUMS <- read.table("Collapsed_Sums.csv", sep = ",", header = TRUE)

#Normalize the data, take UMI and normalize to UMI/Million UMI
row_names <- SUMS[[1]]
SUMS <- SUMS[, -1]
rownames(SUMS) <- row_names
column_sums <- colSums(SUMS)

normalized_df <- SUMS
for (i in seq_along(column_sums)) {
  normalized_df[, i] <- SUMS[, i] * (1000000 / column_sums[i])
}
column_sums_norm <- colSums(normalized_df)

#Double Check Sums=1000000
print(column_sums_norm)

#Reshape Data into Long Format with melt 
normalized_df_melted <- melt(normalized_df, id.vars = normalized_df$rownames ) 

#Merge dataframes based on Name and Input.ID columns
merged <- merge(Bozdech, GENE_ID_KEY, by.x = "Name", by.y = "Input.ID")

# Replace Name with Gene.ID where there is a match
merged$Name <- ifelse(is.na(merged$Gene.ID), merged$Name, merged$Gene.ID)

#Remove Gene.ID Column
merged$Gene.ID <- NULL

#Pearson Correlation Coefficient Calculation
timepoint.cor<-function(CorMe,CorWith,genelist=NULL,method="pearson"){
  if(is.null(genelist)){
    genelist<-intersect(rownames(CorMe),rownames(CorWith))
  }else{
    genelist<-intersect(rownames(CorMe[genelist]),rownames(CorWith[genelist,]))
  }
  cor.value<-cor(CorMe[genelist,],CorWith[genelist,])
  return(cor.value)
}

#Log transform
normalized_df<- log(normalized_df+0.0001)
merged <- unique(merged)
row.names(merged) <- merged$Name
merged$Name <- NULL

####If negative value throws error do the following
negative_value_removal<- "PF3D7_0931900"
merged <- merged[!(row.names(merged) %in% negative_value_removal),]

####
windowRNA<-log(merged)
YOURDATA.COR<-timepoint.cor(CorMe=normalized_df,CorWith=windowRNA)
Best.Match3<-as.data.frame(cbind(names(merged)[apply(YOURDATA.COR,1,which.max)],apply(YOURDATA.COR,1,max)))
colnames(Best.Match3)<-c("Time Point Match","Correlation Coefficient")
View(Best.Match3)

View(YOURDATA.COR)
str(YOURDATA.COR)

# Convert data frame to long format
data_long <- melt(YOURDATA.COR)

data_df <- as.data.frame(YOURDATA.COR)

# Add row names as a column
data_df$Variable <- rownames(YOURDATA.COR)

# Melt the data frame to long format
data_long <- melt(data_df, id.vars = "Variable", variable.name = "TimePoint", value.name = "Value")

#Plot
ggplot(data_long, aes(x = TimePoint, y = Value, color = Variable, group = Variable)) +
  # geom_line() +
  geom_point(size = .8) +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF")) +
  geom_smooth(method = "loess", se = TRUE, level = .95, alpha = .1) +
  labs(x = "Time Point", y = "Correlation", title = "Matrix Plot") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

