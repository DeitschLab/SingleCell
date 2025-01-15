library(ggplot2)
library(ggrepel)

setwd("/home/joe/raid5/DropSeqArchive/Joe/Blant-Altman")

Data <- read.csv("varProportion.csv", header = TRUE)
row.names(Data) <- Data$var

#select the comparison
comparison <- c("qPCR_SAMPLE", "HIVE_SAMPLE") 
Data <- Data[, comparison] 

#calculate average
Data$avg <- rowMeans(Data) 

#calculate difference
Data$diff <- Data[ , 1] - Data[ , 2]

# calculate mean difference
mean_diff <- mean(Data$diff)

# 90% confidence interval 
lower_limit <- mean_diff - 1.96*sd(Data$diff)
upper_limit <- mean_diff + 1.96*sd(Data$diff)

# PlotBA 
ggplot(Data, aes(x = avg, y = difference)) +
  geom_point(size=3) +
  geom_hline(yintercept = mean_difference, color= "black", lwd=1.5) +
  geom_hline(yintercept = lower_limit, color = "red", lwd=1.5, linetype = "dashed") +
  geom_hline(yintercept = upper_limit, color = "red", lwd=1.5, linetype = "dashed") +
  geom_text_repel(aes(label = rownames(Data))) +
  coord_cartesian(ylim = c(-.025, .025), xlim = c(0, .025)) +
  labs(title = paste("Plot of", names(Data)[1], "vs.", names(Data)[2]))+
  theme_bw()
