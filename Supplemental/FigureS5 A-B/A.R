#import top20 genus table
#read corrected genus table(pick top 20 based on sum)
library("readxl")
library(tidyverse)
library(ggplot2)
top20.genus.ligature<- read_xlsx("genus_ligature_top20_barplot.xlsx")
head(top20.genus.ligature)
#transpose the df and keep the sampleIDs as the header
top20.genus.ligature.t <- as.data.frame(t(top20.genus.ligature[,-1]))
colnames(top20.genus.ligature.t) <- top20.genus.ligature$SampleID
head(top20.genus.ligature.t)
top20.genus.ligature.t <- (top20.genus.ligature.t)*100
library(tibble)
top20.genus.ligature.t <- tibble::rownames_to_column(top20.genus.ligature.t, "SampleID")
#merge with metadata
ann<- read.table("annotation_ligature.txt", header = T)
ann <- tibble::rownames_to_column(ann, "SampleID")
top20_full <- merge(ann, top20.genus.ligature.t , by = "SampleID")
#transform data for plotting barplot
top20_full_melt <- reshape2::melt(top20_full)
#genus_reshaped <- reshape2::dcast(top20_full, Sample+Treatment ~ variable, value.var='value') #transform data
#calculate stats
genus_summary <-top20_full_melt %>%
  group_by(Treatment,variable) %>%  # the grouping variable
  summarise(mean_Abun = mean(value),# calculates the mean of each group
            sd_Abun = sd(value), # calculates the standard deviation of each group
            n_Abun = n(),  # calculates the sample size per group
            SE_Abun = sd(value)/sqrt(n()))# calculates the standard error of each group



ann_colors = list(Treatment = c("7a"="#9629cc", "Vehicle"="#FF8C00"))
pdf(file.path(path, "genus(barplot).pdf"), width = 6, height = 5)
ggplot(genus_summary, aes(fill=Treatment, y=mean_Abun, x=variable)) + 
  geom_bar(position="dodge", stat="identity", width = 0.8)+
  geom_errorbar(aes(ymin=mean_Abun,ymax=mean_Abun+SE_Abun),width = .2, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#9629cc", "#FF8C00"))+
  ylab("Mean Relative Abundance(%)") +
  theme_classic()+
  theme(axis.text.x=element_text(size=8, angle = 45, color = "black",vjust=0.5, hjust=1))+
  theme(axis.text.y=element_text(size=11, color = "black"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=11, color = "black"))+
  theme(legend.text=element_text(size=10, color = "black"))+
  theme(legend.title=element_text(size=10, color = "black"))+
  theme(legend.position="top")
dev.off()

