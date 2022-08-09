#import top40 genus table
#read corrected genus table(pick top 20 based on sum)
library("readxl")
top40.genus<- read_xlsx("top40genus.xlsx")
head(top40.genus)
#transpose the df and keep the sampleIDs as the header
top40.genus.t <- as.data.frame(t(top40.genus[,-1]))
colnames(top40.genus.t) <- top40.genus$SampleID
head(top40.genus.t)
#read in annotations
ann<- read.table("annotation.txt", header = T)
#reorder martrix based on annotation
top40_ordered <- top40.genus.t[rownames(ann), ]
#Z score transformation
top40_ordered_zscore <- scale(top40_ordered)
summary(top40_ordered_zscore)
#specify colors for annotation
ann_colors = list(Group = c("Normal/mild"="#000000", "Severe"="#FF0000"))

#sample based heatmap
library(pheatmap)
library(RColorBrewer)
library(viridis)
#png(file.path(path, "Top40genus(sample_based).png"), width=15,height=10.5,units="cm",res=1200)
pdf(file.path(path, "Top40genus_heatmap_no_name.pdf"), width = 6, height =5)
pheatmap(t(top40_ordered_zscore), scale = "row", color = colorRampPalette(c("blue", "white", "red"))(50), 
         cellheight = 7,
         cellwidth = 6,
         fontsize_row=8, 
         show_colnames = F,show_rownames = T, annotation_col = ann, annotation_colors = ann_colors,
         angle_row = "45",  cluster_cols = F)
dev.off()
################################################################################

#end of the script
