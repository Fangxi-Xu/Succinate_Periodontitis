library(readxl)
data_corr <- read_excel("data_corr.xlsx")
#correlation plot
library(ggplot2)
library(ggpubr)
# create plot folder
path <- "./plot"
dir.create(path)

png(file.path(path, "MSD_metabolomics.png"), width=30,height=32,units="cm",res=1200, bg = "transparent")
ggscatter(data_corr, x = "value", y = "Succinate_Acid", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = FALSE, cor.method = "spearman", ylab = "Succinic Acid (Relative Intensity)")  + geom_point(aes(fill=Group,color=Group), size=1) + 
  facet_wrap(vars(variable), scales="free", nrow = 8)+
  scale_fill_manual(values=c("white", "white"), guide = FALSE)+
  scale_color_manual(values = c("Normal/mild"="#000000", "Severe"="#FF0000"))+
  theme(
    axis.text.x=element_text(size=6, face = "bold", color="black"),
    axis.text.y=element_text(size=8, face = "bold", color="black"),
    axis.title.y=element_text(size=12, face = "bold", color="black"),
    strip.text = element_text(size = 10, face = "bold", color="black"))+
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  theme(
    legend.position="top",
    legend.text=element_text(size=10, face = "bold", color="black"))

dev.off()

