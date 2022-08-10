#load phyloseq object
phylo_clean <- readRDS("phylo1.rds")
##################################
library(MicrobiotaProcess)
diffres <- diff_analysis(obj=phylo_clean, #a phyloseq object
                         classgroup="Group",#the factor name in sampledata
                         subclass=NULL,#no subclass compare
                         standard_method=NULL,#the relative abundance of taxonomy will be used
                         mlfun="lda",
                         firstcomfun = "kruskal.test",
                         padjust=NULL,
                         filtermod="pvalue",
                         firstalpha=0.05,
                         strictmod=TRUE,
                         clwilc=TRUE, 
                         subclwilc=FALSE,
                         lda=2)
diffres

plotes <- ggeffectsize(obj=diffres,
                       pointsize=2.5,
                       linecolor="grey",
                       linewidth=1,
                       lineheight=0.4,
                       removeUnknown=TRUE) + scale_color_manual(values=c("Normal/mild"="#000000", "Severe"="#FF0000"))+
          theme(axis.text.x=element_text(size=12, color = "black"))+
          theme(axis.title.x=element_text(size=12, color = "black"))+
          theme(axis.title.y=element_text(size=12, color = "black"))+
          theme(axis.text.y=element_text(size=12, color = "black"))+
          theme(legend.title=element_text(size=12, color = "black"), 
                legend.text=element_text(size=12,  color = "black"),
    axis.text.y=element_text(size=12, color="black"),
    strip.text = element_blank(),
    legend.position="bottom")

pdf(file.path(path, "Lefse_r.pdf"), width=7,height=5)
#png(file.path(path, "Lefse_r.png"), width=12,height=8,units="cm",res=1200, bg = "transparent")
plotes
dev.off()


diffcladeplot <- ggdiffclade(obj=diffres,#diffAnalysisClass Object2
                             alpha=0.4, size=0.8, 
                             skpointsize=0.8,
                             taxlevel=5,
                             cladetext=4,
                             settheme=FALSE, 
                             setColors=FALSE,
                             removeUnknown=TRUE) +
  scale_fill_manual(values=c("Normal/mild"="#000000", "Severe"="#FF0000"))+
  guides(color = guide_legend(keywidth = 0.2,
                              keyheight = 0.6,
                              order = 3, 
                              ncol=1)) + 
  theme(panel.background=element_rect(fill=NA),
        legend.position="right",
        plot.margin=margin(0,0,0,0),
        legend.spacing.y = unit(0.02, "cm"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.box.spacing=unit(0.02,"cm")
        )
pdf(file.path(path, "Lefse_cladogram_r.pdf"), width=8,height=6)
diffcladeplot
dev.off()
