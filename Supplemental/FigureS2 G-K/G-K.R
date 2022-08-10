#load phyloseq object
phylo_geno <- readRDS("phylo2.rds")
################################################################################
#select WT 
phylo_geno_WT<-subset_samples(phylo_geno, Group== "WT"|Group=="WT_Perio")
phylo_geno_WT
phylo_geno_WT_clean<-prune_taxa(taxa_sums(phylo_geno_WT) > 0, phylo_geno_WT)
phylo_geno_WT_clean
#Beta Diveristy
################################################################################
#method:Phylogenetic beta-diversity metrics: Weighted Unifrac
#consider the abundances of different taxa
#PERMANOVA
library(vegan)
metadf <- data.frame(sample_data(phylo_geno_WT_clean))

unifrac.dist <- UniFrac(phylo_geno_WT_clean, 
                        weighted = T, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ Group, data = metadf)
permanova#check significance
swab.wt.uni <- ordinate(phylo_geno_WT_clean, "PCoA", "unifrac", weighted=T)

wt.unifrac <- plot_ordination(phylo_geno_WT_clean, 
                              swab.wt.uni, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_color_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05"))+
  scale_fill_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05"))+
  ggtitle("Weighted UniFrac") + geom_point(size = 4)+
  theme_classic() + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.text.x=element_text(size=14, face="bold", color = "black"))+
  theme(axis.title.x=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.title.y=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.text.y=element_text(size=14, face="bold", color = "black"))+
  theme(legend.title=element_text(size=14, face="bold", color = "black"), 
        legend.text=element_text(size=14, face="bold", color = "black"))+
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  annotate("text", x = -0.06, y = -0.04, size = 5,
           label = "paste(italic(p), \" = 0.07\")", parse = TRUE)
#png(file.path(path, "WT_Weighted_Uni.png"), width=10,height=10,units="cm",res=1200,bg ="transparent")
pdf(file.path(path, "WT_WU.pdf"), width = 4, height =4)
wt.unifrac
dev.off()

################################################################################
#method:bray
#PERMANOVA
dist <- phyloseq::distance(phylo_geno_WT_clean, method = "bray")
perma <- adonis(dist~Group, data = as(sample_data(phylo_geno_WT_clean), "data.frame"), permutations = 999)
#dist <- permutest(betadisper(dist, sample_data(phylo_geno_WT)$Group), permutations = 999)
perma
swab.bray <- ordinate(phylo_geno_WT_clean, "PCoA", "bray")
bray <- plot_ordination(phylo_geno_WT_clean, 
                        swab.bray, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_color_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05"))+
  scale_fill_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05"))+
  ggtitle("Bray Curits") + geom_point(size = 4)+
  theme_classic() + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.text.x=element_text(size=14, face="bold", color = "black"))+
  theme(axis.title.x=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.title.y=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.text.y=element_text(size=14, face="bold", color = "black"))+
  theme(legend.title=element_text(size=14, face="bold", color = "black"), 
        legend.text=element_text(size=14, face="bold", color = "black"))+
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))+
  annotate("text", x = 0.4, y = 0.75, size = 5,
           label = "paste(italic(p), \" = 0.001\")", parse = TRUE)
#png(file.path(path, "WT_Bray_Curits.png"), width=10,height=10,units="cm",res=1200,bg ="transparent")
pdf(file.path(path, "WT_BC.pdf"), width = 4, height =4)
bray
dev.off()
################################################################################
#select KO
phylo_geno_KO<-subset_samples(phylo_geno, Group== "KO"|Group=="KO_Perio")
phylo_geno_KO
phylo_geno_KO_clean<-prune_taxa(taxa_sums(phylo_geno_KO) > 0, phylo_geno_KO)
phylo_geno_KO_clean
#Beta Diveristy
#method:Phylogenetic beta-diversity metrics: Weighted Unifrac
#consider the abundances of different taxa
#PERMANOVA
library(vegan)
metadf <- data.frame(sample_data(phylo_geno_KO_clean))
unifrac.dist <- UniFrac(phylo_geno_KO_clean, 
                        weighted = T, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ Group, data = metadf)
permanova#check significance
swab.wt.uni <- ordinate(phylo_geno_KO_clean, "PCoA", "unifrac", weighted=T)

wt.unifrac <- plot_ordination(phylo_geno_KO_clean, 
                              swab.wt.uni, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_color_manual(values = c("KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  scale_fill_manual(values = c("KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  ggtitle("Weighted UniFrac") + geom_point(size = 4)+
  theme_classic() + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.text.x=element_text(size=14, face="bold", color = "black"))+
  theme(axis.title.x=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.title.y=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.text.y=element_text(size=14, face="bold", color = "black"))+
  theme(legend.title=element_text(size=14, face="bold", color = "black"), 
        legend.text=element_text(size=14, face="bold", color = "black"))+
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))
#pdf(file.path(path, "swab_KO_wt.unifrac_PCOA.pdf"), height = 6, width = 6)
#wt.unifrac
#dev.off()
#png(file.path(path, "KO_Weighted_Uni.png"), width=10,height=10,units="cm",res=1200, bg ="transparent")
pdf(file.path(path, "KO_WU.pdf"), width = 4, height =4)
wt.unifrac
dev.off()
################################################################################
#method:bray
#PERMANOVA
dist <- phyloseq::distance(phylo_geno_KO_clean, method = "bray")
perma <- adonis(dist~Group, data = as(sample_data(phylo_geno_KO_clean), "data.frame"), permutations = 999)
#dist <- permutest(betadisper(dist, sample_data(phylo_geno_KO)$Group), permutations = 999)
perma
swab.bray <- ordinate(phylo_geno_KO_clean, "PCoA", "bray")
bray <- plot_ordination(phylo_geno_KO_clean, 
                        swab.bray, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_color_manual(values = c("KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  scale_fill_manual(values = c("KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  ggtitle("Bray Curits") + geom_point(size = 4)+
  theme_classic() + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.text.x=element_text(size=14, face="bold", color = "black"))+
  theme(axis.title.x=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.title.y=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.text.y=element_text(size=14, face="bold", color = "black"))+
  theme(legend.title=element_text(size=14, face="bold", color = "black"), 
        legend.text=element_text(size=14, face="bold", color = "black"))+
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))
#pdf(file.path(path, "swab_KO_bray_PCOA.pdf"), height = 6, width = 6)
#bray
#dev.off()
#png(file.path(path, "KO_Bray_Curits.png"), width=10,height=10,units="cm",res=1200,bg ="transparent")
pdf(file.path(path, "KO_BC.pdf"), width = 4, height =4)
bray
dev.off()
################################################################################
#subset perio
phylo_geno_perio<-subset_samples(phylo_geno, Group == "WT_Perio"|Group=="KO_Perio")
phylo_geno_perio
phylo_geno_perio_clean <-prune_taxa(taxa_sums(phylo_geno_perio) > 0, phylo_geno_perio)
phylo_geno_perio_clean
################################################################################
#method:Phylogenetic beta-diversity metrics: Weighted Unifrac
#consider the abundances of different taxa
#PERMANOVA
library(vegan)
metadf <- data.frame(sample_data(phylo_geno_perio_clean))

unifrac.dist <- UniFrac(phylo_geno_perio_clean, 
                        weighted = T, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ Group, data = metadf)
permanova#check significance
swab.wt.uni <- ordinate(phylo_geno_perio_clean, "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(phylo_geno_perio_clean, 
                              swab.wt.uni, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_fill_manual(values = c("WT_Perio"="#c71b05", "KO_Perio"="#ed9e0a"))+
  scale_color_manual(values = c("WT_Perio"="#c71b05", "KO_Perio"="#ed9e0a"))+
  ggtitle("Weighted UniFrac") + geom_point(size = 4)+
  theme_classic() + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.text.x=element_text(size=14, face="bold", color = "black"))+
  theme(axis.title.x=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.title.y=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.text.y=element_text(size=14, face="bold", color = "black"))+
  theme(legend.title=element_text(size=14, face="bold", color = "black"), 
        legend.text=element_text(size=14, face="bold", color = "black"))+
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))
 # annotate("text", x = 0.15, y = -0.1, size = 5,
          # label = "paste(italic(p), \" = 0.002\")", parse = TRUE)
#png(file.path(path, "swab_perio_Weighted_Uni.png"), width=10,height=10,units="cm",res=1200,bg ="transparent")
pdf(file.path(path, "swab_perio_Weighted_Uni.pdf"), width=4,height=4)
wt.unifrac
dev.off()
################################################################################
#method:bray
#PERMANOVA
dist <- phyloseq::distance(phylo_geno_perio_clean, method = "bray")
perma <- adonis(dist~Group, data = as(sample_data(phylo_geno_perio_clean), "data.frame"), permutations = 999)
#dist <- permutest(betadisper(dist, sample_data(phylo_geno_perio)$Group), permutations = 999)
perma
#plot
swab.bray <- ordinate(phylo_geno_perio, "PCoA", "bray")
bray <- plot_ordination(phylo_geno_perio, 
                        swab.bray, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_fill_manual(values = c("WT_Perio"="#c71b05", "KO_Perio"="#ed9e0a"))+
  scale_color_manual(values = c("WT_Perio"="#c71b05", "KO_Perio"="#ed9e0a"))+
  ggtitle("Bray Curtis") + geom_point(size = 4)+
  theme_classic() + 
  theme(plot.title = element_text(size = 16, face = "bold", color = "black"))+
  theme(axis.text.x=element_text(size=14, face="bold", color = "black"))+
  theme(axis.title.x=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.title.y=element_text(size=14, face = "bold", color = "black"))+
  theme(axis.text.y=element_text(size=14, face="bold", color = "black"))+
  theme(legend.title=element_text(size=14, face="bold", color = "black"), 
        legend.text=element_text(size=14, face="bold", color = "black"))+
  theme(legend.position="bottom")+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))
  #annotate("text", x = 0.3, y = -1, size = 5,
         #  label = "paste(italic(p), \" = 0.001\")", parse = TRUE)
#png(file.path(path, "swab_perio_Bray_Curits.png"), width=10,height=10,units="cm",res=1200, bg = "transparent")
pdf(file.path(path, "swab_perio_Bray_Curits.pdf"), width=4,height=4)
bray
dev.off()

################################################################################
#phylum level stacked bar-plot
phylo_relative = transform_sample_counts(phylo_geno, function(x) {(x/sum(x))*100} )
glom_phylum <- tax_glom(phylo_relative, taxrank = 'Phylum')
data_phylum <- psmelt(glom_phylum) # create dataframe from phyloseq object
data_phylum$phylum <- as.character(data_phylum$Phylum) #convert to character
#simple way to rename phyla with < 1% abundance
data_phylum$Phylum[data_phylum$Phylum == "p__"] <- "Unclassified"
#calculate statistics
phylum_swab_summary <-data_phylum %>%
  group_by(Group,Phylum) %>%  # the grouping variable
  summarise(mean_Abun = mean(Abundance),# calculates the mean of each group
            sd_Abun = sd(Abundance), # calculates the standard deviation of each group
            n_Abun = n(),  # calculates the sample size per group
            SE_Abun = sd(Abundance)/sqrt(n()))# calculates the standard error of each group

#png(file.path(path, "phylum(stacked).png"), width=14,height=10,units="cm",res=1200,bg ="transparent")
pdf(file.path(path, "phylum(stacked).pdf"), width = 6, height =4)
ggplot(phylum_swab_summary, aes(fill=Phylum, y=mean_Abun, x=Group)) + 
  scale_x_discrete(limits=c("WT", "WT_Perio", "KO", "KO_Perio"))+
  scale_fill_manual(values = c( "#972C8DFF","#FA6B09FF",  "#149BEDFF" ,"#EE0011FF","#9A703EFF",
                               "#15983DFF", "#A1C720FF","#16A08CFF", "#FEC10BFF","#0C5BB0FF",
                               "#BA6222FF","#808080","#EC579AFF"))+
  ylab("Mean Relative Abundance(%)") +
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, angle = 0, color = "black"))+
  theme(axis.text.y=element_text(size=12, color = "black"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16, color = "black"))+
  theme(legend.text=element_text(size=14, color = "black"))+
  theme(legend.title=element_text(size=16, color = "black"))+
  theme(legend.position="right")
#  theme(plot.background = element_rect(fill = "transparent",colour = NA))
dev.off()
####################################################################################################################################################################################
#read corrected genus table for WT samples(pick top 40 based on sum)
library("readxl")
library("tidyverse")
top40_WT<- read_xlsx("genus_WT_top40.xlsx")
head(top40_WT)
top40_WT.t <- as.data.frame(t(top40_WT[,-1]))
colnames(top40_WT.t) <- top40_WT$SampleID
head(top40_WT.t)
#read in annotations
ann<- read.table("WT_annotation.txt", header = T)
ann <- ann %>% remove_rownames %>% column_to_rownames(var="SampleID")
#reorder martrix based on annotation
top40_WT.t.ordered <- top40_WT.t[rownames(ann), ]
#Z score transformation
top40_WT.t.ordered.norm <- scale(top40_WT.t.ordered)
summary(top40_WT.t.ordered.norm )
################################################################################
#specify colors for annotation
ann_colors = list(Group = c('WT'="#1b271a", 'WT_Perio'="#c71b05"))

#sample based heatmap
library(pheatmap)
library(RColorBrewer)
#png(file.path(path, "blankvsperio(sample_based).png"), width=25,height=11,units="cm",res=1200)
col.pal <- RColorBrewer::brewer.pal(9, "GnBu")
pdf(file.path(path, "WT_top40_genus_individual_sample_based_new.pdf"), width = 8, height =10)
pheatmap(t(top40_WT.t.ordered.norm), scale = "row", color = colorRampPalette(c("navy", "white", "firebrick3"))(15), 
         cellheight = 11,
         cellwidth = 10,
         fontsize_row=12, show_colnames = F, annotation_col = ann, annotation_colors = ann_colors,
         cluster_cols = F)
dev.off()
####################################################################################################################################################################################
#the end of the script


