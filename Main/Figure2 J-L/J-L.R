#load phyloseq object
phylo_geno <- readRDS("phylo2.rds")
####################################################################################################################################################################################
#Alpha Diversity
# Setup environment
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools 
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling  
library(rstatix)
#check difference in the number of reads
#summary(sample_sums(phylo_clean))

#calculate alpha diversity
phylo.div <- alpha(phylo_geno, index = "all")
datatable(phylo.div)
# get the metadata out as seprate object
phylo.meta <- meta(phylo_geno)
# Add the rownames as a new colum for easy integration later.
phylo.meta$SampleID <- rownames(phylo.meta)
# Add the rownames to diversity table
phylo.div$SampleID <- rownames(phylo.div)
#merge these two data frames into one
phylo.df <- merge(phylo.div,phylo.meta, by = "SampleID")
#check the tables
colnames(phylo.df)
#calculate PD
library("btools")
alpha_pd<-estimate_pd(phylo_geno)
alpha_pd$SampleID <- rownames(alpha_pd)
phylo.df <- merge(phylo.df,alpha_pd, by = "SampleID")
# check the tables
colnames(phylo.df)
#select columns
phylo.df2 <-phylo.df[, c("Group", "observed", "diversity_shannon", "PD")]
#replace column names
colnames(phylo.df2) <- c("Group", "Observed ASVs", "Shannon", "PD")
# check
colnames(phylo.df2)
#Using Group as id variables
phylo_df_melt <- reshape2::melt(phylo.df2)
#check
head(phylo_df_melt)
#pairwise compare
stat.test <- phylo_df_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ Group) %>% add_significance("p") %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.signif != "ns") #filter non-sig pairs
stat.test

write.csv(stat.test, file = "~/Desktop/Succinate_perio_rerun/R/Figure2-swab-blank_vs_perio/stats/alpha_diversity.csv")

#plot
png(file.path(path, "alpha_div.png"), width=14,height=12,units="cm",res=1200, bg = "transparent")
ggboxplot(phylo_df_melt, x = "Group", y = "value", fill = "Group", color = "white",legend= "right", facet.by = "variable", scales = "free")+
  stat_pvalue_manual(stat.test, label = "p.signif", y.position =c(195, 205, 185, 4.8,4.3,5.6,5.3,16,17,12.5) , bracket.size=0.6, size=8 )+
  geom_point(aes(fill=Group,color=Group), size=2.5)+
  geom_boxplot(aes(fill=Group,color=Group),alpha=0, fatten = 0.8, lwd=0.8)+ 
  scale_x_discrete(limits=c("WT", "WT_Perio", "KO", "KO_Perio"))+
  scale_fill_manual(values=c("white", "white", "white", "white"), guide = FALSE)+
  scale_color_manual(breaks = c("WT", "WT_Perio", "KO", "KO_Perio"),
                     values = c("WT"="#1b271a", "WT_Perio"="#c71b05", "KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  ylab("Alpha Diversity Measure")+
  theme(axis.title.y=element_text(size=16, color = "black"))+
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y=element_text(size=14, color="black"),
    strip.text = element_text(size = 12,  color="black"))+
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face = "bold"))+
  theme(
    legend.position="top",
    legend.text=element_text(size=16, color="black"),
    legend.title = element_blank())+
  theme(plot.background = element_rect(fill = "transparent",colour = NA))

dev.off()
####################################################################################################################################################################################
#Beta Diveristy
#method:Phylogenetic beta-diversity metrics: Weighted Unifrac
#consider the abundances of different taxa
#PERMANOVA
library(vegan)
metadf <- data.frame(sample_data(phylo_geno))

unifrac.dist <- UniFrac(phylo_geno, 
                        weighted = T, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ Group, data = metadf)
permanova#check significance
swab.wt.uni <- ordinate(phylo_geno, "PCoA", "unifrac", weighted=T)

wt.unifrac <- plot_ordination(phylo_geno, 
                              swab.wt.uni, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_fill_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05", "KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  scale_color_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05", "KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  ggtitle("Weighted UniFrac") + geom_point(size = 4)+
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

pdf(file.path(path, "swab_Blankvsperio_wt.unifrac_PCOA.pdf"), height = 4, width = 4)
#png(file.path(path, "Weighted_Uni.png"), width=14,height=10,units="cm",res=1200,bg = "transparent")
wt.unifrac
dev.off()
################################################################################
#method:bray
#PERMANOVA
dist <- phyloseq::distance(phylo_geno, method = "bray")
perma <- adonis(dist~Group, data = as(sample_data(phylo_geno), "data.frame"), permutations = 999)
perma
dist <- permutest(betadisper(dist, sample_data(phylo_geno)$Group), permutations = 999)
perma
swab.bray <- ordinate(phylo_geno, "PCoA", "bray")
bray <- plot_ordination(phylo_geno, 
                        swab.bray, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_fill_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05", "KO"="#006f71", "KO_Perio"="#ed9e0a"))+
  scale_color_manual(values = c("WT"="#1b271a", "WT_Perio"="#c71b05", "KO"="#006f71", "KO_Perio"="#ed9e0a"))+
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

pdf(file.path(path, "swab_Blankvsperio_braycurtis_PCOA.pdf"), height = 4, width = 4)
bray
dev.off()
####################################################################################################################################################################################
#plot genus heatmap
####################################################################################################################################################################################
#read corrected genus table(pick top 20 based on sum)
library("readxl")
top20.genus<- read_xlsx("genus_top20.xlsx")
head(top20.genus)
#transpose the df and keep the sampleIDs as the header
top20.genus.t <- as.data.frame(t(top20.genus[,-1]))
colnames(top20.genus.t) <- top20.genus$SampleID
head(top20.genus.t)
#read in annotations
ann<- read.table("annotation.txt", header = T)
ann <- ann %>% remove_rownames %>% column_to_rownames(var="SampleID")
#reorder martrix based on annotation
top20.genus.t.ordered <- top20.genus.t[rownames(ann), ]
#Z score transformation
top20.genus.t.ordered.norm <- scale(top20.genus.t.ordered)
summary(top20.genus.t.ordered.norm)
################################################################################
#specify colors for annotation
ann_colors = list(Group = c('WT'="#1b271a", 'WT_Perio'="#c71b05", 'KO'="#006f71", 'KO_Perio'="#ed9e0a"))

#sample based heatmap
library(pheatmap)
#png(file.path(path, "blankvsperio(sample_based).png"), width=25,height=11,units="cm",res=1200)
pdf(file.path(path, "top20_genus_individual_sample_based.pdf"), width = 12, height =4)
pheatmap(t(top20.genus.t.ordered.norm), scale = "row", color = colorRampPalette(c("navy", "white", "firebrick3"))(25), 
         cellheight = 11,
         cellwidth = 10,
         fontsize_row=12, show_colnames = F, annotation_col = ann, annotation_colors = ann_colors,
         cluster_cols = F)
dev.off()
####################################################################################################################################################################################
#the end of the script


