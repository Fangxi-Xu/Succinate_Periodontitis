#load phyloseq object
phylo_clean <- read_rds("phylo4.rds")
############################################################################################
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
phylo.div <- alpha(phylo_clean, index = "all")
datatable(phylo.div)
# get the metadata out as seperate object
phylo.meta <- meta(phylo_clean)
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
alpha_pd<-estimate_pd(phylo_clean)
alpha_pd$SampleID <- rownames(alpha_pd)
phylo.df <- merge(phylo.df,alpha_pd, by = "SampleID")
# check the tables
colnames(phylo.df)
#select columns
phylo.df2 <-phylo.df[, c("Treatment", "observed", "diversity_shannon", "PD")]
#replace column names
colnames(phylo.df2) <- c("Group", "Observed ASVs", "Shannon", "PD")
# check
colnames(phylo.df2)
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

#plot
png(file.path(path, "alpha_div.png"), width=14,height=12,units="cm",res=1200, bg = "transparent")
ggboxplot(phylo_df_melt, x = "Group", y = "value", fill = "Group", color = "white",legend= "right", facet.by = "variable", scales = "free")+
  geom_point(aes(fill=Group,color=Group), size=2.5)+
  geom_boxplot(aes(fill=Group,color=Group),alpha=0, fatten = 0.8, lwd=0.8)+ 
  scale_x_discrete(limits=c("Vehicle", "7a"))+
  scale_fill_manual(values=c("white", "white"),guide=FALSE)+
  scale_color_manual(values = c("7a"="#9629cc", "Vehicle"="#db2b10"))+
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

############################################################################################
############################################################################################
##Beta Diversity
#method:Phylogenetic beta-diversity metrics: Weighted Unifrac
#consider the abundances of different taxa
#PERMANOVA
library(vegan)
metadf <- data.frame(sample_data(phylo_clean))

unifrac.dist <- UniFrac(phylo_clean, 
                        weighted = T, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ Treatment, data = metadf)
permanova#check significance
wt.uni <- ordinate(phylo_clean, "PCoA", "unifrac", weighted=T)

wt.unifrac <- plot_ordination(phylo_clean, 
                              wt.uni, color="Treatment") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Treatment, fill=Treatment))+
  scale_fill_manual(values = c("7a"="#9629cc", "Vehicle"="#db2b10"))+
  scale_color_manual(values = c("7a"="#9629cc", "Vehicle"="#db2b10"))+
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

#png(file.path(path, "Ligature_Weighed_Uni.png"), width=10,height=10,units="cm",res=1200,bg = "transparent")
pdf(file.path(path, "7a_UW.pdf"), width = 4, height = 4)
wt.unifrac 
dev.off()

#method:bray
#PERMANOVA
dist <- phyloseq::distance(phylo_clean, method = "bray")
perma <- adonis(dist~Treatment, data = as(sample_data(phylo_clean), "data.frame"), permutations = 999)
perma
dist <- permutest(betadisper(dist, sample_data(phylo_clean)$Treatment), permutations = 999)
perma

ligature.bray <- ordinate(phylo_clean, "PCoA", "bray")
bray <- plot_ordination(phylo_clean, 
                        ligature.bray, color="Treatment") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Treatment, fill=Treatment))+
  scale_fill_manual(values = c("7a"="#9629cc", "Vehicle"="#db2b10"))+
  scale_color_manual(values = c("7a"="#9629cc", "Vehicle"="#db2b10"))+
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
#png(file.path(path, "Ligature_bray_PCOA.png"), width=10,height=10,units="cm",res=1200, bg = "transparent")
pdf(file.path(path, "7a_BC.pdf"), width = 4, height = 4)
bray
dev.off()
################################################################################################
################################################################################################
#phylum level stacked bar-plot
phylo_relative = transform_sample_counts(phylo_clean, function(x) {(x/sum(x))*100} )
glom_phylum <- tax_glom(phylo_relative, taxrank = 'Phylum')
data_phylum <- psmelt(glom_phylum) # create dataframe from phyloseq object
data_phylum$phylum <- as.character(data_phylum$Phylum) #convert to character
#simple way to rename phyla with < 1% abundance
data_phylum$Phylum[data_phylum$Phylum == "p__"] <- "Unclassified"
#calculate statistics
phylum_ligature_summary <-data_phylum %>%
  group_by(Treatment,Phylum) %>%  # the grouping variable
  summarise(mean_Abun = mean(Abundance),# calculates the mean of each group
            sd_Abun = sd(Abundance), # calculates the standard deviation of each group
            n_Abun = n(),  # calculates the sample size per group
            SE_Abun = sd(Abundance)/sqrt(n()))# calculates the standard error of each group
#png(file.path(path, "phylum(stacked).png"), width=12,height=8,units="cm",res=1200)
pdf(file.path(path, "phylum(stacked).pdf"), width = 5, height = 4)
ggplot(phylum_ligature_summary, aes(fill=Phylum, y=mean_Abun, x=Treatment)) + 
  scale_x_discrete(limits=c("Vehicle", "7a"))+
  scale_fill_manual(values = c("#149BEDFF", "#FF0000", "#EC579AFF", "#FEC10BFF","#AA4371",
                               "#15983DFF", "#A1C720FF", "#16A08CFF", "#808080"))+
  ylab("Mean Relative Abundance(%)") +
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.text.x=element_text(size=10, angle = 0, color = "black"))+
  theme(axis.text.y=element_text(size=12, color = "black"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16, color = "black"))+
  theme(legend.text=element_text(size=14, color = "black"))+
  theme(legend.title=element_text(size=16, color = "black"))+
  theme(legend.position="right")
dev.off()

#phylum statistics
#calculate statistics
phylum_melt <- reshape2::melt(data_phylum)
phylum.stat.test <- phylum_melt  %>%
  group_by(Phylum) %>%
  wilcox_test(value ~ Treatment) %>% add_significance("p") %>%
  adjust_pvalue() %>%
  add_significance() #%>%
  #filter(p.signif != "ns") #filter non-sig pairs
phylum.stat.test
#######################################################################################################################################
#import top20 genus table
#read corrected genus table(pick top 20 based on sum)
library("readxl")
top20.genus.ligature<- read_xlsx("top20_genus.xlsx")
head(top20.genus.ligature)
#transpose the df and keep the sampleIDs as the header
top20.genus.ligature.t <- as.data.frame(t(top20.genus.ligature[,-1]))
colnames(top20.genus.ligature.t) <- top20.genus.ligature$SampleID
head(top20.genus.ligature.t)
#read in annotations
ann<- read.table("annotation_ligature.txt", header = T)
#reorder martrix based on annotation
top20_ordered <- top20.genus.ligature.t[rownames(ann), ]
#Z score transformation
top20_ordered_zscore <- scale(top20_ordered)
summary(top20_ordered_zscore)
#specify colors for annotation
ann_colors = list(Treatment = c("7a"="#9629cc", "Vehicle"="#FF8C00"))
#individual sample based heatmap
#heatmap plot for each sample
library(pheatmap)
pdf(file.path(path, "top20_genus_individual_sample_based.pdf"), width = 6, height =4)
pheatmap(t(top20_ordered_zscore), scale = "row", color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         cellwidth = 10, cellheight = 10,
         fontsize_row=8, fontsize_col = 8, annotation_col = ann, annotation_colors = ann_colors,
         angle_col = "45",  cluster_cols = F,
         annotation_legend = F)
dev.off()
#######################################################################################################################################
#######################################################################################################################################
#end of the script

