# create plot folder
path <- "./plot"
dir.create(path)
#load phyloseq object
phylo_clean <- readRDS("phylo1.rds")
##########################################################################################
#Alpha Diversity
# Setup environment
library(microbiome) 
library(phyloseq) 
library(microbiomeutilities) 
library(ggpubr) 
library(DT)
library(data.table) 
library(dplyr)
library(rstatix)

#check difference in the number of reads
summary(sample_sums(phylo_clean))

#calculate alpha diversity
phylo.div <- alpha(phylo_clean, index = "all")
datatable(phylo.div)
# get the metadata out as seprate object
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
  wilcox_test(value ~ Group) %>%
  add_significance() %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.signif != "ns") #filter non-sig pairs
stat.test

#plot
png(file.path(path, "alpha_div.png"), width=13,height=12,units="cm",res=1200, bg = "transparent")
ggboxplot(phylo_df_melt, x = "Group", y = "value", fill = "Group", color = "white",legend= "right", facet.by = "variable", scales = "free")+
  stat_pvalue_manual(stat.test, label = "p.signif", y.position =c(300,5,17) , bracket.size=0.6, size=8 )+
  geom_point(aes(fill=Group,color=Group), size=2.5)+
  geom_boxplot(aes(fill=Group,color=Group),alpha=0, fatten = 1, lwd=0.8)+ 
  scale_fill_manual(values=c("white", "white", "white", "white"))+
  scale_color_manual(values = c("Normal/mild"="#000000", "Severe"="#FF0000"))+
  ylab("Alpha Diversity Measure")
dev.off()
################################################################################################
################################################################################################
#Beta Diveristy
#method:Phylogenetic beta-diversity metrics: Weighted Unifrac
#consider the abundances of different taxa
##PERMANOVA
library(vegan)
metadf <- data.frame(sample_data(phylo_clean))
unifrac.dist <- UniFrac(phylo_clean, 
                        weighted = T, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ Group, data = metadf)
permanova#check significance
wt.uni <- ordinate(phylo_clean, "PCoA", "unifrac", weighted=T)

wt.unifrac<-plot_ordination(phylo_clean, 
                            wt.uni, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_fill_manual(values = c("Normal/mild"="#000000", "Severe"="#FF0000"))+
  scale_color_manual(values = c("Normal/mild"="#000000", "Severe"="#FF0000"))+
  ggtitle("Weighted UniFrac") + geom_point(size = 4)

pdf(file.path(path, "Weighted_Unifrac_PCOA.pdf"),width=4,height=4)
wt.unifrac
dev.off()

#method:bray
dist <- phyloseq::distance(phylo_clean, method = "bray")
perma <- adonis(dist~Group, data = as(sample_data(phylo_clean), "data.frame"), permutations = 999)
perma
#dist <- permutest(betadisper(dist, sample_data(phylo_clean)$Group), permutations = 999)
#perma
bray1 <- ordinate(phylo_clean, "PCoA", "bray")
bray2 <- plot_ordination(phylo_clean, 
                         bray1, color="Group") + 
  stat_ellipse(geom = "polygon", type="norm", size=1, alpha=0.1, aes(color=Group, fill=Group))+
  scale_fill_manual(values = c("Normal/mild"="#000000", "Severe"="#FF0000"))+
  scale_color_manual(values = c("Normal/mild"="#000000", "Severe"="#FF0000"))+
  ggtitle("Bray Curtis") + geom_point(size = 4)

#png(file.path(path, "Bray_Curtis.png"), width=12,height=11,units="cm",res=1200, bg = "transparent")
pdf(file.path(path, "Bray_Curtis.png.pdf"),width=4,height=4)
bray2 
dev.off()
################################################################################################
#phylum level stacked bar-plot
phylo_relative = transform_sample_counts(phylo_clean, function(x) {(x/sum(x))*100} )
glom_phylum <- tax_glom(phylo_relative, taxrank = 'Phylum')
data_phylum <- psmelt(glom_phylum) # create dataframe from phyloseq object
data_phylum$phylum <- as.character(data_phylum$Phylum) #convert to character
#simple way to rename phyla with < 1% abundance
data_phylum$Phylum[data_phylum$Phylum == "p__"] <- "Unclassified"
#calculate statistics
phylum_summary <-data_phylum %>%
  group_by(Group,Phylum) %>%  # the grouping variable
  summarise(mean_Abun = mean(Abundance),# calculates the mean of each group
            sd_Abun = sd(Abundance), # calculates the standard deviation of each group
            n_Abun = n(),  # calculates the sample size per group
            SE_Abun = sd(Abundance)/sqrt(n()))# calculates the standard error of each group
#png(file.path(path, "phylum(stacked).png"), width=13,height=10,units="cm",res=1200,bg = "transparent")
pdf(file.path(path, "phylum(stacked).pdf"), width = 6, height =5)
ggplot(phylum_summary, aes(fill=Phylum, y=mean_Abun, x=Group)) + 
  scale_x_discrete(limits=c("Normal/mild", "Severe"))+
  scale_fill_manual(values = c("#DAA520","#FA6B09FF", "#149BEDFF", "#AA4371",
                               "#15983DFF", "#FEC10BFF","#A1C720FF","#16A08CFF", "#0C5BB0FF", "#EE0011FF", "#9A703EFF","808080"))+
  ylab("Mean Relative Abundance (%)") +
  geom_bar(position="stack", stat="identity")

dev.off()

#phylum statistics
library(rstatix)
#calculate statistics
phylum_melt <- reshape2::melt(data_phylum)
phylum.stat.test <- phylum_melt  %>%
  group_by(Phylum) %>%
  wilcox_test(value ~ Group) %>% add_significance("p") %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.signif != "ns") #filter non-sig pairs
phylum.stat.test

##########################

#end of the script
