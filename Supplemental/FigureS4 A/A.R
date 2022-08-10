phylo_clean <- read_rds("phylo3.rds")
#phylum level stacked bar-plot
phylo_relative = transform_sample_counts(phylo_clean, function(x) {(x/sum(x))*100} )
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
  scale_x_discrete(limits=c("WT_PBS", "WT_Suc", "KO_PBS", "KO_Suc"))+
  scale_fill_manual(values = c( "#972C8DFF","#FA6B09FF",  "#149BEDFF" ,"#EE0011FF","#9A703EFF",
                                "#15983DFF", "#A1C720FF","#16A08CFF", "#FEC10BFF","#0C5BB0FF",
                                "#BA6222FF","#808080","#EC579AFF"))+
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

