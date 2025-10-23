library(data.table)
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)
library(magrittr)
library("dendextend")
library("vegan")
library(ade4)
library('factoextra')
library(grid)
library(gridExtra)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)


#################################################################################
#   Figure 3 Supplement 1 : VANISH taxa (16S)
#################################################################################


##################################################################################
#####   From NCBI Taxonomy: genera belonging to each VANISH family

#f__Paraprevotellaceae = c() # DOES NOT EXIST
f__Succinivibrionaceae = c('Anaerobiospirillum','Succinatimonas','Succinivibrio','Ruminobacter','Succinimonas')
f__Prevotellaceae = c('Hoylesella','Leyella','Marseilla','Massiliprevotella','Metaprevotella','Palleniella','Prevotellamassilia',
                      'Pseudoprevotella','Segatella','Hallella', 'Paraprevotella', 'Prevotella', 'Alloprevotella',
                      'Prevotellaceae_NK3B31_group','Prevotellaceae_UCG-003','Xylanibacter','Palleniella')
p__Spirochaetes = c('Brachyspira','Brevinema','Borrelia','Borreliella','Breznakiella','Brucepastera','Bullifex',
                    'Thermospira','Turneriella', 'Thiospirochaeta',  'Treponema',
                    'Longinema','Leptonema','Leptospira','Leadbettera',
                    'Cristispira', 'Clevelandina',
                    'Gracilinema','Helmutkoenigia','Zuelzera', 'Marispirochaeta','Oceanispirochaeta',
                    'Rectinema','Alkalispirochaeta',
                    'Parasphaerochaeta','Pleomorphochaeta',
                    'Salinispira','Sediminispirochaeta','Spirochaeta','Sphaerochaeta',
                    'Entomospira','Exilispira')



####################################################################################
#####   INDIAN TRIBES: This study

# Relative abundance data
d16<-read.table("Fig2_relab_genus16S.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")

# Exclude outliers 
d16$AK_SK_49<-NULL
d16$AK_SK_10<-NULL
d16$AK_SK_32<-NULL
d16$AK_SR_6<-NULL
d16$AK_SG_17<-NULL
d16$AK_SG_18<-NULL
d16$AK_SR_1<-NULL
d16$AK_SR_4 <- NULL # Streptococcus 
d16$AK_SK_31 <- NULL # oral sample included by mistake


# Subset to genera in VANISH families 
IM_Prevotellaceae_relab <- subset(d16,d16$Genus %in% f__Prevotellaceae)
row.names(IM_Prevotellaceae_relab) <- IM_Prevotellaceae_relab$Genus
IM_Prevotellaceae_relab = IM_Prevotellaceae_relab[,!(names(IM_Prevotellaceae_relab) %in% c("Genus"))]
IM_Prevotellaceae_indsum = colSums(IM_Prevotellaceae_relab) # total of this family per sample

IM_Spirochaetes_relab <- subset(d16,d16$Genus %in% p__Spirochaetes)
row.names(IM_Spirochaetes_relab) <- IM_Spirochaetes_relab$Genus
IM_Spirochaetes_relab = IM_Spirochaetes_relab[,!(names(IM_Spirochaetes_relab) %in% c("Genus"))]
IM_Spirochaetes_indsum = colSums(IM_Spirochaetes_relab) # total of this family per sample

IM_Succinivibrionaceae_relab <- subset(d16,d16$Genus %in% f__Succinivibrionaceae)
row.names(IM_Succinivibrionaceae_relab) <- IM_Succinivibrionaceae_relab$Genus
IM_Succinivibrionaceae_relab = IM_Succinivibrionaceae_relab[,!(names(IM_Succinivibrionaceae_relab) %in% c("Genus"))]
IM_Succinivibrionaceae_indsum = colSums(IM_Succinivibrionaceae_relab) # total of this family per sample


# create data frame to plot
IM_plot <- as.data.frame(cbind(IM_Prevotellaceae_indsum,IM_Succinivibrionaceae_indsum,IM_Spirochaetes_indsum ))
names(IM_plot) <- c("Prevotellaceae",'Succinivibrionaceae','Spirochaetes')

# add metadata
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T)
meta<-subset(meta,meta$sample %in% row.names(IM_plot))
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
#
IM_plot$Tribe = NA
IM_plot$Region = NA
for (row in 1:nrow(IM_plot)){
  samplename = row.names(IM_plot)[row]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  IM_plot[row,4] = tribe
  IM_plot[row,5] = region
}

IM_plot_melt <- reshape2::melt(IM_plot,id.vars=c("Tribe","Region"))





#####################################################################################################
#####   URBAN INDIA: Tandon et al., 2018


tandon <-read.table("Tandon2018_GenusRelab.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss",row.names=1)

# Subset to genera in VANISH families 
tandon_Prevotellaceae_relab <- subset(tandon,row.names(tandon) %in% f__Prevotellaceae)
tandon_Prevotellaceae_relab_t <- as.data.frame(t(tandon_Prevotellaceae_relab)) 
tandon_Prevotellaceae_relab_t = as.data.frame(sapply(tandon_Prevotellaceae_relab_t, as.numeric))
tandon_Prevotellaceae_indsum = rowSums(tandon_Prevotellaceae_relab_t) 

tandon_Succinivibrionaceae_relab <- subset(tandon,row.names(tandon) %in% f__Succinivibrionaceae)
tandon_Succinivibrionaceae_relab_t <- as.data.frame(t(tandon_Succinivibrionaceae_relab)) 
tandon_Succinivibrionaceae_relab_t = as.data.frame(sapply(tandon_Succinivibrionaceae_relab_t, as.numeric))
tandon_Succinivibrionaceae_indsum = rowSums(tandon_Succinivibrionaceae_relab_t) 

tandon_Spirochaetes_relab <- subset(tandon,row.names(tandon) %in% p__Spirochaetes) # this is empty, so:
tandon_Spirochaetes_indsum = rep(0,80)


# create data frame to plot
tandon_plot <- as.data.frame(cbind(tandon_Prevotellaceae_indsum,tandon_Succinivibrionaceae_indsum,tandon_Spirochaetes_indsum ))
names(tandon_plot) <- c("Prevotellaceae",'Succinivibrionaceae','Spirochaetes')
tandon_plot$Tribe = "Urban India" # city, not really a tribe # "Ahmedabad"
tandon_plot$Region = "Urban India" # Indian state
tandon_plot_melt <- reshape2::melt(tandon_plot,id.vars=c("Tribe","Region"))
tandon_plot_melt$value = tandon_plot_melt$value / 100





####################################################################################
#####   URBAN CALIFORNIA: Wastyk et al., 2021


# Relative abundance data labeled as ASVs
d1 <- readRDS("Wastyk2021_ASVrelab.rds") 

# Links ASVs to taxonomy
d2 <- readRDS("Wastyk2021_ASVtaxonomy.rds") 

# Since these data include multiple timepoints (including after a dietary intervention),
# limit to the two baseline samples and take their average rel_ab
t1 <- subset(d1,d1$Timepoint=='1')
t2 <- subset(d1,d1$Timepoint=='2')
# remove the timepoint columns
t1 = t1[,!(names(t1) %in% c("Timepoint"))] 
t2 = t2[,!(names(t2) %in% c("Timepoint"))]
# and average the rest 
baseline <- rbindlist(list(t1,t2))[,lapply(.SD,mean), list(Participant,Group,Group_value)]
baseline_t <- t(baseline)


# Subset to genera in VANISH families 
## Identify the ASVs in each taxonomic group
fefifo_Prevotellaceae_ASVs_temp <- subset(d2,d2$Family %in% c("f__Prevotellaceae","f__[Paraprevotellaceae]"))
fefifo_Prevotellaceae_ASVs <- fefifo_Prevotellaceae_ASVs_temp$rep_ASV_label # get ASV IDs

fefifo_Succinivibrionaceae_temp <- subset(d2,d2$Family == 'f__Succinivibrionaceae')
fefifo_Succinivibrionaceae_ASVs <- fefifo_Succinivibrionaceae_temp$rep_ASV_label # get ASV IDs

fefifo_Spirochaetes_temp  <- subset(d2,d2$Phylum == 'p__Spirochaetes')
fefifo_Spirochaetes_ASVs <- fefifo_Spirochaetes_temp$rep_ASV_label # get ASV IDs

## Subset the rel_ab data to these ASVs
fefifo_Prevotellaceae_relab <- subset(baseline_t,row.names(baseline_t) %in% fefifo_Prevotellaceae_ASVs)
fefifo_Prevotellaceae_relab_t <- as.data.frame(t(fefifo_Prevotellaceae_relab))
fefifo_Prevotellaceae_relab_t = as.data.frame(sapply(fefifo_Prevotellaceae_relab_t, as.numeric))
fefifo_Prevotellaceae_indsum = rowSums(fefifo_Prevotellaceae_relab_t) 
#
fefifo_Succinivibrionaceae_relab <- subset(baseline_t,row.names(baseline_t) %in% fefifo_Succinivibrionaceae_ASVs)
fefifo_Succinivibrionaceae_relab_t <- as.data.frame(t(fefifo_Succinivibrionaceae_relab)) # this loses individual label info, but that's ok
fefifo_Succinivibrionaceae_relab_t = as.data.frame(sapply(fefifo_Succinivibrionaceae_relab_t, as.numeric))
fefifo_Succinivibrionaceae_indsum = rowSums(fefifo_Succinivibrionaceae_relab_t) # total rel_ab per ind of this taxonomic rank 
#
fefifo_Spirochaetes_relab <- subset(baseline_t,row.names(baseline_t) %in% fefifo_Spirochaetes_ASVs)
fefifo_Spirochaetes_relab_t <- as.data.frame(t(fefifo_Spirochaetes_relab)) # this loses individual label info, but that's ok
fefifo_Spirochaetes_relab_t = as.data.frame(sapply(fefifo_Spirochaetes_relab_t, as.numeric))
fefifo_Spirochaetes_indsum = rowSums(fefifo_Spirochaetes_relab_t) # total rel_ab per ind of this taxonomic rank 



# Create df for plot 
fefifo_plot <- as.data.frame(cbind(fefifo_Prevotellaceae_indsum,fefifo_Succinivibrionaceae_indsum,fefifo_Spirochaetes_indsum ))
names(fefifo_plot) <- c("Prevotellaceae",'Succinivibrionaceae','Spirochaetes')
fefifo_plot$Tribe = "Urban California"
fefifo_plot$Region = "Urban California"
fefifo_plot_melt <- reshape2::melt(fefifo_plot,id.vars=c("Tribe","Region"))





####################################################################################
#####   COMBINE POPULATIONS


combine_plot <- rbind(fefifo_plot_melt,IM_plot_melt,tandon_plot_melt)
combine_plot$Region = factor(combine_plot$Region,levels=c("Urban California","Urban India" ,"Northeast Hills","Coast","Deccan Plateau","Trans-Himalayas"))
combine_plot$variable = factor(combine_plot$variable,levels=c("Spirochaetes","Succinivibrionaceae","Prevotellaceae"))
combine_plot$Tribe = factor(combine_plot$Tribe,levels=c("Urban India","Kabui","Warli","Madia","Gondia","Boto","Brokpa","Balti","Purigpa","Urban California"))


### PLOT Figure 3 Supplement 1 
ggplot(combine_plot, aes(x=Region, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme_cowplot() + ylab("Relative abundance") + xlab("") +
  scale_fill_manual(values = c("#525050","#edc42f",  "#b33f25")) + coord_flip() + 
  theme(legend.title=element_blank()) +
  theme(legend.position="none")  + 
  ggtitle("VANISH taxa")

# with legend 
ggplot(combine_plot, aes(x=Region, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme_cowplot() + ylab("Relative abundance") + xlab("") +
  scale_fill_manual(values = c("#525050","#edc42f",  "#b33f25")) + coord_flip() + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme(legend.title=element_blank()) +
  theme(legend.position="top") 
