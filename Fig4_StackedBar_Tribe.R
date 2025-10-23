library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(cowplot)
library("pheatmap")
library(viridis)


path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)


#####################################################################
### Figure 4
#       Stacked bars by tribe 
#            A) 16S
#     Supp 2 C) metagenomics 
#####################################################################



###################################
###  16S  genus-level abundance 

d16<-read.table("Fig2_relab_genus16S.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
d16$AK_SK_31 = NULL # sample included by mistake 

### remove outliers 
d16$AK_SK_49 = NULL
d16$AK_SK_10 = NULL
d16$AK_SK_32 = NULL
d16$AK_SR_6 = NULL
d16$AK_SG_17 = NULL
d16$AK_SG_18 = NULL
d16$AK_SR_1 = NULL
d16$AK_SR_2 = NULL
d16$AK_SR_4 = NULL # STREPTOCOCCUS 

## MELT to three-column format 
d16melt <- reshape::melt(d16, id = c("Genus"))
names(d16melt) <- c("Genus","Sample","Rel_ab") # now equivalent to 'genusRA' for metgenomic
d16melt$Sample <- as.character(d16melt$Sample)

## Add region/tribe metadata
meta16 <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta16 <- meta16 %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))
meta16<-subset(meta16,meta16$sample %in% names(d16))

d16melt$Tribe = NA
d16melt$Region = NA 
for (row in 1:nrow(d16melt)){
  samplename = d16melt[row,2]
  tribe <- as.character(subset(meta16,meta16$sample==samplename)$Tribe)
  region <- as.character(subset(meta16,meta16$sample==samplename)$Region)
  d16melt[row,4] = tribe
  d16melt[row,5] = region
}






##################################################################
#### SELECT TOP GENERA (TO HAVE LABELS/COLORS IN PLOT LEGEND) 
##################################################################

# calculate total abundance per genus, by tribe
totals16 <- d16melt %>% dplyr::group_by(Genus, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Rel_ab))

# pick a rel_ab threshold that yields 15 genera 
common_genera <- subset(totals16,totals16$mean_Rel_ab >= 0.0245) #  a threshold 0.0245 yields 15 genera for SILVA 138 (no outliers). With outliers, need 0.027
length(unique(common_genera$Genus))
top_genera = unique(common_genera$Genus)

# get abundances for only these genera
d16melt_top = subset(d16melt,d16melt$Genus%in%top_genera) 
d16melt_top <- subset(d16melt_top,is.na(d16melt_top$Tribe)==F)



## ADD 'OTHER' COLUMN (genera not abundant enough to get their own label):
# First, sum up the abundances per person of genera that DO get a label
labeled_per_person <- d16melt_top %>%
  dplyr::group_by(Sample)  %>%
  dplyr::summarise(percent_with_genus_label = sum(Rel_ab))
mean(labeled_per_person$percent_with_genus_label )

# Then define other as 1 minus that sum
labeled_per_person$other = 1-labeled_per_person$percent_with_genus_label 

# Add that sum to dmelt16_top in the form of one row per person, where 'genus' = 'other'
dplot <- d16melt_top

for (i in (unique(d16melt_top$Sample))){
  tribe = unique(subset(d16melt_top,d16melt_top$Sample==i)$Tribe)
  region = unique(subset(d16melt_top,d16melt_top$Sample==i)$Region)
  other_number = as.numeric( subset(labeled_per_person,labeled_per_person$Sample==i)[1,3] )
  vector = c('other',i,other_number,tribe,region)
  dplot = rbind(dplot,vector)
}


## Prep df for plotting, including setting the order for the genera to appear in the legend   
temp=c("Segatella","Leyella","Bifidobacterium","Catenibacterium","UCG-002",
       "Megasphaera","Ligilactobacillus","Dialister","Asteroleplasma",
       "Faecalibacterium","Agathobacter","Succinivibrio",
       "Treponema","Escherichia-Shigella","Ruminobacter","other")
dplot$Genus = factor(dplot$Genus,levels=temp)
dplot$Rel_ab = as.numeric(dplot$Rel_ab) # and make rel_ab numeric
dplot$Tribe <- factor(dplot$Tribe, levels = c("Purigpa","Brokpa","Balti","Boto","Kabui","Madia","Gond","Warli"))




## Set up colors for the plot
col_list = c("#E6AB02", #  Segatella
             '#fc8803', # Leyella
             '#b30707',# Bifidobacterium ******
             '#f5051d',# Catenibacterium  ******
             '#cc3835',# Ruminococcaceae UCG-002  ******
             "#f714a0", # Megasphaera ******
             '#f57fc0',# Ligilactobacillus ******
             '#f5a2b1',# Dialister *** more abundant in TH than NEH
             '#f74f60', # Asteroleplasma *****
             "#1F78B4",# Faecalibacterium
             '#29992b',# Agathobacter
             "#A6CEE3", # Succinivibrio
             "#B2DF8A" ,# Treponema 
             '#1b1d87',# Eschericia/Shigella 
             "#3d611e", # Ruminobacter  
             'gray'# other
)


#### Plot FIGURE 4A 
ggplot(dplot, aes(fill=Genus, y=Rel_ab, x=Tribe)) +  theme_classic() +
  geom_bar(position="fill", stat="identity") + ylab("Relative abundance") +
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black",size=10))+
  scale_fill_manual("legend",values = col_list) + 
  xlab("")+
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.4, "cm")) +
  theme(axis.title=element_text(size=10)) +
  theme(axis.text.x = element_text(size=10))  + 
  coord_flip() 




##################################################################################################################
##################################################################################################################


##########################################
###  Metagenomic  genus-level abundance 

d<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 

# Remove outliers 
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2","AK_SR_4")
d<-subset(d,d$sample %in% outliers == FALSE)


#### Convert to table with headers: Sample Genus Rel_ab
# Make the empty table 
samples <- unique(d$sample); length(samples)
genera = unique(d$genus); length(genera)
genusRA <- data.frame(matrix(ncol = 3, nrow = length(samples)*length(genera)))
names(genusRA) <- c("Sample","Genus","Rel_ab")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(d,d$sample==samples[[s]]) # for each sample 
  #sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab)) # DON'T scale RA to sum to 1 per sample here, because unclassified reads will mathematically ignored later
  for (t in 1:length(genera)){ # for each genus 
    count = count+1
    genusRA[count,1] = as.character(samples[[s]])
    genusRA[count,2]= as.character(genera[[t]])
    if ( (genera[[t]] %in% sdf$genus)==TRUE){
      rows <- subset(sdf,sdf$genus == genera[[t]]) 
      genusRA[count,3] = sum(rows$rel_ab) 
    } else (genusRA[count,3]=0)
  }
}

# Add metadata
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))
meta<-subset(meta,meta$sample %in% d$sample)
#
genusRA$Tribe = NA
genusRA$Region = NA 
for (row in 1:nrow(genusRA)){
  samplename = genusRA[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  genusRA[row,4] = tribe
  genusRA[row,5] = region
}

genusRA<-subset(genusRA,genusRA$Genus!="")



##################################################################
#### SELECT TOP GENERA (TO HAVE LABELS/COLORS IN PLOT LEGEND) 
##################################################################

# calculate total abundance per genus, by tribe
totals <- genusRA %>% dplyr::group_by(Genus, Tribe)  %>% dplyr::summarise(mean_Rel_ab = mean(Rel_ab))

# select the most abundant
common_genera <- subset(totals,totals$mean_Rel_ab >= 0.014) 
length(unique(common_genera$Genus))
top_genera = unique(common_genera$Genus)

# Add Ligilactobacillus (and drop UBA1436, the least abundant) 
top_genera = c(top_genera,"Ligilactobacillus")
top_genera = top_genera[ !top_genera == "UBA1436"]

# get abundances for only these genera
d_top = subset(genusRA,genusRA$Genus%in%top_genera) 


## ADD 'OTHER' COLUMN (genera not abundant enough to get their own label):
# First, sum up the abundances per person of genera that DO get a label
labeled_per_person <- d_top %>%
  dplyr::group_by(Sample)  %>%
  dplyr::summarise(percent_with_genus_label = sum(Rel_ab))
mean(labeled_per_person$percent_with_genus_label ) 

# But in metagenomics, each person has a lot of reads unclassified 
classified_per_person <- genusRA %>%
  dplyr::group_by(Sample)  %>%
  dplyr::summarise(percent_with_any_genus_label = sum(Rel_ab))
mean(classified_per_person$percent_with_any_genus_label ) 

# we can effectively ignore the unclassified by defining 'other' as 'total' (ignoring unclassified) minus 'top'
labeled_per_person$other = classified_per_person$percent_with_any_genus_label-labeled_per_person$percent_with_genus_label 
# these will not sum to 100% because the remaining is unclassified

# Add that sum to d_top in the form of one row per person, where 'genus' = 'other'
dplot <- d_top
for (i in (unique(d_top$Sample))){
  tribe = unique(subset(d_top,d_top$Sample==i)$Tribe)
  region = unique(subset(d_top,d_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person,labeled_per_person$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}



## Prep dplot for plotting, including setting the order for the genera to appear in the legend    
temp = c("Segatella","Leyella","Hallella",
         "Bifidobacterium","Ligilactobacillus",
         "Faecalibacterium","Agathobacter","Prevotellamassilia","Succinivibrio",
         "CAG-127","RC9","Treponema_D","CAG-83","Klebsiella","Alistipes",
         "other")
dplot$Genus = factor(dplot$Genus,levels=temp)
dplot$Rel_ab = as.numeric(dplot$Rel_ab) 
dplot$Tribe <- factor(dplot$Tribe,levels = c("Purigpa","Brokpa","Balti","Boto","Kabui","Madia","Gond","Warli"))


## Set up colors for the plot
col_list = c("#E6AB02", #  Segatella 
             '#fc8803',# Leyella
             '#ffbc05',# Hallella
             '#b30707',# Bifidobacterium
             '#f57fc0',# Ligilactobacillus
             "#1F78B4",# Faecalibacterium
             '#29992b',# Agathobacter
             '#ccf24e',# Prevotellamassilia
             "#b5dff5",# Succinivibrio
             '#4a13f0',# CAG-127
             '#45b7f5',# RC9
             "#B2DF8A",# Treponema_D
             "#fcf40d",# CAG-83
             '#717311',# Klebsiella
             "#f3f5ab", # Alistipes 
             'gray'# other
)


## PLOT

ggplot(dplot, aes(fill=Genus, y=Rel_ab, x=Tribe)) +  theme_classic() +
  geom_bar(position="fill", stat="identity") + ylab("Relative abundance") +
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black",size=10))+
  scale_fill_manual("legend",values = col_list) + 
  xlab("")+ 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.4, "cm")) +
  theme(axis.title=element_text(size=10)) +
  theme(axis.text.x = element_text(size=10))  + 
  coord_flip() 
