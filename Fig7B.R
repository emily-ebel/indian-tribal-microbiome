library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(cowplot)
library("pheatmap")
library(viridis)
library(stringr)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)


##################################################
###### Metagenomic abundance data - Indian Microbiome
##################################################

d <-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 

meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta<-subset(meta,meta$sample %in% d$sample)

#############################
#### CONVERT to genus abundance 

samples <- unique(d$sample); length(samples)
genera = unique(d$genus); length(genera)
genusRA <- data.frame(matrix(ncol = 3, nrow = length(samples)*length(genera)))
names(genusRA) <- c("Sample","Genus","Rel_ab")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(d,d$sample==samples[[s]]) # for each sample 
   sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab)) # scale RA 
   
  for (t in 1:length(genera)){ # for each genus 
    count = count+1
    genusRA[count,1] = as.character(samples[[s]])
    genusRA[count,2]= as.character(genera[[t]])
    
    if ( (genera[[t]] %in% sdf$genus)==TRUE){
      rows <- subset(sdf,sdf$genus == genera[[t]]) 
      genusRA[count,3] = sum(rows$rel_ab_adj)  # save the rel_ab for all genomes in the genus 
    } else (genusRA[count,3]=0)
  }
}

genusRA <- subset(genusRA,genusRA$Genus != "")

### ADD TRIBE
genusRA$Tribe = NA
genusRA$Region = NA 
for (row in 1:nrow(genusRA)){
  samplename = genusRA[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  genusRA[row,4] = tribe
  genusRA[row,5] = region
}

genusRA <- genusRA %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
genusRA <- genusRA %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
genusRA <- genusRA %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
genusRA <- genusRA %>% mutate(Region = str_replace(Region, "West", "Coast"))
genusRA <- genusRA %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
genusRA <- genusRA %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))

########
## Add factor variable for group (Bacteroides outlier or not) 

# remove AK_SR_4 first for being a Streptococcus outlier
genusRA <- subset(genusRA,genusRA$Sample != "AK_SR_4")
genusRA$group =NA
outliers = c("AK_SK_49","AK_SK_10","AK_SK_32","AK_SR_6","AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2")
g.out <-subset(genusRA,genusRA$Sample %in% outliers == TRUE)
g.out$group = "Tribal Outlier"
g.non<- subset(genusRA,genusRA$Sample %in% outliers == FALSE )
g.non$group = "Tribal Non-outlier"
genusRA_IM = rbind(g.out,g.non)






##################################################
###### FeFiFo METAGENOMES from Wastyk et al. ########
##################################################

# abundance data
fffd <-read.table("Wastyk2021_metagenomics.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss") 
fffd <- subset(fffd,fffd$breadth>=0.5)

# add columns for genus and species
fffd$genus = NA
fffd$species = NA
# fill in
for (row in 1:dim(fffd)[[1]]){
  classification = fffd[row,14]  
  temp = strsplit(classification,'g__')[[1]][2]
  genus = strsplit(temp,';s__')[[1]][1]
  species = strsplit(temp,';s__')[[1]][2]
  fffd[row,15] = genus
  fffd[row,16] = species
}
fffd <- subset(fffd,fffd$genus != "")


#### Need to update some species/genus names 
## EUBACTERIUM -- new names already present in FFF (unlike humann2, which was fixed in Fig5_stackedbars.R)
## LACTOBACILLUS -- two species (delbrueckii and acidophilus) have current NCBI names already 
#######################################
## PREVOTELLA - needs updating !
### FFF timepoints 1 and 2 (see below) contain the following Prevotella species:
# bivia -- no name change (NCBI)
# disiens -- no name change (NCBI)
# sp003447235 -- no name change (Fig2_S4)
# copri -- update to Segatella (Fig2_S4; NCBI)
# sp000435635 -- update to Segatella (Fig2_S4)
# sp001275135 -- update to Xylanibacter (Fig2_S4)
# sp900545525 -- update to Xylanibacter (Fig2_S4)
# sp900551055 -- update to Xylanibacter (Fig2_S4)

for (row in 1:nrow(fffd)){
  genus = fffd[row,15]
  if (genus == "Prevotella"){
    species = fffd[row,16] 
    if (species =="Prevotella copri") {
      species = "Segatella copri" 
      genus = "Segatella"
    } else if (species =="Prevotella sp000435635"){
      species = "Segatella sp000435635"
      genus = "Segatella"
    } else if  (species =="Prevotella sp001275135"){
      species = "Xylanibacter sp001275135"
      genus = "Xylanibacter"
    } else if  (species =="Prevotella sp900545525"){
      species = "Xylanibacter sp900545525"
      genus = "Xylanibacter"
    } else if  (species =="Prevotella sp900551055"){
      species = "Xylanibacter sp900551055"
      genus = "Xylanibacter"
    } 
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    fffd[row,16] = species
    fffd[row,15] = genus  } }



#############################
#### FILTER for timepoints 
# we only want timepoints 1 and 2 (baseline) -- we'll average them after converting to genus rel_ab

fffmeta <-read.table("Wastyk_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = F,numerals="allow.loss") 
fffmeta_keep = subset(fffmeta,fffmeta$timepoint %in% c('1','2'))
fffd_filt <- subset(fffd,fffd$sample %in% fffmeta_keep$SampleNameCORRECT)


#############################
#### CONVERT to genus abundance 

samples <- unique(fffd_filt$sample); length(samples)
genera = unique(fffd_filt$genus); length(genera)
genusRA_fff <- data.frame(matrix(ncol = 3, nrow = length(samples)*length(genera)))
names(genusRA_fff) <- c("Sample","Genus","Rel_ab")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(fffd_filt,fffd_filt$sample==samples[[s]]) # for each sample 
  sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab)) # scale RA 
  
  for (t in 1:length(genera)){ # for each genus 
    count = count+1
    genusRA_fff[count,1] = as.character(samples[[s]])
    genusRA_fff[count,2]= as.character(genera[[t]])
    
    if ( (genera[[t]] %in% sdf$genus)==TRUE){
      rows <- subset(sdf,sdf$genus == genera[[t]]) 
      genusRA_fff[count,3] = sum(rows$rel_ab_adj)  # save the rel_ab for all genomes in the genus 
    } else (genusRA_fff[count,3]=0)
  }
}


#############################
#### Average across two timepoints for each sample 
# to summarize over timepoint, need to use each sample name to pull the right timepoint
genusRA_fff$timepoint = NA
genusRA_fff$subject = NA
for (row in 1:dim(genusRA_fff)[[1]]){
  sampname = as.character(as.data.frame(genusRA_fff[row,])$Sample)
  subject =as.character(subset(fffmeta_keep,fffmeta_keep$SampleNameCORRECT == sampname)$subject)
  timepoint = as.character(subset(fffmeta_keep,fffmeta_keep$SampleNameCORRECT == sampname)$timepoint)
  genusRA_fff[row,4] = timepoint
  genusRA_fff[row,5] = subject
}

# remove 'sample' column and summarize over timepoint
genusRA_fff_sum <- genusRA_fff[ , -which(names(genusRA_fff) %in% c('Sample'))]
genusRA_fff_sum = genusRA_fff_sum  %>% dplyr::group_by(subject,Genus) %>%  dplyr::summarise(mean_relab = mean(Rel_ab))


## Now we have something similar to genusRA 
genusRA_FFF = genusRA_fff_sum 



####### PREPARE TO COMBINE ACROSS STUDIES

# Add tribe/region/group to FFF
genusRA_FFF$Tribe = "California"
genusRA_FFF$Region = "California"
genusRA_FFF$group = "California"


#### MATCH NAMES
names(genusRA_IM) = c("Sample","Genus","Rel_ab","Tribe", "Region","Group" )
names(genusRA_FFF) = c("Sample","Genus","Rel_ab","Tribe", "Region","Group" )

# How many genera are shared? 
sum(unique(genusRA_IM$Genus) %in% genusRA_FFF$Genus) # 136
length(unique(genusRA_IM$Genus)) # of 203 in IM
length(unique(genusRA_FFF$Genus)) # and 369 in FFF

# COMBINE
genusRA = rbind(genusRA_IM,genusRA_FFF)

# note that genera that are only in one group -- like "1XD42-69" which is only in California --
# don't have any entry (row) in genusRA. This won't matter to plotting unless those genera are very 
# common in one group. So let's just see.



##########################################################
###### SELECT TOP GENERA (TO HAVE LABELS IN STACKED BAR PLOT) 
##########################################################

#  calculate total abundance per genus, by group, so I can make sure to include taxa that meet a minimum abundance in at least 1 group)
totals <- genusRA %>% dplyr::group_by(Genus, Group)  %>% dplyr::summarise(mean_Rel_ab = mean(Rel_ab))

##################
### Select the most abundant genera overall to label in plot 
common_genera <- subset(totals,totals$mean_Rel_ab >= 0.028) 
length(unique(common_genera$Genus))
top_genera = unique(common_genera$Genus) =
top_genera

# get abundances for only these genera
d_top = subset(genusRA,genusRA$Genus%in%top_genera) 


######################
## ADD 'OTHER' COLUMN for  genera not abundant enough to get their own label

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person <- d_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_genus_label = sum(Rel_ab))
mean(labeled_per_person$percent_with_genus_label ) 

# since we scaled the rel_ab at the beginning, any rel_ab unaccounted for belongs to genomes 
# that couldn't be classified to genus, which we'll group as 'other'
classified_per_person <- genusRA %>%
  group_by(Sample)  %>%
  summarise(percent_with_any_genus_label = sum(Rel_ab))
mean(classified_per_person$percent_with_any_genus_label ) 
# define 'other' as 'total'  minus 'top'
labeled_per_person$other = 1-labeled_per_person$percent_with_genus_label 

# Add that sum to d_top in the form of one row per person, where 'genus' = 'other'
dplot <- d_top
for (i in (unique(d_top$Sample))){
  tribe = unique(subset(d_top,d_top$Sample==i)$Tribe)
  region = unique(subset(d_top,d_top$Sample==i)$Region)
  group = unique(subset(d_top,d_top$Sample==i)$Group)
  other_number = as.numeric(subset(labeled_per_person,labeled_per_person$Sample==i)[1,3])
  vector = data.frame(list(i,'other',other_number,tribe,region,group)); names(vector) = names(dplot)
  dplot = rbind(dplot,vector)
}



##########################################
###### BAR PLOT BY GROUP ########
##########################################

### Order top genera 
dplot$Genus = factor(dplot$Genus,levels=c("Segatella","Leyella",
                            "Bacteroides","Phocaeicola","Parabacteroides",
                            "Faecalibacterium","Agathobacter",
                            "Bifidobacterium","Alistipes",
                            "Acetatifactor","Prevotellamassilia",
                            "Blautia_A","other"))


## organize data / adjust class
dplot$Rel_ab = as.numeric(dplot$Rel_ab) 

## Set group labels for plot
unique(dplot$Group)
dplot <- dplot %>% mutate(Group = str_replace(Group, "Tribal Non-outlier", "Indian Tribes"))
dplot <- dplot %>% mutate(Group = str_replace(Group, "Tribal Outlier", "Tribal Outliers"))
unique(dplot$Group)
dplot$Group <- factor(dplot$Group, levels = c("Indian Tribes","California","Tribal Outliers"))


## PLOT
col_list = c("#E6AB02", #  Segatella
             '#fc8803', # Leyella
             '#6d048a', # Bacteroides 
             "#ca91ed", #Phocaeicola, 
             '#9522b5',  #Parabacteroides 
             "#1F78B4",# Faecalibacterium 
             '#29992b',# Agathobacter
             '#b30707',# Bifidobacterium   
             "#f3f5ab", # Alistipes 
             '#66ded9',# Acetatifactor
             '#ccf24e',# "Prevotellamassilia"  
           '#5c1906',# Blautia_A 
             'gray'# other
) 

p<- ggplot(dplot, aes(fill=Genus, y=Rel_ab, x=Group)) +  theme_classic() +
  geom_bar(position="fill", stat="identity") + ylab("Relative abundance") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.9, hjust=1))+
  scale_fill_manual("legend",values = col_list) + xlab("")+
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.4, "cm")) +
  theme(axis.title=element_text(size=8.5)) + #theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) 
p





# save legend (bars look bad...)
ggsave("Fig7B_legend.pdf",  plot=p,width=0.84*3,height=2.02*1.5,units="in")

# screenshot the bars
