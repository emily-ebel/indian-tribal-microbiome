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
### Figure 2 Supplement 3:
#      Stacked bars of metagenomic genera per sample
#####################################################################


##########################################
###  Metagenomic  genus-level abundance 

d<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 

### Add metadata
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta<-subset(meta,meta$sample %in% d$sample)
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))


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
  #sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab)) # DON'T scale RA to sum to 1 per sample here, because unclassified reads will be mathematically ignored later
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
common_genera <- subset(totals,totals$mean_Rel_ab >= 0.017) 
length(unique(common_genera$Genus))
top_genera = unique(common_genera$Genus)

# To highlight outliers, make sure to include Phocaeicola & Streptococcus 
top_genera = c(top_genera,"Streptococcus","Phocaeicola")
# and remove 2 least-abundant to balance
top_genera = top_genera[ !top_genera == 'Treponema_D']
top_genera = top_genera[ !top_genera == 'Alistipes']

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
mean(classified_per_person$percent_with_any_genus_label ) # on average, 49% of the rel_ab gets any genus label

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
         "Bifidobacterium","Streptococcus",
         "Bacteroides","Phocaeicola", "Parabacteroides" ,
         "Faecalibacterium","Agathobacter","Prevotellamassilia","Succinivibrio","Acetatifactor","CAG-127","RC9","other")
dplot$Genus = factor(dplot$Genus,levels=temp)
dplot$Rel_ab = as.numeric(dplot$Rel_ab) #
dplot$Tribe <- factor(dplot$Tribe,
                      levels = c("Kabui","Warli","Gond","Madia","Boto","Balti","Brokpa","Purigpa"))


## Set up colors for the plot
col_list = c("#E6AB02", #  Segatella ]
             '#fc8803',# Leyella  
             '#ffbc05',# Hallella
             '#b30707',# Bifidobacterium    
             "#f714a0", # Streptococcus
             '#6d048a', # Bacteroides ****
             "#ca91ed", #Phocaeicola, ****
             '#9522b5',  #Parabacteroides ****
             "#1F78B4",# Faecalibacterium 
             '#29992b',# Agathobacter
             '#ccf24e',# Prevotellamassilia
            "#b5dff5",# Succinivibrio
              '#66ded9',# Acetatifactor
             '#4a13f0',# CAG-127
             '#45b7f5',# RC9
             'gray'# other
)



#### Plot

s<- ggplot(dplot, aes(fill=Genus, y=Rel_ab, x=Sample)) +  theme_classic() +
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  scale_fill_manual("legend",values = col_list)  + ylab("Relative abundance") + xlab("") +
  theme(legend.title=element_blank()) 
L = s + facet_wrap(~Tribe,scales="free_x",nrow=2) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2


# ragg::agg_png("Fig2_Supp3.png", width = 7*1.65, height = 6*1.25, units = "in", res = 300)
# L
# dev.off()    
