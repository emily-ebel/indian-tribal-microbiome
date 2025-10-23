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
### Figure 2 Supplement 2: 
#      Stacked bars of 16S genera per sample
#####################################################################


###################################
###  16S  genus-level abundance 

d16<-read.table("Fig2_relab_genus16S.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
d16$AK_SK_31 = NULL # sample has no metadata

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
common_genera <- subset(totals16,totals16$mean_Rel_ab >= 0.027) 
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
    "Streptococcus", "Ligilactobacillus",
    "Bacteroides", "Faecalibacterium","Agathobacter", "Succinivibrio",
   "Lachnospiraceae NK4A136 group","Escherichia-Shigella","Treponema","Ruminobacter","other")
dplot$Genus = factor(dplot$Genus,levels=temp)
dplot$Rel_ab = as.numeric(dplot$Rel_ab) # and make rel_ab numeric
dplot$Tribe <- factor(dplot$Tribe,
                      levels = c("Kabui","Warli","Gond","Madia","Boto","Balti","Brokpa","Purigpa"))


## Set up colors for the plot
col_list = c("#E6AB02", #  Segatella
             '#fc8803', # Leyella
             '#b30707',# Bifidobacterium ******
             '#f5051d',# Catenibacterium  ******
             '#cc3835',# Ruminococcaceae UCG-002  ******
             "#f714a0", # Streptococcus ******
             '#f57fc0',# Ligilactobacillus ******
             '#6d048a', # Bacteroides ****,# Bacteroides
             "#1F78B4",# Faecalibacterium
             '#29992b',# Agathobacter
             "#A6CEE3", # Succinivibrio
             '#46dbac',# "Lachnospiraceae NK4A136 group"
             '#1b1d87',# Eschericia/Shigella
             "#B2DF8A",# Treponema
             "#3d611e",# Ruminobacter
             'gray'# other
)



#### Plot

s<- ggplot(dplot, aes(fill=Genus, y=Rel_ab, x=Sample)) +  theme_classic() +
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
  scale_fill_manual("legend",values = col_list)  + ylab("Relative abundance") + xlab("") +
  theme(legend.title=element_blank()) 
s + facet_wrap(~Tribe,scales="free_x",nrow=2) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
