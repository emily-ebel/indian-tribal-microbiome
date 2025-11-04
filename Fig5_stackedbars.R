library('ggplot2')
library(reshape2)
library(stringr)
library(cowplot)
library(grid)
library(gridExtra)
library(data.table)
library(dplyr)
library(RColorBrewer)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)


########################################################
###### Fig 5: FUNCTION STRATIFIED BY TAXA #############
########################################################

### Read in the GO data

d<-read.delim("GOcombined.tsv", header = TRUE, sep = "\t",numerals="allow.loss",colClasses=c('character',rep("numeric",100)))
# Cells contain relative abundance. 
# Columns sum to >1 because stratified, unclassified, and total are all included  ?? 
# Here are the commands that made this table:
# humann2_renorm_table --input  ~/humann2-2.8.2/AK_SG_10_output/AK_SG_10_GO1000_readable.tsv --units relab --output ~/humann2-2.8.2/GO_RA/AK_SG_10_GO_RA.tsv  
# humann2_join_tables --input ~/humann2-2.8.2/GO_RA --output ~/humann2-2.8.2/GO_RA/GOcombined.tsv

# clean up the column names (to sample name only)
names1 <- names(d)[2:length(names(d))]
names2 <- sub('NODUP_HMN_UNMAPPED_TRIM_MARKED_', '', names1)
names3 <- sub('_Abundance.RPKs', '', names2)
names(d)[2:length(names(d))] <- names3
length(names(d)) # 1 is function column, so 100 samples 

# exclude outliers
outliers = c("AK_SK_49","AK_SK_10","AK_SK_32","AK_SR_6","AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2","AK_SR_4")
others = c('AK_SK_24',"AK_SK_19","AK_SK_34") # the first two have very little sequencing data and were excluded from rel_ab data; the third has no metadata
d <- d[ , -which(names(d) %in% c(outliers,others))]
length(names(d)) # still 88 samples ... need to exclude enrichment samples still

# use metadata to exclude enrichment samples 
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- subset(meta,meta$sample %in% names(d))
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau")); meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond")) 
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui")); meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas")); meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast")); meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
fecalkeep = as.character(subset(meta,meta$sample_type=='fecal')$sample)
d <- d[ , which(names(d) %in% c(fecalkeep,'GeneFamily'))]
dsamps <- names(d)[2:length(names(d))] 
meta<-subset(meta,meta$sample %in% dsamps) # now 67 non-outlier samples, as expected



################################################################################
####  pull out one pathway at a time 
### PANEL B: # GO:0004459: [MF] L-lactate dehydrogenase activity
################################################################################

func <- d[d$GeneFamily %like% "GO:0004459", ]  
# the top row of func is OVERALL and equal to all the other rows summed together

#### transpose this one stratified pathway  so rows are samples
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
# drop first row (pathway names) -- add back in a few lines 
df = df0[-1,]
# for some reason all then numbers are characters now, so convert to numeric 
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 

### simplify these GO names -- not to number, b/c they're all the same; but to organism
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified 
# total thus contains the RA of the pathway for each sample, among all the pathways in that sample 
# each sample has some % UNMAPPED and UNGROUPED reads that comprise some of this RA 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names

# melt df to have columns: c(Taxon,Rel_ab, Sample, Tribe)
df$Sample = row.names(df)
dmelt <-  reshape::melt(df, id = c("Sample"))
dmelt <- subset(dmelt,dmelt$variable != "total") # remove 'total' but leave 'unclassified'
dmelt$Tribe = NA; dmelt$Region = NA 

for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  dmelt[row,4] = tribe
  dmelt[row,5] = region
}

names(dmelt) <- c("Sample","Taxon","Pathway_relab", "Tribe", "Region") 


##########################################################
#### Combine individual species into genera ######## 
##########################################################

# Might actually be easiest just to add a 'genus' column to dmelt
dmelt$Taxon = as.character(dmelt$Taxon)
dmelt$Genus = NA
for (row in 1:nrow(dmelt)){
  species = dmelt[row,2]
  genus = strsplit(species, '.s__')[[1]][1]
  dmelt[row,6] = genus
}



################## ############# ########################## ########################## 
#  UPDATE EUBACTERIUM NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  species = dmelt[row,2]
  
  if (genus == "Erysipelotrichaceae_noname"){
    if (species =="Erysipelotrichaceae_noname.s__Eubacterium_biforme"){
      species = "Holdemanella.s__Holdemanella_biformis"
      genus = "Holdemanella" } 
  } else if (genus == "Eubacterium") {
    if (species =="Eubacterium.s__Eubacterium_eligens"){
      species = "Lachnospira.s__Lachnospira_eligens"
      genus = "Lachnospira"
    } else if (species =="Eubacterium.s__Eubacterium_hallii"){
      species = "Anaerobutyricum.s__Anaerobutyricum_hallii"
      genus = "Anaerobutyricum" 
  #  } else if (species =="Eubacterium.s__Eubacterium_limosum"){ # no change in NCBI Taxonomy
  #  } else if (species =="Eubacterium.s__Eubacterium_ramulus"){  # no change in NCBI Taxonomy
    } else if (species =="Eubacterium.s__Eubacterium_rectale"){
      species = "Agathobacter.s__Agathobacter_rectalis"
      genus = "Agathobacter" 
    } else if (species =="Eubacterium.s__Eubacterium_siraeum"){ # new genus name pending. Make sure these special characters aren't an issue ... works at this stage of table, anyway 
      species = "[Eubacterium].s__[Eubacterium]_siraeum"
      genus = "[Eubacterium]"
#    } else if (species =="Eubacterium.s__Eubacterium_ventriosum"){ # no change in NCBI Taxonomy
    } 
  } # end if genus == Eubacterium
    
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
} # end whole loop


################## ############# ########################## ########################## 
#  UPDATE PREVOTELLA NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Prevotella"){
    species = dmelt[row,2] 
    if (species =="Prevotella.s__Prevotella_buccae"){
      species = "Segatella.s__Segatella_buccae" # according to NCBI taxonomy
      genus = "Segatella"
    } else if (species =="Prevotella.s__Prevotella_copri"){
      species = "Segatella.s__Segatella_copri"
      genus = "Segatella"
    }else if (species =="Prevotella.s__Prevotella_stercorea"){
      species = "Leyella.s_Leyella_stercorea"
      genus = "Leyella"
    } else if  (species =="Prevotella.s__Prevotella_timonensis"){
      species = "Hoylesella.s__Hoylesella_timonensis"
      genus = "Hoylesella"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Prevotella loops
}

################## ############# ########################## ########################## 
#  UPDATE LACTOBACILLUS NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lactobacillus"){
    species = dmelt[row,2] 
    if (species =="Lactobacillus.s__Lactobacillus_casei_paracasei"){
      species = "Lacticaseibacillus.s__Lacticaseibacillus_paracasei" # according to chatGPT
      genus = "Lacticaseibacillus"
  #  } else if (species =="Lactobacillus.s__Lactobacillus_delbrueckii"){ # no change 
    }else if (species =="Lactobacillus.s__Lactobacillus_plantarum"){
      species = "Lactiplantibacillus.s_Lactiplantibacillus_plantarum"
      genus = "Lactiplantibacillus"
    }else if (species =="Lactobacillus.s__Lactobacillus_ruminis"){
      species = "Ligilactobacillus.s_Ligilactobacillus_ruminis"
      genus = "Ligilactobacillus"
    } else if  (species =="Lactobacillus.s__Lactobacillus_salivarius"){
      species = "Ligilactobacillus.s__Ligilactobacillus_salivarius"
      genus = "Ligilactobacillus"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Lactobacillus loops
}

################## ############# ########################## ########################## 


## and then create a new table, where each sample has a summed entry for each genus
dmelt_glom = data.frame(matrix(ncol=5,nrow=1))
names(dmelt_glom) = names(dmelt)[1:5]
genera = unique(dmelt$Genus)

## then for each sample, for each genus:
samples = unique(dmelt$Sample)
for (samp in samples){
  samp_dmelt = subset(dmelt,dmelt$Sample==samp)
  tribe = unique(samp_dmelt$Tribe)
  region = unique(samp_dmelt$Region)
  
  for (genus in genera){
    g_samp_dmelt = subset(samp_dmelt, samp_dmelt$Genus ==genus)
    genus_total = sum(g_samp_dmelt$Pathway_relab)
    
    # save a new row
    newrow = c(samp,genus,genus_total,tribe,region)
    dmelt_glom=rbind(dmelt_glom,newrow)
  }
}

dmelt_glom = na.omit(dmelt_glom)





##########################################################
#### SELECT TOP TAXA (for legend) FOR THIS ONE PATHWAY        
##########################################################

dmelt_saved = dmelt
dmelt = dmelt_glom
dmelt$Pathway_relab = as.numeric(dmelt$Pathway_relab)

# calculate total abundance per taxon and tribe across  fecal samples
totals <- dmelt %>% dplyr::group_by(Taxon, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Pathway_relab))

# select most abundant taxa to include in legend
common_taxa<- subset(totals,totals$mean_Rel_ab >=  3.8e-6)  
length(unique(common_taxa$Taxon)) # AIM FOR 8 (including unclassified) 
top_taxa = unique(common_taxa$Taxon)
top_taxa

# get abundances for only these genera
dmelt_top = subset(dmelt,dmelt$Taxon%in%top_taxa)

######################
## ADD 'OTHER' COLUMN (genera that are not abundant enough to get their own label)

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person_top <- dmelt_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Because the relabs of this pathway don't sum to 1 -- it's relab of ALL pathways --
# 'other' is not just 1-x. It's the total -x. So, calculate

#  the total rel_ab for this pathway per person
labeled_per_person_all <- dmelt %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Then define other as 1 minus that sum
labeled_per_person_top$other = labeled_per_person_all$percent_with_taxon_label - labeled_per_person_top$percent_with_taxon_label
#labeled_per_person_top$checksum = labeled_per_person_top$other + labeled_per_person_top$percent_with_taxon_label

# Add 'other' to dmelt_top in the form of one row per person, where 'genus' = 'other'
dplot <- dmelt_top
dplot$Taxon = as.character(dplot$Taxon)
for (i in (unique(dmelt_top$Sample))){
  tribe = unique(subset(dmelt_top,dmelt_top$Sample==i)$Tribe)
  region = unique(subset(dmelt_top,dmelt_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person_top,labeled_per_person_top$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}
unique(dplot$Taxon)
##### at this point, unclassified and other are both there


# order genera by mean abundance
# Common taxa has mean_rel_ab of each taxon by tribe, but DOES NOT INCLUDE 'OTHER'
common_taxa$Taxon = as.character(common_taxa$Taxon) ## THIS does not include 'other'
t <- common_taxa %>%
  group_by(Taxon)  %>%
  summarise(mean_Rel_ab = mean(mean_Rel_ab))
t <- t[order(t$mean_Rel_ab,decreasing=TRUE),]
top_taxa_ordered = t$Taxon
top_taxa_ordered 

# then edit the order so that red (TH-enriched) taxa are on top
L_lactate_dehydrogenase_activity = c("Bifidobacterium","Collinsella","Catenibacterium", "Ligilactobacillus","Agathobacter","Roseburia","Treponema","other","unclassified")
L_lactate_dehydrogenase_activity %in% dplot$Taxon
length(L_lactate_dehydrogenase_activity) == length(unique(dplot$Taxon))
dplot$Taxon = factor(dplot$Taxon,levels=L_lactate_dehydrogenase_activity)


### Master table of colors (each genus has same color across plots)

color_vec = c('#29992b',"#f3f5ab",'#6d048a',"#b30707","#07eafa",
              '#07fa02','#f5051d', '#8878bf','#f09292','#523e0f',
              '#f60cfa', '#1b1d87',"#1F78B4",'#6e3154','#717311',
             '#b792d4','lightpink','#054ffc', '#fc8803','#f57fc0',
              'black','#84b9fa',"#E6AB02",'#d43fb6', "#B2DF8A",
              'darkgray','lightgray',"tan")

genus_vec = c("Agathobacter","Alistipes","Bacteroides","Bifidobacterium","Blautia",
  "Butyrivibrio","Catenibacterium","Citrobacter","Collinsella","Enterobacter",
  "Enterococus","Escherichia","Faecalibacterium","Fusobacterium","Klebsiella",  
  "5_1_57FAA","Lactobacillus","Lactococcus","Leyella","Ligilactobacillus",
  "Odoribacter", "Roseburia","Segatella","Streptococcus","Treponema",
  "other","unclassified","Enterococcus")  
colortable= as.data.frame(cbind(color_vec,genus_vec))
names(colortable)=c("color","genus")  

# check if any top taxa in this pathway are NOT in the color table 
for (i in L_lactate_dehydrogenase_activity){
  if (i %in% genus_vec == FALSE){
    print(i)
  }
}

## adjust class
dplot$Pathway_relab = as.numeric(dplot$Pathway_relab) 
dplot$Tribe <- factor(dplot$Tribe, levels = c("Boto","Balti","Brokpa","Purigpa","Warli", "Gond","Madia","Kabui"))


## subset and order colors for this particular function / set of genera
specific_colortable <- subset(colortable,colortable$genus %in% L_lactate_dehydrogenase_activity)
specific_colortable = specific_colortable[match(L_lactate_dehydrogenase_activity, specific_colortable$genus),]
col_list = specific_colortable$color

## PLOT ! FACET WRAP ALL TRIBES TOGETHER
### TO AUTOMATICALLY CALCULATE THE DISTANCE BETWEEN DATA AND TOP OF PLOT:
temp <- dplot %>% group_by(Sample) %>% summarise(sum_relab = sum(Pathway_relab))
margin = max(temp$sum_relab)*1.05

s <- ggplot(data=dplot, aes(x=reorder(Sample, -Pathway_relab, function(x){ sum(x) }), y=Pathway_relab, fill=Taxon)) +
  geom_col() + theme_classic() +
  ylab("Relative abundance") + theme_cowplot() +
  scale_fill_manual("legend",values = col_list) + 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.35, "cm")) +
  ggtitle("L-lactate dehydrogenase activity (GO:0004459)") +
  theme(plot.title = element_text(size = 12,face="bold")) +
  theme(axis.title=element_text(size=12,face="plain")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, margin),labels =  ~ sprintf(fmt = "%2.e", .)) + # removes space between 0 and x-axis
 theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

L <- s + facet_wrap(~Tribe,scales="free_x",nrow=1) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
L

ggsave("Fig5B.pdf", plot=L,width=4.64*2.2,height=0.95*2.2,units="in")







################################################################################
####  pull out one pathway at a time 
### PANEL C: # galactose metabolic process
################################################################################

func <- d[d$GeneFamily %like% "GO:0006012", ] 

#### transpose this one stratified pathway  so rows are samples
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,] # drop first row (pathway names) -- add back in a few lines 
# for some reason all then numbers are characters now, so convert to numeric 
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 

### simplify these GO names -- not to number, b/c they're all the same; but to organism
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names

# melt df to have columns: c(Taxon,Rel_ab, Sample, Tribe)
df$Sample = row.names(df)
dmelt <-  reshape::melt(df, id = c("Sample"))
dmelt <- subset(dmelt,dmelt$variable != "total") # remove 'total' but leave 'unclassified'
dmelt$Tribe = NA; dmelt$Region = NA 

for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  dmelt[row,4] = tribe
  dmelt[row,5] = region
}

names(dmelt) <- c("Sample","Taxon","Pathway_relab", "Tribe", "Region") 


##########################################################
#### Combine individual species into genera ######## 
##########################################################

# Might actually be easiest just to add a 'genus' column to dmelt
dmelt$Taxon = as.character(dmelt$Taxon)
dmelt$Genus = NA
for (row in 1:nrow(dmelt)){
  species = dmelt[row,2]
  genus = strsplit(species, '.s__')[[1]][1]
  dmelt[row,6] = genus
}


################## ############# ########################## ########################## 
#  UPDATE EUBACTERIUM NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  species = dmelt[row,2]
  
  if (genus == "Erysipelotrichaceae_noname"){
    if (species =="Erysipelotrichaceae_noname.s__Eubacterium_biforme"){
      species = "Holdemanella.s__Holdemanella_biformis"
      genus = "Holdemanella" } 
  } else if (genus == "Eubacterium") {
    if (species =="Eubacterium.s__Eubacterium_eligens"){
      species = "Lachnospira.s__Lachnospira_eligens"
      genus = "Lachnospira"
    } else if (species =="Eubacterium.s__Eubacterium_hallii"){
      species = "Anaerobutyricum.s__Anaerobutyricum_hallii"
      genus = "Anaerobutyricum" 
      #  } else if (species =="Eubacterium.s__Eubacterium_limosum"){ # no change in NCBI Taxonomy
      #  } else if (species =="Eubacterium.s__Eubacterium_ramulus"){  # no change in NCBI Taxonomy
    } else if (species =="Eubacterium.s__Eubacterium_rectale"){
      species = "Agathobacter.s__Agathobacter_rectalis"
      genus = "Agathobacter" 
    } else if (species =="Eubacterium.s__Eubacterium_siraeum"){ # new genus name pending. Make sure these special characters aren't an issue ... works at this stage of table, anyway 
      species = "[Eubacterium].s__[Eubacterium]_siraeum"
      genus = "[Eubacterium]"
      #    } else if (species =="Eubacterium.s__Eubacterium_ventriosum"){ # no change in NCBI Taxonomy
    } 
  } # end if genus == Eubacterium
  
  ## ACTUALLY UPDATE SPECIES/GENUS IN DF
  dmelt[row,2] = species
  dmelt[row,6] = genus
} # end whole loop


################## ############# ########################## ########################## 
#  UPDATE PREVOTELLA NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Prevotella"){
    species = dmelt[row,2] 
    if (species =="Prevotella.s__Prevotella_buccae"){
      species = "Segatella.s__Segatella_buccae" # according to NCBI taxonomy
      genus = "Segatella"
    } else if (species =="Prevotella.s__Prevotella_copri"){
      species = "Segatella.s__Segatella_copri"
      genus = "Segatella"
    }else if (species =="Prevotella.s__Prevotella_stercorea"){
      species = "Leyella.s_Leyella_stercorea"
      genus = "Leyella"
    } else if  (species =="Prevotella.s__Prevotella_timonensis"){
      species = "Hoylesella.s__Hoylesella_timonensis"
      genus = "Hoylesella"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Prevotella loops
}

################## ############# ########################## ########################## 
#  UPDATE LACTOBACILLUS NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lactobacillus"){
    species = dmelt[row,2] 
    if (species =="Lactobacillus.s__Lactobacillus_casei_paracasei"){
      species = "Lacticaseibacillus.s__Lacticaseibacillus_paracasei" # according to chatGPT
      genus = "Lacticaseibacillus"
      #  } else if (species =="Lactobacillus.s__Lactobacillus_delbrueckii"){ # no change 
    }else if (species =="Lactobacillus.s__Lactobacillus_plantarum"){
      species = "Lactiplantibacillus.s_Lactiplantibacillus_plantarum"
      genus = "Lactiplantibacillus"
    }else if (species =="Lactobacillus.s__Lactobacillus_ruminis"){
      species = "Ligilactobacillus.s_Ligilactobacillus_ruminis"
      genus = "Ligilactobacillus"
    } else if  (species =="Lactobacillus.s__Lactobacillus_salivarius"){
      species = "Ligilactobacillus.s__Ligilactobacillus_salivarius"
      genus = "Ligilactobacillus"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Lactobacillus loops
}

################## ############# ########################## ########################## 

## and then create a new table, where each sample has a summed entry for each genus
dmelt_glom = data.frame(matrix(ncol=5,nrow=1))
names(dmelt_glom) = names(dmelt)[1:5]
genera = unique(dmelt$Genus)

## then for each sample, for each genus:
samples = unique(dmelt$Sample)
for (samp in samples){
  samp_dmelt = subset(dmelt,dmelt$Sample==samp)
  tribe = unique(samp_dmelt$Tribe)
  region = unique(samp_dmelt$Region)
  
  for (genus in genera){
    g_samp_dmelt = subset(samp_dmelt, samp_dmelt$Genus ==genus)
    genus_total = sum(g_samp_dmelt$Pathway_relab)
    
    # save a new row
    newrow = c(samp,genus,genus_total,tribe,region)
    dmelt_glom=rbind(dmelt_glom,newrow)
  }
}

dmelt_glom = na.omit(dmelt_glom)



##########################################################
#### SELECT TOP TAXA (for legend) FOR THIS ONE PATHWAY        
##########################################################

dmelt_saved = dmelt
dmelt = dmelt_glom
dmelt$Pathway_relab = as.numeric(dmelt$Pathway_relab)

# calculate total abundance per taxon and tribe across  fecal samples
totals <- dmelt %>% dplyr::group_by(Taxon, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Pathway_relab))

# select most abundant taxa to include in legend
common_taxa<- subset(totals,totals$mean_Rel_ab >= 2.0e-5)  
length(unique(common_taxa$Taxon)) # AIM FOR 8 (including unclassified) 
top_taxa = unique(common_taxa$Taxon)
top_taxa

# get abundances for only these genera
dmelt_top = subset(dmelt,dmelt$Taxon%in%top_taxa)

######################
## ADD 'OTHER' COLUMN (genera that are not abundant enough to get their own label)

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person_top <- dmelt_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

#  the total rel_ab for this pathway per person
labeled_per_person_all <- dmelt %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Then define other as 1 minus that sum
labeled_per_person_top$other = labeled_per_person_all$percent_with_taxon_label - labeled_per_person_top$percent_with_taxon_label

# Add 'other' to dmelt_top in the form of one row per person, where 'genus' = 'other'
dplot <- dmelt_top
dplot$Taxon = as.character(dplot$Taxon)
for (i in (unique(dmelt_top$Sample))){
  tribe = unique(subset(dmelt_top,dmelt_top$Sample==i)$Tribe)
  region = unique(subset(dmelt_top,dmelt_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person_top,labeled_per_person_top$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}
unique(dplot$Taxon)
##### at this point, unclassified and other are both there


# order genera by mean abundance
# Common taxa has mean_rel_ab of each taxon by tribe, but DOES NOT INCLUDE 'OTHER'
common_taxa$Taxon = as.character(common_taxa$Taxon) ## THIS does not include 'other'
t <- common_taxa %>%
  group_by(Taxon)  %>%
  summarise(mean_Rel_ab = mean(mean_Rel_ab))
t <- t[order(t$mean_Rel_ab,decreasing=TRUE),]
top_taxa_ordered = t$Taxon
top_taxa_ordered 

# then edit the order so that red (TH-enriched) taxa are on top
galactose_metabolic_process = c("Bifidobacterium","Catenibacterium","Collinsella","Faecalibacterium","Segatella","Agathobacter","Escherichia","other","unclassified")
galactose_metabolic_process  %in% dplot$Taxon
length(galactose_metabolic_process) == length(unique(dplot$Taxon))
dplot$Taxon = factor(dplot$Taxon,levels=galactose_metabolic_process)

# check if any top taxa in this pathway are NOT in the color table 
for (i in galactose_metabolic_process){
  if (i %in% genus_vec == FALSE){
    print(i)
  }
}

## adjust class
dplot$Pathway_relab = as.numeric(dplot$Pathway_relab) 
dplot$Tribe <- factor(dplot$Tribe, levels = c("Boto","Balti","Brokpa","Purigpa","Warli", "Gond","Madia","Kabui"))


## subset and order colors for this particular function / set of genera
specific_colortable <- subset(colortable,colortable$genus %in% galactose_metabolic_process)
specific_colortable = specific_colortable[match(galactose_metabolic_process, specific_colortable$genus),]

col_list = specific_colortable$color

## PLOT ! FACET WRAP ALL TRIBES TOGETHER
### TO AUTOMATICALLY CALCULATE THE DISTANCE BETWEEN DATA AND TOP OF PLOT:
temp <- dplot %>% group_by(Sample) %>% summarise(sum_relab = sum(Pathway_relab))
margin = max(temp$sum_relab)*1.05

s <- ggplot(data=dplot, aes(x=reorder(Sample, -Pathway_relab, function(x){ sum(x) }), y=Pathway_relab, fill=Taxon)) +
  geom_col() + theme_classic() +
  ylab("Relative abundance") + theme_cowplot() +
  scale_fill_manual("legend",values = col_list) + 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.35, "cm")) +
  ggtitle("Galactose metabolic process (GO:0006012)") +
  theme(plot.title = element_text(size = 12,face="bold")) +
  theme(axis.title=element_text(size=12,face="plain")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, margin),labels =  ~ sprintf(fmt = "%2.e", .)) + # removes space between 0 and x-axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

L <- s + facet_wrap(~Tribe,scales="free_x",nrow=1) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
L

ggsave("Fig5C.pdf", plot=L,width=4.64*2.2,height=0.95*2.2,units="in")


################################################################################
####  pull out one pathway at a time 
### PANEL E: # glucose catabolic process
################################################################################

func <- d[d$GeneFamily %like% "GO:0006007", ] # glucose catabolic process 

#### transpose this one stratified pathway  so rows are samples
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,] # drop first row (pathway names) -- add back in a few lines 
# for some reason all then numbers are characters now, so convert to numeric 
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 

### simplify these GO names -- not to number, b/c they're all the same; but to organism
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names

# melt df to have columns: c(Taxon,Rel_ab, Sample, Tribe)
df$Sample = row.names(df)
dmelt <-  reshape::melt(df, id = c("Sample"))
dmelt <- subset(dmelt,dmelt$variable != "total") # remove 'total' but leave 'unclassified'
dmelt$Tribe = NA; dmelt$Region = NA 

for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  dmelt[row,4] = tribe
  dmelt[row,5] = region
}

names(dmelt) <- c("Sample","Taxon","Pathway_relab", "Tribe", "Region") 


##########################################################
#### Combine individual species into genera ######## 
##########################################################

# Might actually be easiest just to add a 'genus' column to dmelt
dmelt$Taxon = as.character(dmelt$Taxon)
dmelt$Genus = NA
for (row in 1:nrow(dmelt)){
  species = dmelt[row,2]
  genus = strsplit(species, '.s__')[[1]][1]
  dmelt[row,6] = genus
}


################## ############# ########################## ########################## 
#  UPDATE EUBACTERIUM NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  species = dmelt[row,2]
  
  if (genus == "Erysipelotrichaceae_noname"){
    if (species =="Erysipelotrichaceae_noname.s__Eubacterium_biforme"){
      species = "Holdemanella.s__Holdemanella_biformis"
      genus = "Holdemanella" } 
  } else if (genus == "Eubacterium") {
    if (species =="Eubacterium.s__Eubacterium_eligens"){
      species = "Lachnospira.s__Lachnospira_eligens"
      genus = "Lachnospira"
    } else if (species =="Eubacterium.s__Eubacterium_hallii"){
      species = "Anaerobutyricum.s__Anaerobutyricum_hallii"
      genus = "Anaerobutyricum" 
      #  } else if (species =="Eubacterium.s__Eubacterium_limosum"){ # no change in NCBI Taxonomy
      #  } else if (species =="Eubacterium.s__Eubacterium_ramulus"){  # no change in NCBI Taxonomy
    } else if (species =="Eubacterium.s__Eubacterium_rectale"){
      species = "Agathobacter.s__Agathobacter_rectalis"
      genus = "Agathobacter" 
    } else if (species =="Eubacterium.s__Eubacterium_siraeum"){ # new genus name pending. Make sure these special characters aren't an issue ... works at this stage of table, anyway 
      species = "[Eubacterium].s__[Eubacterium]_siraeum"
      genus = "[Eubacterium]"
      #    } else if (species =="Eubacterium.s__Eubacterium_ventriosum"){ # no change in NCBI Taxonomy
    } 
  } # end if genus == Eubacterium
  
  ## ACTUALLY UPDATE SPECIES/GENUS IN DF
  dmelt[row,2] = species
  dmelt[row,6] = genus
} # end whole loop


################## ############# ########################## ########################## 
#  UPDATE PREVOTELLA NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Prevotella"){
    species = dmelt[row,2] 
    if (species =="Prevotella.s__Prevotella_buccae"){
      species = "Segatella.s__Segatella_buccae" # according to NCBI taxonomy
      genus = "Segatella"
    } else if (species =="Prevotella.s__Prevotella_copri"){
      species = "Segatella.s__Segatella_copri"
      genus = "Segatella"
    }else if (species =="Prevotella.s__Prevotella_stercorea"){
      species = "Leyella.s_Leyella_stercorea"
      genus = "Leyella"
    } else if  (species =="Prevotella.s__Prevotella_timonensis"){
      species = "Hoylesella.s__Hoylesella_timonensis"
      genus = "Hoylesella"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Prevotella loops
}

################## ############# ########################## ########################## 
#  UPDATE LACTOBACILLUS NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lactobacillus"){
    species = dmelt[row,2] 
    if (species =="Lactobacillus.s__Lactobacillus_casei_paracasei"){
      species = "Lacticaseibacillus.s__Lacticaseibacillus_paracasei" # according to chatGPT
      genus = "Lacticaseibacillus"
      #  } else if (species =="Lactobacillus.s__Lactobacillus_delbrueckii"){ # no change 
    }else if (species =="Lactobacillus.s__Lactobacillus_plantarum"){
      species = "Lactiplantibacillus.s_Lactiplantibacillus_plantarum"
      genus = "Lactiplantibacillus"
    }else if (species =="Lactobacillus.s__Lactobacillus_ruminis"){
      species = "Ligilactobacillus.s_Ligilactobacillus_ruminis"
      genus = "Ligilactobacillus"
    } else if  (species =="Lactobacillus.s__Lactobacillus_salivarius"){
      species = "Ligilactobacillus.s__Ligilactobacillus_salivarius"
      genus = "Ligilactobacillus"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Lactobacillus loops
}

################## ############# ########################## ########################## 

## and then create a new table, where each sample has a summed entry for each genus
dmelt_glom = data.frame(matrix(ncol=5,nrow=1))
names(dmelt_glom) = names(dmelt)[1:5]
genera = unique(dmelt$Genus)

## then for each sample, for each genus:
samples = unique(dmelt$Sample)
for (samp in samples){
  samp_dmelt = subset(dmelt,dmelt$Sample==samp)
  tribe = unique(samp_dmelt$Tribe)
  region = unique(samp_dmelt$Region)
  
  for (genus in genera){
    g_samp_dmelt = subset(samp_dmelt, samp_dmelt$Genus ==genus)
    genus_total = sum(g_samp_dmelt$Pathway_relab)
    
    # save a new row
    newrow = c(samp,genus,genus_total,tribe,region)
    dmelt_glom=rbind(dmelt_glom,newrow)
  }
}

dmelt_glom = na.omit(dmelt_glom)



##########################################################
#### SELECT TOP TAXA (for legend) FOR THIS ONE PATHWAY        
##########################################################

dmelt_saved = dmelt
dmelt = dmelt_glom
dmelt$Pathway_relab = as.numeric(dmelt$Pathway_relab)

# calculate total abundance per taxon and tribe across  fecal samples
totals <- dmelt %>% dplyr::group_by(Taxon, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Pathway_relab))

# select most abundant taxa to include in legend
common_taxa<- subset(totals,totals$mean_Rel_ab >= 5.0e-6)  
length(unique(common_taxa$Taxon)) # AIM FOR 8 (including unclassified) 
top_taxa = unique(common_taxa$Taxon)
top_taxa

# get abundances for only these genera
dmelt_top = subset(dmelt,dmelt$Taxon%in%top_taxa)

######################
## ADD 'OTHER' COLUMN (genera that are not abundant enough to get their own label)

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person_top <- dmelt_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

#  the total rel_ab for this pathway per person
labeled_per_person_all <- dmelt %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Then define other as 1 minus that sum
labeled_per_person_top$other = labeled_per_person_all$percent_with_taxon_label - labeled_per_person_top$percent_with_taxon_label

# Add 'other' to dmelt_top in the form of one row per person, where 'genus' = 'other'
dplot <- dmelt_top
dplot$Taxon = as.character(dplot$Taxon)
for (i in (unique(dmelt_top$Sample))){
  tribe = unique(subset(dmelt_top,dmelt_top$Sample==i)$Tribe)
  region = unique(subset(dmelt_top,dmelt_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person_top,labeled_per_person_top$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}
unique(dplot$Taxon)
##### at this point, unclassified and other are both there


# order genera by mean abundance
# Common taxa has mean_rel_ab of each taxon by tribe, but DOES NOT INCLUDE 'OTHER'
common_taxa$Taxon = as.character(common_taxa$Taxon) ## THIS does not include 'other'
t <- common_taxa %>%
  group_by(Taxon)  %>%
  summarise(mean_Rel_ab = mean(mean_Rel_ab))
t <- t[order(t$mean_Rel_ab,decreasing=TRUE),]
top_taxa_ordered = t$Taxon
top_taxa_ordered 

# then edit the order so that red (TH-enriched) taxa are on top
glucose_catabolic_process = c("Catenibacterium","Faecalibacterium","Agathobacter","Leyella","Alistipes","Escherichia","Treponema","other","unclassified"   )
glucose_catabolic_process %in% dplot$Taxon                                    
dplot$Taxon = factor(dplot$Taxon,levels=glucose_catabolic_process)

# check if any top taxa in this pathway are NOT in the color table 
for (i in glucose_catabolic_process){
  if (i %in% genus_vec == FALSE){
    print(i)
  }
}

## adjust class
dplot$Pathway_relab = as.numeric(dplot$Pathway_relab) 
dplot$Tribe <- factor(dplot$Tribe, levels = c("Boto","Balti","Brokpa","Purigpa","Warli", "Gond","Madia","Kabui"))


## subset and order colors for this particular function / set of genera
specific_colortable <- subset(colortable,colortable$genus %in% glucose_catabolic_process)
specific_colortable = specific_colortable[match(glucose_catabolic_process, specific_colortable$genus),]

col_list = specific_colortable$color

## PLOT ! FACET WRAP ALL TRIBES TOGETHER
### TO AUTOMATICALLY CALCULATE THE DISTANCE BETWEEN DATA AND TOP OF PLOT:
temp <- dplot %>% group_by(Sample) %>% summarise(sum_relab = sum(Pathway_relab))
margin = max(temp$sum_relab)*1.05

s <- ggplot(data=dplot, aes(x=reorder(Sample, -Pathway_relab, function(x){ sum(x) }), y=Pathway_relab, fill=Taxon)) +
  geom_col() + theme_classic() +
  ylab("Relative abundance") + theme_cowplot() +
  scale_fill_manual("legend",values = col_list) + 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.35, "cm")) +
  ggtitle("Glucose catabolic process (GO:0006007)")+
  theme(plot.title = element_text(size = 12,face="bold")) +
  theme(axis.title=element_text(size=12,face="plain")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, margin),labels =  ~ sprintf(fmt = "%2.e", .)) + # removes space between 0 and x-axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

L <- s + facet_wrap(~Tribe,scales="free_x",nrow=1) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
L

ggsave("Fig5E.pdf", plot=L,width=4.64*2.2,height=0.95*2.2,units="in")



################################################################################
####  pull out one pathway at a time 
### SUPP 1 PANEL C: # lactose metabolic process
################################################################################

func <- d[d$GeneFamily %like% "GO:0005988", ] # lactose metabolic process

#### transpose this one stratified pathway  so rows are samples
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,] # drop first row (pathway names) -- add back in a few lines 
# for some reason all then numbers are characters now, so convert to numeric 
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 

### simplify these GO names -- not to number, b/c they're all the same; but to organism
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names

# melt df to have columns: c(Taxon,Rel_ab, Sample, Tribe)
df$Sample = row.names(df)
dmelt <-  reshape::melt(df, id = c("Sample"))
dmelt <- subset(dmelt,dmelt$variable != "total") # remove 'total' but leave 'unclassified'
dmelt$Tribe = NA; dmelt$Region = NA 

for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  dmelt[row,4] = tribe
  dmelt[row,5] = region
}

names(dmelt) <- c("Sample","Taxon","Pathway_relab", "Tribe", "Region") 


##########################################################
#### Combine individual species into genera ######## 
##########################################################

dmelt$Taxon = as.character(dmelt$Taxon)
dmelt$Genus = NA
for (row in 1:nrow(dmelt)){
  species = dmelt[row,2]
  genus = strsplit(species, '.s__')[[1]][1]
  dmelt[row,6] = genus
}


################## ############# ########################## ########################## 
#  UPDATE EUBACTERIUM NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  species = dmelt[row,2]
  
  if (genus == "Erysipelotrichaceae_noname"){
    if (species =="Erysipelotrichaceae_noname.s__Eubacterium_biforme"){
      species = "Holdemanella.s__Holdemanella_biformis"
      genus = "Holdemanella" } 
  } else if (genus == "Eubacterium") {
    if (species =="Eubacterium.s__Eubacterium_eligens"){
      species = "Lachnospira.s__Lachnospira_eligens"
      genus = "Lachnospira"
    } else if (species =="Eubacterium.s__Eubacterium_hallii"){
      species = "Anaerobutyricum.s__Anaerobutyricum_hallii"
      genus = "Anaerobutyricum" 
      #  } else if (species =="Eubacterium.s__Eubacterium_limosum"){ # no change in NCBI Taxonomy
      #  } else if (species =="Eubacterium.s__Eubacterium_ramulus"){  # no change in NCBI Taxonomy
    } else if (species =="Eubacterium.s__Eubacterium_rectale"){
      species = "Agathobacter.s__Agathobacter_rectalis"
      genus = "Agathobacter" 
    } else if (species =="Eubacterium.s__Eubacterium_siraeum"){ # new genus name pending. Make sure these special characters aren't an issue ... works at this stage of table, anyway 
      species = "[Eubacterium].s__[Eubacterium]_siraeum"
      genus = "[Eubacterium]"
      #    } else if (species =="Eubacterium.s__Eubacterium_ventriosum"){ # no change in NCBI Taxonomy
    } 
  } # end if genus == Eubacterium
  
  ## ACTUALLY UPDATE SPECIES/GENUS IN DF
  dmelt[row,2] = species
  dmelt[row,6] = genus
} # end whole loop


################## ############# ########################## ########################## 
#  UPDATE PREVOTELLA NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Prevotella"){
    species = dmelt[row,2] 
    if (species =="Prevotella.s__Prevotella_buccae"){
      species = "Segatella.s__Segatella_buccae" # according to NCBI taxonomy
      genus = "Segatella"
    } else if (species =="Prevotella.s__Prevotella_copri"){
      species = "Segatella.s__Segatella_copri"
      genus = "Segatella"
    }else if (species =="Prevotella.s__Prevotella_stercorea"){
      species = "Leyella.s_Leyella_stercorea"
      genus = "Leyella"
    } else if  (species =="Prevotella.s__Prevotella_timonensis"){
      species = "Hoylesella.s__Hoylesella_timonensis"
      genus = "Hoylesella"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Prevotella loops
}

################## ############# ########################## ########################## 
#  UPDATE LACTOBACILLUS NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lactobacillus"){
    species = dmelt[row,2] 
    if (species =="Lactobacillus.s__Lactobacillus_casei_paracasei"){
      species = "Lacticaseibacillus.s__Lacticaseibacillus_paracasei" # according to chatGPT
      genus = "Lacticaseibacillus"
      #  } else if (species =="Lactobacillus.s__Lactobacillus_delbrueckii"){ # no change 
    }else if (species =="Lactobacillus.s__Lactobacillus_plantarum"){
      species = "Lactiplantibacillus.s_Lactiplantibacillus_plantarum"
      genus = "Lactiplantibacillus"
    }else if (species =="Lactobacillus.s__Lactobacillus_ruminis"){
      species = "Ligilactobacillus.s_Ligilactobacillus_ruminis"
      genus = "Ligilactobacillus"
    } else if  (species =="Lactobacillus.s__Lactobacillus_salivarius"){
      species = "Ligilactobacillus.s__Ligilactobacillus_salivarius"
      genus = "Ligilactobacillus"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Lactobacillus loops
}

################## ############# ########################## ########################## 

## and then create a new table, where each sample has a summed entry for each genus
dmelt_glom = data.frame(matrix(ncol=5,nrow=1))
names(dmelt_glom) = names(dmelt)[1:5]
genera = unique(dmelt$Genus)

## then for each sample, for each genus:
samples = unique(dmelt$Sample)
for (samp in samples){
  samp_dmelt = subset(dmelt,dmelt$Sample==samp)
  tribe = unique(samp_dmelt$Tribe)
  region = unique(samp_dmelt$Region)
  
  for (genus in genera){
    g_samp_dmelt = subset(samp_dmelt, samp_dmelt$Genus ==genus)
    genus_total = sum(g_samp_dmelt$Pathway_relab)
    
    # save a new row
    newrow = c(samp,genus,genus_total,tribe,region)
    dmelt_glom=rbind(dmelt_glom,newrow)
  }
}

dmelt_glom = na.omit(dmelt_glom)



##########################################################
#### SELECT TOP TAXA (for legend) FOR THIS ONE PATHWAY        
##########################################################

dmelt_saved = dmelt
dmelt = dmelt_glom
dmelt$Pathway_relab = as.numeric(dmelt$Pathway_relab)

# calculate total abundance per taxon and tribe across  fecal samples
totals <- dmelt %>% dplyr::group_by(Taxon, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Pathway_relab))

# select most abundant taxa to include in legend
common_taxa<- subset(totals,totals$mean_Rel_ab >= 3.4e-6)  
length(unique(common_taxa$Taxon)) # AIM FOR 8 (including unclassified) 
top_taxa = unique(common_taxa$Taxon)
top_taxa

# get abundances for only these genera
dmelt_top = subset(dmelt,dmelt$Taxon%in%top_taxa)

######################
## ADD 'OTHER' COLUMN (genera that are not abundant enough to get their own label)

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person_top <- dmelt_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

#  the total rel_ab for this pathway per person
labeled_per_person_all <- dmelt %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Then define other as 1 minus that sum
labeled_per_person_top$other = labeled_per_person_all$percent_with_taxon_label - labeled_per_person_top$percent_with_taxon_label

# Add 'other' to dmelt_top in the form of one row per person, where 'genus' = 'other'
dplot <- dmelt_top
dplot$Taxon = as.character(dplot$Taxon)
for (i in (unique(dmelt_top$Sample))){
  tribe = unique(subset(dmelt_top,dmelt_top$Sample==i)$Tribe)
  region = unique(subset(dmelt_top,dmelt_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person_top,labeled_per_person_top$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}
unique(dplot$Taxon)
##### at this point, unclassified and other are both there


# order genera by mean abundance
# Common taxa has mean_rel_ab of each taxon by tribe, but DOES NOT INCLUDE 'OTHER'
common_taxa$Taxon = as.character(common_taxa$Taxon) ## THIS does not include 'other'
t <- common_taxa %>%
  group_by(Taxon)  %>%
  summarise(mean_Rel_ab = mean(mean_Rel_ab))
t <- t[order(t$mean_Rel_ab,decreasing=TRUE),]
top_taxa_ordered = t$Taxon
top_taxa_ordered 

# then edit the order so that red (TH-enriched) taxa are on top
lactose_metabolic = c("Collinsella", "Catenibacterium" ,"Ligilactobacillus", "Streptococcus","Faecalibacterium","Agathobacter","Butyrivibrio", "other", "unclassified")
lactose_metabolic %in% dplot$Taxon
length(lactose_metabolic) == length(unique(dplot$Taxon))
dplot$Taxon = factor(dplot$Taxon,levels=lactose_metabolic)

# check if any top taxa in this pathway are NOT in the color table 
for (i in lactose_metabolic){
  if (i %in% genus_vec == FALSE){
    print(i)
  }
}

## adjust class
dplot$Pathway_relab = as.numeric(dplot$Pathway_relab) 
dplot$Tribe <- factor(dplot$Tribe, levels = c("Boto","Balti","Brokpa","Purigpa","Warli", "Gond","Madia","Kabui"))


## subset and order colors for this particular function / set of genera
specific_colortable <- subset(colortable,colortable$genus %in% lactose_metabolic)
specific_colortable = specific_colortable[match(lactose_metabolic, specific_colortable$genus),]
col_list = specific_colortable$color

## PLOT ! FACET WRAP ALL TRIBES TOGETHER
### TO AUTOMATICALLY CALCULATE THE DISTANCE BETWEEN DATA AND TOP OF PLOT:
temp <- dplot %>% group_by(Sample) %>% summarise(sum_relab = sum(Pathway_relab))
margin = max(temp$sum_relab)*1.05

s <- ggplot(data=dplot, aes(x=reorder(Sample, -Pathway_relab, function(x){ sum(x) }), y=Pathway_relab, fill=Taxon)) +
  geom_col() + theme_classic() +
  ylab("Relative abundance") + theme_cowplot() +
  scale_fill_manual("legend",values = col_list) + 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.35, "cm")) +
  ggtitle("Lactose metabolic process (GO:0005988)")+
  theme(plot.title = element_text(size = 12,face="bold")) +
  theme(axis.title=element_text(size=12,face="plain")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, margin),labels =  ~ sprintf(fmt = "%2.e", .)) + # removes space between 0 and x-axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

L <- s + facet_wrap(~Tribe,scales="free_x",nrow=1) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
L

ggsave("Fig5_S1C.pdf", plot=L,width=4.64*2.2,height=0.95*2.2,units="in")




################################################################################
####  pull out one pathway at a time 
### SUPP 1 PANEL D: # glucose transmembrane transporter
################################################################################

func <- d[d$GeneFamily %like% "GO:0005355", ] # glucose transmembrane transporter 

#### transpose this one stratified pathway  so rows are samples
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,] # drop first row (pathway names) -- add back in a few lines 
# for some reason all then numbers are characters now, so convert to numeric 
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 

### simplify these GO names -- not to number, b/c they're all the same; but to organism
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names

# melt df to have columns: c(Taxon,Rel_ab, Sample, Tribe)
df$Sample = row.names(df)
dmelt <-  reshape::melt(df, id = c("Sample"))
dmelt <- subset(dmelt,dmelt$variable != "total") # remove 'total' but leave 'unclassified'
dmelt$Tribe = NA; dmelt$Region = NA 

for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  dmelt[row,4] = tribe
  dmelt[row,5] = region
}

names(dmelt) <- c("Sample","Taxon","Pathway_relab", "Tribe", "Region") 


##########################################################
#### Combine individual species into genera ######## 
##########################################################

dmelt$Taxon = as.character(dmelt$Taxon)
dmelt$Genus = NA
for (row in 1:nrow(dmelt)){
  species = dmelt[row,2]
  genus = strsplit(species, '.s__')[[1]][1]
  dmelt[row,6] = genus
}


################## ############# ########################## ########################## 
#  UPDATE EUBACTERIUM NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  species = dmelt[row,2]
  
  if (genus == "Erysipelotrichaceae_noname"){
    if (species =="Erysipelotrichaceae_noname.s__Eubacterium_biforme"){
      species = "Holdemanella.s__Holdemanella_biformis"
      genus = "Holdemanella" } 
  } else if (genus == "Eubacterium") {
    if (species =="Eubacterium.s__Eubacterium_eligens"){
      species = "Lachnospira.s__Lachnospira_eligens"
      genus = "Lachnospira"
    } else if (species =="Eubacterium.s__Eubacterium_hallii"){
      species = "Anaerobutyricum.s__Anaerobutyricum_hallii"
      genus = "Anaerobutyricum" 
      #  } else if (species =="Eubacterium.s__Eubacterium_limosum"){ # no change in NCBI Taxonomy
      #  } else if (species =="Eubacterium.s__Eubacterium_ramulus"){  # no change in NCBI Taxonomy
    } else if (species =="Eubacterium.s__Eubacterium_rectale"){
      species = "Agathobacter.s__Agathobacter_rectalis"
      genus = "Agathobacter" 
    } else if (species =="Eubacterium.s__Eubacterium_siraeum"){ # new genus name pending. Make sure these special characters aren't an issue ... works at this stage of table, anyway 
      species = "[Eubacterium].s__[Eubacterium]_siraeum"
      genus = "[Eubacterium]"
      #    } else if (species =="Eubacterium.s__Eubacterium_ventriosum"){ # no change in NCBI Taxonomy
    } 
  } # end if genus == Eubacterium
  
  ## ACTUALLY UPDATE SPECIES/GENUS IN DF
  dmelt[row,2] = species
  dmelt[row,6] = genus
} # end whole loop


################## ############# ########################## ########################## 
#  UPDATE PREVOTELLA NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Prevotella"){
    species = dmelt[row,2] 
    if (species =="Prevotella.s__Prevotella_buccae"){
      species = "Segatella.s__Segatella_buccae" # according to NCBI taxonomy
      genus = "Segatella"
    } else if (species =="Prevotella.s__Prevotella_copri"){
      species = "Segatella.s__Segatella_copri"
      genus = "Segatella"
    }else if (species =="Prevotella.s__Prevotella_stercorea"){
      species = "Leyella.s_Leyella_stercorea"
      genus = "Leyella"
    } else if  (species =="Prevotella.s__Prevotella_timonensis"){
      species = "Hoylesella.s__Hoylesella_timonensis"
      genus = "Hoylesella"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Prevotella loops
}

################## ############# ########################## ########################## 
#  UPDATE LACTOBACILLUS NAMES
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lactobacillus"){
    species = dmelt[row,2] 
    if (species =="Lactobacillus.s__Lactobacillus_casei_paracasei"){
      species = "Lacticaseibacillus.s__Lacticaseibacillus_paracasei" # according to chatGPT
      genus = "Lacticaseibacillus"
      #  } else if (species =="Lactobacillus.s__Lactobacillus_delbrueckii"){ # no change 
    }else if (species =="Lactobacillus.s__Lactobacillus_plantarum"){
      species = "Lactiplantibacillus.s_Lactiplantibacillus_plantarum"
      genus = "Lactiplantibacillus"
    }else if (species =="Lactobacillus.s__Lactobacillus_ruminis"){
      species = "Ligilactobacillus.s_Ligilactobacillus_ruminis"
      genus = "Ligilactobacillus"
    } else if  (species =="Lactobacillus.s__Lactobacillus_salivarius"){
      species = "Ligilactobacillus.s__Ligilactobacillus_salivarius"
      genus = "Ligilactobacillus"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Lactobacillus loops
}

################## ############# ########################## ########################## 

## and then create a new table, where each sample has a summed entry for each genus
dmelt_glom = data.frame(matrix(ncol=5,nrow=1))
names(dmelt_glom) = names(dmelt)[1:5]
genera = unique(dmelt$Genus)

## then for each sample, for each genus:
samples = unique(dmelt$Sample)
for (samp in samples){
  samp_dmelt = subset(dmelt,dmelt$Sample==samp)
  tribe = unique(samp_dmelt$Tribe)
  region = unique(samp_dmelt$Region)
  
  for (genus in genera){
    g_samp_dmelt = subset(samp_dmelt, samp_dmelt$Genus ==genus)
    genus_total = sum(g_samp_dmelt$Pathway_relab)
    
    # save a new row
    newrow = c(samp,genus,genus_total,tribe,region)
    dmelt_glom=rbind(dmelt_glom,newrow)
  }
}

dmelt_glom = na.omit(dmelt_glom)



##########################################################
#### SELECT TOP TAXA (for legend) FOR THIS ONE PATHWAY        
##########################################################

dmelt = dmelt_glom
dmelt$Pathway_relab = as.numeric(dmelt$Pathway_relab)

# calculate total abundance per taxon and tribe across  fecal samples
totals <- dmelt %>% dplyr::group_by(Taxon, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Pathway_relab))

# select most abundant taxa to include in legend
common_taxa<- subset(totals,totals$mean_Rel_ab >= 3e-7)  
length(unique(common_taxa$Taxon)) # AIM FOR 8 (including unclassified) 
top_taxa = unique(common_taxa$Taxon)
top_taxa

# get abundances for only these genera
dmelt_top = subset(dmelt,dmelt$Taxon%in%top_taxa)

######################
## ADD 'OTHER' COLUMN (genera that are not abundant enough to get their own label)

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person_top <- dmelt_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

#  the total rel_ab for this pathway per person
labeled_per_person_all <- dmelt %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Then define other as 1 minus that sum
labeled_per_person_top$other = labeled_per_person_all$percent_with_taxon_label - labeled_per_person_top$percent_with_taxon_label

# Add 'other' to dmelt_top in the form of one row per person, where 'genus' = 'other'
dplot <- dmelt_top
dplot$Taxon = as.character(dplot$Taxon)
for (i in (unique(dmelt_top$Sample))){
  tribe = unique(subset(dmelt_top,dmelt_top$Sample==i)$Tribe)
  region = unique(subset(dmelt_top,dmelt_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person_top,labeled_per_person_top$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}
unique(dplot$Taxon)
##### at this point, unclassified and other are both there


# order genera by mean abundance
# Common taxa has mean_rel_ab of each taxon by tribe, but DOES NOT INCLUDE 'OTHER'
common_taxa$Taxon = as.character(common_taxa$Taxon) ## THIS does not include 'other'
t <- common_taxa %>%
  group_by(Taxon)  %>%
  summarise(mean_Rel_ab = mean(mean_Rel_ab))
t <- t[order(t$mean_Rel_ab,decreasing=TRUE),]
top_taxa_ordered = t$Taxon
top_taxa_ordered 

# then edit the order so that red (TH-enriched) taxa are on top (NA for this function)
glucose_transporter = c("Escherichia","Bacteroides","Klebsiella","Enterobacter","Alistipes", "other","unclassified") # more minor ones aren't visible
glucose_transporter %in% dplot$Taxon                                    
dplot$Taxon = factor(dplot$Taxon,levels=glucose_transporter)

# check if any top taxa in this pathway are NOT in the color table 
for (i in glucose_transporter){
  if (i %in% genus_vec == FALSE){
    print(i)
  }
}

## adjust class
dplot$Pathway_relab = as.numeric(dplot$Pathway_relab) 
dplot$Tribe <- factor(dplot$Tribe, levels = c("Boto","Balti","Brokpa","Purigpa","Warli", "Gond","Madia","Kabui"))


## subset and order colors for this particular function / set of genera
specific_colortable <- subset(colortable,colortable$genus %in% glucose_transporter)
specific_colortable = specific_colortable[match(glucose_transporter, specific_colortable$genus),]
col_list = specific_colortable$color

## PLOT ! FACET WRAP ALL TRIBES TOGETHER
### TO AUTOMATICALLY CALCULATE THE DISTANCE BETWEEN DATA AND TOP OF PLOT:
temp <- dplot %>% group_by(Sample) %>% summarise(sum_relab = sum(Pathway_relab))
margin = max(temp$sum_relab)*1.05

s <- ggplot(data=dplot, aes(x=reorder(Sample, -Pathway_relab, function(x){ sum(x) }), y=Pathway_relab, fill=Taxon)) +
  geom_col() + theme_classic() +
  ylab("Relative abundance") + theme_cowplot() +
  scale_fill_manual("legend",values = col_list) + 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.35, "cm")) +
  ggtitle("Glucose transmembrane transporter (GO:0005355)")+
  theme(plot.title = element_text(size = 12,face="bold")) +
  theme(axis.title=element_text(size=12,face="plain")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, margin),labels =  ~ sprintf(fmt = "%2.e", .)) + # removes space between 0 and x-axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

L <- s + facet_wrap(~Tribe,scales="free_x",nrow=1) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
L

ggsave("Fig5_S1D.pdf", plot=L,width=4.64*2.2,height=0.95*2.2,units="in")





########################################################################################



###########################################################################################################################
#### METACYC functional abundance data    ## manually remove the '#' from the header first!! 
##############################################################################################################################

# the total pathway here != constituent parts , unlike GO

d<-read.delim("metaCyc_100.tsv", header = TRUE, sep = "\t",numerals="allow.loss",colClasses=c('character',rep("numeric",100)))
# simplify the column names (which are samples )
names1 <- names(d)[2:length(names(d))]
names2 <- sub('NODUP_HMN_UNMAPPED_TRIM_MARKED_', '', names1)
names3 <- sub('_Abundance', '', names2)
names(d)[2:length(names(d))] <- names3

# exclude outliers 
outliers = c("AK_SK_49","AK_SK_10","AK_SK_32","AK_SR_6","AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2","AK_SR_4")
others = c('AK_SK_24',"AK_SK_19","AK_SK_34") # the first two have very little sequencing data and were excluded from rel_ab data; the third has no metadata
d <- d[ , -which(names(d) %in% c(outliers,others))]

# use metadata to exclude enrichment samples 
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- subset(meta,meta$sample %in% names(d))
fecalkeep = as.character(subset(meta,meta$sample_type=='fecal')$sample)
d <- d[ , which(names(d) %in% c(fecalkeep,'Pathway'))] # FOR METACYC
dsamps <- names(d)[2:length(names(d))] 
meta<-subset(meta,meta$sample %in% dsamps)


################################################################################
####  pull out one pathway at a time 
### MAIN PANEL D: # leloir
################################################################################

func1 <- d[d$Pathway %like% "PWY-6317", ] # Leloir I
func2 <- d[d$Pathway %like% "PWY66-422", ] # Leloir V
# combine both together 
df0 <- as.data.frame(t(func1),stringsAsFactors = FALSE)
df_func1 = df0[-1,]
df_func1[] <- lapply(df_func1, type.convert, as.is = TRUE, dec='.')
names(df_func1) = func1$Pathway 
new_names <-str_split_fixed(names(df_func1), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified
new_names[length(new_names)]  = 'unclassified'
names(df_func1) <- new_names
origtotal <- df_func1$total
sum <- rowSums(df_func1)
sum <- sum-df_func1$total
df_func1$total<-sum
#
df0 <- as.data.frame(t(func2),stringsAsFactors = FALSE)
df_func2 = df0[-1,]
df_func2[] <- lapply(df_func2, type.convert, as.is = TRUE, dec='.')
names(df_func2) = func2$Pathway 
new_names <-str_split_fixed(names(df_func2), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified
new_names[length(new_names)]  = 'unclassified'
names(df_func2) <- new_names
origtotal <- df_func2$total
sum <- rowSums(df_func2)
sum <- sum-df_func2$total
df_func2$total<-sum
### Need to add together these 2 Leloir abundances 
df <- data.frame(matrix(ncol = 1, nrow = nrow(df_func2))); names(df)= "total"
df$total = df_func1$total + df_func2$total
row.names(df) = row.names(df_func2)
new_col_names = c("total")
for (c in 2:ncol(df_func1)){
  col_name = names(df_func1)[c]
  if (col_name %in% names(df_func2)){ # if the col is shared in both df, add together 
    new_col = df_func1[,c] + df_func2[,c]
  } else {new_col = df_func1[,c]}
  new_col_names = c(new_col_names,col_name)
  df <- cbind(df,new_col)
}
names(df) <- new_col_names
# also add cols only in df_func2
for (c in 2:ncol(df_func2)){
  col_name = names(df_func2)[c]
  if (col_name %in% names(df_func1)){ # if the col is shared with the other df, it's already taken care of 
    "do nothing"
  } else {
    new_col = df_func2[,c]
    new_col_names = c(new_col_names,col_name)
    df <- cbind(df,new_col)
  }
}
names(df) <- new_col_names



#### melt df to have columns: c(Taxon,Rel_ab, Sample, Tribe)

df$Sample = row.names(df)
dmelt <-  reshape::melt(df, id = c("Sample"))
dmelt <- subset(dmelt,dmelt$variable != "total") # remove 'total' but leave 'unclassified'
dmelt$Tribe = NA; dmelt$Region = NA 

for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  dmelt[row,4] = tribe
  dmelt[row,5] = region
}

names(dmelt) <- c("Sample","Taxon","Pathway_relab", "Tribe", "Region") 

##########################################################
#### Combine individual species into genera ######## 
##########################################################

dmelt$Taxon = as.character(dmelt$Taxon)
dmelt$Genus = NA
for (row in 1:nrow(dmelt)){
  species = dmelt[row,2]
  genus = strsplit(species, '.s__')[[1]][1]
  dmelt[row,6] = genus
}


################## ############# ########################## ########################## 
#  UPDATE EUBACTERIUM 
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  species = dmelt[row,2]
  
  if (genus == "Erysipelotrichaceae_noname"){
    if (species =="Erysipelotrichaceae_noname.s__Eubacterium_biforme"){
      species = "Holdemanella.s__Holdemanella_biformis"
      genus = "Holdemanella" } 
  } else if (genus == "Eubacterium") {
    if (species =="Eubacterium.s__Eubacterium_eligens"){
      species = "Lachnospira.s__Lachnospira_eligens"
      genus = "Lachnospira"
    } else if (species =="Eubacterium.s__Eubacterium_hallii"){
      species = "Anaerobutyricum.s__Anaerobutyricum_hallii"
      genus = "Anaerobutyricum" 
      #  } else if (species =="Eubacterium.s__Eubacterium_limosum"){ # no change in NCBI Taxonomy
      #  } else if (species =="Eubacterium.s__Eubacterium_ramulus"){  # no change in NCBI Taxonomy
    } else if (species =="Eubacterium.s__Eubacterium_rectale"){
      species = "Agathobacter.s__Agathobacter_rectalis"
      genus = "Agathobacter" 
    } else if (species =="Eubacterium.s__Eubacterium_siraeum"){ # new genus name pending. Make sure these special characters aren't an issue ... works at this stage of table, anyway 
      species = "[Eubacterium].s__[Eubacterium]_siraeum"
      genus = "[Eubacterium]"
      #    } else if (species =="Eubacterium.s__Eubacterium_ventriosum"){ # no change in NCBI Taxonomy
    } 
  } # end if genus == Eubacterium
  
  ## ACTUALLY UPDATE SPECIES/GENUS IN DF
  dmelt[row,2] = species
  dmelt[row,6] = genus
} # end whole loop


################## ############# ########################## ########################## 
#  UPDATE PREVOTELLA
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Prevotella"){
    species = dmelt[row,2] 
    if (species =="Prevotella.s__Prevotella_buccae"){
      species = "Segatella.s__Segatella_buccae" # according to NCBI taxonomy
      genus = "Segatella"
    } else if (species =="Prevotella.s__Prevotella_copri"){
      species = "Segatella.s__Segatella_copri"
      genus = "Segatella"
    }else if (species =="Prevotella.s__Prevotella_stercorea"){
      species = "Leyella.s_Leyella_stercorea"
      genus = "Leyella"
    } else if  (species =="Prevotella.s__Prevotella_timonensis"){
      species = "Hoylesella.s__Hoylesella_timonensis"
      genus = "Hoylesella"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Prevotella loops
}

################## ############# ########################## ########################## 
#  UPDATE LACTOBACILLUS
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lactobacillus"){
    species = dmelt[row,2] 
    if (species =="Lactobacillus.s__Lactobacillus_casei_paracasei"){
      species = "Lacticaseibacillus.s__Lacticaseibacillus_paracasei" # according to chatGPT
      genus = "Lacticaseibacillus"
      #  } else if (species =="Lactobacillus.s__Lactobacillus_delbrueckii"){ # no change 
    }else if (species =="Lactobacillus.s__Lactobacillus_plantarum"){
      species = "Lactiplantibacillus.s_Lactiplantibacillus_plantarum"
      genus = "Lactiplantibacillus"
    }else if (species =="Lactobacillus.s__Lactobacillus_ruminis"){
      species = "Ligilactobacillus.s_Ligilactobacillus_ruminis"
      genus = "Ligilactobacillus"
    } else if  (species =="Lactobacillus.s__Lactobacillus_salivarius"){
      species = "Ligilactobacillus.s__Ligilactobacillus_salivarius"
      genus = "Ligilactobacillus"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Lactobacillus loops
}


################## ############# ########################## ########################## 
#  UPDATE "Lachnospiraceae_bacterium_5_1_57FAA"
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lachnospiraceae_noname"){
    species = dmelt[row,2] 
    if (species =="Lachnospiraceae_noname.s__Lachnospiraceae_bacterium_5_1_57FAA"){
      #  species = "Lachnospiraceae_bacterium_5_1_57FAA" 
      species = "5_1_57FAA" 
      #  genus = "Lachnospiraceae_5_1_57FAA"
      genus = "5_1_57FAA" 
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } 
}


################## ############# ########################## ########################## 

##  create a new table, where each sample has a summed entry for each genus
dmelt_glom = data.frame(matrix(ncol=5,nrow=1))
names(dmelt_glom) = names(dmelt)[1:5]
genera = unique(dmelt$Genus)

## then for each sample, for each genus:
samples = unique(dmelt$Sample)
for (samp in samples){
  samp_dmelt = subset(dmelt,dmelt$Sample==samp)
  tribe = unique(samp_dmelt$Tribe)
  region = unique(samp_dmelt$Region)
  
  for (genus in genera){
    g_samp_dmelt = subset(samp_dmelt, samp_dmelt$Genus ==genus)
    genus_total = sum(g_samp_dmelt$Pathway_relab)
    
    # save a new row
    newrow = c(samp,genus,genus_total,tribe,region)
    dmelt_glom=rbind(dmelt_glom,newrow)
  }
}

dmelt_glom = na.omit(dmelt_glom)

dmelt_glom  <- dmelt_glom  %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))


##########################################################
#### SELECT TOP TAXA (for legend) FOR THIS ONE PATHWAY        
##########################################################

dmelt_saved = dmelt
dmelt = dmelt_glom
dmelt$Pathway_relab = as.numeric(dmelt$Pathway_relab)

# calculate total abundance per taxon and tribe across  fecal samples
totals <- dmelt %>% dplyr::group_by(Taxon, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Pathway_relab))

# select most abundant taxa to include in legend
common_taxa<- subset(totals,totals$mean_Rel_ab >= 8.0e-6)  
length(unique(common_taxa$Taxon)) # AIM FOR 8 (including unclassified) 
top_taxa = unique(common_taxa$Taxon)
top_taxa

# get abundances for only these genera
dmelt_top = subset(dmelt,dmelt$Taxon%in%top_taxa)

######################
## ADD 'OTHER' COLUMN (genera that are not abundant enough to get their own label)

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person_top <- dmelt_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Because the relabs of this pathway don't sum to 1 -- it's relab of ALL pathways --
# 'other' is not just 1-x. It's the total -x. So, calculate

#  the total rel_ab for this pathway per person
labeled_per_person_all <- dmelt %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Then define other as 1 minus that sum
labeled_per_person_top$other = labeled_per_person_all$percent_with_taxon_label - labeled_per_person_top$percent_with_taxon_label
#labeled_per_person_top$checksum = labeled_per_person_top$other + labeled_per_person_top$percent_with_taxon_label

# Add 'other' to dmelt_top in the form of one row per person, where 'genus' = 'other'
dplot <- dmelt_top
dplot$Taxon = as.character(dplot$Taxon)
for (i in (unique(dmelt_top$Sample))){
  tribe = unique(subset(dmelt_top,dmelt_top$Sample==i)$Tribe)
  region = unique(subset(dmelt_top,dmelt_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person_top,labeled_per_person_top$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}
unique(dplot$Taxon)
##### at this point, unclassified and other are both there


# order genera by mean abundance
# Common taxa has mean_rel_ab of each taxon by tribe, but DOES NOT INCLUDE 'OTHER'
common_taxa$Taxon = as.character(common_taxa$Taxon) ## THIS does not include 'other'
t <- common_taxa %>%
  group_by(Taxon)  %>%
  summarise(mean_Rel_ab = mean(mean_Rel_ab))
t <- t[order(t$mean_Rel_ab,decreasing=TRUE),]
top_taxa_ordered = t$Taxon
top_taxa_ordered 

# then edit the order so that red (TH-enriched) taxa are on top
leloir_combined_metacyc = c('Faecalibacterium',"5_1_57FAA",'Agathobacter','Fusobacterium','Treponema','Escherichia',"Blautia", 'other','unclassified')
leloir_combined_metacyc  %in% dplot$Taxon
dplot$Taxon = factor(dplot$Taxon,levels=leloir_combined_metacyc)


## adjust class
dplot$Pathway_relab = as.numeric(dplot$Pathway_relab) 
dplot$Tribe <- factor(dplot$Tribe, levels = c("Boto","Balti","Brokpa","Purigpa","Warli", "Gond","Madia","Kabui"))


## subset and order colors for this particular function / set of genera
specific_colortable <- subset(colortable,colortable$genus %in% leloir_combined_metacyc)
specific_colortable = specific_colortable[match(leloir_combined_metacyc, specific_colortable$genus),]
col_list = specific_colortable$color

## PLOT ! FACET WRAP ALL TRIBES TOGETHER
### TO AUTOMATICALLY CALCULATE THE DISTANCE BETWEEN DATA AND TOP OF PLOT:
temp <- dplot %>% group_by(Sample) %>% summarise(sum_relab = sum(Pathway_relab))
margin = max(temp$sum_relab)*1.05

s <- ggplot(data=dplot, aes(x=reorder(Sample, -Pathway_relab, function(x){ sum(x) }), y=Pathway_relab, fill=Taxon)) +
  geom_col() + theme_classic() +
  ylab("Relative abundance") + theme_cowplot() +
  scale_fill_manual("legend",values = col_list) + 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.35, "cm")) +
  ggtitle("Galactose degradation Leloir (PWY-6317+PWY-66-422)")+
  theme(plot.title = element_text(size = 12,face="bold")) +
  theme(axis.title=element_text(size=12,face="plain")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, margin),labels =  ~ sprintf(fmt = "%2.e", .)) + # removes space between 0 and x-axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

L <- s + facet_wrap(~Tribe,scales="free_x",nrow=1) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
L

ggsave("Fig5D.pdf", plot=L,width=4.64*2.2,height=0.95*2.2,units="in")





################################################################################
####  pull out one pathway at a time 
### SUPP PANEL B: # LACTOSECAT-PWY
################################################################################

func <- d[d$Pathway %like% "LACTOSECAT-PWY", ] # Lactose and galactose degradation I
#### transpose so rows are samples
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
# drop first row (pathway names) -- add back in a few lines 
df = df0[-1,]
# for some reason all then numbers are characters now, so convert to numeric 
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$Pathway 

### simplify these metacyc names  to organism
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names 
# ## FOR METACYC, CALCULATE A STRATIFIED TOTAL (since all strats+unclassified < total, b/c of pathway definitions)
origtotal <- df$total
sum <- rowSums(df)
sum <- sum-df$total
plot(origtotal,sum)
df$total<-sum


#### melt df to have columns: c(Taxon,Rel_ab, Sample, Tribe)

df$Sample = row.names(df)
dmelt <-  reshape::melt(df, id = c("Sample"))
dmelt <- subset(dmelt,dmelt$variable != "total") # remove 'total' but leave 'unclassified'
dmelt$Tribe = NA; dmelt$Region = NA 

for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  dmelt[row,4] = tribe
  dmelt[row,5] = region
}

names(dmelt) <- c("Sample","Taxon","Pathway_relab", "Tribe", "Region") 

##########################################################
#### Combine individual species into genera ######## 
##########################################################

dmelt$Taxon = as.character(dmelt$Taxon)
dmelt$Genus = NA
for (row in 1:nrow(dmelt)){
  species = dmelt[row,2]
  genus = strsplit(species, '.s__')[[1]][1]
  dmelt[row,6] = genus
}


################## ############# ########################## ########################## 
#  UPDATE EUBACTERIUM 
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  species = dmelt[row,2]
  
  if (genus == "Erysipelotrichaceae_noname"){
    if (species =="Erysipelotrichaceae_noname.s__Eubacterium_biforme"){
      species = "Holdemanella.s__Holdemanella_biformis"
      genus = "Holdemanella" } 
  } else if (genus == "Eubacterium") {
    if (species =="Eubacterium.s__Eubacterium_eligens"){
      species = "Lachnospira.s__Lachnospira_eligens"
      genus = "Lachnospira"
    } else if (species =="Eubacterium.s__Eubacterium_hallii"){
      species = "Anaerobutyricum.s__Anaerobutyricum_hallii"
      genus = "Anaerobutyricum" 
      #  } else if (species =="Eubacterium.s__Eubacterium_limosum"){ # no change in NCBI Taxonomy
      #  } else if (species =="Eubacterium.s__Eubacterium_ramulus"){  # no change in NCBI Taxonomy
    } else if (species =="Eubacterium.s__Eubacterium_rectale"){
      species = "Agathobacter.s__Agathobacter_rectalis"
      genus = "Agathobacter" 
    } else if (species =="Eubacterium.s__Eubacterium_siraeum"){ # new genus name pending. Make sure these special characters aren't an issue ... works at this stage of table, anyway 
      species = "[Eubacterium].s__[Eubacterium]_siraeum"
      genus = "[Eubacterium]"
      #    } else if (species =="Eubacterium.s__Eubacterium_ventriosum"){ # no change in NCBI Taxonomy
    } 
  } # end if genus == Eubacterium
  
  ## ACTUALLY UPDATE SPECIES/GENUS IN DF
  dmelt[row,2] = species
  dmelt[row,6] = genus
} # end whole loop


################## ############# ########################## ########################## 
#  UPDATE PREVOTELLA
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Prevotella"){
    species = dmelt[row,2] 
    if (species =="Prevotella.s__Prevotella_buccae"){
      species = "Segatella.s__Segatella_buccae" # according to NCBI taxonomy
      genus = "Segatella"
    } else if (species =="Prevotella.s__Prevotella_copri"){
      species = "Segatella.s__Segatella_copri"
      genus = "Segatella"
    }else if (species =="Prevotella.s__Prevotella_stercorea"){
      species = "Leyella.s_Leyella_stercorea"
      genus = "Leyella"
    } else if  (species =="Prevotella.s__Prevotella_timonensis"){
      species = "Hoylesella.s__Hoylesella_timonensis"
      genus = "Hoylesella"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Prevotella loops
}

################## ############# ########################## ########################## 
#  UPDATE LACTOBACILLUS
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lactobacillus"){
    species = dmelt[row,2] 
    if (species =="Lactobacillus.s__Lactobacillus_casei_paracasei"){
      species = "Lacticaseibacillus.s__Lacticaseibacillus_paracasei" # according to chatGPT
      genus = "Lacticaseibacillus"
      #  } else if (species =="Lactobacillus.s__Lactobacillus_delbrueckii"){ # no change 
    }else if (species =="Lactobacillus.s__Lactobacillus_plantarum"){
      species = "Lactiplantibacillus.s_Lactiplantibacillus_plantarum"
      genus = "Lactiplantibacillus"
    }else if (species =="Lactobacillus.s__Lactobacillus_ruminis"){
      species = "Ligilactobacillus.s_Ligilactobacillus_ruminis"
      genus = "Ligilactobacillus"
    } else if  (species =="Lactobacillus.s__Lactobacillus_salivarius"){
      species = "Ligilactobacillus.s__Ligilactobacillus_salivarius"
      genus = "Ligilactobacillus"
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } # end if genus==Lactobacillus loops
}


################## ############# ########################## ########################## 
#  UPDATE "Lachnospiraceae_bacterium_5_1_57FAA"
################## ############# ########################## ########################## 

for (row in 1:nrow(dmelt)){
  genus = dmelt[row,6]
  if (genus == "Lachnospiraceae_noname"){
    species = dmelt[row,2] 
    if (species =="Lachnospiraceae_noname.s__Lachnospiraceae_bacterium_5_1_57FAA"){
    #  species = "Lachnospiraceae_bacterium_5_1_57FAA" 
      species = "5_1_57FAA" 
    #  genus = "Lachnospiraceae_5_1_57FAA"
      genus = "5_1_57FAA" 
    } # end if species loops
    ## ACTUALLY UPDATE SPECIES/GENUS IN DF
    dmelt[row,2] = species
    dmelt[row,6] = genus
  } 
}


################## ############# ########################## ########################## 

##  create a new table, where each sample has a summed entry for each genus
dmelt_glom = data.frame(matrix(ncol=5,nrow=1))
names(dmelt_glom) = names(dmelt)[1:5]
genera = unique(dmelt$Genus)

## then for each sample, for each genus:
samples = unique(dmelt$Sample)
for (samp in samples){
  samp_dmelt = subset(dmelt,dmelt$Sample==samp)
  tribe = unique(samp_dmelt$Tribe)
  region = unique(samp_dmelt$Region)
  
  for (genus in genera){
    g_samp_dmelt = subset(samp_dmelt, samp_dmelt$Genus ==genus)
    genus_total = sum(g_samp_dmelt$Pathway_relab)
    
    # save a new row
    newrow = c(samp,genus,genus_total,tribe,region)
    dmelt_glom=rbind(dmelt_glom,newrow)
  }
}

dmelt_glom = na.omit(dmelt_glom)

dmelt_glom  <- dmelt_glom  %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))


##########################################################
#### SELECT TOP TAXA (for legend) FOR THIS ONE PATHWAY        
##########################################################

dmelt_saved = dmelt
dmelt = dmelt_glom
dmelt$Pathway_relab = as.numeric(dmelt$Pathway_relab)

# calculate total abundance per taxon and tribe across  fecal samples
totals <- dmelt %>% dplyr::group_by(Taxon, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Pathway_relab))

# select most abundant taxa to include in legend
common_taxa<- subset(totals,totals$mean_Rel_ab >=  3.0e-7)  
length(unique(common_taxa$Taxon)) # AIM FOR 8 (including unclassified) 
top_taxa = unique(common_taxa$Taxon)
top_taxa

# get abundances for only these genera
dmelt_top = subset(dmelt,dmelt$Taxon%in%top_taxa)

######################
## ADD 'OTHER' COLUMN (genera that are not abundant enough to get their own label)

# First, sum up the abundances per person of genera that DO get a label
labeled_per_person_top <- dmelt_top %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Because the relabs of this pathway don't sum to 1 -- it's relab of ALL pathways --
# 'other' is not just 1-x. It's the total -x. So, calculate

#  the total rel_ab for this pathway per person
labeled_per_person_all <- dmelt %>%
  group_by(Sample)  %>%
  summarise(percent_with_taxon_label = sum(Pathway_relab))

# Then define other as 1 minus that sum
labeled_per_person_top$other = labeled_per_person_all$percent_with_taxon_label - labeled_per_person_top$percent_with_taxon_label
#labeled_per_person_top$checksum = labeled_per_person_top$other + labeled_per_person_top$percent_with_taxon_label

# Add 'other' to dmelt_top in the form of one row per person, where 'genus' = 'other'
dplot <- dmelt_top
dplot$Taxon = as.character(dplot$Taxon)
for (i in (unique(dmelt_top$Sample))){
  tribe = unique(subset(dmelt_top,dmelt_top$Sample==i)$Tribe)
  region = unique(subset(dmelt_top,dmelt_top$Sample==i)$Region)
  other_number = as.numeric(subset(labeled_per_person_top,labeled_per_person_top$Sample==i)[1,3])
  vector = c(i,'other',other_number,tribe,region)
  dplot = rbind(dplot,vector)
}
unique(dplot$Taxon)
##### at this point, unclassified and other are both there


# order genera by mean abundance
# Common taxa has mean_rel_ab of each taxon by tribe, but DOES NOT INCLUDE 'OTHER'
common_taxa$Taxon = as.character(common_taxa$Taxon) ## THIS does not include 'other'
t <- common_taxa %>%
  group_by(Taxon)  %>%
  summarise(mean_Rel_ab = mean(mean_Rel_ab))
t <- t[order(t$mean_Rel_ab,decreasing=TRUE),]
top_taxa_ordered = t$Taxon
top_taxa_ordered 

# then edit the order so that red (TH-enriched) taxa are on top
lactose_galactose = c("Ligilactobacillus", "Streptococcus","Lactobacillus","Enterococcus" , "other", "unclassified")
lactose_galactose %in% dplot$Taxon
length(lactose_galactose) == length(unique(dplot$Taxon))
dplot$Taxon = factor(dplot$Taxon,levels=lactose_galactose)


## adjust class
dplot$Pathway_relab = as.numeric(dplot$Pathway_relab) 
dplot$Tribe <- factor(dplot$Tribe, levels = c("Boto","Balti","Brokpa","Purigpa","Warli", "Gond","Madia","Kabui"))


## subset and order colors for this particular function / set of genera
specific_colortable <- subset(colortable,colortable$genus %in% lactose_galactose)
specific_colortable = specific_colortable[match(lactose_galactose, specific_colortable$genus),]
col_list = specific_colortable$color

## PLOT ! FACET WRAP ALL TRIBES TOGETHER
### TO AUTOMATICALLY CALCULATE THE DISTANCE BETWEEN DATA AND TOP OF PLOT:
temp <- dplot %>% group_by(Sample) %>% summarise(sum_relab = sum(Pathway_relab))
margin = max(temp$sum_relab)*1.05

s <- ggplot(data=dplot, aes(x=reorder(Sample, -Pathway_relab, function(x){ sum(x) }), y=Pathway_relab, fill=Taxon)) +
  geom_col() + theme_classic() +
  ylab("Relative abundance") + theme_cowplot() +
  scale_fill_manual("legend",values = col_list) + 
  theme(legend.title=element_blank()) + theme(legend.key.size = unit(0.35, "cm")) +
  ggtitle("Lactose and galactose degradation (LACTOSECAT-PWY)")+
  theme(plot.title = element_text(size = 12,face="bold")) +
  theme(axis.title=element_text(size=12,face="plain")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, margin),labels =  ~ sprintf(fmt = "%2.e", .)) + # removes space between 0 and x-axis
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) 

L <- s + facet_wrap(~Tribe,scales="free_x",nrow=1) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
L

ggsave("Fig5_S1B.pdf", plot=L,width=4.64*2.2,height=0.95*2.2,units="in")

