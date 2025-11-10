library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(cowplot)
library("pheatmap")
library(viridis)
library(tidyr)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)


#################################################################################
#   Figure 4 Supplement 2 
#################################################################################


##################################################################################
#####   A: boxplots of Segatella, Bifidobacterium, and Ligilactobacillus by Region

#########
# 16S

### Read in abundance data 
d16 <-read.table("Fig2_relab_genus16S.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
d16 <- subset(d16,select=-c(AK_SK_31)) # has no metadata

## Remove outliers 
d16$AK_SK_49 = NULL
d16$AK_SK_10 = NULL
d16$AK_SK_32 = NULL
d16$AK_SR_6 = NULL
d16$AK_SG_17 = NULL
d16$AK_SG_18 = NULL
d16$AK_SR_1 = NULL
d16$AK_SR_2 = NULL
d16$AK_SR_4 = NULL # high Streptococcus

## Melt data
d16melt <- reshape::melt(d16, id = c("Genus"))
names(d16melt) <- c("Genus","Sample","Rel_ab") 
d16melt$Sample <- as.character(d16melt$Sample)

## Add tribe/region
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta<-subset(meta,meta$sample %in% d16melt$Sample)
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))
#
d16melt$Tribe = NA
d16melt$Region = NA 
for (row in 1:nrow(d16melt)){
  samplename = d16melt[row,2]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  d16melt[row,4] = tribe
  d16melt[row,5] = region
}

## Subset target genera &  Prepare df for plotting
target = c("Segatella","Bifidobacterium","Ligilactobacillus")
S42_16 = subset(d16melt,d16melt$Genus %in% target)
S42_16$Genus = factor(S42_16$Genus,levels=c("Segatella","Bifidobacterium","Ligilactobacillus"))
S42_16$Tribe = factor(S42_16$Tribe,levels=c("Balti","Brokpa","Purigpa","Boto","Warli","Madia","Gond","Kabui"))
S42_16$Region = factor(S42_16$Region,levels=c("Trans-Himalayas","Coast","Deccan Plateau","Northeast Hills"))

########################################################
###   Plot Figure 4 Supplement 2A (top, 16S)
#ragg::agg_png("Fig4_S2BA_16S.png", width = 4.5*1, height = 2.5*1, units = "in", res = 300)
ggplot(S42_16, aes(fill=Genus, y=Rel_ab, x=Region)) +  theme_classic() +
  geom_boxplot() +ylab("Relative abundance") + ggtitle("16S") +
  theme(legend.position="top")
#dev.off()    




##############
# metagenomics

d<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss") 

## Remove outliers 
outliers =c('AK_SK_49','AK_SK_10',"AK_SK_32","AK_SR_6",'AK_SG_17','AK_SG_18','AK_SR_1',"AK_SR_2","AK_SR_4")
d <- subset(d,(d$sample %in% outliers) == F)

## Instead of melting, make empty table with headers: Sample Genus Rel_ab
samples <- unique(d$sample)
genera = unique(d$genus)
genusRA <- data.frame(matrix(ncol = 3, nrow = length(samples)*length(genera)))
names(genusRA) <- c("Sample","Genus","Rel_ab")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(d,d$sample==samples[[s]]) # for each sample 
  sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab)) # scale RA to sum to 1 (aka ignore unclassified reads)
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

genusRA <- subset(genusRA,genusRA$Genus!='')


## Add tribe/region
meta2 <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T)  # 1 fewer sample than 16S b/c AK_SG_19 failed metagenomic sequence
meta2<-subset(meta2,meta2$sample %in% genusRA$Sample)
meta2 <- meta2 %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta2 <- meta2 %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta2 <- meta2 %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta2 <- meta2 %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta2 <- meta2 %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta2 <- meta2 %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))
#
genusRA$Tribe = NA
genusRA$Region = NA 
for (row in 1:nrow(genusRA)){
  samplename = genusRA[row,1]
  tribe <- as.character(subset(meta2,meta2$sample==samplename)$Tribe)
  region <- as.character(subset(meta2,meta2$sample==samplename)$Region)
  genusRA[row,4] = tribe
  genusRA[row,5] = region
}

## Subset target genera &  Prepare df for plotting
S42_meta = subset(genusRA, genusRA$Genus %in% target)
S42_meta$Genus = factor(S42_meta$Genus,levels=c("Segatella","Bifidobacterium","Ligilactobacillus"))
S42_meta$Tribe = factor(S42_meta$Tribe,levels=c("Balti","Brokpa","Purigpa","Boto","Warli","Madia","Gond","Kabui"))
S42_meta$Region = factor(S42_meta$Region,levels=c("Trans-Himalayas","Coast","Deccan Plateau","Northeast Hills"))

####### Plot Figure 4 Supplement 2A (bottom, metagenomics)
#ragg::agg_png("Fig4_S2A_meta.png", width = 4.5*1, height = 2.5*1, units = "in", res = 300)
ggplot(S42_meta, aes(fill=Genus, y=Rel_ab, x=Region)) +  theme_classic() +
  geom_boxplot() +ylab("Relative abundance") + ggtitle("Metagenomics")  +
  theme(legend.position="top")
#dev.off()    



##################################################################################
#####   B: Scatterplots of Segatella X Bifidobacterium per tribe

# Summarize each data type by genus and tribe 
totals16 <- d16melt %>% dplyr::group_by(Genus, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Rel_ab))
totalsmeta <- genusRA %>% dplyr::group_by(Genus, Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Rel_ab))

### Subset Segatella and Bifido to calculate correlation
temp16Bifido <- subset(totals16,totals16$Genus %in% c("Bifidobacterium") )
temp16Sega <- subset(totals16,totals16$Genus %in% c("Segatella") )
summary(lm(temp16Bifido$mean_Rel_ab~temp16Sega$mean_Rel_ab)) # R2 = 0.674, p = 0.007656
#
tempMetaBifido <- subset(totalsmeta,totalsmeta$Genus %in% c("Bifidobacterium") )
tempMetaSega <- subset(totalsmeta,totalsmeta$Genus %in% c("Segatella") )
summary(lm(tempMetaBifido$mean_Rel_ab~tempMetaSega$mean_Rel_ab)) # R2 = 0.1045, p = 0.226


### Organize to plot
# 16S
plotdata16= as.data.frame(cbind(temp16Bifido$mean_Rel_ab,temp16Sega$mean_Rel_ab,temp16Sega$Tribe))
names(plotdata16) =c("Bifidobacterium","Segatella","Tribe")
plotdata16$Tribe=factor(plotdata16$Tribe,levels = c("Warli","Kabui","Gond","Madia","Balti","Boto","Brokpa","Purigpa"))
plotdata16$data.type = "16S"
# metagenomics
plotdataM= as.data.frame(cbind(tempMetaBifido$mean_Rel_ab,tempMetaSega$mean_Rel_ab,tempMetaSega$Tribe))
names(plotdataM) =c("Bifidobacterium","Segatella","Tribe")
plotdataM$Tribe=factor(plotdataM$Tribe,levels = c("Warli","Kabui","Gond","Madia","Balti","Boto","Brokpa","Purigpa"))
plotdataM$data.type = "metagenomics"
# combine
plotdata = rbind(plotdata16,plotdataM)
plotdata$Bifidobacterium=as.numeric(plotdata$Bifidobacterium)
plotdata$Segatella=as.numeric(plotdata$Segatella)



########################################################
###   Plot Figure 4 Supplement 2B 
x=ggplot(plotdata, aes(x=Bifidobacterium, y=Segatella)) +
  geom_point(aes(color=Tribe),size=2.5) + theme_cowplot() +
  theme(legend.position = "top") +
  scale_colour_manual(values = c("#b33f25","#edc42f","#2f3cb5","#acb2e6","#7CAE00","#a5e602","#cde09d","#4b6901") )
M <- x + facet_wrap(~data.type,scales="free_y",nrow=2) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2


ggsave("Fig4_S2B.png",plot=M,width=3.9*1.2,height=4.7*1.2,units="in")



########################################################
### Code for Panel C is in "Fig4_StackedBar_Tribe.R"





#################################################################################
#   Correlations between taxa types in abundances of main taxa
#################################################################################

# how many genera are shared between data types? 
genera_16S = unique(d16melt$Genus)
genera_M = unique(genusRA$Genus)
#
which(genera_16S %in% genera_M)
shared = genera_16S[which(genera_16S %in% genera_M)]
length(shared); shared

# let's look at the abundant genera in Figure 4A
fig4a = c("Segatella","Leyella","Bifidobacterium","Catenibacterium","UCG-002",
             "Megasphaera","Ligilactobacillus","Dialister","Asteroleplasma","Faecalibacterium",
             "Agathobacter","Succinivibrio","Treponema","Escherichia-Shigella","Ruminobacter")

# some of these genera AREN'T in the metagenomic data 
fig4a[-which(fig4a %in% genera_M)]  # UCG-002 Asteroleplasma Treponema Eschericia-Shigella Ruminobacter 

# metagenomics does have "Treponema_D" and "Escherichia" -- rename to match 16S 
# for the other taxa, we'll add to genusRA at 0 abundance after summarizing shortly 
genusRA2 <- genusRA %>% mutate(Genus = str_replace(Genus, "Treponema_D", "Treponema"))
genusRA2 <- genusRA2 %>% mutate(Genus = str_replace(Genus, "Escherichia", "Escherichia-Shigella"))

plotM = subset(genusRA2,genusRA2$Genus %in% fig4a)
plot16 = subset(d16melt,d16melt$Genus %in% fig4a)
length(unique(plotM$Genus))
length(unique(plot16$Genus))

### get the mean rel_ab of each of these taxa by tribe and data type
mean16 <- plot16 %>% dplyr::group_by(Genus,Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Rel_ab))
meanmeta <- plotM %>% dplyr::group_by(Genus,Tribe) %>% dplyr::summarise(mean_Rel_ab = mean(Rel_ab))

# drop the 3 taxa from 16S that don't have equivalents in metagenomics 
mean16<-subset(mean16,mean16$Genus %in% meanmeta$Genus)

# now we can plot abundance correlations!
# if we combine both datatypes into one df
mean16$data.type = "16S"
meanmeta$data.type = "meta"
plotdata <- rbind(mean16,meanmeta)


# Pivot wider so that each Genus has both 16S and metag values
df_wide <- plotdata %>%
  pivot_wider(
    names_from = data.type,
    values_from = mean_Rel_ab )

# Define the colors 
genus_vec = c("Segatella","Leyella","Bifidobacterium","Catenibacterium",
                        "Megasphaera","Ligilactobacillus","Dialister","Faecalibacterium",
                        "Agathobacter","Succinivibrio","Treponema","Escherichia-Shigella")
              

color_vec = c("#E6AB02",'#fc8803',"#b30707",'#f5051d',
              '#d43fb6',  '#f57fc0','#f09292',"#1F78B4",
  '#29992b','#84b9fa',"#B2DF8A",'#1b1d87')
colortable= as.data.frame(cbind(color_vec,genus_vec))
names(colortable)=c("color","Genus")  

# combine with df_wide
df_wide2 <- df_wide %>%
  left_join(colortable, by = "Genus")

# Then make scatterplot faceted by Tribe
z = ggplot(df_wide, aes(x = `16S`, y = meta)) +
  geom_point(aes(color = Genus), size = 2, alpha = 0.8) +
  facet_wrap(~ Tribe, scales = "free") +
  scale_color_manual(values = setNames(df_wide2$color, df_wide2$Genus)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_bw() +
  labs( x = "16S Mean Relative Abundance",
    y = "Metagenomic Mean Relative Abundance" )
z

# Looks good -- save!
ggsave("Fig4_SX.pdf", plot=z,width=8,height=6,units="in")


# calculate linear models for each tribe
tribes = unique(meanmeta$Tribe)
genera = unique(df_wide2$Genus)

for (t in 1:length(tribes)) {
    temp = subset(df_wide,df_wide$Tribe == tribes[t])
    print(tribes[t])
    print(summary(lm(temp$meta~temp$`16S`)))
  }

# Tribe R2 p-val
# Balti  0.9202 5.09e-07 
# Boto  0.9347 1.85e-07 
# Brokpa  0.9467 6.69e-08 
# Gond 0.9046  1.257e-06
# Kabui 0.9303  2.577e-07
# Madia 0.8501  1.227e-05
# Purigpa  0.9519 3.989e-08
# Warli 0.9139  7.463e-07

