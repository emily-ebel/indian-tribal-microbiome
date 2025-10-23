library(reshape2)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library("ggsignif")   

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)

#####################################################################
### Figure 2 Panels D, E,F
### Need genus-level abundance for both data types + outlier status


###################################
###  16S  genus-level abundance 

d <-read.table("Fig2_relab_genus16S.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
row.names(d) <- d$Genus
d <- d[ , !(names(d) %in% c("AK_SK_31"))] # included by mistake
d16melt <- reshape2::melt(d, id = c("Genus"))
names(d16melt) <- c("Genus","Sample","Rel_ab") 
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2") # 8 high-Bacteroides 
d16melt$Sample <- as.character(d16melt$Sample)
d16melt$datatype = "16S"

# add outlier status 
outliers1<-subset(d16melt,d16melt$Sample %in% outliers) 
outliers1$outlier = 'yes'
others1 <- subset(d16melt,(d16melt$Sample %in% outliers) == FALSE )
others1$outlier = 'no'
d16melt2 <- rbind(outliers1, others1)
names(d16melt2)
d16melt2 <- d16melt2[,c(2,1,3,4,5)] # new df


##########################################
###  Metagenomic taxon abundance 

d<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss") 

### Remove enrichment samples 
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta<-subset(meta,meta$sample %in% d$sample)

# Make empty table with headers: Sample Genus Rel_ab
samples <- unique(d$sample); length(samples)
genera = unique(d$genus); length(genera)
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
genusRA$datatype = "metagen."

# add Bacteroides outlier status (N=8 from PCA, see methods) 
outliers =c('AK_SK_49','AK_SK_10',"AK_SK_32","AK_SR_6",'AK_SG_17','AK_SG_18','AK_SR_1',"AK_SR_2")
outliers1<-subset(genusRA,genusRA$Sample %in% outliers) 
outliers1$outlier = 'yes'
others1 <- subset(genusRA,(genusRA$Sample %in% outliers) == FALSE )
others1$outlier = 'no'
genusRA2 <- rbind(outliers1, others1)
names(genusRA2)


################################
# Combine across both data types

genusboth = rbind(d16melt2,genusRA2)



#####################################################################
### FIGURE 2D
#####################################################################

### Pull out Bacteroides and relatives

bact <- subset(genusboth,genusboth$Genus %in% c("Bacteroides","Parabacteroides","Phocaeicola"))
# combine 
bact2 <- bact %>% dplyr::group_by(Sample,datatype,outlier)  %>% dplyr::summarise(bactsum = sum(Rel_ab))
bact2$outlier = factor(bact2$outlier,levels=c("yes","no"))

bact_by_outlier <- ggplot(bact2, aes(x=outlier, y=bactsum)) + 
  geom_boxplot(aes(fill=datatype)) +
  ylab("Bacteroides* relative abundance") + xlab("Outlier") +
  scale_fill_manual(values = c('#ffffff',"#a3a3a2")) + 
  theme_cowplot() + ylim(0,0.7)+  theme(legend.position = "top") +  theme(legend.title=element_blank()) +
  theme(plot.margin = margin(10,10,10,10, "pt")) 
bact_by_outlier

#ggsave("Fig2D.pdf", plot=bact_by_outlier,width=3,height=4,units="in")


# stats
outlier16s<-subset(bact2,bact2$outlier=='yes' & bact2$datatype=="16S") 
outliermeta<-subset(bact2,bact2$outlier=='yes' & bact2$datatype!="16S") 
others16s <- subset(bact2,bact2$outlier == 'no' &  bact2$datatype=="16S")
othersmeta <- subset(bact2,bact2$outlier == 'no' &  bact2$datatype!="16S")

t.test(outlier16s$bactsum,others16s$bactsum) # p = 0.00587
t.test(outliermeta$bactsum,othersmeta$bactsum) # p = 0.000759




#####################################################################
### FIGURE 2E
#####################################################################

### Need alpha-diversity for both data types (RSV/species) + outlier status

################
#### 16S
d <-read.table("Fig2_relab_RSV.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
row.names(d) <- d$RSV
d <- subset(d,select=-c(RSV))
d <- d[ , !(names(d) %in% c("AK_SK_31"))] # drop sample with no metadata

alpha <- data.frame(matrix(ncol = 2, nrow = 0))
for (col in 1:ncol(d)){
  data = d[,col]
  N = sum(data>0)
  subject = names(d)[col]
  newrow = c(subject,N)
  alpha = rbind(alpha,newrow)  }
colnames(alpha) <- c('subject', 'RSVs') # 
alpha$RSVs = as.numeric(alpha$RSVs)

# add outlier status  USING SAME 8 HIGH-BACTEROIDES OUTLIERS FOR BOTH METHODS  
outliers =c('AK_SK_49','AK_SK_10',"AK_SK_32","AK_SR_6",'AK_SG_17','AK_SG_18','AK_SR_1',"AK_SR_2")
outliers1<-subset(alpha,alpha$subject %in% outliers) 
outliers1$outlier = 'yes'
others1 <- subset(alpha,(alpha$subject %in% outliers) == FALSE )
others1$outlier = 'no'
t.test(outliers1$RSVs,others1$RSVs) # p = 0.00238
alpha_16S <- rbind(outliers1, others1)

####################
### metagenomics  

d<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss") 

## Make  table with headers: Sample SPECIES Rel_ab
samples <- unique(d$sample); length(samples)
species_list = unique(d$genome)
length(species_list)
meta_species <- data.frame(matrix(ncol = 3, nrow = length(samples)*length(species_list)))
names(meta_species) <- c("Sample","Species","Rel_ab")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(d,d$sample==samples[[s]]) # for each sample 
  sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab))  # scale rel_ab so each sample sums to 1 (i.e. ignore unmapped/unclassified reads)
  for (t in 1:length(species_list)){ # for each genus 
    count = count+1
    meta_species[count,1] = as.character(samples[[s]])
    meta_species[count,2]= as.character(species_list[[t]])
    if ( (species_list[[t]] %in% sdf$genome)==TRUE){
      rows <- subset(sdf,sdf$genome == species_list[[t]]) 
      meta_species[count,3] = sum(rows$rel_ab_adj)  # save the rel_ab for all genomes in the genus 
    } else (meta_species[count,3]=0)
  }
}

# Remove entries where rel_ab is 0
meta_species2 <- subset(meta_species,meta_species$Rel_ab != 0)

# Count species per sample 
Nspecies <- meta_species2 %>% dplyr::group_by(Sample)  %>% dplyr::summarise(n = n())


# add outlier status 
outliers =c('AK_SK_49','AK_SK_10',"AK_SK_32","AK_SR_6",'AK_SG_17','AK_SG_18','AK_SR_1',"AK_SR_2")
outliers1<-subset(Nspecies,Nspecies$Sample %in% outliers) 
outliers1$outlier = 'yes'
others1 <- subset(Nspecies,(Nspecies$Sample %in% outliers) == FALSE )
others1$outlier = 'no'
t.test(outliers1$n,others1$n) # p = 0.975
alpha_meta <- rbind(outliers1, others1)


### Combine both data types 
alpha_16S$datatype = "16S"
alpha_meta$datatype = "metagenomics"
names(alpha_16S) = c("sample","taxa","outlier","datatype")
names(alpha_meta) = c("sample","taxa","outlier","datatype")
alpha_both = rbind(alpha_16S,alpha_meta)
alpha_both$outlier = factor(alpha_both$outlier,levels=c("yes","no"))

otu_by_outlier <- ggplot(alpha_both, aes(x=datatype, y=taxa)) + 
  geom_boxplot(aes(fill=outlier)) +
  ylab("Taxa per subject") + xlab("") + ylim(0,300) +
  scale_fill_manual(values = c("#a3a3a2",'#ffffff')) + 
  theme_cowplot() + theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  theme(plot.margin = margin(10,10,-10,10, "pt")) 
otu_by_outlier

#ggsave("Fig2E.pdf", plot=otu_by_outlier,width=3,height=4,units="in")


#####################################################################
### FIGURE 2F
#####################################################################

## see how the one old genus (Prevotella) shakes out into the new genera
prev = subset(genusboth,genusboth$Genus %in% c("Prevotella","Segatella","Leyella","Hallella","Xylanibacter"))
prev2 <- prev %>% group_by(Genus,datatype,outlier)  %>% summarise(prevsum = sum(Rel_ab))
# --> it's mostly Segatella (mean 28% 16S, 11% metagenomics). Next highest is Leyella (mean ~1.5% both data types)

## Make the plot with Segatella
prev = subset(genusboth,genusboth$Genus %in% c("Segatella"))
prev2 <- prev %>% dplyr::group_by(Sample,datatype,outlier)  %>% dplyr::summarise(prevsum = sum(Rel_ab))
prev2$outlier = factor(prev2$outlier,levels=c("yes","no"))

## add prev data to bact df
bact2$prevsum = prev2$prevsum 

## plot
scatter2 <- ggplot(bact2, aes(x=bactsum, y=prevsum, fill=datatype, col=datatype, shape=outlier)) + 
  geom_point(size=4) +
  ylab("Segatella relative abundance") + xlab("Bacteroides* relative abundance") +
  theme_cowplot() + 
  scale_color_manual(values = c('black',"black")) + 
  scale_fill_manual(values = alpha(c('#ffffff',"#a3a3a2"),0.7)) + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  theme(plot.margin = margin(10,10,10,10, "pt"))  + 
  scale_shape_manual(values=c(23,21))  
scatter2

#ggsave("Fig2F.pdf", plot=scatter2,width=4.5,height=4,units="in")

#### LINEAR MODEL
temp <- subset(bact2,bact2$datatype=="16S")
summary(lm(temp$prevsum~temp$bactsum)) # p = 4.245e-06, R = 0 .24

temp <- subset(bact2,bact2$datatype!="16S") 
summary(lm(temp$prevsum~temp$bactsum)) # p = 0.006, R = 0 .0859


