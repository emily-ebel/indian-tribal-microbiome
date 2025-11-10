library(reshape2)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
library("ggsignif") 
library(vegan)
library(randomcoloR)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)

# Function for plotting overlapping histograms
makeTransparent<-function(someColor, alpha=100)
{  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                            blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}




##################################################################################
#######    16S Alpha Diversity (RSVs)     ########################################
##################################################################################

### Read in abundance data 
d <-read.table("Fig2_relab_RSV.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
row.names(d) <- d$RSV
d <- subset(d,select=-c(RSV,AK_SK_31)) # has no metadata

### Create df "alpha" to contain number of taxa per sample for the 16S data
alpha <- data.frame(matrix(ncol = 2, nrow = 0))

for (col in 1:ncol(d)){
  data = d[,col]
  N = sum(data>0)
  subject = names(d)[col]
  newrow = c(subject,N)
  alpha = rbind(alpha,newrow)  }

colnames(alpha) <- c('subject', 'RSVs') #  RDP, non-agglomerated ONLY
alpha$RSVs = as.numeric(alpha$RSVs)



### Get sequencing depth from a file with samples (rows) by RSVs (columns) with depth in cells  
readcount_16S <- read.table("Fig2_16S_depth.csv", header = F, sep = " ",stringsAsFactors = F,numerals="allow.loss") 
goodnames <- c("AK_SG_10","AK_SG_11","AK_SG_12","AK_SG_15","AK_SG_16","AK_SG_17","AK_SG_18","AK_SG_19","AK_SG_2","AK_SG_25",
               "AK_SG_26","AK_SG_4","AK_SG_5","AK_SG_7","AK_SG_9",
               
               "AK_SK_9","AK_SK_14","AK_SK_18","AK_SK_2","AK_SK_26","AK_SK_27","AK_SK_29","AK_SK_31","AK_SK_36","AK_SK_39",
               "AK_SK_47","AK_SK_47.2","AK_SK_49","AK_SW_1",
               
               "AK_SW_11","AK_SW_12","AK_SW_14","AK_SW_2","AK_SW_3","AK_SW_4","AK_SW_5","AK_SW_7","AK_SW_8","AK_SW_9","AK_SG_20",
               "AK_SG_22","AK_SG_28","AK_SG_29","AK_SG_8","AK_SK_1","AK_SK_10","AK_SK_11","AK_SK_12",
               
               "AK_SK_13","AK_SK_15","AK_SK_17","AK_SK_19","AK_SK_25","AK_SK_30","AK_SK_32","AK_SK_33","AK_SK_35","AK_SK_37",
               "AK_SK_38","AK_SK_40","AK_SK_44","AK_SK_45","AK_SK_48","AK_SK_49.2","AK_SK_5","AK_SK_8","AK_SR_1","AK_SR_10",
               "AK_SR_2.2","AK_SR_3","AK_SR_4","AK_SR_5","AK_SR_6","AK_SR_7","AK_SR_8","AK_SR_9")
row.names(readcount_16S) = goodnames
readcount_16S$V1 <- NULL
depth16 = as.data.frame(rowSums(readcount_16S)); names(depth16)='total_reads'
depth16 = subset(depth16,row.names(depth16) != "AK_SK_31")
mean(depth16$total_reads)
depth16$sample = row.names(depth16)

## add the depth to "alpha"
alpha <- alpha[order(alpha$subject),]
depth16S <- depth16[order(depth16$sample),]
alpha$subject == depth16S$sample
alpha$depth = depth16S$total_reads


################################################
### Plot RSVs by depth (Figure 2 supplement 1a)

ggplot(alpha, aes(x=depth, y=RSVs)) + 
  geom_point(colour="black",pch=21, size=3) + theme_cowplot() + ylab("RSVs") + xlab("High-quality reads")+
  theme(legend.position = "none") +
  geom_smooth(method = "lm", se = FALSE,color="black",linewidth=.5) +xlim(0,102000)+ggtitle('16S')

# statistics - linear model
summary(lm(alpha$RSVs~alpha$depth)) # R2=0.3773, p =2.767e-09




#### Add tribe and region from metadata to "alpha"

meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T)
meta<-subset(meta,meta$sample %in% alpha$subject)

alpha$tribe = NA
alpha$region = NA 
for (row in 1:nrow(alpha)){
  samplename = alpha[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  alpha[row,4] = tribe
  
  if (is.na(tribe)==F) {
    if (tribe == "Warli"){
      alpha[row,5] = "Coast"
    } else if (tribe == "Kabui"){
      alpha[row,5] = "Northeast Hills"
    } else if (tribe == "Gondia" | tribe == "Madia"){
      alpha[row,5] = "Deccan Plateau"
    } else {alpha[row,5]="Trans-Himalayas"} 
  } else {alpha[row,5]="Trans-Himalayas"}  

} 

################################################
### 16S rarefaction curve 

# Need matrix where cells are reads, rows are samples, and columns are RSVs --> readcount_16 has the reads (for 3007 RSVs)
rarecurve(readcount_16S,label=F,step=200,xlab="16S reads",ylab="RSVs detected",col=randomColor(count = 76))
# --> most of these curves do plateau 




############################################################
## Plot 16S alpha diversity by Region -- Figure 2A
############################################################

alpha$region <- factor(alpha$region, levels = c("Coast","Deccan Plateau","Northeast Hills","Trans-Himalayas"))
col_list<-c("#b33f25", "#2f3cb5",   "#edc42f","#2a8022")

c <- ggplot(alpha, aes(x=region, y=RSVs, group=region)) + geom_boxplot(aes(fill=region)) +
  theme_cowplot() + ylab("RSV count") + theme(legend.position = "none") + xlab("") +
  scale_fill_manual(values = col_list) +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) +
  theme(plot.margin = margin(10,10,-10,10, "pt"))

# add bars for significant pairwise comparisons
c +  geom_signif(comparisons = list(c("Trans-Himalayas","Northeast Hills"),  c("Trans-Himalayas","Deccan Plateau"),c("Trans-Himalayas","Coast")),
                        map_signif_level = T,  y_position = c(280, 304,328), textsize=0, # add asterisks in ppt
                        test='t.test') 

# statistics - one-way ANOVA 
res.aov <- aov(RSVs ~ region, data = alpha)  
summary(res.aov) # p = 0.001
p2 <- as.data.frame(TukeyHSD(res.aov)[1]) # TH vs all are significant; nothing else 





###############################################################
## Plot 16S alpha diversity by Tribe -- Figure 1 Supplement 1D
###############################################################

alpha2 <-alpha %>% mutate(tribe = str_replace(tribe, "Gondia", "Gond"))
alpha2$tribe <- factor(alpha2$tribe, levels = rev(c("Purigpa","Balti","Brokpa","Boto","Madia","Gond","Kabui","Warli")))
alpha2<-subset(alpha2,is.na(alpha2$tribe)==F)
ggplot(alpha2, aes(x=tribe, y=RSVs, group=tribe)) + geom_boxplot(aes(fill=tribe)) +
  theme_cowplot() + ylab("RSV count") + theme(legend.position = "none") + xlab("") +
  scale_fill_manual(values = c('#b33f25',"#edc42f","#2f3cb5",'#acb2e6',"#7CAE00","#a5e602","#cde09d","#4b6901")) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) + ggtitle('16S')

# statistics - one-way ANOVA 
res.aov <- aov(RSVs ~ tribe, data = alpha2) # one-way ANOVA was performed 
summary(res.aov) # p = 0.008
p2 <- as.data.frame(TukeyHSD(res.aov)[1]) # only purigpa-Warli significant 





##################################################################################
#######           Metagenomic Alpha Diversity 
#######          (Species Representative Genomes)     
##################################################################################


### Read in abundance data (at genome level)
d<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss") 

### Read in metadata
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T)
meta<-subset(meta,meta$sample %in% d$sample)


#### Create Species Abundance Table with headers:  Sample Species Rel_ab
samples <- unique(d$sample); length(samples) # start with empty tables
species_list = unique(d$genome)
meta_species <- data.frame(matrix(ncol = 3, nrow = length(samples)*length(species_list)))
names(meta_species) <- c("Sample","Species","Rel_ab")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(d,d$sample==samples[[s]]) # for each sample 
  
  for (t in 1:length(species_list)){ # for each genus 
    count = count+1
    meta_species[count,1] = as.character(samples[[s]])
    meta_species[count,2]= as.character(species_list[[t]])
    
    if ( (species_list[[t]] %in% sdf$genome)==TRUE){
      rows <- subset(sdf,sdf$genome == species_list[[t]]) 
      meta_species[count,3] = sum(rows$rel_ab)  # save the rel_ab for all genomes in the genus 
    } else (meta_species[count,3]=0)
  }
}

# Remove entries where species is not specified
meta_species <- subset(meta_species,meta_species$Species!='')

# To summarize alpha diversity, remove entries where rel_ab is 0
meta_species2 <- subset(meta_species,meta_species$Rel_ab != 0)

# Then create table with [samples, N_species]
Nspecies <- meta_species2 %>% dplyr::group_by(Sample)  %>% dplyr::summarise(n = n())
hist(Nspecies$n,xlab="genomes")



############ Quantify species detection ~ depth
metadepth <- as.data.frame(cbind(d$sample,d$pairs))
metadepth <- unique(metadepth)
names(metadepth) <- c("sample","readpairs")
metadepth$readpairs = as.numeric(metadepth$readpairs)

# add the depth to alpha
Nspecies <- Nspecies[order(Nspecies$Sample),]
metadepth <- metadepth[order(metadepth$sample),]
Nspecies$Sample == metadepth$sample
Nspecies$depth = metadepth$readpairs


# CONVERT 'pairs' to 'Gb'
Nspecies$bp = Nspecies$depth*280 # 2x140 sequencing
Nspecies$Gbp = Nspecies$bp /1000000000

### plot Figure 2 Supplement 1B 
Nspecies$dummy = 1
x=ggplot(Nspecies, aes(x=Gbp, y=n)) + 
  geom_point(aes(colour=factor(dummy), 
                 fill = factor(dummy)), shape=21, size = 1) + 
  scale_fill_manual(values=c('orange')) + 
  scale_colour_manual(values=c('orange'))+
  theme_cowplot() + ylab("Species") + xlab("Gbp")+
  theme(legend.position = "none") +
  geom_smooth(method = "lm", se = FALSE,color='orange') +ggtitle('Metagenomics')
x

#ggsave("Fig2_S1B.pdf", plot=x,height=1.79*1.8,width=2.13*1.8,units="in")

summary(lm(Nspecies$n~Nspecies$Gbp)) #Adjusted R-squared:  0.305 , 1.675e-07
# this linear regression is skewed by the outliers AK_SR_1 and AK_SR_2 and, to a lesser extent, AK_SK_32


################################################
### metagenomic rarefaction curve 

# Need matrix where cells are reads, rows are samples, and columns are GENERA (fewer columns/more manageable than species)
samples <- unique(d$sample)
taxa <- unique(d$genome)
df <- data.frame(matrix(ncol = length(taxa), nrow = length(samples)))
names(df) <- taxa; rownames(df) <- samples
# 
for (s in 1:length(samples)){
  sdf <- subset(d,d$sample==samples[[s]])
   for (t in 1:length(taxa)){
    if ( (taxa[[t]] %in% sdf$genome)==TRUE){
      row <- subset(sdf,sdf$genome == taxa[[t]])
      df[s,t]= row$reads # save the READS from this row/genome/species 
    } else (df[s,t]=0)
  }
} 


## Examine range of total depth per sample
totals <- rowSums(df)
min_depth <- min(totals) # 182,810 
max_depth <- max(totals) # 37,901,839

# choose 100 points log-spaced to sample early growth finely
n_points <- 100
depths <- unique(round(exp(seq(log(1), log(min_depth), length.out = n_points))))
depths <- depths[depths >= 1]

# compute expected richness for each sample at these depths
raref_mat <- t(sapply(1:nrow(df), function(i) {
  rarefy(df[i, ], sample = depths)
}))
rownames(raref_mat) <- rownames(df)
colnames(raref_mat) <- depths

# plot: one curve per sample
matplot(depths, t(raref_mat), type = "l", lty = 1, xlab = "Metagenomic Read Pairs (2x140 bp)", ylab = "Expected species",xlim=c(0,3000))

# --> plateau happens before 3000 read pairs for all curves, which is <1 Mb of data 
# so while more sequencing always has potential to add more species, 
# given the species/genomes observed in this dataset, we did enough sequencing to detect those 
# but this curve makes less sense for metagenomic data, because 1 read pair is never enough for a species --
# the minimum is 4000

################################################################## 
## Figure 2 Supplement 1C
## (Histogram of taxa per subject for both data types)
################################################################## 

hist(Nspecies$n,breaks=seq(0,300,20),col=makeTransparent('black',alpha=180),ylim=c(0,35),xlab="",main="",ylab="")
hist(alpha$RSVs,add=T,breaks=seq(0,300,20),col=makeTransparent('white',alpha=100))
legend(100, 34, legend=c("Metagenomic SRGs", "16S RSVs"), ce = 0.9, fill = c(makeTransparent('black',alpha=180),col=makeTransparent('white',alpha=100)) )
title(ylab="Samples", line=2.3, cex.lab=1.2)
title(xlab="Taxa detected", line=2.3, cex.lab=1.1)



################################################################## 
## FIGURE 1 SUPPLEMENT 1 E  ### 
#   metagenomic alpha diversity (SRGs detected) by tribe
################################################################## 

alphaMetaSpecies <- data.frame(matrix(ncol = 2, nrow = 0))
samples <- unique(meta_species$Sample)

for (s in samples){
  data = subset(meta_species,meta_species$Sample==s)
  N = sum(data$Rel_ab>0)
  newrow = c(s,N)
  alphaMetaSpecies = rbind(alphaMetaSpecies,newrow)  }
colnames(alphaMetaSpecies) <- c('subject', 'Species')
alphaMetaSpecies$Species = as.numeric(alphaMetaSpecies$Species)

# Add tribe/region metadata 
alphaMetaSpecies$tribe = NA
alphaMetaSpecies$region = NA
for (row in 1:nrow(alphaMetaSpecies)){
  samplename = alphaMetaSpecies[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  if (tribe == "Gondia"){
    tribe = "Gond"
  }
  alphaMetaSpecies[row,3] = tribe
  if (tribe == "Warli"){
    alphaMetaSpecies[row,4] = "Coast"
  } else if (tribe == "Kabui"){
    alphaMetaSpecies[row,4] = "Northeast Hills"
  } else if (tribe == "Gond" | tribe == "Madia"){
    alphaMetaSpecies[row,4] = "Deccan Plateau"
  } else {alphaMetaSpecies[row,4]="Trans-Himalayas"}
}

alphaMetaSpecies$tribe = factor(alphaMetaSpecies$tribe,levels=c("Warli","Kabui","Gond","Madia","Boto","Brokpa","Balti","Purigpa"))

col_list = c('#b33f25',"#edc42f","#2f3cb5",'#acb2e6',"#7CAE00","#a5e602","#cde09d","#4b6901") 

a <- ggplot(alphaMetaSpecies, aes(x=tribe, y=Species, group=tribe)) + geom_boxplot(aes(fill=tribe)) +
  theme_cowplot() + ylab("Species detected") + theme(legend.position = "none") + xlab("") +
  scale_fill_manual(values = col_list) + ggtitle("Metagenomics")
b =a + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
b

#ggsave("Fig2_S1E.pdf", plot=b,width=1.85*1.92,height=1.85*2.21,units="in")

# statistics - one-way ANOVA 
res.aov <- aov(Species ~ tribe, data = alphaMetaSpecies) # one-way ANOVA was performed 
summary(res.aov) # p = 0.214
#p2 <- as.data.frame(TukeyHSD(res.aov)[1]) # nothing significant  







