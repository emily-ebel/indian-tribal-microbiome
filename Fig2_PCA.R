library("dendextend")
library('ggplot2')
library("vegan")
library(ade4)
library('factoextra')
library(dplyr)
library(plyr)
library(stringr)
library(cowplot)
library(grid)
library(gridExtra)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]



#############################################################
# PCA of 16S RSV relative abundance
  # by region, Figure 2B
  # by tribe, Figure 2 Supplement 1E
#############################################################

## 16S relative abundance data
d16 <-read.table("Fig2_relab_RSV.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
names(d16)[1] = "RSV" 
d16$AK_SK_31 <- NULL; # this sample has no metameata
row.names(d16) <- d16$RSV
d16 <- subset(d16,select=-c(RSV))
d16 <- as.data.frame(t(d16))


## Metadata 
meta16 <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T)
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta16 <- meta16 %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta16 <- meta16 %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))

## Match metadata samples to 16S samples & order
meta16<-subset(meta16,meta16$sample %in% row.names(d16))
d16 <- subset(d16,row.names(d16) %in% meta16$sample)  
d16 <- d16[ order(row.names(d16)), ]
meta16 <- meta16[order(meta16$sample),]


# arcsin transform relative abundance data, as proposed by Morgan et al. for Maaslin and other linear models
df.arcsine = asin(sqrt(d16)) 

# calculate distance and center 
bc_dist.arcsine = vegan::vegdist(df.arcsine, method = "bray") 
dfc = as.data.frame(scale(bc_dist.arcsine,center=T))  

# perform PCA
pca_bcdist16 <- dudi.pca(dfc,nf=5,scale=T,center=T,scannf=F)

# extract sample x-y coordinates and eigenvalues 
res.ind16 <- get_pca_ind(pca_bcdist16)
coords16 <- res.ind16$coord 
eig.val16 <- get_eigenvalue(pca_bcdist16)
head(eig.val16) 

# construct data frame for plot 
dimsToplot16 <- as.data.frame(cbind(coords16$Dim.1,coords16$Dim.2,meta16$Region,as.character(meta16$Tribe)))
names(dimsToplot16) <- c("PC1","PC2","Region","Tribe")
dimsToplot16$PC1 <- as.numeric(dimsToplot16$PC1); dimsToplot16$PC2 <- as.numeric(dimsToplot16$PC2)
row.names(dimsToplot16) = meta16$sample


######
# PLOT WITH REGIONAL HULLS
#####
col_list<-c("#b33f25", "#2f3cb5",   "#edc42f","#2a8022")
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2") # 8 high-Bacteroides 
dimsToplot16_NOOUTLIERS = subset(dimsToplot16,row.names(dimsToplot16) %in% outliers == FALSE)
hulls_16S_region <- ddply(dimsToplot16_NOOUTLIERS, "Region", find_hull)

ggplot(data = dimsToplot16, aes(y = PC1, x = -PC2, colour=Region, fill = Region)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)  + 
  theme_cowplot() + ggtitle("16S RSVs")+ theme(plot.title = element_text(size = 15, face = "plain"))+
  ylab("PC1 (34.6%)") + xlab("PC2 (15.4%)")+ theme(legend.title=element_blank()) +
  geom_polygon(data = hulls_16S_region, alpha = 0.25) +  theme(legend.position = "none") +
  theme(legend.text=element_text(size=11)) + 
  scale_x_continuous(breaks = c(-5,0,5)  ) 

# calculate significance by permutation (excluding outliers)
d16_nooutliers = subset(d16,(row.names(d16) %in% outliers)==F)
meta16_nooutliers = subset(meta16,(meta16$sample %in% outliers)==F)
df.arcsine_nooutliers = asin(sqrt(d16_nooutliers)) # arcsine transform rel abs
bc_dist.arcsine_nooutliers = vegan::vegdist(df.arcsine_nooutliers, method = "bray") # calculate distance
adonis2(bc_dist.arcsine_nooutliers ~  meta16_nooutliers$Region, by='terms', perm=100001) # p = 1e-5 , R2 = 0.157 excluding 8 high-Bacteroides outliers


######
# PLOT WITH TRIBAL HULLS
#####
col_list<-c("#cde09d","#7CAE00","#a5e602",'#acb2e6',"#edc42f","#2f3cb5","#4b6901",'#b33f25')
hulls_16S_tribe <- ddply(dimsToplot16_NOOUTLIERS, "Tribe", find_hull)

ggplot(data = dimsToplot16, aes(y = PC1, x = -PC2, colour=Tribe, fill = Tribe)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)   + ggtitle("16S RSVs")+
  theme_cowplot() +  
  ylab("PC1 (34.9%)") + xlab("PC2 (15.0%)")+
  geom_polygon(data = hulls_16S_tribe, alpha = 0.25) 

# calculate significance by permutation (excluding outliers)
adonis2(bc_dist.arcsine_nooutliers ~ meta16_nooutliers$Tribe, by='terms', perm=100001) # p = 1e-5 , R2 = 0.2198









#############################################################
# PCA of metagenomic species relative abundance
  # by region, Figure 2C
  # by tribe, Figure 2 Supplement 1G
#############################################################


## metagenomic relative abundance data
d<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 



# convert to df where rows = samples, columns = species-representative genomes, and cells contain relative abundance [0-1]
samples <- unique(d$sample)
taxa <- unique(d$genome)
df <- data.frame(matrix(ncol = length(taxa), nrow = length(samples))); names(df) <- taxa; rownames(df) <- samples
for (s in 1:length(samples)){
  sdf <- subset(d,d$sample==samples[[s]])
   sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab))  # scale rel_ab so each sample sums to 1 (i.e. ignore unmapped/unclassified reads)
  for (t in 1:length(taxa)){
    if ( (taxa[[t]] %in% sdf$genome)==TRUE){
      row <- subset(sdf,sdf$genome == taxa[[t]])
      df[s,t]= row$rel_ab_adj 
    } else (df[s,t]=0)
  }
} 


## Metadata
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))

## Match metadata samples to metagenomic samples & order
meta<-subset(meta,meta$sample %in% row.names(df))
df <- subset(df,row.names(df) %in% meta$sample)
df <- df[ order(row.names(df)), ]
meta <- meta[order(meta$sample),]


# arcsin transform relative abundance data, as proposed by Morgan et al. for Maaslin and other linear models
df.arcsine = asin(sqrt(df)) 

# calculate distance and center 
bc_dist.arcsine = vegan::vegdist(df.arcsine, method = "bray") # calculate distance 
dfc = as.data.frame(scale(bc_dist.arcsine,center=T))  # center

# perform PCA
pca_bcdist <- dudi.pca(dfc,nf=5,scale=T,center=T,scannf=F)

# extract sample x-y coordinates and eigenvalues 
res.ind <- get_pca_ind(pca_bcdist)
coords <- res.ind$coord # Get the coordinate (x-y) values for individuals
eig.val <- get_eigenvalue(pca_bcdist)
head(eig.val)  # PCA1 33.227 %, PC2 13.033 %, 10/21/25

# construct data frame for plot 
dimsToplot <- as.data.frame(cbind(coords$Dim.1,coords$Dim.2,meta$Region,as.character(meta$Tribe)))
names(dimsToplot) <- c("PC1","PC2","Region","Tribe")
dimsToplot$PC1 <- as.numeric(dimsToplot$PC1); dimsToplot$PC2 <- as.numeric(dimsToplot$PC2)
row.names(dimsToplot) = meta$sample




######
# PLOT WITH REGIONAL HULLS, excluding outliers
#####

col_list<-c("#b33f25", "#2f3cb5",   "#edc42f","#2a8022")
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2") # 8 high-Bacteroides 
dimsToplot_NOOUTLIERS = subset(dimsToplot,row.names(dimsToplot) %in% outliers == FALSE)
hulls_meta_region <- ddply(dimsToplot_NOOUTLIERS, "Region", find_hull)
#p=
ggplot(data = dimsToplot, aes(y = PC1, x = PC2, colour=Region, fill = Region)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)  + 
  theme_cowplot() +  theme(legend.position = "none") + theme(legend.title=element_blank()) +
  ggtitle("Metagenomic species") +  theme(plot.title = element_text(size = 15, face = "plain"))+
  ylab("PC1 (33.2%)") + xlab("PC2 (13.0%)")+
  geom_polygon(data = hulls_meta_region, alpha = 0.25) + 
  scale_x_continuous(breaks = c(-5,0,5)  )

#ggsave("Fig2C.pdf", plot=p,height=1.94*2,width=1.2*2,units="in") 



# calculate significance by permutation (excluding outliers)
df_nooutliers = subset(df,(row.names(df) %in% outliers)==F)
meta_nooutliers = subset(meta,(meta$sample %in% outliers)==F)
df.arcsine_nooutliers = asin(sqrt(df_nooutliers)) # arcsine transform rel abs
bc_dist.arcsine_nooutliers = vegan::vegdist(df.arcsine_nooutliers, method = "bray") # calculate distance
adonis2(bc_dist.arcsine_nooutliers ~  meta_nooutliers$Region, by='terms', perm=100001) # p < 1e-5, R2 = 0.128, 10/21/25


######
# PLOT WITH TRIBAL HULLS, excluding outliers
#####
col_list<-c("#cde09d","#7CAE00","#a5e602",'#acb2e6',"#edc42f","#2f3cb5","#4b6901",'#b33f25')
hulls_meta_tribe <- ddply(dimsToplot_NOOUTLIERS, "Tribe", find_hull)
q=ggplot(data = dimsToplot, aes(y = PC1, x = PC2, colour=Tribe, fill = Tribe)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)  + 
  theme_cowplot() +  theme(legend.position = "none") + theme(legend.title=element_blank()) +
  ylab("PC1 (33.2%)") + xlab("PC2 (13.0%)")+ 
  geom_polygon(data = hulls_meta_tribe, alpha = 0.25) 

# calculate significance by permutation (excluding outliers)
adonis2(bc_dist.arcsine_nooutliers ~  meta_nooutliers$Tribe, by='terms', perm=100001) # p < 1e-5 , R2 = 0.191, 10/21/25

ggsave("FigS2G.pdf",
       plot=q,height=1.94*2,width=1.4*2,units="in")
