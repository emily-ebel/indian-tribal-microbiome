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
head(eig.val16)  # 34.6, 15.4

# construct data frame for plot 
dimsToplot16 <- as.data.frame(cbind(coords16$Dim.1,coords16$Dim.2,meta16$Region,as.character(meta16$Tribe)))
names(dimsToplot16) <- c("PC1","PC2","Region","Tribe")
dimsToplot16$PC1 <- as.numeric(dimsToplot16$PC1); dimsToplot16$PC2 <- as.numeric(dimsToplot16$PC2)
row.names(dimsToplot16) = meta16$sample


######
# PLOT WITH REGIONAL HULLS -- Figure 2B
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
#
adonis2(bc_dist.arcsine_nooutliers ~  meta16_nooutliers$Region, by='terms', perm=100001)
# p = 1e-5 , R2 = 0.157 excluding 8 high-Bacteroides outliers
# verified 11/4/25
                                                                                          


######
# PLOT WITH TRIBAL HULLS -- Fig 2 Supp 1F
#####
col_list<-c("#cde09d","#7CAE00","#a5e602",'#acb2e6',"#edc42f","#2f3cb5","#4b6901",'#b33f25')
hulls_16S_tribe <- ddply(dimsToplot16_NOOUTLIERS, "Tribe", find_hull)

u = ggplot(data = dimsToplot16, aes(y = PC1, x = -PC2, colour=Tribe, fill = Tribe)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)   + ggtitle("16S RSVs")+
  theme_cowplot() +  
  ylab("PC1 (34.6%)") + xlab("PC2 (15.4%)")+
  geom_polygon(data = hulls_16S_tribe, alpha = 0.25)  +
  theme( legend.justification = "top" )
u
# calculate significance by permutation (excluding outliers)
adonis2(bc_dist.arcsine_nooutliers ~ meta16_nooutliers$Tribe, by='terms', perm=10001) 
# p = 1e-5 , R2 = 0.2198, verified 110425

ggsave("Fig2_S3B.pdf", plot=u,height=1.94*2,width=1.4*2.7,units="in") 



################################################
#  are TH tribes significantly different from each other?
# 16S
################################################

# subset metadata 
meta_th <- subset(meta16, meta16$Region == "Trans-Himalayas") # N = 34, including outliers 

# subset relab data
d16_th <- subset(d16,row.names(d16) %in% meta_th$sample) # N = 34, including outliers 

# arcsin transform relative abundance data
d16_th.arcsine = asin(sqrt(d16_th)) 

# calculate distance and center 
bc_dist.d16_th.arcsine = vegan::vegdist(d16_th.arcsine, method = "bray") 
dfc = as.data.frame(scale(bc_dist.d16_th.arcsine,center=T))  

# perform PCA
pca_bcdist16 <- dudi.pca(dfc,nf=5,scale=T,center=T,scannf=F)

# extract sample x-y coordinates and eigenvalues 
res.ind16 <- get_pca_ind(pca_bcdist16)
coords16 <- res.ind16$coord 
eig.val16 <- get_eigenvalue(pca_bcdist16)
head(eig.val16)  # PC1 41.4%, PC2 13.3%, ratio = 3.11

# construct data frame for plot 
dimsToplot16 <- as.data.frame(cbind(coords16$Dim.1,coords16$Dim.2,as.character(meta_th$Tribe)))
names(dimsToplot16) <- c("PC1","PC2","Tribe")
dimsToplot16$PC1 <- as.numeric(dimsToplot16$PC1); dimsToplot16$PC2 <- as.numeric(dimsToplot16$PC2)
row.names(dimsToplot16) = meta_th$sample


# find the hull WITHOUT the high-Bacteroides outlier
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2") # 8 high-Bacteroides 
dimsToplot16_NOOUTLIERS = subset(dimsToplot16,row.names(dimsToplot16) %in% outliers == FALSE)
hulls_16S_tribe <- ddply(dimsToplot16_NOOUTLIERS, "Tribe", find_hull)

col_list<-c("#cde09d","#7CAE00","#a5e602","#4b6901")

x = ggplot(data = dimsToplot16, aes(y = PC1, x = -PC2, colour=Tribe, fill = Tribe)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)   + ggtitle("16S RSVs")+
  theme_cowplot() +  
  ylab("PC1 (41.4%)") + xlab("PC2 (13.3%)")+ 
  geom_polygon(data = hulls_16S_tribe, alpha = 0.25) +
  theme(
  #  legend.position = "top"
    legend.justification = "top"
  )
x
#ggsave("Fig2_S1H.pdf", plot=x,height=1.94*2,width=1.4*2.7,units="in") 

# STATS - calculate significance by permutation (excluding outliers)
d16_th_nooutliers = subset(d16_th.arcsine,(row.names(d16_th.arcsine) %in% outliers)==F)
meta_th_nooutliers = subset(meta_th,(meta_th$sample %in% outliers)==F)
#df.arcsine_th_nooutliers = asin(sqrt(d16_th_nooutliers)) # arcsine transform rel abs
bc_dist_th_nooutliers = vegan::vegdist(d16_th_nooutliers, method = "bray") # calculate distance
adonis2(bc_dist_th_nooutliers ~  meta_th_nooutliers$Tribe, by='terms', perm=10001) # p = 0.0180 , R2 = 0.134 excluding 8 high-Bacteroides outliers




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
# PLOT WITH REGIONAL HULLS, excluding outliers (Figure 2C)
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
#
adonis2(bc_dist.arcsine_nooutliers ~  meta_nooutliers$Tribe, by='terms', perm=10001) 
# p < 1e-5 , R2 = 0.191, 10/21/25

adonis2(bc_dist.arcsine_nooutliers ~  meta_nooutliers$Region, by='terms', perm=10001) 
# p = 1 e -4, R2 = 0.13, 11/10/25



######
# PLOT WITH TRIBAL HULLS, excluding outliers -- Figure 2 Supplement 2G
#####
col_list<-c("#cde09d","#7CAE00","#a5e602",'#acb2e6',"#edc42f","#2f3cb5","#4b6901",'#b33f25')
hulls_meta_tribe <- ddply(dimsToplot_NOOUTLIERS, "Tribe", find_hull)
q=ggplot(data = dimsToplot, aes(y = PC1, x = PC2, colour=Tribe, fill = Tribe)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)  + 
  theme_cowplot() +  #theme(legend.position = "none") + theme(legend.title=element_blank()) +
  ylab("PC1 (33.2%)") + xlab("PC2 (13.0%)")+    ggtitle("Metagen. Species")+
  geom_polygon(data = hulls_meta_tribe, alpha = 0.25) +
  theme(legend.justification = "top" )
q

#ggsave("FigS2G.pdf", plot=q,height=1.94*2,width=1.4*2.7,units="in") 


################################################
#  are TH tribes significantly different from each other?
# metagenomics
################################################

# subset metadata
meta_th <- subset(meta, meta$Region == "Trans-Himalayas")

# subset relab data
df_th <- subset(df,row.names(df) %in% meta_th$sample)
meta_th <- subset(meta_th, meta_th$sample %in% row.names(df_th)) # still includes outliers

# arcsin transform relative abundance data
df_th.arcsine = asin(sqrt(df_th)) 

# calculate distance and center 
bc_dist.df_th = vegan::vegdist(df_th.arcsine, method = "bray") 
dfc = as.data.frame(scale(bc_dist.df_th,center=T))  

# perform PCA
pca_meta_th <- dudi.pca(dfc,nf=5,scale=T,center=T,scannf=F)

# extract sample x-y coordinates and eigenvalues 
res.indf <- get_pca_ind(pca_meta_th)
coords16 <- res.indf$coord 
eig.val16 <- get_eigenvalue(pca_meta_th)
head(eig.val16)  # PC1 31.9%, PC2 13.1%, verified 110425


# construct data frame for plot 
dimsToplotmeta <- as.data.frame(cbind(coords16$Dim.1,coords16$Dim.2,meta_th$Region,as.character(meta_th$Tribe)))
names(dimsToplotmeta) <- c("PC1","PC2","Region","Tribe")
dimsToplotmeta$PC1 <- as.numeric(dimsToplotmeta$PC1); dimsToplotmeta$PC2 <- as.numeric(dimsToplotmeta$PC2)
row.names(dimsToplotmeta) = meta_th$sample


# find the hull WITHOUT the high-Bacteroides outlier
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2") 
dimsToplotmeta_NOOUTLIERS = subset(dimsToplotmeta,row.names(dimsToplotmeta) %in% outliers == FALSE)
hulls_16S_tribe <- ddply(dimsToplotmeta_NOOUTLIERS, "Tribe", find_hull)


col_list<-c("#cde09d","#7CAE00","#a5e602","#4b6901")

x = ggplot(data = dimsToplotmeta, aes(y = PC1, x = -PC2, colour=Tribe, fill = Tribe)) +
  geom_point() + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)   + ggtitle("Metagen. Species")+
  theme_cowplot() +  
  ylab("PC1 (31.9%)") + xlab("PC2 (13.1%)")+ 
  geom_polygon(data = hulls_16S_tribe, alpha = 0.25) +
  theme(legend.justification = "top")
x
ggsave("Fig2_S1I.pdf", plot=x,height=1.94*2,width=1.4*2.7,units="in") 

# STATS - calculate significance by permutation (excluding outliers)
df_th_nooutliers = subset(df_th,(row.names(df_th) %in% outliers)==F)
meta_th_nooutliers = subset(meta_th,(meta_th$sample %in% outliers)==F)
df.arcsine_th_nooutliers = asin(sqrt(df_th_nooutliers)) # arcsine transform rel abs
bc_dist_th_nooutliers = vegan::vegdist(df.arcsine_th_nooutliers, method = "bray") # calculate distance
adonis2(bc_dist_th_nooutliers ~  meta_th_nooutliers$Tribe, by='terms', perm=10001)
# p = 0.1862 , R2 = 0.118
