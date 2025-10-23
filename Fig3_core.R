library(data.table)
library(dplyr)
library(grid)
library("pheatmap")
library(gridExtra)
library(stringr)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)

#############################################################
#   Figure 3: Heatmap of core microbiome (16S)
#############################################################

## 16S relative abundance data
d <-read.table("Fig2_relab_genus16S.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
d$AK_SK_31 <- NULL; # this sample has no metadata (including tribe)

### Remove  outliers before defining the core
d$AK_SK_49<-NULL
d$AK_SK_10<-NULL
d$AK_SK_32<-NULL
d$AK_SR_6<-NULL
d$AK_SG_17<-NULL
d$AK_SG_18<-NULL
d$AK_SR_1<-NULL
d$AK_SR_4 <- NULL # Streptococcus 


## MELT so df matches the three-column format 
dmelt <- reshape2::melt(d, id = c("Genus"))
names(dmelt) <- c("Genus","Sample","Rel_ab") 
dmelt$Sample <- as.character(dmelt$Sample)
length(unique(dmelt$Sample))



####################
### Add tribe from metadata

meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta<-subset(meta,meta$sample %in% dmelt$Sample)
#
dmelt$tribe = NA
for (row in 1:nrow(dmelt)){
  samplename = dmelt[row,2]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  dmelt[row,4] = tribe
}
dmelt <- dmelt %>% mutate(tribe = str_replace(tribe, "Gondia", "Gond"))





############## ############## 
############## ############## 
############## ############## 

# For each taxon,
# calculate its prevalence at each RA threshold overall
# and store in a df

############## ############## 
############## ############## 
############## ############## 


header = c("Taxon","any","point001","point0025","point005","point0075","point01","point025","point05","point075","point1","Group")

taxa = unique(dmelt$Genus)
thresh = c(0.1,0.075,0.05,0.025,0.01,0.0075,0.005,0.0025,0.001,0)
thresh = rev(thresh)

regions =c("Deccan Plateau","Trans-Himalayas","Warli","Kabui")


# beginning of df to store results 
results = header

for (t in taxa){
  
  # FOR EACH TAXON
  results_taxon = rep(NA,9) # empty vector to store prevalence for this taxon (at various rel_ab thresholds)
  dmeltT <- subset(dmelt,dmelt$Genus==t) # get rel_ab data for this taxon only
  
  # for each threshold -- ALL 
  for (r in 1:length(thresh)){
    meet_thresh = subset(dmeltT,dmeltT$Rel_ab > thresh[r]) # get rel_ab data that meet the threshold 
    prevalence = nrow(meet_thresh)/length(unique(dmeltT$Sample)) # calculate prevalence
    results_taxon[r] = prevalence # store 
  }  # end r loop
  # save row for this taxon 
  row_taxon = c(t,results_taxon,"All")
  results = rbind(results,row_taxon)
  
  ## NOW DO EACH THRESHOLD FOR EACH REGION
  
  for (tr in regions){
    results_taxon_tribe = rep(NA,9) # new empty row to fill in 
    
    if (tr=="Deccan Plateau"){
      dmeltT_TRIBE = subset(dmeltT,dmeltT$tribe == "Madia" | dmeltT$tribe == "Gondia" ) 
    }
    else if (tr=="Trans-Himalayas"){
      dmeltT_TRIBE = subset(dmeltT,dmeltT$tribe == "Brokpa" | dmeltT$tribe == "Balti"  | dmeltT$tribe == "Boto"  | dmeltT$tribe == "Purigpa" ) 
    }
    else if (tr=="Warli"){
      dmeltT_TRIBE = subset(dmeltT,dmeltT$tribe == "Warli")
    }
    else if (tr=="Kabui"){
      dmeltT_TRIBE = subset(dmeltT,dmeltT$tribe == "Kabui") 
    }
    
    for (r in 1:length(thresh)){
      meet_thresh = subset(dmeltT_TRIBE,dmeltT_TRIBE$Rel_ab > thresh[r]) # get rel_ab data that meet the threshold 
      prevalence = nrow(meet_thresh)/length(unique(dmeltT_TRIBE$Sample)) # calculate prevalence
      results_taxon_tribe[r] = prevalence # store 
    }  # end r loop
    # save row for this taxon 
    row_taxon_tribe = c(t,results_taxon_tribe,tr)
    results = rbind(results,row_taxon_tribe)
  } # end tr loop
  
} # end t loop


# clean up results df 
results = as.data.frame(results)
names(results) = results[1,] 
results2 = results[-1,] # drop the redundant header



############## ############## 
############## ############## 

# COMBINE ALL REGIONS ON ONE PLOT 

############## ############## 
############## ############## 

### First, need to establish the order of the core taxa: 
plotTaxa <- c()
for (region in c("All","Warli","Deccan Plateau", "Trans-Himalayas","Kabui" )) {  # For each group in order of core size
  # subset the group
  temp <- subset(results2,results2$Group==region)
  # define its core taxa ('any' threshold)
  temp2 <- subset(temp, temp$any >= 1) #set the minimum prevalence 
  # sort by '0.1' threshold 
  temp2 = temp2[with(temp2, order(point1,decreasing = T)),]
  # the taxa are now in the right order for this group
  temp_taxa = temp2$Taxon
  for (i in temp_taxa) { if ((i %in% plotTaxa)==FALSE){ plotTaxa = c(plotTaxa,i) } }
} # end region loop  

plotTaxa # 31  taxa are at 100% prevalence in at least 1 region. This is the order of the taxa on the y-axis.



### Now, want to make each region's heatmap with the taxa in this order 

# melt results2 so that Taxon is unique 
rm <- reshape2::melt(results2,id.vars=c('Taxon','Group'))

# remove some of these thresholds to make the figure narrower
rm <- subset(rm,rm$variable != 'point0025')
rm <- subset(rm,rm$variable != 'point0075')
rm <- subset(rm,rm$variable != 'point025')
rm <- subset(rm,rm$variable != 'point075')

# informative colnames and types
names(rm) <- c("Taxon","Group","Threshold","Prevalence")
rm$Prevalence <- as.numeric(rm$Prevalence)


####### MAKE A PANEL FOR EACH GROUP
colfunc<-colorRampPalette(c("#3285B9","#AEDEA4","#FDFEBE","#FCAA5E","#DA281F"))
col=(colfunc(10))

#
rmg_All <- subset(rm,rm$Group == "All")
rg_All <-reshape2::dcast(rmg_All, Taxon~Threshold) # Uses Prevalence as value column:
row.names(rg_All) <- rg_All$Taxon
rg_All <- subset(rg_All,select=-c(Taxon))
rgt_All <- subset(rg_All,row.names(rg_All) %in% plotTaxa) ## Subset to desired taxa and order 
rgt_All <- rgt_All[match(plotTaxa, row.names(rgt_All)),]

## fix names
rgtm_All <- as.matrix(rgt_All)
#
heatmap_All <- pheatmap(rgtm_All, angle_col = 45,col=col,scale="none",cluster_rows=FALSE,cluster_cols = F,cellheight=11,cellwidth=9,
                        labels_col=c(">0",'0.001', "0.005","0.01",'0.05',"0.1" ), labels_row=rep('',length(plotTaxa)),
                        fontsize = 7, fontsize_row = 7, fontsize_col = 6, legend=F)
# to examine / identify taxa:
View(rgt_All)

#
#
rmg_Warli <- subset(rm,rm$Group == "Warli")
rg_Warli <-reshape2::dcast(rmg_Warli, Taxon~Threshold) # Uses Prevalence as value column:
row.names(rg_Warli) <- rg_Warli$Taxon
rg_Warli <- subset(rg_Warli,select=-c(Taxon))
rgt_Warli <- subset(rg_Warli,row.names(rg_Warli) %in% plotTaxa) ## Subset to desired taxa and order 
rgt_Warli <- rgt_Warli[match(plotTaxa, row.names(rgt_Warli)),]
rgtm_Warli <- as.matrix(rgt_Warli)
heatmap_Warli <- pheatmap(rgtm_Warli, angle_col = 45,col=col,scale="none",cluster_rows=FALSE,cluster_cols = F,cellheight=11,cellwidth=9,
                          labels_col=c(">0",'0.001', "0.005","0.01",'0.05',"0.1" ), labels_row=rep('',length(plotTaxa)),
                          fontsize = 7, fontsize_row = 7, fontsize_col = 6, legend=F,
                          plot.margin=unit(c(0,0,0,0), "cm")) 
#
#
rmg_TH <- subset(rm,rm$Group == "Trans-Himalayas")
rg_TH <-reshape2::dcast(rmg_TH, Taxon~Threshold) # Uses Prevalence as value column:
row.names(rg_TH) <- rg_TH$Taxon
rg_TH <- subset(rg_TH,select=-c(Taxon))
rgt_TH <- subset(rg_TH,row.names(rg_TH) %in% plotTaxa) ## Subset to desired taxa and order 
rgt_TH <- rgt_TH[match(plotTaxa, row.names(rgt_TH)),]
rgtm_TH <- as.matrix(rgt_TH)
heatmap_TH <- pheatmap(rgtm_TH, angle_col = 45,col=col,scale="none",cluster_rows=FALSE,cluster_cols = F,cellheight=11,cellwidth=9,
                       labels_col=c(">0",'0.001', "0.005","0.01",'0.05',"0.1" ), labels_row=rep('',length(plotTaxa)),
                       fontsize = 7, fontsize_row = 7, fontsize_col = 6, legend=F) 
#
#
rmg_DP <- subset(rm,rm$Group == "Deccan Plateau")
rg_DP <-reshape2::dcast(rmg_DP, Taxon~Threshold) # Uses Prevalence as value column:
row.names(rg_DP) <- rg_DP$Taxon
rg_DP <- subset(rg_DP,select=-c(Taxon))
rgt_DP <- subset(rg_DP,row.names(rg_DP) %in% plotTaxa) ## Subset to desired taxa and order 
rgt_DP <- rgt_DP[match(plotTaxa, row.names(rgt_DP)),]
rgtm_DP <- as.matrix(rgt_DP)
heatmap_DP <- pheatmap(rgtm_DP, angle_col = 45,col=col,scale="none",cluster_rows=FALSE,cluster_cols = F,cellheight=11,cellwidth=9,
                       labels_col=c(">0",'0.001', "0.005","0.01",'0.05',"0.1" ), labels_row=rep('',length(plotTaxa)),
                       fontsize = 7, fontsize_row = 7, fontsize_col = 6, legend=F) 
#
#
rmg_Kabui <- subset(rm,rm$Group == "Kabui")
rg_Kabui <-reshape2::dcast(rmg_Kabui, Taxon~Threshold) # Uses Prevalence as value column:
row.names(rg_Kabui) <- rg_Kabui$Taxon
rg_Kabui <- subset(rg_Kabui,select=-c(Taxon))
rgt_Kabui <- subset(rg_Kabui,row.names(rg_Kabui) %in% plotTaxa) ## Subset to desired taxa and order 
rgt_Kabui <- rgt_Kabui[match(plotTaxa, row.names(rgt_Kabui)),]
rgtm_Kabui <- as.matrix(rgt_Kabui)
heatmap_Kabui <- pheatmap(rgtm_Kabui, angle_col = 45,col=col,scale="none",cluster_rows=FALSE,cluster_cols = F,cellheight=11,cellwidth=9,
                          labels_col=c(">0",'0.001', "0.005","0.01",'0.05',"0.1" ),labels_row=rep('',length(plotTaxa)),
                          fontsize = 7, fontsize_row = 7, fontsize_col = 6, legend=F) 
#


#### Arrange all five heatmaps in one image
grid.arrange(heatmap_All[[4]],heatmap_Warli[[4]], heatmap_DP[[4]], heatmap_TH[[4]], heatmap_Kabui[[4]], nrow = 1)

# To get the row labels, replot the last region with labels 
for_labels <- pheatmap(rgtm_Kabui, angle_col = 45,col=col,scale="none",cluster_rows=FALSE,cluster_cols = F,cellheight=11,cellwidth=9,
                          labels_col=c(">0",'0.001', "0.005","0.01",'0.05',"0.1" ),
                          fontsize = 7, fontsize_row = 7, fontsize_col = 6, legend=F) 

### TO DRAW THE BLACK BOXES ...
### VIEW e.g. 'rgtm_TH'
### to see WHICH TAXA ARE ACTUALLY CORE IN EACH GROUP 