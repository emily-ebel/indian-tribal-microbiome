library(Maaslin2)
library(readr)
library(tidyverse)
library(data.table)
library(stringi)
library(pheatmap)
library(viridis)
library(tibble)
library(tidyr)
library(dplyr)

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)


#####################################################################
##
##  Run Maaslin model; plot Figure 4B and Figure 4 Supplement 1
##
#####################################################################



########################################
## 16S genus-level abundance data -- needs some modification before Maaslin input
# 
# d<-read.table("Fig2_relab_genus16S.csv", header = TRUE, sep = " ",stringsAsFactors = F,numerals="allow.loss")
# d2 = as.data.frame(t(d)); names(d2) = d2[1,]; d2 = d2[-1,] # transpose
# 
# # one issue: unacceptable characters are in genus names
# namestofix = names(d2)
# namestofix =  stri_replace_all_fixed(namestofix,'[','')
# namestofix =  stri_replace_all_fixed(namestofix,']','')
# namestofix =  stri_replace_all_fixed(namestofix,' ','_')
# d3 = d2
# names(d3) = namestofix
# 
# # another issue: numbers are class character
# d4 <- mutate_all(d3, function(x) as.numeric(as.character(x)))
# 
# # Save the updated df to file 
# d4_out <- cbind(ID = rownames(d4), d4)
# write.table(d4_out,file="Fig4_relab_genus16S_Maaslin.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

df_IM_data = read.table(file             = "Fig4_relab_genus16S_Maaslin.tsv", # table written above with write.table
                        header           = TRUE,
                        sep              = "\t", 
                        row.names        = 1,
                        stringsAsFactors = FALSE) # 68 samples, 166 genera, same as Fig2_relab_genus16S.csv


# Exclude outliers
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2","AK_SR_4") # 8 high-Bacteroides + Streptococcus 
df_IM_data <- subset(df_IM_data,(row.names(df_IM_data)%in%outliers==FALSE)) # leaves 68 non-outliers 
df_IM_data[1:5, 1:5]


########################################
## Metadata table 
   # customized for Maaslin
   # already excludes outlier samples 

df_IM_metadata = read.table(file             = "Maaslin2_metadata_011625.tsv", 
                            header           = TRUE, 
                            sep              = "\t", 
                            row.names        = 1,
                            stringsAsFactors = FALSE)
df_IM_metadata[1:5, ]


########################################
##  RUN THE MODEL (or read in previous results)
# fit_data_IM = Maaslin2(input_data     = df_IM_data, 
#                        input_metadata = df_IM_metadata, 
#                        min_prevalence = 0,
#                        normalization  = "NONE",
#                        output         = "Maaslin_output_IM", 
#                        fixed_effects  = c("Region"),
#                        reference      = c("Region,TransHimalayas"))
# save the results 
#saveRDS(fit_data_IM, file = "fit_data_IM.rds")

# read in previous results
fit_data_IM  <- readRDS("fit_data_IM.rds")

########################################
##  Parse results

all_results<-fit_data_IM$results 

sig_results <- subset(all_results,all_results$qval <= 0.05) 
length(unique(sig_results$feature)) # 58 genera significantly different between TH and at least one other region (q≤0.05)

# taxa more abundant in TH (negative coefficient)
sig_neg_results <- subset(sig_results,sig_results$coef<0) 
length(unique(sig_neg_results$feature)) # 42/58 significantly more abundant in TH
sort(unique(sig_neg_results$feature)) 

# taxa significantly more abundant in TH than all 3 regions 
sub1 <- subset(sig_neg_results,sig_neg_results$qval <0.05) 
e <- sub1 %>% count(feature)
e2 <- subset(e,e$n>=3)
nrow(e2) # 18 
sub2 <- subset(sig_neg_results,sig_neg_results$feature%in%e2$feature)
unique(sub2$feature)
plot_taxa = unique(sub2$feature) # 18 taxa plotted in Figure 4B


#################################################################
# Plot yes/no literature heatmap of these 18 taxa for Fig 4B (right)

# this handmade table summarizes Table S3 
litdf<-read.table("dairy_associations_figure_022825.csv", header = TRUE, sep = ",",stringsAsFactors = F) 
row.names(litdf) <- litdf$Taxon
colnames(litdf) <- c("Taxon","Fermented dairy","Grows on lactose","Raw milk","Cattle rumen","Cattle feces","Cattle reproductive")
litdf <- select(litdf, -c(Taxon)) 
litmat <- as.matrix(sapply(litdf, as.numeric))  

pheatmap(litmat, angle_col = 45,scale="none",col=viridis(3),cluster_rows=FALSE,cluster_cols = F,cellheight=14,cellwidth=14,
         labels_row =  row.names(litmat), fontsize = 9, fontsize_row = 9, fontsize_col = 9)



#################################################################
# Plot rel_ab heatmap of these 18 taxa for Fig 4B (left)

## Melt abundance data into 3 columns: Sample, Genus, Rel_ab
d16melt <- df_IM_data %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Genus",
    values_to = "Rel_ab"
  )

# Shorten one genome name so it fits on the final plot
#d16melt <- d16melt %>% mutate(Genus = str_replace(Genus, "Lachnospiraceae_NK4A136_group", "Lachnospiraceae_NK4A136"))

## Add tribe and region from metadata
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))
meta<-subset(meta,meta$sample %in% d16melt$Sample)
#
d16melt$tribe = NA
d16melt$region = NA
for (row in 1:nrow(d16melt)){
  samplename = as.character(d16melt[row,1])
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  d16melt[row,4] = tribe
  d16melt[row,5]= region
}


### Subset rel_ab data to the 18 taxa plotted in Fig 4B
dplot <- subset(d16melt,d16melt$Genus %in% plot_taxa)

## Summarize by region 
dplot2 <- dplot %>% group_by(region, Genus)  %>% summarise(mean = mean(Rel_ab))

## Further prepare df for plotting 
dplot2$mean <- as.numeric(dplot2$mean) # rel_ab must be numeric
dplot3 <-reshape2::dcast(dplot2, region~Genus) # use mean as value (instead of separate column)
dplot4 <- as.data.frame(t(dplot3))
names(dplot4) = dplot4[1,]
dplot4 <- subset(dplot4,dplot4$Coast != "Coast") # remove row 1 that repeats header
dplot4$Coast = as.numeric(dplot4$Coast) # make values numeric again
dplot4$`Deccan Plateau` = as.numeric(dplot4$`Deccan Plateau`) # make values numeric again
dplot4$`Northeast Hills` = as.numeric(dplot4$`Northeast Hills`) # make values numeric again
dplot4$`Trans-Himalayas` = as.numeric(dplot4$`Trans-Himalayas`) # make values numeric again
dplot4 <- dplot4[ order(dplot4$`Trans-Himalayas`,decreasing = T), ] # sort by abundance in TH  
row.names(dplot4)[7] = "Parolsenella / Libanicoccus"  # edit taxon name

# order regions 
d5 <- as.data.frame(cbind(dplot4$`Northeast Hills`,dplot4$`Deccan Plateau`,dplot4$Coast,dplot4$`Trans-Himalayas`))
names(d5) = c("Northeast Hills","Deccan Plateau","Coast","Trans-Himalayas")
dm <- as.matrix(d5)



####### MAKE THE HEATMAP

# Position the breaks at the quantiles of the data, so each color will represent an equal proportion of the data
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
dm_num = sort(as.vector(dm))
dm_num_min = dm_num[1:61]
mat_breaks <- quantile_breaks(dm_num_min, n = 8)
myb <- c(mat_breaks,0.01,0.2,0.5,0.1)


pheatmap(dm, angle_col = 45,scale="none",col=viridis(length(myb) - 4 ),breaks=mat_breaks,cluster_rows=FALSE,cluster_cols = F,cellheight=14,cellwidth=14,
         labels_col=c("Northeast Hills","Deccan Plateau","Coast","Trans-Himalayas"),
         labels_row =  row.names(dplot4), fontsize = 9, fontsize_row = 9, fontsize_col = 9,legend_breaks=myb)

# save this image as a PDF; make legend manually in stacked_bar_figure.pptx (in 'figure 4 taxa and TH')