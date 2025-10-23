library(dplyr)
library(RColorBrewer)
library("gplots")
require(Cairo)


path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)



#########################################
##  Figure 1C (Daily foods)  ##
#########################################

# input and format data
d <-read.csv("Fig1_Diet_table_daily.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
row.names(d) = d$Food
d <- subset(d,select=-c(Food,Sumbyfood,Sumbygroup))
df <- d %>% mutate_if(is.integer, as.numeric)
df <- as.matrix(df)

# group the foods into types in order, as coded by the color strip
my_group <- c(1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4,5,6,7,7,7,8) 
my_group = sort(my_group,decreasing=T)
colSide <- brewer.pal(9, "Set3")[my_group]

# set heatmap colors
coul <- colorRampPalette(brewer.pal(8, "Purples"))(25) 

## create plot (exported to PDF) 
CairoPDF(file = "daily_heatmap.pdf", width = 5.96, height = 5.08)
heatmap(df,scale="none",Rowv = NA,col=coul,cexRow=1.25,RowSideColors=colSide)
dev.off()

## create legend scale: **SCREENSHOT** the scale and add to figure separately (because PDF export alters the colors)
heatmap.2(as.matrix(df), scale = "none", Rowv = NA, dendrogram = "column", col = coul,
          cexRow = 1.25, RowSideColors = colSide, key = TRUE, key.title = "Value",  
          trace = "none", density.info = "none")




#########################################
##  Figure 1 Supplement 1 (all foods)  ##
#########################################

# input and format data
d <-read.csv("Fig1_Diet_table_all.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
row.names(d) = d$Food
d <- subset(d,select=-c(Food,Sumbyfood,Sumbygroup))
df <- d %>% mutate_if(is.integer, as.numeric)
df <- as.matrix(df)

# group the foods into types in order, as coded by the color strip
my_group <- c(1,1,1,1,1,1,1,1,1,1,1,1, #spices 
              2,2,2,2,2,2,2, #grains
              3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, # fruits and veg 
              4,4,4,4,4,4,4, # dairy
              5,5, # non-dairy tea
              6, # fermented foods 
              7,7,7,7,7,7,7,7, # Animal protein 
              8,8,8, # alcohol 
              9,9,9,9,9,9,9,9,9,9,9, # Legumes 
              10,10,10) # other 
my_group = sort(my_group,decreasing=T)
colSide <- brewer.pal(10, "Set3")[my_group]

# set heatmap colors
coul <- colorRampPalette(brewer.pal(8, "Purples"))(25) 

## create plot (exported to PDF) 
CairoPDF(file = "allfoods_heatmap.pdf", width = 7, height = 24)
heatmap(df,scale="none",Rowv = NA,col=coul,RowSideColors=colSide,cexRow=0.7,margins=c(8,8))
dev.off()

## create legend scale: **SCREENSHOT** the scale and add to figure separately (because PDF export alters the colors)
heatmap.2(as.matrix(df), scale = "none", Rowv = NA, dendrogram = "column", col = coul,
          cexRow = 1.25, RowSideColors = colSide, key = TRUE, key.title = "Value",  
          trace = "none", density.info = "none")