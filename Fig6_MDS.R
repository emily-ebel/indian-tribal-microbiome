

path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)


##### ##### ##### ##### ##### ##### ##### ##### #####
###   Figure 6A 
###   MDS OF BIFIDOBACTERIUM ADOLESCENTIS SNVS   
##### ##### ##### ##### ##### ##### ##### ##### #####

# SNVs came from aligning (with nucmer) each MAG/genome to the ATCC 15703 reference, producing .delta files; then
# delta-filter -r -q -1 ref_qry.delta > ref_qry.filter
# show-snps -Clr ref_qry.filter > ref_qry.snps

# the SNVs were further filtered and arranged into the following table using 
#    Fig6_mummer_allele_table_032525.ipynb      
# and then subsetted to the 132 genomes used for plotting 

snps<-read.table("Fig6_Bif-adol-SNVs.txt",header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss")
snps2 <- as.data.frame(t(snps))
names(snps2) = snps2[1,]
snps2 <- snps2[-1, ]
snps4 <- data.frame(apply(snps2, 2, function(x) as.numeric(as.character(x))))
row.names(snps4) = row.names(snps2)

# need to make sure all of these columns (SNV positions) are phylogenetically informative
# so -- count how many 0s, 1s, and NAs are in each columns (position)
counts <- data.frame(
  column = names(snps4),
  n0 = colSums(snps4 == 0, na.rm = TRUE),
  n1 = colSums(snps4 == 1, na.rm = TRUE),
  nNA = colSums(is.na(snps4))
)

# this shows that there are positions that are all 0/NA, or all 0/NA except a single 1
# these cols are uninformative and should be removed
temp <- subset(counts,counts$n1<2)
snps5 <- snps4[, !names(snps4) %in% temp$column]

# now we have a set of 104,289 informative SNVs across 132 genomes (plus ATCC 15703)


##########################################
##### Read in Metadata (country, region, type) #######
##########################################

names <-read.table("Fig6_Bif-adol_metadata.txt", header = TRUE, sep = "\t",stringsAsFactors = F) 

### make sure row.names(snps5) are in the same order as names (which contains metadata)
snps5 <- snps5[ order(row.names(snps5)), ]
names <- names[ order(names$R_name), ]
row.names(snps5)==names$R_name


#####################
##### DO MDS #######
#####################

## Calculate distance between genomes using SNP matrix
distfull = as.matrix(dist(snps5))

## perform MDS using cmdscale
mds.arcsine <- cmdscale(distfull)

## organize for plotting:
x <- mds.arcsine[, 1]
y <- mds.arcsine[, 2]



##################### #######
##### Add metadata and filter 
##################### #######

plotdf <- as.data.frame(cbind(x,y,names$Country,names$Region,names$Type,names$Plot_geo))
names(plotdf) <- c("x","y","Country","Region","Type","Plot_geo")
plotdf$x = as.numeric(plotdf$x); plotdf$y = as.numeric(plotdf$y)



### Set levels and colors for Location legend
plotdf$Plot_geo = factor(plotdf$Plot_geo, levels=c("India","Nepal","Mongolia","China", "USA","Europe","El Salvador","Hadza"))
col_list = c("#53B400","#b1f576","#C49A00","#f0ce54","#f59ae0","#f50083","#00B6EB","#8f71f5") 
plotdf$Type = factor(plotdf$Type,levels=c("MAG","isolate"))


### PLOT 
x = ggplot(data = plotdf, aes(y = y, x = -x, colour=Plot_geo, fill = Plot_geo,shape=Type)) +
  geom_point(size=3) + geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  theme_cowplot() + theme(legend.title=element_blank()) +
  scale_shape_manual(values=c(16,9)) +
  xlab("MDS1") + ylab("MDS2") + #theme(legend.position = "none") +
  theme(legend.text=element_text(size=12)) + 
  theme(plot.margin = margin(10,10,10,10, "pt")) # top, right, bottom, left +
x

#ggsave("Fig6A.pdf", plot=x,width=3.05*1.7,height=2.47*1.7,units="in")

