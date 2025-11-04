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


############################################################################################
## GOAL: 
#  ONE DATA FRAME WITH COLUMNS c("Sample","Tribe","Region","Dairy (category)", "function (relative abundance")  
#. TO FACET FOR FIGURE 5A
############################################################################################



#### FIRST, GATHER THE GO PATHWAYS:
  # 3 GO pathways in Figure 5, plus 2 in Figure 5 Supplement 1

####################################################################
####      GO functional abundance data  
#####################################################################



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


##################################################################################################################
####     PULL OUT THE 5 GO PATHWAYS, one at a time 
####################################################################################################################

####  GO:0004459: [MF] L-lactate dehydrogenase activity
# create df for this function
func <- d[d$GeneFamily %like% "GO:0004459", ]  
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,]
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 


new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified. 'total' thus contains the RA of the pathway for each sample, among all the pathways in that sample 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names
df$Sample = row.names(df)
lactate_dehydrogenase = as.data.frame(as.matrix(cbind(row.names(df),df$total)))  
names(lactate_dehydrogenase) = c("Sample","lactate_dehydrogenase")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

func <- d[d$GeneFamily %like% "GO:0006012", ] # galactose metabolic process
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,]
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified. total thus contains the RA of the pathway for each sample, among all the pathways in that sample 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names
df$Sample = row.names(df)
galactose_GO = as.data.frame(as.matrix(cbind(row.names(df),df$total)))  
names(galactose_GO) = c("Sample","galactose_GO")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

func <- d[d$GeneFamily %like% "GO:0006007", ] # glucose catabolic process 
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,]
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified. total thus contains the RA of the pathway for each sample, among all the pathways in that sample 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names
df$Sample = row.names(df)
glucose_GO = as.data.frame(as.matrix(cbind(row.names(df),df$total)))  
names(glucose_GO) = c("Sample","glucose_GO")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

func <- d[d$GeneFamily %like% "GO:0005355", ] # glucose transmembrane transporter 
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,]
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified. total thus contains the RA of the pathway for each sample, among all the pathways in that sample 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names
df$Sample = row.names(df)
glucose_transport = as.data.frame(as.matrix(cbind(row.names(df),df$total)))  
names(glucose_transport) = c("Sample","glucose_transport")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

func <- d[d$GeneFamily %like% "GO:0005988", ] # lactose metabolic processs
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,]
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$GeneFamily ### FOR GO 
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified, i.e., the based GO:00XXXX UNstratified. total thus contains the RA of the pathway for each sample, among all the pathways in that sample 
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names
df$Sample = row.names(df)
lactose_metabolic = as.data.frame(as.matrix(cbind(row.names(df),df$total)))  
names(lactose_metabolic) = c("Sample","lactose_metabolic")




##################################################################################################################
####      Metacyc functional abundance data  
###################################################################################################################

#### THEN, GATHER THE TWO METACYC PATHWAYS ###

# file created with these commands:
# humann2_renorm_table --input  ~/humann2-2.8.2/AK_SK_10_output/NODUP_HMN_UNMAPPED_TRIM_MARKED_AK_SK_10_pathabundance.tsv --units relab --output ~/humann2-2.8.2/pathRA/AK_SK_10_pathabundance_RA.tsv
# humann2_join_tables --input ~/humann2-2.8.2/pathRA --output ~/humann2-2.8.2/pathRA/metaCyc_100.tsv

d<-read.delim("metaCyc_100.tsv", header = TRUE, sep = "\t",numerals="allow.loss",colClasses=c('character',rep("numeric",100)))
names1 <- names(d)[2:length(names(d))]
names2 <- sub('NODUP_HMN_UNMAPPED_TRIM_MARKED_', '', names1)
names3 <- sub('_Abundance', '', names2)
names(d)[2:length(names(d))] <- names3
d <- d[ , -which(names(d) %in% c(outliers,others))]
#
# use metadata to exclude enrichment samples 
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- subset(meta,meta$sample %in% names(d))
fecalkeep = as.character(subset(meta,meta$sample_type=='fecal')$sample)
d <- d[ , which(names(d) %in% c(fecalkeep,'Pathway'))] # FOR METACYC
dsamps <- names(d)[2:length(names(d))] 
meta<-subset(meta,meta$sample %in% dsamps) # 67 samples
 


##################################################################################################################
####     PULL OUT TOTALS OF 2 METACYC PATHWAYS 
####################################################################################################################


### PATHWAY ONE: Lactose and galactose degradation I
func <- d[d$Pathway %like% "LACTOSECAT-PWY", ] # Lactose and galactose degradation I
df0 <- as.data.frame(t(func),stringsAsFactors = FALSE)
df = df0[-1,]
df[] <- lapply(df, type.convert, as.is = TRUE, dec='.')
names(df) = func$Pathway 
new_names <-str_split_fixed(names(df), fixed("g__"), 2)[, 2]
new_names[1] = 'total' # first column is the sum of all taxa + unclassified
new_names[length(new_names)]  = 'unclassified'
names(df) <- new_names 
# ## FOR METACYC, CALCULATE A STRATIFIED TOTAL (since all strats+unclassified < total, b/c of pathway definitions)
origtotal <- df$total
sum <- rowSums(df)
sum <- sum-df$total
df$total<-sum
#### SAVE
df$Sample = row.names(df)
lactose_galactose_MC = as.data.frame(as.matrix(cbind(row.names(df),df$total)))  
names(lactose_galactose_MC) = c("Sample","lactose_galactose_MC")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

### PATHWAY TWO: LELOIR

#### combine the two
func1 <- d[d$Pathway %like% "PWY-6317", ] # Leloir I
func2 <- d[d$Pathway %like% "PWY66-422", ] # Leloir V
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
### SAVE
df$Sample = row.names(df)
leloir = as.data.frame(as.matrix(cbind(row.names(df),df$total)))  
names(leloir) = c("Sample","leloir")






##################################################################################################################
####    COMBINE PATHWAYS FOR PLOTS
####################################################################################################################

# put all function dfs in the same (sample) order
lactate <- lactate[order(lactate$Sample),]
galactose_GO <- galactose_GO[order(galactose_GO$Sample),]
leloir <- leloir[order(leloir$Sample),]
glucose_GO <- glucose_GO[order(glucose_GO$Sample),]
lactose_metabolic <- lactose_metabolic[order(lactose_metabolic$Sample),]
glucose_transport  <- glucose_transport[order(glucose_transport$Sample),]
lactose_galactose_MC<- lactose_galactose_MC[order(lactose_galactose_MC$Sample),]




##################################################################################################################
####   BOXPLOTS OF FUNCTIONS BY DAIRY STATUS
## FIGURE 5A and FIGURE 5 SUPPLEMENT 1A
####################################################################################################################

### Figure 5A 

multifunction = as.data.frame(as.matrix(cbind(lactate$Sample,lactate$lactate,galactose_GO$galactose_GO,leloir$leloir,glucose_GO$glucose_GO)))              
names(multifunction) = c("Sample","L-lactate dehydrogenase","Galactose Metabolic Process","Galactose Degradation (Leloir)","Glucose Catabolic Process")


### Add tribe, region, and dairy status
multifunction$Tribe = NA
multifunction$Region = NA
multifunction$Dairy = NA

for (row in 1:nrow(multifunction)){
  samplename = multifunction[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  multifunction[row,6] = tribe
  multifunction[row,7] = region
  if (tribe == "Warli"){multifunction[row,8]="Weekly (Warli)"}
  if (region == "North"){multifunction[row,8]="Daily (Trans-Himalayas)"}
  if (region %in% c("Central","North-East")){multifunction[row,8]="Never (Gond, Madia, Kabui)"}
}

### boxplot
mfmelt <-  reshape::melt(multifunction, id = c("Sample","Tribe","Region","Dairy"))
mfmelt$value=as.numeric(mfmelt$value)
mfmelt$Dairy= factor(mfmelt$Dairy,levels=c("Daily (Trans-Himalayas)","Weekly (Warli)","Never (Gond, Madia, Kabui)"))

## Order the pathways 
mfmelt$variable =  factor(mfmelt$variable,levels=c("L-lactate dehydrogenase","Galactose Metabolic Process","Galactose Degradation (Leloir)","Glucose Catabolic Process"))

# plot
s<- ggplot(mfmelt, aes(x=variable, y=value, fill=Dairy)) +
  geom_boxplot(outlier.size=1) + theme_cowplot() +  ylab("Relative abundance") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(strip.text.x = element_text(size = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.4)),labels =  ~ sprintf(fmt = "%2.e", .)) +   #  before 030625
  theme(axis.text.y = element_text(size = 14)) +
  # scale_fill_manual(values = c("#2a8022","#6ce362", "white" ))
  scale_fill_manual(values = c("#8a1108","#db4c42", "#fcc9c5" ))
totals <- s + facet_wrap(~variable,scales="free",nrow=1, labeller = label_wrap_gen(width=10)) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
totals 

#ggsave("Fig5A.pdf", plot=totals,width=4.64*2.2*1.17,height=0.95*2.3*1.17,units="in")

######## Stats
#
## L-lactate dehydrogenase
#
path1 <- subset(mfmelt,mfmelt$variable=="L-lactate dehydrogenase")
daily <- subset(path1,path1$Dairy == "Daily (Trans-Himalayas)")
weekly <- subset(path1,path1$Dairy == "Weekly (Warli)")
never <- subset(path1,path1$Dairy == "Never (Gond, Madia, Kabui)")
t.test(daily$value,weekly$value) # 1.139e-05
t.test(daily$value,never$value) # 7.158e-07
t.test(never$value,weekly$value) # 0.7826
#
## Galactose Metabolic Process
path1 <- subset(mfmelt,mfmelt$variable=="Galactose Metabolic Process")
daily <- subset(path1,path1$Dairy == "Daily (Trans-Himalayas)")
weekly <- subset(path1,path1$Dairy == "Weekly (Warli)")
never <- subset(path1,path1$Dairy == "Never (Gond, Madia, Kabui)")
t.test(daily$value,weekly$value) #  0.01122
t.test(daily$value,never$value) #  0.003807
t.test(never$value,weekly$value) #  0.4371
#
## Galactose Degradation (Leloir)
path1 <- subset(mfmelt,mfmelt$variable=="Galactose Degradation (Leloir)")
daily <- subset(path1,path1$Dairy == "Daily (Trans-Himalayas)")
weekly <- subset(path1,path1$Dairy == "Weekly (Warli)")
never <- subset(path1,path1$Dairy == "Never (Gond, Madia, Kabui)")
t.test(daily$value,weekly$value) #  0.09187
t.test(daily$value,never$value) #   0.0002352
t.test(never$value,weekly$value) #  0.589
#
## Glucose Catabolic Process
path1 <- subset(mfmelt,mfmelt$variable=="Glucose Catabolic Process")
daily <- subset(path1,path1$Dairy == "Daily (Trans-Himalayas)")
weekly <- subset(path1,path1$Dairy == "Weekly (Warli)")
never <- subset(path1,path1$Dairy == "Never (Gond, Madia, Kabui)")
t.test(daily$value,weekly$value) #  p-value = .44
t.test(daily$value,never$value) # p = 0.006897 
t.test(never$value,weekly$value) #  0.1789







###############################
### Figure 5 Supplement 1A 

multifunction = as.data.frame(as.matrix(cbind(lactose_galactose_MC$Sample,lactose_galactose_MC$lactose_galactose_MC,lactose_metabolic$lactose_metabolic,glucose_transport$glucose_transport)))              
names(multifunction) = c("Sample","Lactose and Galactose Degradation","Lactose Metabolic Process","Glucose Transmembrane Transporter Activity")

### Add tribe, region, and dairy status
multifunction$Tribe = NA
multifunction$Region = NA
multifunction$Dairy = NA

for (row in 1:nrow(multifunction)){
  samplename = multifunction[row,1]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  multifunction[row,5] = tribe
  multifunction[row,6] = region
  if (tribe == "Warli"){multifunction[row,7]="Weekly (Warli)"}
  if (region == "North"){multifunction[row,7]="Daily (Trans-Himalayas)"}
  if (region %in% c("Central","North-East")){multifunction[row,7]="Never (Gond, Madia, Kabui)"}
}

### boxplot
mfmelt <-  reshape::melt(multifunction, id = c("Sample","Tribe","Region","Dairy"))
mfmelt$value=as.numeric(mfmelt$value)
mfmelt$Dairy= factor(mfmelt$Dairy,levels=c("Daily (Trans-Himalayas)","Weekly (Warli)","Never (Gond, Madia, Kabui)"))

## Order the pathways 
mfmelt$variable =  factor(mfmelt$variable,levels=c("Lactose and Galactose Degradation","Lactose Metabolic Process","Glucose Transmembrane Transporter Activity"))

s<- ggplot(mfmelt, aes(x=variable, y=value, fill=Dairy)) +
  geom_boxplot(outlier.size=1) + theme_cowplot() +  ylab("Relative abundance") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(strip.text.x = element_text(size = 12)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.4)),labels =  ~ sprintf(fmt = "%2.e", .)) +   #  before 030625
  theme(axis.text.y = element_text(size = 14)) +
  # scale_fill_manual(values = c("#2a8022","#6ce362", "white" ))
  scale_fill_manual(values = c("#8a1108","#db4c42", "#fcc9c5" ))
totals <- s + facet_wrap(~variable,scales="free",nrow=1, labeller = label_wrap_gen(width=10)) # see https://stats.stackexchange.com/questions/24806/dropping-unused-levels-in-facets-with-ggplot2
totals 

#ggsave("Fig5_S1A.pdf", plot=totals,width=4.64*2.2*1.17*.81,height=0.95*2.3*1.17,units="in")


######## Stats
#
## Lactose and Galactose Degradation
path1 <- subset(mfmelt,mfmelt$variable=="Lactose and Galactose Degradation")
daily <- subset(path1,path1$Dairy == "Daily (Trans-Himalayas)")
weekly <- subset(path1,path1$Dairy == "Weekly (Warli)")
never <- subset(path1,path1$Dairy == "Never (Gond, Madia, Kabui)")
t.test(daily$value,weekly$value) #  0.07036
t.test(daily$value,never$value) #  0.006129
t.test(never$value,weekly$value) #  0.3444
#
## Glucose Transmembrane Transporter Activity
path1 <- subset(mfmelt,mfmelt$variable=="Glucose Transmembrane Transporter Activity")
daily <- subset(path1,path1$Dairy == "Daily (Trans-Himalayas)")
weekly <- subset(path1,path1$Dairy == "Weekly (Warli)")
never <- subset(path1,path1$Dairy == "Never (Gond, Madia, Kabui)")
t.test(daily$value,weekly$value) #  p-value = . 0.005897
t.test(daily$value,never$value) # p =0.002119
t.test(never$value,weekly$value) # 0.7717
#
# Lactose Metabolic Process
path1 <- subset(mfmelt,mfmelt$variable=="Lactose Metabolic Process")
daily <- subset(path1,path1$Dairy == "Daily (Trans-Himalayas)")
weekly <- subset(path1,path1$Dairy == "Weekly (Warli)")
never <- subset(path1,path1$Dairy == "Never (Gond, Madia, Kabui)")
t.test(daily$value,weekly$value) #  p-value = 0.93
t.test(daily$value,never$value) # p = 0.40
t.test(never$value,weekly$value) # 0.64


##################################################################################################################