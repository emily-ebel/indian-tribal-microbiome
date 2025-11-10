## This script reads in dbCAN output (overview) for genomes from this study, FeFiFo (Wastyk 2021),
## and the B. adolescentis genomes in Table S5

## It computes a table with the  # of genes per CAZy family per Mb for each genome 

## And plots Figure 7C, 7D, and 7A 

library(stringr)
library(data.table)
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)
library(plyr)
library(ade4)
library(factoextra)
library(vegan)
library("dendextend")


# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 
#    INPUT DATA 
# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 

# Collect dbCAN output files
# Indian Tribal Microbiome + FeFiFo (California)
files1 <- list.files(path="/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/dbcan/090523", pattern="*overview.txt", full.names=TRUE, recursive=FALSE) # IM + FFF
# Main set of B. adolescentis 
files2 <- list.files(path="/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/dbcan/082324", pattern="*overview.txt", full.names=TRUE, recursive=FALSE) # Bif adol
# Other B. adolescentis from Hadza/Nepal
files3 <- list.files(path="/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/dbcan/101824_partial", pattern="*overview.txt", full.names=TRUE, recursive=FALSE) # Hadza/Nepal
#combine
files4 <- c(files1,files2,files3)



# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 
#    MAKE TABLE OF CAZY FAMS x GENOME, where cells contain  # of genes
# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 

#           GH1	GH2 	...	GH1000 
#genome1
#genome2
#...
#genomeN 

# so we need a list of all the genomes and all the CAZy families 

## genomes can be obtained from file list:
genomes = c()
for (f in files4) {
  g = strsplit(f,'dbcan/')[[1]][2] 
  g2 = strsplit(g,'.overview.txt')[[1]][1]
  genomes = c(genomes, g2)
} 

## for the CAZy families, go through all the files once
allcazy <- c()

for (f in files4) {
  
  d <- read.table(f, sep = "\t",stringsAsFactors = FALSE,numerals="allow.loss",header=FALSE)
  names(d)= c('GeneID','EC','HMMER','eCAMI','DIAMOND','NTools')
  
  # limit to rows with HHMER results only
  d<-subset(d,d$HMMER!='-')
  hcazy = d$HMMER; hcazy <- hcazy[!hcazy == '-']
  
  # PARSE HMMER RESULTS - remove the parentheses and everything between them (domains)
  hcazy2 <- gsub("\\s*\\([^\\)]+\\)","",as.character(hcazy))
  
  #  remove duplicates for now (we're just identifying all the CAZy fams, not counting genes yet)
  hcazy2 <- unique(hcazy2)
  
  # also split anything wiht '+' into 2 
  hcazy3 <- grep("\\+", hcazy2,value=TRUE) 
  hcazy4 <- c()
  for (i in hcazy3){
    both = strsplit(i,"\\+")
    hcazy4 <- c(hcazy4, both[[1]][1])
    hcazy4 <- c(hcazy4, both[[1]][2])
  }
  hcazy4 <- unique(hcazy4) 
  
  # so final HHMER vector will have hcazy4, plus things from hcazy2 with no '+'
  hcazy5 <- unique(c(grep("\\+", hcazy2,value=TRUE,invert=TRUE),hcazy4)) 
  
  allcazy <- unique(c(allcazy,hcazy5))
}

# sort all cazymes and remove ones that are just EC numbers 
allcazy2 <- sort(grep("\\.", allcazy,value=TRUE,invert=TRUE))
# --> 460 CAZy fams in 2570 genomes



###### Make the empty table
df <- data.frame(matrix(ncol = length(allcazy2), nrow = length(genomes)))
colnames(df) <- allcazy2; rownames(df) <- genomes

##### Go back through the files (genomes) and fill in the table.
## THIS TAKES SEVERAL MINUTES
for (f in files4) {
  
  # get genome name
  g = strsplit(f,'dbcan/')[[1]][2]
  g2 = strsplit(g,'.overview.txt')[[1]][1] # genome name 
  d <- read.table(f, sep = "\t",stringsAsFactors = FALSE,numerals="allow.loss",header=FALSE)
  names(d)= c('GeneID','EC','HMMER','eCAMI','DIAMOND','NTools')
  
  d<-subset(d,d$HMMER!='-')
  # clean up HMMER parentheses 
  d$HMMER <- gsub("\\s*\\([^\\)]+\\)","",as.character(d$HMMER))
  
  # Add extra columns to resolve + -- there can be up to 6 in B. adol
  d$HMMER2 = '-'
  d$HMMER3 = '-'
  d$HMMER4 = '-'
  d$HMMER5 = '-'
  d$HMMER6 = '-'
  d$HMMER7 = '-'
  d$HMMER8 = '-'
  
  #### ROW BY ROW, SPLIT '+'
  ## for each row of d,
  for (r in 1:dim(d)[[1]]){
    row = d[r,]
    if (row$HMMER %like% "\\+") { # if it has a +
      both = strsplit(row$HMMER,"\\+")# split it up
      #  save in seprate columns 
      d[r,3] = both[[1]][[1]] # HMMER
      d[r,7] = both[[1]][[2]] # HMMER2
      if (length(both[[1]])>=3) { d[r,8] = both[[1]][[3]] } # HMMER3
      if (length(both[[1]])>=4) { d[r,9] = both[[1]][[4]] } # HMMER4
      if (length(both[[1]])>=5) { d[r,10] = both[[1]][[5]] } # HMMER3
      if (length(both[[1]])>=6) { d[r,11] = both[[1]][[6]] } # HMMER4
      if (length(both[[1]])>=7) { d[r,12] = both[[1]][[7]] } # HMMER7
      if (length(both[[1]])>=8) { d[r,13] = both[[1]][[8]] } # HMMER8
      if (length(both[[1]])>8) { print(row) } 
    }
  } # end for row of d loop 
  
  ## NOW FOR EACH CAZYME, COUNT GENE MATCHES
  for (c in 1:length(allcazy2)){
    targetcazy = allcazy2[c]
    matches <- subset(d,d$HMMER==targetcazy | d$HMMER2==targetcazy | d$HMMER3 ==targetcazy  | d$HMMER4 ==targetcazy | d$HMMER5 ==targetcazy | d$HMMER6 ==targetcazy | d$HMMER7 ==targetcazy | d$HMMER8 ==targetcazy)
    # fill it in
    df[g2,targetcazy] = length(which(matches==targetcazy)) # using 'dim' skips duplicates in the same gene 
  }
}

# remove prefix folders from genome names
oldnames<-row.names(df)
newnames <- oldnames
newnames <- gsub('032525/','',newnames)
newnames <- gsub('090523/','',newnames)
newnames <- gsub('082324/','',newnames)
newnames <- gsub('101824_partial/','',newnames)
row.names(df) <- newnames



###############  Remove CAzymes in families not of interest (Justin: focus on degradation of carbohydrates, not building))
df2 <- df %>% select(-contains('cohesin'))
df2 <- df2 %>% select(-contains('SLH'))
df2 <- df2 %>% select(-contains('AA'))
df2 <- df2 %>% select(-contains('CE'))

# write it to file to save re-computing later
write.table(df2,file="/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/CAZymes_by_genome_103025.csv",sep=",",col.names=TRUE,row.names=TRUE)




# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 
#    CONVERT TO CAZYMES PER MB OF GENOME
# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 

# note: compare to CAZymes_by_genome_PERMB_111924.csv (what's actually used in paper)
# to see if certain families should be excluded, as in CAZymes_figure_update_100923.R

#### read in CAZymes_by_genome file made above
path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)
cazy <-read.table("CAZymes_by_genome_103025.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss",row.names=1) 

#### get genome sizes 
imd <-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 
fffd <-read.table("Wastyk2021_metagenomics.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 
fffd <- subset(fffd,fffd$breadth>=0.5)

### not all the genomes run through dbcan are actually analyzed for the paper (due to low breadth, for example)
### so, filter 'cazy' only for genomes that passed 

# get the genome names used in paper 
im_genomes = as.character(unique(imd$genome))
fff_genomes <- as.character(unique(fffd$genome))

# FFF genome names need to have suffixes added
fff_genomes2 <- c()
for (i in 1:length(fff_genomes)){
  oldname = fff_genomes[[i]]
  if (substr(oldname, 1, 3)=="GUT"){
    newname=paste(oldname,'.fna.fa',sep='')
  } else if (substr(oldname, 1, 4)=="META"){
    newname=paste(oldname,'.fa',sep='')} else {print(oldname)} 
  fff_genomes2 <- c(fff_genomes2,newname)
}

# combine
goodgenomes <- unique(c(im_genomes,fff_genomes2))

# subset cazy data to these genomes only 
cazy2 <- subset(cazy,row.names(cazy) %in% goodgenomes) # for 1708 genomes


###### Now, do the calculation:
  # for each row, find the genome size in Mb and divide the row by it 

cazy3 <- cazy2 # copy to version that will be  normalized by length 
yeslength = c()
nolength = c()

for (row in 1:nrow(cazy2)) {
  
  ### FIND THE GENOME SIZE IN ONE OR THE OTHER DATASETS
  genome1 = row.names(cazy2)[row]
  genomesize = 0 
  genome2= ""
  if (genome1 %in% imd$genome){
    # get the length
    imdrow = subset(imd,imd$genome == genome1)
    genomesize = unique(imdrow$length)
  } else { 
    
    # need to get the name in fff format
    if (substring(genome1,1,3) == "GUT"){ 
      genome2 = strsplit(genome1,'[.]')[[1]][[1]] } # need to drop '.fna.fa'
    if (substring(genome1,1,4) == "META"){ 
      genome2 = substring(genome1,1,nchar(genome1)-3)} # need to drop '.fa'
    
    if (genome2 %in% fffd$genome){
      # get the length
      fffdrow = subset(fffd,fffd$genome == genome2)
      genomesize = unique(fffdrow$length)
    }  } # end else 
  
  if(genomesize == 0){ # just to check we have everything
    nolength = c(nolength,genome1)} else {yeslength = c(yeslength,genome1)}
  
  ### CONVERT TO MB AND NORMALIZE THE ROW 
  genomesizeMB = genomesize/1000000
  cazy3[row,] =  cazy3[row,]/genomesizeMB
}

## save 
write.table(cazy3,file="CAZymes_by_genome_PERMB_103025.csv",sep=",",col.names=TRUE,row.names=TRUE)







# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 
#    PLOT FIGURE 7C
# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 




################
### Want to find average CAZy / Mb for each SAMPLE, given CAZy / Mb per genome & genome relative abundance
## 
## Divide this into two steps:
##  1. make genome x sample table, where cells are rel_ab
##  2. calculate a CAZy profile for each sample



########################################################################
#### STEP 1: MAKE GENOME X SAMPLE TABLES FOR TRIBES AND CALIFORNIANS 
########################################################################

#####################################
## IM abundance data 
path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)
dim<-read.table("Fig2_relab_SRG.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 

## IM metadata 
#meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = F) 
#meta<-subset(meta,meta$sample %in% dim$sample)

#### convert to genome abundance table with columns: Sample, Genome, Rel_ab
samples <- unique(dim$sample); length(samples)
genomes = unique(dim$genome); length(genomes)
genomeRA <- data.frame(matrix(ncol = 3, nrow = length(samples)*length(genomes)))
names(genomeRA) <- c("Sample","Genome","Rel_ab")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(dim,dim$sample==samples[[s]]) # for each sample 
  sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab)) # scale RA 
  
  for (t in 1:length(genomes)){ 
    count = count+1
    genomeRA[count,1] = as.character(samples[[s]])
    genomeRA[count,2]= as.character(genomes[[t]])
    
    if ( (genomes[[t]] %in% sdf$genome)==TRUE){
      rows <- subset(sdf,sdf$genome == genomes[[t]]) 
      genomeRA[count,3] = sum(rows$rel_ab_adj)   
    } else (genomeRA[count,3]=0)
  }
}







##################################
## FeFiFo abundance data
f <-read.table("Wastyk2021_metagenomics.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 
f <- subset(f,f$breadth>=0.5) # remove genomes where < 50% is covered

# remove timepoints 6 and 7
f <- subset(f,f$timepoint %in% c(1,2) )

#### convert to genome abundance table with columns: Sample, Genome, Rel_ab, Timepoint
samples <- unique(f$sample); length(samples)
genomes = unique(f$genome); length(genomes)
genomeRA_fff <- data.frame(matrix(ncol = 5, nrow = length(samples)*length(genomes)))
names(genomeRA_fff) <- c("Sample","Genome","Rel_ab", "Timepoint","Subject")

# Fill it 
count=0
for (s in 1:length(samples)){
  sdf <- subset(f,f$sample==samples[[s]]) # for each sample 
  sdf$rel_ab_adj = sdf$rel_ab*(1/sum(sdf$rel_ab)) # scale RA 
  
  for (t in 1:length(genomes)){ 
    count = count+1
    genomeRA_fff[count,1] = as.character(samples[[s]])
    genomeRA_fff[count,2]= as.character(genomes[[t]])
    
    if ( (genomes[[t]] %in% sdf$genome)==TRUE){
      rows <- subset(sdf,sdf$genome == genomes[[t]]) 
      genomeRA_fff[count,3] = sum(rows$rel_ab_adj)   
    } else (genomeRA_fff[count,3]=0)
    
    # add timepoint for this genome too
    genomeRA_fff[count,4] = rows$timepoint
    # and subject, since the Sample ID includes timepoint
    genomeRA_fff[count,5] = rows$subject
  }
}

## take average rel_ab across two baseline timepoints 
genomeRA_fff_2 <- genomeRA_fff %>% dplyr::group_by(Subject,Genome)  %>% dplyr::summarise(Rel_ab_mean = mean(Rel_ab))

## now, make sure genome names are in the same format for both datasets:
  # IM has .fna and .fa; FFF does not
  # so add the suffixes to FFF

genomeRA_fff_2$fixed_genome_name = NA
genomeRA_fff_2$Genome <- as.character(genomeRA_fff_2$Genome)
for (i in 1:nrow(genomeRA_fff_2)){
  oldname = genomeRA_fff_2[i,2]
  if (substr(oldname, 1, 3)=="GUT"){
    newname=paste(oldname,'.fna.fa',sep='')
  } else if (substr(oldname, 1, 4)=="META"){
    newname=paste(oldname,'.fa',sep='')
    } else {print(oldname)} 
  genomeRA_fff_2[i,4] = newname
}
# cleanup
genomeRA_fff_2$Genome = genomeRA_fff_2$fixed_genome_name
genomeRA_fff_2 <- select(genomeRA_fff_2, -c(fixed_genome_name))

# save for now
saveRDS(genomeRA_fff_2, file = "genomeRA_fff_2.rds")
saveRDS(genomeRA, file = "genomeRA_IM.rds")

# compare 
which(unique(genomeRA_fff_2$Genome) %in% genomeRA$Genome)  # only 20 genomes shared b/w datasets

## combine both datasets 
names(genomeRA_fff_2) = names(genomeRA)
genomeRA_fff_2$Sample = as.character(genomeRA_fff_2$Sample)
d = rbind(genomeRA_fff_2,genomeRA)


#########################
## convert into a dataframe where row = sample, column = genomes, and cells == RA 
## this will add 0s for genomes not observed in one study that were observed in the other

# Make the empty df 
samples <- unique(d$Sample); taxa <- unique(d$Genome)
df <- data.frame(matrix(ncol = length(taxa), nrow = length(samples))); names(df) <- taxa; rownames(df) <- samples
####### Fill it 
for (s in 1:length(samples)){
  # get genome abundance data for this sample only
  sdf <- subset(d,d$Sample==samples[[s]])
  sdf$rel_ab_adj = sdf$Rel_ab*(1/sum(sdf$Rel_ab))   # scale RA so it sums to 1 in this sample -- redundant (but harmless) since we already did this 
  for (t in 1:length(taxa)){
    if ( (taxa[[t]] %in% sdf$Genome)==TRUE){
      row <- subset(sdf,sdf$Genome == taxa[[t]])      
      df[s,t]= row$rel_ab_adj 
    } else (df[s,t]=0) # ## add rel_ab = 0 for genomes not observed in one study that were observed in the other
  }
}

# save for now
saveRDS(df, file = "df_FFF_IM.rds")




########################################################################
## STEP 2: MAKE CAZY PROFILE FOR EACH SAMPLE
########################################################################

#       samp1	samp2	...	sampN
# GH1
# GH2
# ..
# GHN 



# using the genome rel_ab per sample from this script (saveRDS(df, file = "df_FFF_IM.rds"))
path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)
df = readRDS(file = "df_FFF_IM.rds")

# and the CAZys per genome per MB from this scrip (write.table(cazy3,file="CAZymes_by_genome_PERMB_103025.csv",sep=",",col.names=TRUE,row.names=TRUE))
cazy = read.table("CAZymes_by_genome_PERMB_103025.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE,numerals="allow.loss") 
cazy = cazy[,!(colSums(cazy) == 0) ] # remove 13 cazy families with 0 abundance 


# make the empty df
profile <- data.frame(matrix(ncol = dim(df)[[1]], nrow = dim(cazy)[[2]]))
names(profile) <- rownames(df)
row.names(profile) <- names(cazy)

## FILL IT 
  # For IM: takes ~17 seconds for 1 sample, so all 77 samples should take 28 min 
  # For FeFiFo: takes much longer because sequencing is deeper/more genomes are detected
  # as of 103025 it's a couple hours total

for (samp in 3:nrow(df)){  ### ! ! ! ! ! 
  print(samp)
    
  # GET THE GENOMES PRESENT IN THIS SAMPLE 
  samprow <- df[samp,] # genome relative abundance per sample -- a subset of df
  samp_ra_0<-t(samprow)
  samp_ra <- subset(samp_ra_0,samp_ra_0[,1] > 0) # every row is a genome; 1 column for the sample 
  # for samp = 1, there are 38 genomes 
  
  # GET THE CAZY FAMILIES PRESENT IN THIS SAMPLE 
  samp_cazy0 <- subset(cazy,row.names(cazy) %in% row.names(samp_ra)) # subset cazy to genomes present in sample
  samp_cazy <- samp_cazy0[, names(samp_cazy0) %in% c(names(which(colSums(samp_cazy0) > 0)))] # subset again to families present in these genomes
  
  # FOR EACH CAZY FAMILY [IN THIS SAMPLE]  
  for (cazfam in 1:ncol(samp_cazy)){ # cazy is the genomes x genes per family table 
    
    # START TAB -- this is specific to CAZy family, and will be summed over genomes
    cazytab = 0
    
    # FOR EACH GENOME IN THE SAMPLE:
    for (g in 1:nrow(samp_ra)){ # g is a genome index 
      # MULTIPLY GENOME RA IN SAMPLE * NUMBER OF FAMILY MEMBERS <-- THIS IS VALUE IN CELLS OF OUTPUT TABLE 
      # first get number of CAZyme family members in this bug
      genome = row.names(samp_ra)[[g]] # g is an index, genome is genome name
      
      bugsrow <- samp_cazy[genome,] # one row for this genome; columns are CAZy families 
      num_members = bugsrow[,cazfam] # this is the number of genes in the CAZy family in this bug 
      
      # num_members is NA if the genome has no CAZymes in this family -- change to 0
      if (is.na(num_members) == T){
        num_members = 0  }
      
      ## then find RA of genome in sample 
      ra <- as.numeric(subset(samp_ra,row.names(samp_ra)==genome))
      ra_times_genenum = ra*num_members
      
      # ADD TO TAB
      cazytab = cazytab+ra_times_genenum
    } # end genome loop 
    
    # SAVE TAB FOR THIS CAZY FAMILY/SAMPLE 
    # need to determine appropriate column in 'profile', which has more cazy families than 'samp_cazy'
    cazytab_familyname = names(samp_cazy)[cazfam]
    cazytab_familyname_index <- row.names(profile) %>% equals(cazytab_familyname) %>% which
    profile[cazytab_familyname_index,samp] = cazytab # samp is the column
  } # end cazy loop
} # end sample loop


# the NAs in 'profile' are for CAZY families that aren't in a given sample -- make them 0
profile2 <- profile %>% replace(is.na(.), 0)


hist(colSums(profile2)) # mean CAZy density across all samples is 40 genes / Mb, normally distributed 

hist(rowSums(profile2)/ncol(profile2)) # a few CAZy families are very abundant--  >4 copies per Mb 
                                                                                  # (Mean of all samples, both studies) 
## Save! 
#write.table(profile2,file="CAZymes_by_sample_IMplusFFF_103025.csv",sep=",",col.names=TRUE,row.names=TRUE)






#### #### #### #### #### #### #### #### #### #### 
####  PLOT TOTAL CAZY DENSITY BY POPULATION 
####            FIGURE 7C
#### #### #### #### #### #### #### #### #### #### 


profile <-read.table("CAZymes_by_sample_IMplusFFF_103025.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 
fffs <- names(profile)[1:37]
ims <- names(profile)[38:112]

# outliers 
outliers = c("AK_SK_49","AK_SK_10","AK_SK_32","AK_SR_6","AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2") #AK_SK_24 and AK_SK_19 have no/little data , AK_SK_34 is an oral sample; AK_SR_4 has 60+% Streptococcus
others = c('AK_SK_24',"AK_SK_19","AK_SK_34","AK_SR_4")
ims <- subset(ims,ims %in% c(outliers,others)==F)
#
prof_fff <- profile[, fffs]
prof_im <- profile[, ims]
prof_outliers<-profile[,outliers]
#
im_density = colSums(prof_im)
fff_density = colSums(prof_fff)
outlier_density= colSums(prof_outliers)


#########
# Divide IM densities into tribes/populations:
# Ultimately aiming for:    Sample, Density, Region, Tribe 
#########
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- subset(meta,meta$sample %in% names(profile))
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau")); meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond")) 
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui")); meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas")); meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast")); meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))

#
im_d <- as.data.frame(im_density); names(im_d) = "Density"
im_d$Sample = row.names(im_d)
im_d$Tribe = NA; im_d$Region = NA 
for (row in 1:nrow(im_d)){
  samplename = im_d[row,2]
  tribe <- as.character(subset(meta,meta$sample==samplename)$Tribe)
  region <- as.character(subset(meta,meta$sample==samplename)$Region)
  im_d[row,3] = tribe
  im_d[row,4] = region
}
#
out_d <- as.data.frame(outlier_density); names(out_d) = "Density"
out_d$Sample = row.names(out_d)
out_d$Tribe = 'Tribal outlier'; out_d$Region = 'Tribal outlier' 
#
fff_d <- as.data.frame(fff_density); names(fff_d) = "Density"
fff_d$Sample = row.names(fff_d)
fff_d$Tribe = 'California'; fff_d$Region = 'California' 
#
### Combine
plot_df = rbind(im_d,fff_d,out_d)
plot_df$Tribe = factor(plot_df$Tribe, levels = c("Warli","Gond","Madia","Kabui","Boto","Brokpa","Balti","Purigpa","Tribal outlier","California"))
plot_df$Region = factor(plot_df$Region,levels =c("Coast","Deccan Plateau","Northeast Hills","Trans-Himalayas","California","Tribal outlier"))

col_list<-c("#b33f25", "#2f3cb5",   "#edc42f","#2a8022","gray","purple")
density = ggplot(plot_df, aes(x=Region, y=Density)) + 
  geom_violin(aes(fill=Region),draw_quantiles = 0.5) +  ylab("CAZymes per Mb") + xlab("") +
  theme_cowplot() +
  scale_fill_manual(values = col_list)  + ylim(20,65) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(legend.position = "none") 
density
ggsave("Fig7A.pdf", plot=density,width=1.4*2.1,height=1.9*2.25,units="in")



### STATS 
out = subset(plot_df,plot_df$Region=="Tribal outlier") # tribal outliers 
im_all = subset(plot_df,plot_df$Region != "California"); im_nonOUT = subset(im_all,im_all$Region != "Tribal outlier")  # tribal NOT outlier 
CA = subset(plot_df,plot_df$Region == "California") # CA
#
t.test(im_nonOUT$Density,out$Density) #p-value = 0.0394 --> TRIBES VS OUTLIER
t.test(im_nonOUT$Density,CA$Density) #p-value = 1.29e-07 --> TRIBES VS CALIFORNIA
#
### STATS - variation among regions (not outlier, not CA)
data <- subset(plot_df,plot_df$Region!="Tribal outlier")
data <- subset(data,data$Region!="California")
res.aov <- aov(Density ~ Region, data = data) # one-way ANOVA was performed 
summary(res.aov) # p =  0.654





#### #### #### #### #### #### #### #### #### #### 
####  PLOT CAZY PCA
####            FIGURE 7A
#### #### #### #### #### #### #### #### #### #### 


## Read in CAZy density per Mb per sample
path = "/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/Revision_October_2025/Github/"
setwd(path)
cazy <-read.table("CAZymes_by_sample_IMplusFFF_103025.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 
cazyt <- as.data.frame(t(cazy))
cazy2 <- cazyt

# remove the streptococcus outlier 
cazy2 <- subset(cazy2,row.names(cazy2) != "AK_SR_4")

### organize metadata
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T)  # lowercase tribe is from Matt's v2; uppercase is from Dattatray/Abhijit
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau")) # replace regions with better names 
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))
# Subset metadata to match cazy IM data
meta <- subset(meta,meta$sample %in% row.names(cazy2))

### ORDER DF AND META BY SAMPLE so the metadata corresponds appropriately
cdf <- cazy2[ order(row.names(cazy2)), ]
meta <- meta[order(meta$sample),]
## variables for California
row.names(cdf)
regions <- c(as.character(meta$Region),rep("California",37))
tribes <- c(as.character(meta$Tribe),rep("California",37))

# notate the Bacteroides outliers
outliers = c("AK_SK_49","AK_SK_10","AK_SK_32","AK_SR_6","AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2") 
outlier_vector = rep(NA,111) # same length as N samples in cdf -- will add to metadata later
count=0
for (i in row.names(cdf)){
  count=count+1
  if (i %in% outliers){
    outlier_vector[count] = 'yes'
  } else {outlier_vector[count] = 'no'}
} 


#######################################################
####  PERFORM PCA
#######################################################

# MAKE SURE THESE FILES HAVE SAMPLES IN THE SAME ORDER (should be done above)
#cdf <- cazy2[ order(row.names(cazy2)), ] 
#meta <- meta[order(meta$sample),] 

bc_dist = vegan::vegdist(cdf, method = "bray") # calculate distance 
dfc = as.data.frame(scale(bc_dist,center=T))  # center
as.matrix(dfc)[1:5, 1:5]
#
pca_bcdist <- dudi.pca(dfc,nf=5,scale=T,center=T,scannf=F)
#
res.ind <- get_pca_ind(pca_bcdist)
coords <- res.ind$coord # Get the coordinate (x-y) values for indivdiuals
eig.val <- get_eigenvalue(pca_bcdist)
head(eig.val) # 37.2% and 31.6% 
dimsToplot <- as.data.frame(cbind(coords$Dim.1,coords$Dim.2,regions,tribes))
names(dimsToplot) <- c("PC1","PC2","Region","Tribe")
dimsToplot$PC1 <- as.numeric(dimsToplot$PC1); dimsToplot$PC2 <- as.numeric(dimsToplot$PC2)
row.names(dimsToplot) = row.names(cdf)
dimsToplot$Region = factor(dimsToplot$Region,levels=c("Coast","Deccan Plateau","Northeast Hills","Trans-Himalayas","California"))
dimsToplot$outlier = outlier_vector



####
#### PLOT WITHOUT HULL INITIALLY 
####

col_list<-c("#b33f25", "#2f3cb5",   "#edc42f","#2a8022","gray")
ggplot(data = dimsToplot, aes(x = PC1, y = PC2, colour=Region, fill = Region)) +
  geom_point(aes(shape=outlier)) + geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)  + 
  theme_cowplot() +
  xlab("PC1 (37.2%)") + ylab("PC2 (31.6%)")
#

# to compare Indians without Californians
temp <- subset(dimsToplot,dimsToplot$Region!="California")
ind <- subset(cdf,row.names(cdf) %in% row.names(temp))
bc_dist_ind = vegan::vegdist(ind, method = "bray") # calculate distance 
adonis2(bc_dist_ind ~ temp$Region, by='terms', perm=10001) # p = 0.0018  , R2 = 0.104
adonis2(bc_dist_ind ~ temp$Tribe, by='terms', perm=10001) # p = 0.007798   , R2 = 0.158

# to compare CA vs india:
country = c(rep("India",74),rep("USA",37))
adonis2(bc_dist ~ country, by='terms', perm=10001) # p = 9.998e-05 , R2 = 0.192




####
#### PLOT WITH HULL (THAT EXCLUDES OUTLIERS) 
####

dimsToplot$sample = row.names(dimsToplot)

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
dimsToplot_NOOUTLIERS = subset(dimsToplot,row.names(dimsToplot) %in% outliers == FALSE)
hulls_metag_region <- ddply(dimsToplot_NOOUTLIERS, "Region", find_hull)
#
fig7a = ggplot(data = dimsToplot, aes(x = -PC1, y = PC2, colour=Region, fill = Region)) + 
  geom_point(aes(shape=outlier))  +
  geom_hline(yintercept = 0, lty = 2) +geom_vline(xintercept = 0, lty = 2) +
  scale_color_manual(values = col_list) + 
  scale_fill_manual(values = col_list)  + 
  theme_cowplot()  +  theme(legend.position = "right") +
  ylab("PC1 (37.2%)") + xlab("PC2 (31.6%)") + 
  geom_polygon(data = hulls_metag_region, alpha = 0.25)  +
  theme(legend.text=element_text(size=10),legend.title=element_text(size=11))  +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)))
fig7a

ggsave("Fig7A.pdf", plot=fig7a,height=2.05*2,width=2.1*2,units="in")






#### #### #### #### #### #### #### #### #### #### 
#### COMPARE CAZYS ACROSS POPULATIONS
####            FIGURE 7D
# banana
#### #### #### #### #### #### #### #### #### #### 


profile <-read.table("CAZymes_by_sample_IMplusFFF_103025.csv", header = TRUE, sep = ",",stringsAsFactors = TRUE,numerals="allow.loss") 

## Remove Bacteroides outliers + Strep outlier for this analysis 
outliers = c("AK_SK_49","AK_SK_10","AK_SK_32","AK_SR_6","AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2", "AK_SR_4")
profile <- select(profile, -all_of(outliers))

#### metadata 
meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- subset(meta,meta$sample %in% names(profile))
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau")); meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui")); meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas")); meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast")); meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))

##### DIVIDE INDIVIDUALS INTO STUDY 
FFFsamps = c("X8000","X8001","X8002","X8003","X8004","X8006","X8007", "X8008","X8009","X8010","X8011","X8012","X8013","X8014",     
             "X8016","X8017","X8018","X8020","X8021","X8022","X8023","X8024","X8025","X8026","X8027","X8028","X8029","X8030","X8032","X8033","X8034","X8035","X8036","X8037","X8038","X8039","X8041")
prof_fff <- select(profile, all_of(FFFsamps))
prof_IM <- select(profile, -all_of(FFFsamps))



##############################################################
### CALCULATE MEAN DENSITY OF EACH FAMILY IN EACH POPULATION
### AND PERFORM WILCOXIN (Mann-Whitney U) TEST
##############################################################

### FIND AVERAGE CAZYME VALUE PER GROUP #AND COMBINE
fff <- rowSums(prof_fff)/ncol(prof_fff)
im <- rowSums(prof_IM)/ncol(prof_IM) 
prof = as.data.frame(cbind(fff,im))
names(prof) = c('fff_mean','im_mean')

# Mann-Whitney U is non-parametric:  
# for each CAZY family and population, there is a distribution of densities. Does the distribution vary by population?

prof$wilcox.pval = NA
for (row in 1:dim(prof)[[1]]){
  # for each family
  fam = row.names(prof)[row]
  
  # get the densities
  imd = as.numeric(prof_IM[row,])
  fffd = as.numeric(prof_fff[row,])
  
  # assemble into df with indicator variable for population 
  imdf = as.data.frame(cbind(imd,rep(0,length(imd))))
  names(imdf)=c("fam.density",'population')
  #
  fffdf = as.data.frame(cbind(fffd,rep(1,length(fffd))))
  names(fffdf)=c("fam.density",'population')
  # 
  df = rbind(imdf,fffdf)
  
  # run the test and add the p-value to the 'prof' data frame
  t = wilcox.test(fam.density ~ population, data=df) 
  wilcox.pval = t[[3]]
  prof$wilcox.pval[row] = wilcox.pval
  
}

# remove fams that are 0 in both
prof <- subset(prof,is.na(prof$wilcox.pval)==F)

# subset to families with significant adj p-vals
prof$wilcox.pval_adj = p.adjust(prof$wilcox.pval, method="BH")
profsig <- subset(prof,prof$wilcox.pval_adj<=0.05)

# 300/406 =  74% of CAZy families differ significantly in density b/w CA and Indian tribes


##################################################################
### CALCULATE THE RATIO (FFF / IM ) 
### AND FILTER LOW-MEAN FAMILIES, WHICH COULD HAVE INFLATED RATIOS 
##################################################################

profsig$ratio = profsig$im/profsig$fff

## Use a threshold filter of mean 0.1 instances/Mb 
profsig$max = pmax(profsig$fff_mean,profsig$im_mean)
profsigfilt <- subset(profsig,profsig$max >=0.1)
# now down to 77 non-rare families 


### SEPARATE INTO ENRICHMENT GROUP TO ADJUST INFINITY SIGN AND X-INDEX
cali_enrichment <- subset(profsigfilt,profsigfilt$ratio<=1)
indian_enrichment <- subset(profsigfilt,profsigfilt$ratio>1)
# take the inverse of the Californian ratio so it also goes from 1 to +Inf
cali_enrichment$ratioi = cali_enrichment$fff/cali_enrichment$im 
cali_enrichment <- cali_enrichment[ , -which(names(cali_enrichment) %in% c('ratio'))]
names(cali_enrichment)[6] = "ratio"
# 
hist(indian_enrichment$ratio)
hist(cali_enrichment$ratio)


## substitute infinities, take log, then turn Cali negative for now (for indexing)
cali_enrichment$ratio <- as.numeric( gsub("Inf", "10000", cali_enrichment$ratio))
cali_enrichment$log.ratio <- log(cali_enrichment$ratio)
indian_enrichment$ratio <- as.numeric( gsub("Inf", "10000", indian_enrichment$ratio))
indian_enrichment$log.ratio <- log(indian_enrichment$ratio)
cali_enrichment$log.ratio <- cali_enrichment$log.ratio  * -1 

## add column for cazy families for plotting
indian_enrichment$family = row.names(indian_enrichment)
cali_enrichment$family = row.names(cali_enrichment)

## combine and sort
all_enrichment = rbind(cali_enrichment,indian_enrichment)
all_enrichment <- all_enrichment[order(all_enrichment$log.ratio),]


## turn into volcano plot
# 1. give each family an x-index based on its y-value, negative and positive centered at y=0 
all_enrichment$x_index = seq(1,dim(all_enrichment)[[1]],by=1)
# change this so India is on left
all_enrichment$x_index =77-all_enrichment$x_index 

## 2. for negative x-values, turn their y-values positive
all_enrichment$log.ratio = abs(all_enrichment$log.ratio)


## PLOT 
all_enrichment$family=factor(all_enrichment$family,levels=all_enrichment$family)

p =ggplot(all_enrichment, aes(x=x_index,y=log.ratio, fill=family)) + 
  geom_point(show.legend = FALSE, size = 0.5) + theme_cowplot(font_size=12) +
  geom_hline(yintercept=log(1), color = "black") +
  geom_vline(xintercept=36.5, linetype="dashed", color = "black") +# where enrichment switches b/w pops
  theme(axis.text.x=element_blank()) +
  ylab("log(ratio)") + xlab("") + 
  scale_x_discrete(expand = expansion(add = 3))
p  

### Save and label/annotate in pptx
#ggsave("Fig7D_left.pdf",   plot=p,width=1.55*1.75,height=1.42*1.75,units="in")



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
### Same thing for TH versus other tribes
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

# prof_IM contains subject-level data for the Indian tribes 
# separate TH from others 
th_samps <- subset(meta,meta$Region == "Trans-Himalayas")
prof_th <- select(profile, th_samps$sample)
prof_other <- select(profile, -th_samps$sample)
prof_other <- select(prof_other, -FFFsamps)

### FIND AVERAGE CAZYME VALUE PER GROUP AND COMBINE
th <- rowSums(prof_th)/ncol(prof_th)
non <- rowSums(prof_other)/ncol(prof_other) 
prof = as.data.frame(cbind(th,non))
names(prof) = c('TH','others')


### PERFORM MANN-WHITNEY U TEST 

prof$wilcox.pval = NA
for (row in 1:dim(prof)[[1]]){
  # for each family
  
  # get the densities
  thd = as.numeric(prof_th[row,])
  otherd = as.numeric(prof_other[row,])
  
  # assemble into df with indicator variable for population 
  thdf = as.data.frame(cbind(thd,rep(0,length(thd))))
  names(thdf)=c("fam.density",'population')
  #
  otherdf = as.data.frame(cbind(otherd,rep(1,length(otherd))))
  names(otherdf)=c("fam.density",'population')
  # 
  df = rbind(thdf,otherdf)
  
  # run the test and add the p-value to the 'prof' data frame
  t = wilcox.test(fam.density ~ population, data=df) 
  wilcox.pval = t[[3]]
  prof$wilcox.pval[row] = wilcox.pval
  
}

# remove fams that are 0 in both
prof <- subset(prof,is.na(prof$wilcox.pval)==F)


# Remove non-significant families
prof$wilcox.pval_adj = p.adjust(prof$wilcox.pval, method="BH")
profsig <- subset(prof,prof$wilcox.pval_adj<=0.05)
# here, only 51/362 fams (14%) are significantly different -- another indicator of big diff b/w IM and Cali


## calculate the ratio 
profsig$ratio = profsig$TH/profsig$others


## filter out families at less than threshold density 
profsig$max = pmax(profsig$TH,profsig$others)
profsigfilt <- subset(profsig,profsig$max >=0.1) 
# only 17 families left 




### SEPARATE INTO ENRICHMENT GROUP TO ADJUST INFINITY SIGN AND X-INDEX
TH_enrichment <- subset(profsigfilt,profsigfilt$ratio>1)
others_enrichment <- subset(profsigfilt,profsigfilt$ratio<=1)
# take the inverse of the 'others' ratio so it also goes from 1 to +Inf (instead of 0 to 1)
others_enrichment$ratioi = others_enrichment$others/others_enrichment$TH 
others_enrichment <- others_enrichment[ , -which(names(others_enrichment) %in% c('ratio'))]
names(others_enrichment)[6] = "ratio"
# 
hist(others_enrichment$ratio) # 1-2
hist(TH_enrichment$ratio) # up to 30X !
#

## substitute infinities, take log, then turn TH negative for now (for indexing)
TH_enrichment$log.ratio <- log(TH_enrichment$ratio)
others_enrichment$log.ratio <- log(others_enrichment$ratio)
TH_enrichment$log.ratio <- TH_enrichment$log.ratio  * -1 


## add column for cazy families for plotting
others_enrichment$family = row.names(others_enrichment)
TH_enrichment$family = row.names(TH_enrichment)

## combine and sort
all_enrichment = rbind(TH_enrichment,others_enrichment)
all_enrichment <- all_enrichment[order(all_enrichment$log.ratio),]


# # Assign xval
all_enrichment$x_index = seq(1,dim(all_enrichment)[[1]],by=1)


## Flip the negative ratios positive
all_enrichment$log.ratio = abs(all_enrichment$log.ratio)


## PLOT 
all_enrichment$family=factor(all_enrichment$family,levels=all_enrichment$family)

p =ggplot(all_enrichment, aes(x=x_index,y=log.ratio, fill=family)) + 
  geom_point(show.legend = FALSE, size = 0.5) + theme_cowplot(font_size=12) +
  geom_hline(yintercept=log(1), color = "black") +
  geom_vline(xintercept=11.5, linetype="dashed", color = "black") +# where enrichment switches b/w pops
  theme(axis.text.x=element_blank()) +
  ylab("log(ratio)") + xlab("") + 
  scale_x_discrete(expand = expansion(add = 3))
p  


### Save and label/annotate in Figure 7.pptx
#ggsave("Fig7D_right.pdf", plot=p,width=1.55*1.75,height=1.42*1.75,units="in")