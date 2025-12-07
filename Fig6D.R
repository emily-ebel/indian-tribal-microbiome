library(stringr)
library(data.table)
library(ggplot2)
library(cowplot)

############################################################################
## This script reads in dbCAN output (overview.txt files) for multiple genomes
##  and makes a table with the # of genes per CAZy family per genome,
##   which is used for the enrichment tests in Figure 6D
############################################################################

## List of dbCAN output files for 133 B. adolescentis genomes
files4 <- list.files(path="/Users/emily/Documents/SonnenburgPostdoc/IndianMicrobiome/dbcan/032525", pattern="*overview.txt", full.names=TRUE, recursive=FALSE) # Bif adol


## we want to to make this table: 
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

# --> 81 CAZy fams in these 133 B. adolescentis genomes 




###### NOW, WE CAN MAKE AN EMPTY TABLE OF  GENOMES BY CAZYMES ### 
# where each cell is the number of cazyme genes for that family/genome

df <- data.frame(matrix(ncol = length(allcazy2), nrow = length(genomes)))
colnames(df) <- allcazy2; rownames(df) <- genomes

## Go back through the files (genomes) and fill in the table.
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


# remove the directory prefix from the genome names 
oldnames<-row.names(df)
newnames <- oldnames
newnames <- gsub('032525/','',newnames)
row.names(df) <- newnames



# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 
## What CAZy families are enriched in the Indian tribal B. adolescentis? (and Mongolian MAGs)
# @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ 

all_samples = row.names(df) 


######## divide genomes into categories

IM = c("Bifidobacterium-adolescentis_AK_SK_11_Brokpa.fa","Bifidobacterium-adolescentis_AK_SK_1_Brokpa.fa",  
  "Bifidobacterium-adolescentis_AK_SK_12_Brokpa.fa","Bifidobacterium-adolescentis_AK_SK_13_Brokpa.fa",
  "Bifidobacterium-adolescentis_AK_SK_2_Purigpa.fa","Bifidobacterium-adolescentis_AK_SK_32_Balti.fa" ,
  "Bifidobacterium-adolescentis_AK_SK_33_Brokpa.fa","Bifidobacterium-adolescentis_AK_SK_35_Balti.fa", 
  "Bifidobacterium-adolescentis_AK_SK_36_Boto.fa" ,"Bifidobacterium-adolescentis_AK_SK_44_Balti.fa" ,
       "Bifidobacterium-adolescentis_AK_SK_25_Purigpa.fa","Bifidobacterium-adolescentis_AK_SK_37_Brokpa.fa",
  "Bifidobacterium-adolescentis_AK_SK_48_Purigpa.fa","Bifidobacterium-adolescentis_AK_SK_8_Balti.fa",
  "Bifidobacterium-adolescentis_AK_SK_9_Boto.fa" ,
       "Bifidobacterium-adolescentis_AK_SK_45_Purigpa.fa", "Bifidobacterium-adolescentis_AK_SK_5_Balti.fa",
       "Bifidobacterium-adolescentis_RS12_A_Warli.fa","Bifidobacterium-adolescentis_RS5_A_Warli.fa",
  "Bifidobacterium-adolescentis_RS12_B_Warli.fa")
  
Mongolia = row.names(df)[71:90]
Hadza = row.names(df)[94:101]
Nepal = row.names(df)[102:112]

industrialized <- all_samples[ (all_samples %in% IM == FALSE)]
industrialized <- industrialized[ (industrialized %in% Hadza == FALSE)]
industrialized <- industrialized[ (industrialized %in% Nepal == FALSE)]
industrialized <- industrialized[ (industrialized %in% Mongolia == FALSE)] # leaves 74 'industrialized'
# for El Salvadoran samples, we don't know if they're industrialized or not, so remove them
# and also remove ATCC
drop = c(industrialized[37:46],"ATCC-15703_dwnld031625_sbmtrJapan_GCA_000010425.fna")
industrialized <- industrialized[ (industrialized %in% drop == FALSE)] # leaves 63 'industrialized'


######## divide cazy df into categories
indust_df <- subset(df,row.names(df)%in%industrialized)  
IM_df <- subset(df,row.names(df)%in%IM)
Mongolia_df=subset(df,row.names(df)%in%Mongolia)




# ### Do isolates have more CAZys than MAGs? --> yes
# notmags = c(all_samples[1:18],all_samples[38],all_samples[49:55],all_samples[66:69],all_samples[132])
# notmags_df <- subset(df2,row.names(df2)%in%notmags)
# notmags_totals = as.data.frame(rowSums(notmags_df)); names(notmags_totals)='num_cazy_genes'
# b= seq(from=40,to=120,by=5)
# hist(notmags_totals$num_cazy_genes,xlim=c(40,120),breaks=b)
# #
# allmags = all_samples[!all_samples %in% notmags]
# allmags_df <- subset(df2,row.names(df2)%in%allmags)
# allmags_totals = as.data.frame(rowSums(allmags_df)); names(allmags_totals)='num_cazy_genes'
# hist(allmags_totals$num_cazy_genes,xlim=c(40,120),col="lightpink",breaks=b)
# hist(notmags_totals$num_cazy_genes,xlim=c(40,120),breaks=b,add=T)
# t.test(allmags_totals$num_cazy_genes,notmags_totals$num_cazy_genes) # 76 vs 99,  p-value < 2.2e-16




##################################################################################
### which CAZY fams are enrichmd in Indian tribal B. adolescentis vs. industrialized? 
##################################################################################

# save the pvals
pval_ = c()

count = 0
countsig=0

for (c in 1:length(names(df))){ # for each cazyme family
  non_IM_count = indust_df[,c]
  IM_count = IM_df[,c]
  
  if (mean(non_IM_count) != mean(IM_count)) {
      count = count+1
      x=t.test(non_IM_count,IM_count,alternative='less')
      pval_ = c(pval_,x[[3]])

      if (x[[3]]<=0.0006){ # # Bonferroni here is on 81 CAZy fams, so 0.0006
        countsig = countsig+1
        print(names(df)[c]) # print the significant tests
        print(x)
        print('')
      } } else { pval_ = c(pval_,1)}
}


# work with pval
p_IM = as.data.frame(cbind(pval_,names(df)))
p_IM$adj = p.adjust(pval_, method="BH")

# for IM vs 63 industrialized:

# CBM23 p = 0.0002555
# 0.2539683 (indust) vs 0.8500000 (IM)
# adj_p = 0.0103

# GH95 p 0.0002308
# 0.06349206 vs 0.55000000 
# adj_p = 0.0103



##################################################################################
### which CAZY fams are enriched in MONGOLIAN B. adolescentis vs. industrialized? 
##################################################################################

# save the pvals
pval_ = c()

count = 0

for (c in 1:length(names(df))){ # for each cazyme family
  non_IM_count = indust_df[,c]
  Mongolia_count = Mongolia_df[,c]
  
  if (mean(non_IM_count) != mean(Mongolia_count)) {
    count = count+1
    x=t.test(non_IM_count,Mongolia_count,alternative='less')
    pval_ = c(pval_,x[[3]])
    
    if (x[[3]]<=0.0006){ # # Bonferroni here is on 81 CAZy fams, so 0.0006
      countsig = countsig+1
      print(names(df)[c]) # print the significant tests
      print(x)
      print('')
    } } else { pval_ = c(pval_,1)}
}

# work with pval
p_Mongol = as.data.frame(cbind(pval_,names(df)))
p_Mongol$adj = p.adjust(pval_, method="BH")

# for Mongolia vs 63 industrialized:

# CBM23 p = 0.0004308
#0.2539683 (indust) vs 1.1500000 (Mongolia)
# adj_p = 0.0103

# CE9 p = 0.00164373893040849
# adj_p = 0.04438095

# GH42 p = 	0.00158342450015461
# adj_p = 0.04438095



#########################################################
####  Make a table of copynumbers for  barplot
#########################################################

## GH95

## Get copy nubmers by group 
length(which(IM_df$GH95 ==0))
length(which(IM_df$GH95 ==1))
length(IM_df$GH95)
#
length(which(Mongol_df$GH95 ==0))
length(which(Mongol_df$GH95 ==1))
length(Mongol_df$GH95)
#
length(which(indust_df$GH95 ==0))
length(which(indust_df$GH95 ==1))
length(indust_df$GH95)


### Now put them in a table 
population <- c(rep("India" , 4) , rep("Mongolia" , 4) , rep("Industrialized" , 4) )
copynum <- rep(c("0", "1","2","3+") , 3)
gene = rep("GH95",12)
value <- c(9/20,11/20,0,0,     18/20,2/20,0,0,     59/63,4/63,0,0) 
gh95 <- data.frame(population,copynum,gene,value)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## CBM23

length(which(IM_df$CBM23 ==0))
length(which(IM_df$CBM23 ==1))
length(which(IM_df$CBM23 ==2))
length(IM_df$CBM23)
#
length(which(Mongol_df$CBM23 ==0))
length(which(Mongol_df$CBM23 ==1))
length(which(Mongol_df$CBM23 ==2))
length(Mongol_df$CBM23)
#
length(which(indust_df$CBM23 ==0))
length(which(indust_df$CBM23 ==1))
length(which(indust_df$CBM23 ==2))
length(indust_df$CBM23)


population <- c(rep("India" , 4) , rep("Mongolia" , 4) , rep("Industrialized" , 4) )
copynum <- rep(c("0", "1","2","3+") , 3)
gene = rep("CBM23",12)
value <- c(5/20,13/20,2/20,0,  8/20,1/20,11/20,0,  55/63,0/63,8/63,0) ### automate this banana 
CBM23 <- data.frame(population,copynum,gene,value)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

 ## GH42
length(which(IM_df$GH42 ==0))
length(which(IM_df$GH42 ==1))
length(which(IM_df$GH42 ==2))
length(which(IM_df$GH42 >=3))
length(IM_df$GH42)
#
length(which(Mongol_df$GH42 ==0))
length(which(Mongol_df$GH42 ==1))
length(which(Mongol_df$GH42 ==2))
length(which(Mongol_df$GH42 >=3))
length(Mongol_df$GH42)
#
#
length(which(indust_df$GH42 ==0))
length(which(indust_df$GH42 ==1))
length(which(indust_df$GH42 ==2))
length(which(indust_df$GH42 >=3))
length(indust_df$GH42)
#
#
population <- c(rep("India" , 4) , rep("Mongolia" , 4) , rep("Industrialized" , 4) )
copynum <- rep(c("0", "1","2","3+" ) , 3)
gene = rep("GH42",12)
value <- c(1/20,5/20,7/20,7/20,    0/20,0/20,0/20,20/20,  0,0/63,6/63,57/63)  
GH42 <- data.frame(population,copynum,gene,value)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## CE9

length(which(IM_df$CE9 ==0))
length(which(IM_df$CE9 ==1))
length(which(IM_df$CE9 ==2))
length(IM_df$CE9)
#
length(which(Mongol_df$CE9 ==0))
length(which(Mongol_df$CE9 ==1))
length(which(Mongol_df$CE9 ==2))
length(Mongol_df$CE9)
#
length(which(indust_df$CE9 ==0))
length(which(indust_df$CE9 ==1))
length(which(indust_df$CE9 ==2))
length(indust_df$CE9)

population <- c(rep("India" , 4) , rep("Mongolia" , 4) , rep("Industrialized" , 4) )
copynum <- rep(c("0", "1","2","3+") , 3)
gene = rep("CE9",12)
value <- c(19/20,1/20,0/20,0,  5/20,15/20,0/20,0,  39/63,24/63,0,0) ### automate this banana
CE9 <- data.frame(population,copynum,gene,value)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


### Combine and plot 
data <-rbind(gh95,CBM23,GH42,CE9)
data$population=factor(data$population,levels=c("India","Mongolia","Industrialized"))
data$gene = factor(data$gene,levels=c("GH95","CBM23","GH42","CE9"))
#
#
col_list=c('forestgreen','goldenrod' ,'gray')
f <- ggplot(data, aes(x=copynum, y=value, fill=population)) + 
  geom_bar(position="dodge", stat="identity",width=0.7) +theme_cowplot()+
  theme(legend.title=element_blank()) +
  ylab("Fraction B. adolescentis genomes")+xlab("Domain copy number")+  scale_fill_manual(values = col_list) 
g = f + facet_wrap(~gene,scales="free",nrow=2) 
g

ggsave("Fig6D.pdf", plot=g,width=2.38*2.25,height=2.15*2.25,units="in")
