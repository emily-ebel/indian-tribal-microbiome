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
##  Maaslin with covariates (sex, age, BMI) with 10-25% missing data
##  
##
#####################################################################

df_IM_data = read.table(file             = "Fig4_relab_genus16S_Maaslin.tsv", # Fig2_relab_genus16S.csv formatted for Maaslin 
                        header           = TRUE,
                        sep              = "\t", 
                        row.names        = 1,
                        stringsAsFactors = FALSE) # 68 samples, 166 genera, same as Fig2_relab_genus16S.csv


# Exclude outliers
outliers = c("AK_SK_49","AK_SK_10",'AK_SK_32','AK_SR_6',"AK_SG_17","AK_SG_18","AK_SR_1","AK_SR_2","AK_SR_4") # 8 high-Bacteroides + Streptococcus 
df_IM_data <- subset(df_IM_data,(row.names(df_IM_data)%in%outliers==FALSE)) # leaves 68 non-outliers 


########################################
## Format metadata for Maaslin

meta <-read.table("sample_metadata.csv", header = TRUE, sep = ",",stringsAsFactors = T) 
meta <- subset(meta,meta$sample %in% row.names(df_IM_data))
meta <- subset(meta,meta$sample != "AK_SK_31") # this sample has no metadata
#
meta <- meta %>% mutate(Region = str_replace(Region, "Central", "Deccan Plateau"))
meta <- meta %>% mutate(Region = str_replace(Region, "North-East", "Kabui"))
meta <- meta %>% mutate(Region = str_replace(Region, "North", "Trans-Himalayas"))
meta <- meta %>% mutate(Region = str_replace(Region, "West", "Coast"))
meta <- meta %>% mutate(Region = str_replace(Region, "Kabui", "Northeast Hills"))
meta <- meta %>% mutate(Tribe = str_replace(Tribe, "Gondia", "Gond"))

# Now there are 67 samples 

# Reduce to: Sample, Tribe, Region, Age, Sex, BMI
meta_maaslin <- meta[, c("sample", "Tribe", "Region","sex","Age","BMI")]
meta_maaslin <- meta_maaslin %>% mutate(sex = str_replace(sex, "female", "F"))
meta_maaslin <- meta_maaslin %>% mutate(sex = str_replace(sex, "male", "M"))
meta_maaslin$sex[meta_maaslin$sex == "missing"] <- NA 
row.names(meta_maaslin) = meta_maaslin$sample
meta_maaslin[1:5, ]


###################################################
## Start by looking at each covariate by region
###################################################


####### SEX / GENDER ###########
#################################

counts <- meta_maaslin %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_NA     = sum(is.na(sex)),
    n_F   = sum(sex == "F", na.rm = TRUE),
    n_M   = sum(sex == "M", na.rm = TRUE)
  )
counts2=counts
counts2$missing = counts2$n_NA / (counts2$n_NA + counts2$n_M + counts2$n_F) # Around 10% of subjects have NA for sex
counts2$male = counts2$n_M / (counts2$n_M + counts2$n_F) # 
counts2

####### --> Trans-Himalayan subjects have  more missing gender data / fewer males ####### 

# Test significance with ANOVA by region, melting the data first
counts_long  <- counts %>% 
  pivot_longer(
    cols = n_NA:n_M,
    names_to = "Type",
    values_to = "Count"
  )
fit <- aov(Count ~ Region + Type, data = counts_long)
summary(fit) # p = 0.01 (region) and 0.01 (type)
# So yes, the gender data vary by region


##################
#### Is sex related to 16S gut composition?

##### Run the model, ignoring missing data for now 
meta_sex <- subset(meta_maaslin,is.na(meta_maaslin$sex)==F)
data_sex <- subset(df_IM_data, row.names(df_IM_data) %in% meta_sex$sample)
meta_sex$sex = as.factor(meta_sex$sex)

fit_data_1 <- Maaslin2(
  data_sex, meta_sex, 'SexOnly',
  fixed_effects = c("sex"),
  standardize = FALSE)

sex_results<-fit_data_1$results 
subset(sex_results,sex_results$qval <= 0.05) 

# --> only one taxon, UCG.002, is associated with sex (less abundant in males)


# Since TH tribes have fewer males, based on this model alone, they might be expected to have more UCG.002
# But if there aren't really fewer males....then UCG.002 might be more abundant in TH for another reason. Like diet.
# Modeling non-missing data only is misleading 
# (also, if sex is treated as three categories, nothing is significant); # meta_maaslin$sex = as.factor(meta_maaslin$sex) + #   reference = c("sex", "M"),


##### Run the model again, this time imputing missing sex data 
counts2 
# most likely, the 1 DP NA is F and the 6 TH NAs are M ( based on line 673, "5 men and 5 women" per tribe)
# fix this in counts_long
temp0 <- subset(meta_maaslin,is.na(meta_maaslin$sex)==F)
temp <- subset(meta_maaslin,is.na(meta_maaslin$sex)==T)
temp[1,4]= "F"
temp[2:7,4]= "M"
meta_imputed = rbind(temp0,temp)
#
fit_data_3 <- Maaslin2(
  df_IM_data, meta_imputed, 'SexOnly',
  fixed_effects = c("sex"),
  standardize = FALSE)

#""There are no associations to plot!""

# So, I think this means that modeling sex will be misleading, not only due to omitting data, 
# but by artificially (for the purposes of this study) conflating sex and Trans-Himalayas




##########    AGE    ###########
#################################

# Does age vary by region? --> not significantly
counts <- meta_maaslin %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
      n_NA     = sum(is.na(Age)),
      n_obs    = sum(!is.na(Age)),
      mean_Age = mean(Age, na.rm = TRUE),
      sd_Age   = sd(Age, na.rm = TRUE),
      min_Age  = min(Age, na.rm = TRUE),
      max_Age  = max(Age, na.rm = TRUE)
    )

boxplot(Age ~ Region, data = meta_age,
        main = "Age by Region",
        xlab = "Region", ylab = "Age")
anova_fit <- aov(Age ~ Region, data = meta_age)
summary(anova_fit) # p = 0.21, not significant when missing data are ignored


##################
#### Is age related to 16S gut composition?
meta_age <- subset(meta_maaslin, is.na(meta_maaslin$Age)==F)
data_age <- subset(df_IM_data, row.names(df_IM_data) %in% row.names(meta_age))
fit_data_2 <- Maaslin2(
  data_age, meta_age, 'AgeOnly',
  fixed_effects = c("Age"),
  standardize = FALSE)

# "There are no associations to plot!"

## Age isn't related to microbiome composition in the available data. 
# So, trying to model it would reduce sample size/power for no clear reason.
# Unlike sex, age can't be imputed in this data set.




##########    BMI    ###########
#################################

# BMI has a lot of missing missing data
counts <- meta_maaslin %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    n_NA     = sum(is.na(BMI)),
    n_obs    = sum(!is.na(BMI))
  )
counts

# No Kabui samples have BMI, so they'll need to be removed to summarize the rest
meta_temp <- subset(meta_maaslin,meta_maaslin$Tribe != "Kabui")

# Look at how BMI varies among other regions
boxplot(BMI ~ Region, data = meta_temp,
        main = "BMI by Region",
        xlab = "Region", ylab = "BMI")
anova_fit <- aov(BMI ~ Region, data = meta_temp)
summary(anova_fit) # p = 0.00024 -- significant variation across regions


##################
#### Is BMI related to 16S gut composition?
meta_BMI <- subset(meta_maaslin,is.na(meta_maaslin$BMI)==F)
data_BMI <- subset(df_IM_data, row.names(df_IM_data) %in% meta_BMI$sample)

fit_data_4 <- Maaslin2(
  data_BMI, meta_BMI, 'BMIOnly',
  fixed_effects = c("BMI"),
  standardize = FALSE)

BMI_results<-fit_data_4$results 
sig_BMI_results <- subset(BMI_results,BMI_results$qval <= 0.05) 

# ## 
# Bifidobacterium*
# Lactobacillus*
# Megasphaera*
# Eubacterium_eligens_group
# Ligilactobacillus*
# UCG.008
# * = associated with higher BMI

## Do the TH tribes have higher BMI? 
# Who is missing BMI?
summary_df <- meta_maaslin %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    mean_BMI = mean(BMI, na.rm = TRUE),
    sd_BMI   = sd(BMI, na.rm = TRUE),
    n_obs    = sum(!is.na(BMI)),
    n_NA     = sum(is.na(BMI))
  )
# BMI reporting ranges from 0% to 100%

# from available data, TH are substantially heavier (mean 24.7 vs 19.7)





###################################################
## Do covariates (sex, age, BMI) vary by region?





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
#fit_data_IM  <- readRDS("fit_data_IM.rds")





# EXAMPLE 
# fit_data <- Maaslin2(
#   taxo, meta_example, 'demo_output',
#   fixed_effects = c('diagnosis', 'dysbiosisnonIBD','dysbiosisUC','dysbiosisCD', 'antibiotics', 'age'), # character or itegerclas
#   random_effects = c('site', 'subject'), # character
#   reference = "diagnosis,nonIBD",
#   standardize = FALSE)


##########################################
# WHAT METADATA VARIABLES NEED TO BE CORRECTED FOR? NEW 10/22/25


#######

# Actual model with SEX, AGE, and BMI covariates



### OVERALL,
# Gender, Run, and BMI are all associated with Region; so I fear "correcting" for those will
# remove any real signal from Region. 

# Let's try with BMI just to see.
IM_fit_data <- Maaslin2(
  taxo2, meta16_dt, 'region_plus_BMI',
  fixed_effects = c("Region","BMI"),
  # random_effects = c('sample'),
  reference = c("Region,Trans-Himalayas"),
  standardize = FALSE)

# to isolate Region from BMI covariate:
maaslin2_all_results<-IM_fit_data$results 
maaslin2_results<-maaslin2_all_results %>% filter(metadata == 'Region') # Discard covariate associations
maaslin2_results$qval<-p.adjust(maaslin2_results$pval, method = 'BH') # re-calculate the q-values 
sig_results <- subset(maaslin2_results,maaslin2_results$qval < 0.1)
# --> now, there are no results from Kabui because none of them have BMI data
# so this doesn't seem like it's going to work 
# maybe just run the test with no covariates and note that TH people are heavier + more female


#######

# Actual model with SEX and AGE covariates




#######

# Actual model with SEX covariate




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