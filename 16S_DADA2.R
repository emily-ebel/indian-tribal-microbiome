
# source: https://f1000researchdata.s3.amazonaws.com/manuscripts/10545/a565a995-e2d3-4069-a0e1-52ed93a99b14_8986_-_susan_holmes_v2.pdf?doi=10.12688/f1000research.8986.2&numberOfBrowsableCollections=77&numberOfBrowsableInstitutionalCollections=4&numberOfBrowsableGateways=4

### PREP WORKSPACE ###
 
setwd("/home/ebel/user_data/16S/IndianMicrobiome")

library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn") # to install these, need to use Bioconductor (see above)
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
set.seed(101)


#######################################################################
### Load FASTQ paths 
#######################################################################
rawseq_path = "/home/ebel/user_data/16S/IndianMicrobiome/all_fastq"
filt_path = "/home/ebel/user_data/16S/IndianMicrobiome/all_fastq/trimmed_reads" # this path must EXIST -- will be filled later
fns <- sort(list.files(rawseq_path, full.names = TRUE))
fnFs <- fns[grepl("_R1", fns)]
fnRs <- fns[grepl("_R2", fns)]


#######################################################################
### Inspect data, trim, and write to filt_path
#######################################################################

ii <- sample(length(fnFs), 3) # randomly pick 3 samples 
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) } # plot function from dada2
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

# for the IM samples, quality is consistently over 30 across the whole 250 bp read 
# but, the first 10 bp are lower-quality, so will still trim those  ("likely to contain pathological errors")
# --> actually trim 20, based on after pics

# We combine these trimming parameters with standard filtering parameters, the most important 
# being the enforcement of a maximum of 2 expected errors per-read. Trimming and filtering is 
# performed on paired reads jointly, i.e. both reads must pass the filter for the pair to pass.

filtFs <- file.path(filt_path, basename(fnFs)) # R1 
filtRs <- file.path(filt_path, basename(fnRs)) # R2

for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=20, truncLen=c(250, 250),  
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)
}


## check the data by resetting some variables
rawseq_path = filt_path 
fns <- sort(list.files(rawseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]
ii <- sample(length(fnFs), 3) # randomly pick 3 samples 
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) } # plot function from dada2
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }





#######################################################################
# Build model to identify errors and correct them 
#######################################################################

# Import filtered sequence data and simultaneously dereplicate.
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
#
# Name the resulting derep-class objects by their sample name.
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sam.names
names(derepFs) <- sam.names
names(derepRs) <- sam.names


# Infer parameters for a model to distinguish sequencing errors from substitution errors.
# Using a subset of the data
# This takes forever without multithreading, so multithread! Also using 5 fewer samples. 
ddF <- dada(derepFs[1:40], err=NULL, selfConsist=TRUE, multithread=TRUE) # 
ddR <- dada(derepRs[1:40], err=NULL, selfConsist=TRUE,  multithread=TRUE) # convergence after ___ round (7 rounds tutorial data), taking about 5 minutes 
# takes 9 rounds

# Inspect the fit between the observed (points) and fitted (lines) error rates
plotErrors(ddF); plotErrors(ddR)

# Now do the sequence inference, which removes sequencing errors.
# This step can be pooled (more accurate, slower) or unpooled. 
# Pooling is intractable on the scale of tens of millions of reads. 
# Here we have 5 million reads, but had issues pooling, so don't pool 

dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=F, multithread = TRUE)  # start 4:57, end by 5:12
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=F, multithread = TRUE)  # start 5:12

# Now merge together the inferred forward and reverse sequences, removing 
# paired sequences that do not perfectly overlap as a final control against residual errors.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# these steps take a while, so write object to file 
saveRDS(mergers, file = "IM_mergers_110823.rds")


#######################################################################
# Construct sequence table and remove chimeras 
#######################################################################

mergers <- readRDS("IM_mergers_110823.rds")
seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
dim(seqtab.all)

# remove chimeric sequences 
seqtab <- removeBimeraDenovo(seqtab.all)
dim(seqtab)

# Variant numbers go from 10704 -> 3007, so removes 72% of variants 



#######################################################################
# Assign taxonomy
#######################################################################

ref_fasta <- "/home/ebel/user_data/16S/dada2_Holmes_tutorial/data/silva_nr99_v138.2_toGenus_trainset_copy.fa.gz"  # SILVA update -- 138.2 vs 138 -- hopefully containing Segatella
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
taxtab2 <- as.data.frame(taxtab)
length(unique(taxtab2$Genus)) 
length(unique(taxtab2$Species))  


#######################################################################
# Construct phylogenetic tree
#######################################################################

# First, perform multiple alignment using DECIPHER 
seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree -- just the sequences themselves 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# Then build a tree using phangorn. NJ is used as a starting point 
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order



#######################################################################
# Combine data into phyloseq object
#######################################################################

# First, get and clean up metadata 
# the sample IDs need to correspond to:
rownames(seqtab)

mimarks_path <- "/home/ebel/user_data/16S/IndianMicrobiome/IM_16S_metadata_R.csv" # 
samdf <- read.csv(mimarks_path, header=TRUE)
all(rownames(seqtab) %in% samdf$sample_id) # TRUE
rownames(samdf) <- samdf$sample_id


# The full suite of data for this study – the sample-by-sequence feature table,
# the sample metadata, the sequence taxonomies, and the phylogenetic tree – 
# can now be combined into a single object.
# ps <- phyloseq(tax_table(taxtab), sample_data(samdf),
#                otu_table(seqtab.copy, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

ps <- phyloseq(tax_table(taxtab), sample_data(samdf),
               otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(treeNJ))




#######################################################################
# Prevalence filtering
#######################################################################
library("phyloseq")
library("gridExtra")
ps

# If you expect the set of all organisms from all samples to be well-represented 
# in the available taxonomic reference database,it is reasonable or even 
# advisable to filter taxonomic features for which a high-rank taxonomy could 
# not be assigned. Such ambiguous features in this setting are almost always 
# sequence artifacts that don’t exist in nature. It should be obvious that such 
# a filter is not appropriate for samples from poorly characterized or novel specimens.

# To begin filtering, create a table of read counts for each Phylum present
rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL)
# remove NAs (without phylum)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))



# Compute prevalence of each feature & store as data.frame
# Prevalence is defined as the number of samples in which a taxon appears at least once.
prevdf = apply(X = otu_table(ps0),
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame and examine 
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))
phylumcount = plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
phylumcount

# Choose phyla to filter, based on lowest prevalence 
# In this case, phyla with a total prevlanece (col2) of < 10
filterPhyla = c("Thermoplasmatota","Patescibacteria","Methanobacteriota","Synergistota") # SILVA 138.2
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1
# reexamine 
prevdf = apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps1),
                    tax_table(ps1))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

## how many genera are left now? 
length(unique(prevdf$Genus)) # 245 SILVA 138.2


## explore the relationship of prevalence and total read count
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter yintercept, which will be a threshold
  geom_hline(yintercept = 0.026, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# In this case, I"m deciding to keep taxa that are detected in 2 or more samples (2.6% or greater)

# Based on these results, choose a prevalence threshold & execute filter
prevalenceThreshold = 0.026 * nsamples(ps0)
prevalenceThreshold = 2
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)

# how many genera left now?
prevdf = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps2),
                    tax_table(ps2))
length(unique(prevdf$Genus)) # 151 RDP, 210 GTDP, 167 SILVA 138.2



# get a new taxtab now
taxtab_filt = as.data.frame(tax_table(ps2))
n <- taxtab_filt %>%
  group_by(Genus) %>%
  dplyr::summarise(rsv_n = n())
n <- subset(n,is.na(n$Genus)==F)
names(n)=c("Genus","16S variants")

## WRITE THIS TAX TAB TO FILE
write.table(taxtab_filt,file="IM_16S_relab_SILVA-138.2_non-aglommerated_010725.csv")



#######################################################################
# Agglomerate taxa
#######################################################################

# Combine all features that descend from the same genus.
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)


# If taxonomy is not available, can do tree-based agglomeration by defining
# a tree height (similar to OTU clustering)
#h1 = 0.4
#ps4 = tip_glom(ps2, h = h1)

# plot the original tree + 2 types of agglomerated trees
multiPlotTitleTextSize = 8
#p2tree = plot_tree(ps2, method = "treeonly", ladderize = "left",
#title = "Before Agglomeration") +  theme(plot.title = element_text(size = multiPlotTitleTextSize))

p3tree = plot_tree(ps3, method = "treeonly", 
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree
#p4tree = plot_tree(ps4, method = "treeonly",
#                   ladderize = "left", title = "By Height") +
#  theme(plot.title = element_text(size = multiPlotTitleTextSize))
# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)




#######################################################################
# Abundance value transformation 
#######################################################################

# transformation of count data [to relative abundance] is necessary to account for differences in library size, 
# variance, scale, etc. # transform_sample_counts() function provides a flexible way, 
#by providing (as the second argument) a function that defines the transformation

# In this example, we'll first define a custom plot function that allows us
# to compare abundance scale/distribution before and after transformation 

plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Region",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# The actual transformation here is simply relative abundance. Define 
# this simple transformation function within phyloseq's transform_sample_counts


#### accurate variable name
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})

### try getting rel_ab table here
percentages <- psmelt(ps3ra)
write.table(percentages,file="IM_16S_relab_SILVA-138.2_aglommerated_010725.csv")



#plot the abundance values before and after transformation.
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2, plotBefore, plotAfter)

# plot one order if desired
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)

