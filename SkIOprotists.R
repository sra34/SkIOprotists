## Run DADA2 analysis, following DADA2 pipeline tutorial (1.12) - https://benjjneb.github.io/dada2/tutorial.html and workflow from AstrobioMike - https://astrobiomike.github.io/amplicon/dada2_workflow_ex

# Install and load library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
library(dada2)

# Define the path to the folder which has all the .fastq files (you can keep them zipped!)
path <- "C:/Users/sean_/Desktop/18S metagenomic samples/SkIOprotists"
list.files(path)

# Partition forward and reverse reads into separate strings to view quality profiles of the reads
forward <- sort(list.files(pattern="_R1_001.fastq", full.names = TRUE))
reverse <- sort(list.files(pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names based on naming format (date_replicate_plate position)
sample.names <- sapply(strsplit(basename(forward), "_"), `[`, 1)

# Plot and observe the quality profiles of forward and reverse reads (allows you to cutoff reads at certain threshold). In this case, we plot the first 10 forward and reverse reads.
plotQualityProfile(forward[1:10])
plotQualityProfile(reverse[1:10])

# Place filtered samples into a new folder directory called "filtered"
filtered_forward <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtered_reverse <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Trim forward and reverse sequences based on specifications and prior observations of quality profiles. The trimLeft function is used to manually trim primers that have not yet been removed.
out <- filterAndTrim(forward, filtered_forward, reverse, filtered_reverse, truncLen=c(200,240), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, matchIDs=TRUE, trimLeft = c(20, 18), multithread=FALSE) # On Windows set multithread=FALSE

# Learn the error rates of the forward and reverse reads testing the first 100M bp (default)
errF <- learnErrors(filtered_forward, multithread=FALSE)
errR <- learnErrors(filtered_reverse, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplication of samples to combine identical sequences and give abundances of sequences 
derepFs <- derepFastq(filtered_forward, verbose=TRUE)
derepRs <- derepFastq(filtered_reverse, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Apply core sample inference algorithm to the deprelicated data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# Construct a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove sequence lengths much shorter or larger than the region of interest
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(365,390)]

# Remove chimeras and check to see how many chimeras you have (% of reads)
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

# Track number of reads that have made it through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxonomy using PR2 database .fasta file (install most recent version)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/sean_/Desktop/18S metagenomic samples/SkIOprotists/pr2_version_4.11.1_dada2.fasta.gz", multithread=TRUE, taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"))

# Extract all the files needed for handoff to phyloseq (e.g. ASV count and taxonomy tables)
asv_seqs <- colnames(seqtab.nochim)
write.table(asv_seqs, "ASVs_rawseqs.tsv", sep="\t", quote=F, col.names=NA)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)


## Load libraries needed for phyloseq analysis
library(phyloseq)
library(ggplot2)
library (viridis)
library(dplyr)
library(devtools)
library(ranacapa)
library(vegan)
library(plyr)


# Read in your count, taxonomy, and raw sequences tables. Create data frame that will be used to create Supplemental Table 2 with filtered ASV counts, annotations for each ASV, and raw seqeunces. This data frame will be curated for easy phyloseq input.  
count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=NULL, check.names=F, sep="\t")
tax_tab18S <- read.table("ASVs_taxonomy_061119.tsv", header=T, row.names=NULL, check.names=F, sep="\t") # Taxonomy file was manually edited to add missing species level annotation (e.g. Chaetoceros sp.) if genus was annotated properly (e.g. Chaetoceros).
rawseq <- read.table("ASVs_rawseqs.tsv", header=T, row.names=NULL, check.names=F, sep="\t") # Add sequence information for Supplemental Table 2
merge <- cbind(count_tab,tax_tab18S,rawseq)
colnames(merge)[1]<-"OTU.ID"
merge$Var.98 <- NULL
merge$Var.107 <- NULL

# Number of sequences per sample and mean
colsum<-apply(merge[2:97],2,sum); colsum
mean(colsum)

# Remove metazoan sequences, expected Syndiniales parasites of metazoans (Group IV), and "NA" sequences at the Supergroup level
mergenew <- subset(merge, !grepl("Metazoa",merge$Division)); dim(mergenew)
mergenew <- subset(mergenew, !grepl("Dino-Group-IV", mergenew$Order)); dim(mergenew)
mergenew <- mergenew[ !is.na(mergenew$Supergroup),]; dim(mergenew)

# Rearrange large data frame for Supplemental Table 2 (ASV counts, annotation, and sequences)
merge_table <- mergenew[, c(1, 98:105, 106, 17:97, 2:16)]
colnames(merge_table)[c(1:106)]<-c("ASV ID","Kingdom", "Supergroup", "Division", "Class", "Order", "Family",	"Genus", "Species", "Sequence", "3/16/17 A","3/16/17 B","3/16/17 C","3/22/17 A","3/22/17 B","3/22/17 C","3/29/17 A","3/29/17 B","3/29/17 C","4/12/17 A","4/12/17 B","4/12/17 C","4/20/17 A","4/20/17 B","4/20/17 C","5/3/17 A","5/3/17 B","5/3/17 C","6/7/17 A","6/7/17 B","6/7/17 C","6/15/17 A","6/15/17 B","6/15/17 C","6/22/17 A","6/22/17 B","6/22/17 C","7/6/17 A","7/6/17 B","7/6/17 C","7/21/17 A","7/21/17 B","7/21/17 C","7/27/17 A","7/27/17 B","7/27/17 C","8/3/17 A","8/3/17 B","8/3/17 C","8/10/17 A","8/10/17 B","8/10/17 C","8/16/17 A","8/16/17 B","8/16/17 C","8/23/17 A","8/23/17 B","8/23/17 C","8/30/17 A","8/30/17 B","9/6/17 A","9/6/17 B","9/6/17 C","9/20/17 A","9/20/17 B","9/20/17 C","9/26/17 A","9/26/17 B","9/26/17 C","10/4/17 A","10/4/17 B","10/4/17 C","10/11/17 A","10/11/17 B","10/18/17 A","10/18/17 B","10/18/17 C","11/8/17 A","11/8/17 B","11/8/17 C","11/16/17 A","11/16/17 B","11/16/17 C","11/21/17 A","11/21/17 B","12/8/17 A","12/8/17 B","12/8/17 C","12/14/17 A","12/14/17 B","12/14/17 C","1/18/18 A","1/18/18 B","1/18/18 C","1/24/18 A","1/24/18 B","1/24/18 C","2/1/18 A","2/1/18 B","2/1/18 C","2/8/18 A","2/8/18 B","2/8/18 C", "2/21/18 A", "2/21/18 B", "2/21/18 C")
rowsum<-apply(merge_table[11:106],1,sum) # Remove global singletons 
merge_table = merge_table[ rowsum>1, ]
merge_table[, c(11:106)] <- lapply( merge_table[ ,c(11:106)], function(x) x/sum(x, na.rm=TRUE) )
write.csv(merge_table, file="merge_table.csv", row.names=FALSE)


# Prepare count and tax table to be imported into phyloseq object (this may seem like extra work - you can directly input DADA2 files into phyloseq, but I wanted to create a master table first)
count_tab2 <- mergenew[, c(1:97)]
write.table(count_tab2, file="count_tab2.txt", quote=FALSE, sep="\t", row.names=FALSE)
tax_tab2 <- mergenew[, c(1, 98:105)]
write.table(tax_tab2, file="tax_tab2.txt", quote=FALSE, sep="\t", row.names=FALSE)

# Import new filtered count and tax files and create phyloseq object
tax_tab <- as.matrix(read.table("tax_tab2.txt", header=T, row.names=1, check.names=F, sep="\t"))
count_tab <- read.table("count_tab2.txt", header=T, row.names=1, check.names=F, sep="\t")
sample_info_tab <- read.table("Sampleinfo.txt", header=T, row.names=1, check.names=F, sep="\t") # Add metadata file that corresponds to samples
ps <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE),tax_table(tax_tab), sample_data(sample_info_tab)) # Merge everything together - we can also add a phylogenetic tree file if we want here

# Plot rarefaction curves and facet them by sampling month - Supplemental Figure 1
rare <- ggrare(ps, step = 10, plot = TRUE, parallel = FALSE, se = FALSE)
rare$data$Month <- factor(rare$data$Month, levels = c("Mar", "Apr", "May", "Jun","Jul", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb"))
rare + theme(legend.position = "none") + facet_wrap(Month ~ ., ncol=3) +theme_bw()+theme(legend.position = "none")

# Plot alpha Shannon diversity for total 18S community and only Syndiniales over time - Figure 2 (a, b)
alphadiv <- plot_richness(ps, x = "Date", measures=c("Shannon"), color = "Temp") + geom_boxplot(stat = "boxplot", lwd =1) + scale_colour_viridis() + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12)) + theme(axis.text.y = element_text(size=12))+ theme(legend.text = element_text(colour="black", size=12))
alphadiv$data$Date = factor(alphadiv$data$Date, levels = c(31617, 32217, 32917, 41217, 42017, 50317,60717, 61517, 62217, 70617, 72117, 72717, 80317, 81017, 81617, 82317, 83017, 90617, 92017, 92617, 100417, 101117, 101817, 110817, 111617, 112117, 120717, 121417, 11618, 12418, 20118, 20818, 22118))
alphadiv$layers <- alphadiv$layers[-1]
plot(alphadiv)

ps_Syn <- subset_taxa(ps, Class=="Syndiniales") # Subset for only Syndiniales sequences
alphadiv_Syn <- plot_richness(ps_Syn, x = "Date", measures=c("Shannon"), color = "Temp") + geom_boxplot(stat = "boxplot", lwd =1) + scale_colour_viridis()+ theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12)) + theme(axis.text.y = element_text(size=12))+theme(legend.text = element_text(colour="black", size=12))
alphadiv_Syn$data$Date = factor(alphadiv_Syn$data$Date, levels = c(31617, 32217, 32917, 41217, 42017, 50317,60717, 61517, 62217, 70617, 72117, 72717, 80317, 81017, 81617, 82317, 83017, 90617, 92017, 92617, 100417, 101117, 101817, 110817, 111617, 112117, 120717, 121417, 11618, 12418, 20118, 20818, 22118))
alphadiv_Syn$layers <- alphadiv_Syn$layers[-1]
plot(alphadiv_Syn)

# NMDS bray curtis for total 18S community and only Syndiniales - Figure 2(c, d)
ps = filter_taxa(ps, function (x) {sum(x) > 1}, prune=TRUE) # Phyloseq object will now have global singletons removed for subequent analyses
ps_rel  <- transform_sample_counts(ps, function(x) x / sum(x) ) # Convert absolute counts to relative abundances
ord_total <- ordinate(ps_rel, "NMDS", "bray")
p = plot_ordination(ps_rel, ord_total, color="Temp") + guides(color=FALSE)
p + theme_bw() + theme(text = element_text(size = 16)) + geom_point(aes(fill=Temp),size = 5, shape = 21, colour = "black") + scale_fill_viridis() +theme(legend.text = element_text(colour="black", size=12))

ps_rel_Syn <- subset_taxa(ps_rel, Class=="Syndiniales")
ord_Syn <- ordinate(ps_rel_Syn, "NMDS", "bray")
p = plot_ordination(ps_rel_Syn, ord_Syn, color="Temp") + guides(color=FALSE)
p +  theme_bw() + theme(text = element_text(size = 16)) + geom_point(aes(fill=Temp),size = 5, shape = 21, colour = "black") + scale_fill_viridis() +theme(legend.text = element_text(colour="black", size=12))


# Plot faceted taxonomy bar plots at class level - Figure 1a
x1 <- tax_glom(ps, taxrank = 'Class') # Agglomerate taxa based on level
x2 <- transform_sample_counts(x1, function(x) x/sum(x)) # Transform to relative abundance
x3 <- psmelt(x2) # Melt the data for ggplot data frame
x3 <- x3 %>% # Calculate mean and sd of replicate abundance values
  dplyr::group_by(Class, Date) %>%
  dplyr::summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup

x3$Class <- as.character(x3$Class) # Convert to character
x3$Class[x3$Mean < 0.05]= "Other"
write.table(x3, file="18S_abundances.txt", quote=FALSE, sep="\t") # Unable to get error bars for the "other" category. Write table for accurate mean and sd of "other group" and than reimport new .txt file which has proper mean and sd for "other". 
x5 <- read.table("18S_abundances_new.txt", header=T, row.names = 1, check.names=F, sep = "\t") # This new table will be very similar to the one you exported, just have proper mean and sd for the "other" group and the same values for the other major class level groups.
p <- ggplot(data=x5, aes(x=Date, y=Mean, fill=Class)) # Plot your bars, adjust settings as necessary
p$data$Date <- factor(p$data$Date, levels = c(31617, 32217, 32917, 41217, 42017, 50317,60717, 61517, 62217, 70617, 72117, 72717, 80317, 81017, 81617, 82317, 83017, 90617, 92017, 92617, 100417, 101117, 101817, 110817, 111617, 112117, 120717, 121417, 11618, 12418, 20118, 20818, 22118))
p$data$Class <- factor(p$data$Class, levels = c("Other", "Filosa-Thecofilosea", "Spirotrichea", "Cryptophyceae","Dinophyceae", "Mamiellophyceae", "Syndiniales", "Bacillariophyta"))
p + geom_bar(stat="identity", position = "stack")+scale_y_continuous(expand = c(0, 0), breaks=seq(0,0.6,0.2), limits=c(0, 0.6)) + geom_hline(yintercept=0) + facet_wrap(~Class, ncol=2) + geom_errorbar(aes(ymax = Mean+SD, ymin = Mean-SD, width=0.5)) + scale_fill_manual(values= c("#757575","#EFE343","#A0CEEF","#DB4141","#7fbf7b","#0273B2","#762a83","#1b7837"))+ theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5, size=10)) + theme(axis.text.y=element_text(size=10))+ theme(legend.position="bottom") + guides(fill=guide_legend(nrow=3)) + theme(axis.title.x=element_blank())+ labs(y = "Relative Abundance")

# Plot faceted taxonomy bar plots Within Syndiniales (Order level) - Figure 1b
ps_Syn_bar <- subset_taxa(ps, Class=="Syndiniales")
x1_Syn <- tax_glom(ps_Syn_bar, taxrank = 'Order')
x2_Syn <- transform_sample_counts(x1_Syn, function(x) x/sum(x))
x3_Syn <- psmelt(x2_Syn)
x3_Syn <- x3_Syn %>% # Calculate mean and sd of replicate abundance values
  dplyr::group_by(Order, Date) %>%
  dplyr::summarize(Mean = mean(Abundance), SD = sd(Abundance)) %>%
  ungroup
x3_Syn$Order <- as.character(x3_Syn$Order)
x3_Syn$Order[x3_Syn$Mean < 0.05]= "Other"
p <- ggplot(data=x3_Syn, aes(x=Date, y=Mean, fill=Order)) # Plot your bars, adjust settings as necessary
p$data$Date <- factor(p$data$Date, levels = c(31617, 32217, 32917, 41217, 42017, 50317,60717, 61517, 62217, 70617, 72117, 72717, 80317, 81017, 81617, 82317, 83017, 90617, 92017, 92617, 100417, 101117, 101817, 110817, 111617, 112117, 120717, 121417, 11618, 12418, 20118, 20818, 22118))
p$data$Order <- factor(p$data$Order, levels = c("Other", "Dino-Group-I", "Dino-Group-II", "Dino-Group-III"))
p + geom_bar (aes(), stat="identity", position="stack")+scale_y_continuous(expand = c(0, 0), breaks=seq(0,1,0.2), limits=c(0, 1)) + geom_hline(yintercept=0)+ facet_wrap(~Order, ncol=1) + geom_errorbar(aes(ymax = Mean+SD, ymin = Mean-SD, width=0.5))+ scale_fill_manual(values= c("#757575","#762a83","#9970ab","#c2a5cf")) + theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5, size=10)) + theme(axis.text.y=element_text(size=10))+ theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2)) + theme(axis.title.x=element_blank())+ labs(y="Relative Abundance")

 
# Prepare and execute CCA analysis (get rid of sample day 090617 that does not have nutrient data) - Supplemental Figure 2a
ps_subset <- subset_samples(ps, sample_names(ps) != "090617-A-E3" & sample_names(ps) != "090617-B-E4" & sample_names(ps) != "090617-C-E5")
metadata <- as(sample_data(ps_subset), "data.frame") # Subset sample metadata to a data frame
metadata[, 5:15] <- log1p((metadata[5:15])) # Log transform data that will be added to CCA
sample_data(ps_subset) <- metadata # Upload the data frame back to the phyloseq object used for CCA
ps2  <- transform_sample_counts(ps_subset, function(x) x / sum(x) )
ps2_ord1 <- ordinate(physeq = ps2, method = "CCA", distance = "bray", formula = ~ 1)
ps2_ord2 <- ordinate(physeq = ps2, method = "CCA", distance = "bray", formula = ~ Temp + Salinity + Chlorophyll + Nitrate + Ammonium + Phosphate + Silicate + Solar + POC + PON + DO)
ordi <- ordistep(ps2_ord1, perm.max=999, steps = 50, scope=formula(ps2_ord2)) # Step through factors to see which are significant for the ordination compared to just the intercept
cca_plot <- plot_ordination(physeq = ps2, ordination = ps2_ord2, axes = c(1,2))  + xlim(-3,3)+ ylim(-3,4) + theme_bw() # Edit to include significant factors via ordistep and run the ordination
cca_plot$layers <- cca_plot$layers[-1]
cca_plot$data$Month <- factor(cca_plot$data$Month, levels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb"))
arrowmat <- vegan::scores(ps2_ord2, display = "bp") # Include factors as biplot arrows on the plot
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = CCA1, yend = CCA2, x = 0, y = 0, shape = NULL, color = NULL, label = labels) # Maps the arrows on the CCA plot
label_map <- aes(x = 1.3 * CCA1, y = 1.3 * CCA2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))
cca_plot + geom_point(aes(fill=Month),size = 5, shape = 21, colour = "black") + scale_fill_viridis(discrete=TRUE, option="plasma", direction = -1) + geom_segment(mapping = arrow_map, size = 1, data = arrowdf, color = "black", arrow = arrowhead) + geom_text(mapping = label_map, size = 4, color="black", data = arrowdf, show.legend = FALSE)



# Prepare and execute CCA analysis of only Syndiniales (get rid of sample day 090617 that does not have nutrient data) - Supplemental Figure 2b
ps_Syn_CCA <- subset_taxa(ps2, Class=="Syndiniales") # Subset to only have Syndiniales - this phyloseq object will have the log-transformed metadata from the prior total 18S CCA
ps_Syn_ord1 <- ordinate(physeq = ps_Syn_CCA, method = "CCA", distance = "bray", formula = ~ 1)
ps_Syn_ord2 <- ordinate(physeq = ps_Syn_CCA, method = "CCA", distance = "bray", formula = ~ Temp + Salinity + Chlorophyll + Nitrate + Ammonium + Phosphate + Silicate + Solar + POC + PON + DO)
ordi_Syn <- ordistep(ps_Syn_ord1, perm.max=999, steps = 50, scope=formula(ps_Syn_ord2)) # Step through factors to see which are significant for the ordination
cca_Syn_plot <- plot_ordination(physeq = ps_Syn_CCA, ordination = ps_Syn_ord2, axes = c(1,2)) + xlim(-4,3)+ ylim(-5,3) + theme_bw()
cca_Syn_plot$layers <- cca_Syn_plot$layers[-1]
cca_Syn_plot$data$Month <- factor(cca_Syn_plot$data$Month, levels = c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec", "Jan", "Feb"))
arrowmat <- vegan::scores(ps_Syn_ord2, display = "bp")
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = CCA1, yend = CCA2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
label_map <- aes(x = 1.3 * CCA1, y = 1.3 * CCA2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))
cca_Syn_plot + geom_point(aes(fill=Month),size = 5, shape = 21, colour = "black") + scale_fill_viridis(discrete=TRUE, option="plasma", direction = -1)+ geom_segment(mapping = arrow_map, size = 1, data = arrowdf, color = "black", arrow = arrowhead) + geom_text(mapping = label_map, size = 4, color="black", data = arrowdf, show.legend = FALSE)


## Prepare for CoNet network analysis with original merged ASV and taxonomy data frame - Used for Figure 3 and Supplemental Figure 3. 
merge_Conet <- mergenew[ !is.na(mergenew$Genus),]; dim(merge_Conet) # Remove taxa not identified to genus level
rowsum<-apply(merge_Conet[2:97],1,sum) # Remove global singletons 
merge_Conet = merge_Conet[ rowsum>1, ]
df <- merge_Conet[c(1:97)] # Get samples only file
names(merge_Conet) # Get key or attributes file
key<-merge_Conet[c(1, 98:105)]
head(key)
dfnew <- df[rowMeans(df!=0)>0.50,] # ASVs have to be found in greater than half of all sample days 

# Filter key and attach to data frame for CoNet network analyses (total year data)
tmp<-as.data.frame(dfnew$`OTU.ID`); colnames(tmp)<-"OTU.ID"
key_filtered<-join(tmp, key)
merge_final <- cbind(dfnew,key_filtered)
dfnew2<- merge_final[, c(1:66, 70:97,102:106)] # Filter out 9/6/17 from the dataset (incomplete metadata for this day) and only filter to class level taxonomy and below
colnames(dfnew2)[c(1:94)]<-c("OTU.ID","11618A","11618B","11618C","12418A","12418B","12418C","20118A","20118B","20118C","20818A","20818B","20818C","22118A","22118B","22118C","31617A","31617B","31617C","32217A","32217B","32217C","32917A","32917B","32917C","41217A","41217B","41217C","42017A","42017B","42017C","50317A","50317B","50317C","60717A","60717B","60717C","61517A","61517B","61517C","62217A","62217B","62217C","70617A","70617B","70617C","72117A","72117B","72117C","72717A","72717B","72717C","80317A","80317B","80317C","81017A","81017B","81017C","81617A","81617B","81617C","82317A","82317B","82317C","83017A","83017B","92017A","92017B","92017C","92617A","92617B","92617C","100417A","100417B","100417C","101117A","101117B","101817A","101817B","101817C","110817A","110817B","110817C","111617A","111617B","111617C","112117A","112117B","120717A","120717B","120717C","121417A","121417B","121417C")
dfnew2=dfnew2[c(1:99)]


dfnew2$Class <- sub("^", "c_", dfnew2$Class) # Add identifier before each taxonomy level to format for CoNet program
dfnew2$Order <- sub("^", "o_", dfnew2$Order)
dfnew2$Family <- sub("^", "f_", dfnew2$Family)
dfnew2$Genus <- sub("^", "g_", dfnew2$Genus)
dfnew2$Species <- sub("^", "s_", dfnew2$Species)
dfnew2$taxonomy <- paste(dfnew2$Class, dfnew2$Order, dfnew2$Family, dfnew2$Genus, dfnew2$Species, sep = "; ") # Only need up to class level, makes new column called "taxonomy" that combines annotation 
dfnew2 <-dfnew2[, c(1:94,100)] # This table will be ready for easy input into CoNet
write.table(dfnew2, file="18S_fornetwork.txt", quote=FALSE, sep="\t", row.names=FALSE)

## Some more tips before using CoNet:
# Manually change labels of Syndiniales ASV to SSV to rapidly distinguish between Syndiniales and other protist taxa (e.g. for filtering overall network)
# Slightly modify this .txt file to be readable in CoNet, which is a plugin for Cytoscape (see Faust and Raes 2016). CoNet also has some options to run networks in R, though I have not tested this. 
# Add "# Constructed from biom file" as the first line in the .txt file. Add "#" before OTU.ID. Then upload file to CoNet - you will have the option to normalize ASV abundance table in CoNet.
# Upload .txt file of metadata with matching sample names
# I will provide the files needed to run the network in CoNet (18S_fornetwork.txt and Sampleinfo_network.txt) and the actual cytoscape file with raw network details (18S_network.cys). Network information (edges and nodes) can be exported as .csv files and sorted to assess number of negative/positive edges between nodes (ASVs). 


## See Supplemental Tables 3 and 4 for network information (edges and nodes) that can be filtered and used to create circos plots
## Import network information from CoNet into R to make circos plots 
library("circlize")


# Positive edges only for total 18S network - read in .csv file that is constructed from the network - Supplemental Figure 3a
pos_edges <- read.csv("Positive_edges.csv", header=T, row.names = 1, check.names=F) # Number of positive edges between groups from the overall 18S network
pos_mat <- data.matrix(pos_edges, rownames.force = NA)
grid.col = c(ConThree = "#0C9899", Other = "#777578", Abiotic = "#CFCACB",
             Prymn = "#EFE326", Mamiell = "#2E749F", MAST = "#E76896", Crypto = "#D94143", Diatoms = "#1A793F", Dino = "#92B076", Ciliate = "#A0CBED", Syndiniales = "#772B8C") # Labels and colors for your nodes
tiff("18S_positive.tiff", width = 8, height = 8, units = 'in', res = 600) # Useful to export the large network file 
chordDiagram(pos_mat, annotationTrack="grid", preAllocateTracks = 1, grid.col=grid.col, transparency = 0.2, order = c("ConThree", "Abiotic", "Prymn", "Mamiell", "MAST", "Crypto","Ciliate","Syndiniales" , "Dino", "Other","Diatoms")) # Main function to control grids, transparency of the connections, label positions around the circle, etc. 
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) { # Orients the labels on the outside of the network for better visualization
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off() # Figure file will show up in working directory folder

# Negative edges only for total 18S network - read in .csv file constructed from the network - Supplemental Figure 3b
neg_edges <- read.csv("Negative_edges.csv", header=T, row.names = 1, check.names=F) # Number of negative edges between groups from the overall 18S network
neg_mat <- data.matrix(neg_edges, rownames.force = NA)
grid.col = c(ConThree = "#0C9899", Other = "#777578", Abiotic = "#CFCACB",
             Prymn = "#EFE326", Mamiell = "#2E749F", MAST = "#E76896", Crypto = "#D94143", Diatoms = "#1A793F", Dino = "#92B076", Ciliate = "#A0CBED", Syndiniales = "#772B8C")
tiff("18S_negative.tiff", width = 8, height = 8, units = 'in', res = 600)
chordDiagram(neg_mat, annotationTrack="grid", preAllocateTracks = 1, grid.col=grid.col, transparency = 0.2, order = c("ConThree", "Abiotic", "Prymn", "Mamiell", "MAST", "Crypto","Ciliate","Syndiniales" , "Dino", "Other","Diatoms"))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off() # Same steps, just with negative edges from the 18S network


# Negative and positive edges for Syndiniales ASVs only - Figure 3
all_edges <- read.csv("Syndiniales_edges.csv", header=T, row.names = 1, check.names=F) # Number of positive (first column) and negative (second column) edges between Syndiniales and protist groups
all_mat <- data.matrix(all_edges, rownames.force = NA)
rownames(all_mat)=rownames(all_mat)
colnames(all_mat)=colnames(all_mat)
df = expand.grid(rownames(all_mat)[1:11], colnames(all_mat)[1:1])
df1 = df
df1$value = c(15,21,6,31,16,0,2,1,0,4,34) # Make sure the positive edge values are in the correct order
df2 = df
df2$value = c(0,-64,-54,-14,-29,-23,-14,-10,-10,-2,-3) # Correct the order of the data; to denote negative from positive edges, make the negative edge values negative
df = rbind(df1, df2) # Recombine these two accurate data frames with positive and negative values
grid.col = structure(1:13, names = c(rownames(all_mat)[1:11], colnames(all_mat)[1:2])) 
grid.col = c(ConThree = "#0C9899", Other = "#777578", Abiotic = "#CFCACB", Prymn = "#EFE326", Mamiell = "#2E749F", MAST = "#E76896", Crypto = "#D94143", Diatoms = "#1A793F", Dino = "#92B076", Ciliate = "#A0CBED", Syndiniales = "#772B8C")
tiff("Syndiniales.tiff", width = 8, height = 8, units = 'in', res = 600)
chordDiagram(df, col = ifelse(df$value > 0, "#003E57", "#DD6624"), annotationTrack="grid", preAllocateTracks = 1, grid.col = grid.col, transparency = 0.2, link.rank=rank(df[[3]])) # If the value is positive, it will assign it as blue; negative will be orange
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

## Export data for PLS models - need average relative abundance over time for the major class level groups (PLS analysis will yield VIP values used in Supplemental Table 1)
mergedps = merge_samples(ps, group = "Date", fun = mean) # Merge samples based on date
ps1  <- transform_sample_counts(mergedps, function(x) x / sum(x) ) # Relative abundance 
ps_pls <- subset_taxa(ps1, Class=="Filosa-Thecofilosea" | Class== "Cryptophyceae" | Class=="Spirotrichea" | Class=="Mamiellophyceae" | Class=="Dinophyceae" | Class=="Bacillariophyta" | Class=="Syndiniales") # Include only the major groups at class level
x1_pls <- tax_glom(ps_pls, taxrank = 'Class') # Agglomerate to class
OTU1 = as(otu_table(x1_pls), "matrix") # Convert ASV table to matrix
if(taxa_are_rows(x1_pls)){OTU1 <- t(OTU1)} # Switch rows and columns
OTUdf.ps_rel = as.data.frame(OTU1) # Convert to data frame
colnames(OTUdf.ps_rel)[c(1:7)]<-c("Mamiellophyceae", "Syndiniales", "Dinophyceae", "Bacillariophyta", "Cryptophyceae", "Spirotrichea", "Filosa-Thecofilosea") # Label rows and columns
rownames(OTUdf.ps_rel)[c(1:33)] <- c(31617, 32217, 32917, 41217, 42017, 50317,60717, 61517, 62217, 70617, 72117, 72717, 80317, 81017, 81617, 82317, 83017, 90617, 92017, 92617, 100417, 101117, 101817, 110817, 111617, 112117, 120717, 121417, 11618, 12418, 20118, 20818, 22118)
write.csv(OTUdf.ps_rel, file="18S_pls.csv", quote=FALSE, row.names = FALSE) # This file will be uploaded for pls analysis, get rid of date column by rows.names = False

# Run pls analysis with relative abundances
library (pls)
library (caret)
library (plsVarSel)

# Read in .csv file of environmental factors and change to matrix format; log-transform data
factors <- read.csv('factors_pls.csv', header = TRUE)
factorslog <- log1p(factors)

# Read in 18S relative abundances, which correspond to the factors and change to correct format
all_euks <- read.csv('18S_pls.csv', header = TRUE)
all_euks = as.matrix(all_euks)
Syn_only <- all_euks[, c(2)]
Syn_only = as.matrix(Syn_only)
Dino_only <- all_euks[, c(3)]
Dino_only = as.matrix(Dino_only)
Mamiell_only <- all_euks[, c(1)]
Mamiell_only = as.matrix(Mamiell_only)
Diatom_only <- all_euks[, c(4)]
Diatom_only = as.matrix(Diatom_only)
Crypto_only <- all_euks[, c(5)]
Crypto_only = as.matrix(Crypto_only)
Ciliate_only <- all_euks[, c(6)]
Ciliate_only = as.matrix(Ciliate_only)
Filosa_only <- all_euks[, c(7)]
Filosa_only <- as.matrix(Filosa_only)

# Run PLSR model to identify appropriate # of components needed in the model - name model (plsr_model); can also use validation = LOO (Leave on out); add na.action = na.omit (omit na values in model - good for missing data); can also add jackknife = TRUE
require(pls)
set.seed (10)
plsr_model <- plsr(Syn_only~., na.action = na.omit, data = factorslog, scale = TRUE, validation = "LOO", jackknife = TRUE)
summary (plsr_model) 


# Validate plot base on MSEP or variance to identify # of components (take the fist local minimum to avoid over-fitting)
validationplot(plsr_model, val.type = "RMSEP")
plot(RMSEP(plsr_model), legendpos = "topright")

# Calculate variable importance of the projection (VIP), which tells you how important each variable is to the model projection that involves protist group abundance (VIP > 1, means variables were very important to model)
vip <- VIP(plsr_model, opt.comp = 3) # Specify optimum components for VIP (all groups were computed with a 3-component model)
barplot(vip, ylab = "Variable Importance of Projection (VIP)", axes=TRUE)
summary (vip)
vip # This will give you values for each variable in the model

## Re-run the above plsr function using each protist group (e.g. Ciliate_only, Filosa_only, etc.)

## Export transformed otu count table for ASV-ASV pairing graphs (Figure 5 and Supplemental Figures 4 and 5); will give you relative abundance of all ASVs in the filtered dataset
ASV_abundance_rel = as(otu_table(ps_rel), "matrix") # Use the original phyloseq object transformed to relative abundance (ps_rel)
if(taxa_are_rows(ps_rel)){ASV_abundance_rel <- t(ASV_abundance_rel)}
ASV_abundance_rel = as.data.frame(ASV_abundance_rel)
write.csv(ASV_abundance_rel, file="ASV_abundance_rel.csv", quote=FALSE, row.names=TRUE)

