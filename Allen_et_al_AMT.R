# Load required packages
library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(reshape2)
library(gridExtra)
library(zoo)
library(tibble)


# Set working directory and path to raw sequences
path <- ("~/")
list.files(path)


# Extract sample names (used in filtering and trimming)
f.names = as.vector(list.files(path, pattern = "_R1_001.fastq.gz", 
                               full.names = F))
r.names = as.vector(list.files(path, pattern = "_R2_001.fastq.gz", 
                               full.names = F))


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# Plot and save quality profiles
qpf <- plotQualityProfile(fnFs[1:2])
qpr <- plotQualityProfile(fnRs[1:2])
ggsave("amt_dna_rna_quality_F.jpeg", qpf, width = 60, height = 40, units = "cm", device = "jpeg")
ggsave("amt_dna_rna_quality_R.jpeg", qpr, width = 60, height = 40, units = "cm", device = "jpeg")


# Filter and trim reads (quality and primer removal)
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(f.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(r.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,210), trimLeft=c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


# Calculate and plot the error model
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plot.errF <- plotErrors(errF, nominalQ = TRUE) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot.errR <- plotErrors(errR, nominalQ = TRUE) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("amt_dna_rna_error_F.jpeg", plot.errF, width = 60, height = 40, units = "cm", device = "jpeg")
ggsave("amt_dna_rna_error_R.jpeg", plot.errR, width = 60, height = 40, units = "cm", device = "jpeg")


# Dereplicate sequences for computational efficiency
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- f.names
names(derepRs) <- r.names


# Run DADA2 algorithm on dereplicated sequences to resolve unique sequences at single-nucleotide resolution
dadaFs <- dada(derepFs, err <- errF, multithread = TRUE)
dadaRs <- dada(derepRs, err <- errR, multithread = TRUE)
dadaFs[[1]]
dadaRs[[1]]


# Merge paired-end reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)


# Generate table of sequences  
seqtab <- makeSequenceTable(mergers)


# Plot sequence length distribution 
ex <- as.data.frame(table(nchar(getSequences(seqtab))))
ex.plot <- ggplot(ex, aes(Var1, Freq)) +
  geom_bar(stat = "identity", aes(fill = Var1)) +
  ylab("Frequency") +
  xlab("Merged Sequence Length") +
  scale_x_discrete() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")
ggsave("amt_dna_rna_length_plot.jpeg", ex.plot, width = 15, height = 7.5, units = "cm", device = "jpeg")


# Remove chimeric sequences 
seqtab.ex.chi <- removeBimeraDenovo(seqtab, method = "consensus", verbose = T)


# Track sequence loss throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), 
               rowSums(seqtab), rowSums(seqtab.ex.chi))
colnames(track) <- c("input", "filtered", "denoised", "merged", 
                     "tabled", "no chim")
rownames(track) <- f.names
head(track)


# Assign taxonomy with minimum bootstrap confidence of 80% (SILVA v132)
silva.taxa <- assignTaxonomy(seqtab.ex.chi, "silva_nr_v132_train_set.fa.gz", minBoot = 80, outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)


# Assign species level taxonomy (SILVA v132) 
silva.taxa.species <- addSpecies(silva.taxa$tax, "silva_species_assignment_v132.fa.gz")


# Upload associated metadata (available in github repository)
amt.meta <- read.csv("amt_dna_rna_mapping.csv", header = T)


# Setup sample names to correspond with sequence data
sample.names
amt.meta$sample_id
amt.meta$seq_id <- substring(amt.meta$sample_id, 3)


# Apply sample names to metadata rows (for downstream phyloseq recognition)
row.names(amt.meta) <- amt.meta$seq_id
amt.meta <- as.data.frame(amt.meta)


# Apply sample names to sequence table
row.names(seqtab.ex.chi) = sample.names


# Construct integrated phyloseq object
amt.phy <- phyloseq(tax_table(silva.taxa.species), otu_table(seqtab.ex.chi, taxa_are_rows = FALSE), sample_data(amt.meta))


# Create vector of unique ASV names
dim(otu_table(amt.phy))
dim(tax_table(amt.phy))
a.vec <- as.vector(1:11047) # Value reflects total unique ASVs
a.nam <- cbind("asv_", a.vec)
a.nam <- as.data.frame(a.nam)
asv.names <- paste0(a.nam$V1, a.nam$a.vec)
asv.names <- as.data.frame(asv.names)
taxa_names(amt.phy) <- asv.names$asv.names


# Restructure tax_table to append column with most detailed taxonomic classification concatenated with uniquie ASV number
bc.t <- t(as.data.frame(tax_table(amt.phy)))
bc.t[bc.t==""] <- NA
bc.fill <- na.locf(bc.t, na.rm = TRUE)
t.bc.fill <- as.data.frame(t(bc.fill))
head(t.bc.fill)
rnc.bc <- rownames_to_column(t.bc.fill, "ASV")
rnc.bc$taxa_ASV <- paste(rnc.bc$Species,rnc.bc$ASV)
safe.bc <- as.data.frame(tax_table(amt.phy))
safe.bc$taxa_ASV <- paste(rnc.bc$taxa_ASV)
bc.tax <- tax_table(safe.bc)
colnames(bc.tax) <- colnames(safe.bc)
rownames(bc.tax) <- rownames(safe.bc)


# Check all other columns remain identical and update phyloseq object with the updated tax_table 
identical(bc.tax[1:11047,1:7], tax_table(amt.phy)) # Should be true
tax_table(amt.phy) <- bc.tax
View(as.data.frame(tax_table(amt.phy)))


# Remove sequences which are not assigned as Kingdom Bacteria, then remove sequences assigned as chloroplasts or mitochondria
amt.phy.x <- split_proteo_silva(amt.phy)
amt.phy.x <- subset_taxa(amt.phy.x, Kingdom == "Bacteria")
rowSums(otu_table(amt.phy.x))
amt.phy.x = subset_taxa(amt.phy.x, (Order!="Chloroplast") | is.na(Order))
dim(tax_table(amt.phy.x))
amt.phy.x = subset_taxa(amt.phy.x, (Family!="Mitochondria") | is.na(Family))
rowSums(otu_table(amt.phy.x))


# Rarefy sequences to an even depth across all samples. This results in the removal of Station 5 1% PAR gDNA and consequently Station 5 1% PAR cDNA. Samples unrelated to the transect are also removed.
set.seed(711)
amt.phy.rare <- rarefy_even_depth(amt.phy.x, sample.size = 15973, trimOTUs = TRUE) 
amt.phy.rare <- subset_samples(amt.phy.rare, seq_id != 67)
amt.phy.rare <- subset_samples(amt.phy.rare, station != "E1") 

  

