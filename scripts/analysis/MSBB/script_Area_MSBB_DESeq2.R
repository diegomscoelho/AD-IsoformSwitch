#######################
# Necessary libraries #
#######################

library(DESeq2)
library(dplyr)
library(tximport)

kallisto_path = "PATH-TO/reprocessing/kallisto/MSBB/" # Set path
# Import data files
files = paste0(list.files(kallisto_path, full.names=T),"/abundance.tsv")

# Reference
tx2gene <- read.table("../refs/tr2g.All.txt.gz", header = T, stringsAsFactors = F)


## Import metadata
# https://www.synapse.org/#!Synapse:syn6100548 - MSBB_RNAseq_covariates.csv
MSBB_metadata <- read.csv("../refs/MSBB_RNAseq_covariates.csv", header = T)
# https://www.synapse.org/#!Synapse:syn6101474 - MSBB_individual_metadata.csv (deprecated)
MSBB_clinical <- read.csv("../refs/MSBB_clinical.csv.gz", header = T)
metadata <- merge(MSBB_metadata, MSBB_clinical, by="individualIdentifier")
metadata <- metadata %>% select(SampleID = sampleIdentifier, Diagnosis = NP.1, Area = BrodmannArea, RIN, Sex) %>%
  filter(Diagnosis %in% c(1,2)) %>% arrange(SampleID) %>% distinct()

# Change values to nomeclature
metadata$Diagnosis <- plyr::mapvalues(metadata$Diagnosis, c(1,2,3,4), c("Control","AD","prAD","poAD"))
metadata$SampleID <- as.character(metadata$SampleID)

# Select Areas
areas <- levels(metadata$Area)

lapply(areas, function(area){

myDesign = metadata %>% filter(Area %in% area) %>% select(SampleID, Diagnosis) %>% distinct()

# Files to choose
samples <- grep(files, pattern= paste(paste0("/",myDesign$SampleID,"/"), collapse="|"), value = T)
valid <- stringr::str_split_fixed(samples, "/", 5)[,4]

myDesign <- myDesign %>% filter(SampleID %in% valid)
rownames(myDesign) <- myDesign$SampleID
myDesign$Diagnosis <- as.factor(myDesign$Diagnosis)

head(samples)
head(myDesign)

tail(samples)
tail(myDesign)

# Import kallisto counts gene-leve through tximport
kallistoQuant <- tximport(samples, type = "kallisto", tx2gene = tx2gene[,c(1,2)])
# Import kallisto counts transcript-level
kallistoQuant.tx <- tximport(samples, type = "kallisto", txOut = T)

# Gene-level DESeq2 analysis
dds <- DESeqDataSetFromTximport(kallistoQuant, myDesign, ~ Sex + RIN + Diagnosis)
dds$Diagnosis <- relevel(dds$Diagnosis, ref = "Control") # Reference

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
dres <- DESeq2::results(dds)
summary(subset(dres, padj < 0.01 & abs(log2FoldChange) > log2(1.3)))
# Save analysis
saveRDS(dres, file= paste0("MSBB.",area,".DESeq2.dres.rds"))
saveRDS(dds, file= paste0("MSBB.",area,".DESeq2.dds.rds"))

# Transcript-level DESeq2 analysis
dds.tx <- DESeqDataSetFromTximport(kallistoQuant.tx, myDesign,  ~ Sex + RIN + Diagnosis)
dds.tx$Diagnosis <- relevel(dds.tx$Diagnosis, ref = "Control") # Reference

keep.tx <- rowSums(counts(dds.tx)) >= 3
dds.tx <- dds.tx[keep.tx,]
dds.tx <- DESeq(dds.tx)
dres.tx <- DESeq2::results(dds.tx)
summary(subset(dres.tx, padj < 0.01 & abs(log2FoldChange) > log2(1.3)))
# Save analysis
saveRDS(dres.tx, file= paste0("MSBB.",area,".DESeq2.dres.tx.rds"))
saveRDS(dds, file= paste0("MSBB.",area,".DESeq2.dds.tx.rds"))

print(paste("End of",area,"analysis"))

})



