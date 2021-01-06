#######################
# Necessary libraries #
#######################

library(DESeq2)
library(dplyr)
library(tximport)

### Import Data
kallisto_path = "PATH-TO/reprocessing/kallisto/ROSMAP/" # Set path
# Set directory with kallisto tsv files separated by samples
files = paste0(list.files(kallisto_path, full.names=T), "/abundance.tsv")
tx2gene <- read.table("../refs/tr2g.All.txt.gz", header = T, stringsAsFactors = F)

## Import metadata
# https://www.synapse.org/#!Synapse:syn3382527 - ROSMAP_IDkey.csv
ROSMap_metadata <- read.csv("../refs/ROSMAP_IDkey.csv", stringsAsFactors=F, header = T)
# https://www.synapse.org/#!Synapse:syn3191087 - ROSMAP_Clinical.csv
ROSMap_clinical <- read.csv("../refs/ROSMAP_Clinical.csv", stringsAsFactors=F, header = T)

All_metadata <- merge(ROSMap_metadata, ROSMap_clinical, by = "projid")
All_metadata <- All_metadata %>% select(SampleID = rnaseq_id, Diagnosis = cogdx, RIN, Sex) %>% filter(Diagnosis %in% c(1,4)) %>% filter(SampleID != "") %>% arrange(SampleID) %>% distinct()

# Rename
All_metadata$Diagnosis <- plyr::mapvalues(All_metadata$Diagnosis, c(1,4), c("Control", "AD"))
rownames(All_metadata) <- All_metadata$SampleID
All_metadata$Diagnosis <- as.factor(All_metadata$Diagnosis)

# Samples with AD or Control
files = grep(files, pattern = paste0(All_metadata$SampleID, collapse="|"), value = T)

# Tximport data
kallistoQuant <- tximport(files, type = "kallisto", tx2gene = tx2gene[,c(1,2)])
kallistoQuant.tx <- tximport(files, type = "kallisto", txOut= T)

# DESeq2 gene-level analysis
dds <- DESeqDataSetFromTximport(kallistoQuant, All_metadata, , ~ Sex + RIN + Diagnosis)
dds$Diagnosis <- relevel(dds$Diagnosis, ref = "Control") # Reference

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
dres <- DESeq2::results(dds)
summary(subset(dres, padj < 0.01 & abs(log2FoldChange) > log2(1.3)))
# Save files
saveRDS(dres, file= "ROSMAP.DESeq2.dres.rds")
saveRDS(dds, file= "ROSMAP.DESeq2.dds.rds")

# DESeq2 transcript-level analysis
dds.tx <- DESeqDataSetFromTximport(kallistoQuant.tx, All_metadata, ~ Sex + RIN + Diagnosis)
dds.tx$Diagnosis <- relevel(dds.tx$Diagnosis, ref = "Control") # Reference

keep.tx <- rowSums(counts(dds.tx)) >= 3
dds.tx <- dds.tx[keep.tx,]
dds.tx <- DESeq(dds.tx)
dres.tx <- DESeq2::results(dds.tx)
summary(subset(dres.tx, padj < 0.01 & abs(log2FoldChange) > log2(1.3)))
# Save files
saveRDS(dres.tx, file= "ROSMAP.DESeq2.dres.tx.rds")
saveRDS(dds, file= "ROSMAP.DESeq2.dds.tx.rds")





