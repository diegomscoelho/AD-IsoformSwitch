#######################
# Necessary libraries #
#######################

library(DESeq2)
library(tximport)
library(dplyr)

### Import Data
kallisto_path = "PATH-TO/reprocessing/kallisto/MAYO_TCX/" # Set path
## Import metadata
MayoTCX_metadata <- read.csv("../refs/MayoRNAseq_RNAseq_TCX_covariates.csv", stringsAsFactors=F, header = T)
myDesign <- MayoTCX_metadata %>% select(SampleID, Diagnosis, RIN, Sex) %>% filter(Diagnosis %in% c("Control", "AD")) %>% arrange(SampleID)
rownames(myDesign) = myDesign$SampleID
myDesign$Diagnosis <- as.factor(myDesign$Diagnosis)

## Import kallisto files
files = paste0(list.files(kallisto_path, full.names=T), "/abundance.tsv")
files = grep(paste0(paste0("/",myDesign$SampleID), collapse="|"), files, value=T)
tx2gene <- read.table("../refs/tr2g.All.txt.gz", header = T, stringsAsFactors = F)
colnames(tx2gene) <- c("TXNAME", "GENEID", "GENE_NAME")
kallistoQuant <- tximport(files, type = "kallisto", tx2gene = tx2gene[,c(1,2)])
kallistoQuant.tx <- tximport(files, type = "kallisto", txOut = T)

# Create gene-level DESeq2 object
dds <- DESeqDataSetFromTximport(kallistoQuant, myDesign, ~ Sex + RIN + Diagnosis)
dds$Diagnosis <- relevel(dds$Diagnosis, ref = "Control") # Reference
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
dres <- DESeq2::results(dds)

summary(subset(dres, padj < 0.01 & abs(log2FoldChange) > log2(1.3)))

saveRDS(dres, file= "Mayo.TCX.DESeq2.dres.rds")
saveRDS(dds, file= "Mayo.TCX.DESeq2.dds.rds")

# Create transcript-level DESeq2 object
dds.tx <- DESeqDataSetFromTximport(kallistoQuant.tx, myDesign, ~ Sex + RIN + Diagnosis)
dds.tx$Diagnosis <- relevel(dds.tx$Diagnosis, ref = "Control") # Reference

keep.tx <- rowSums(counts(dds.tx)) >= 3
dds.tx <- dds.tx[keep.tx,]
dds.tx <- DESeq(dds.tx)
dres.tx <- DESeq2::results(dds.tx)

summary(subset(dres.tx, padj < 0.01 & abs(log2FoldChange) > log2(1.3)))

saveRDS(dres.tx, file= "Mayo.TCX.DESeq2.dres.tx.rds")
saveRDS(dds.tx, file= "Mayo.TCX.DESeq2.dds.tx.rds")

