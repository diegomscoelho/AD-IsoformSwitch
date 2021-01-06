#######################
# Necessary libraries #
#######################

library(DESeq2)
library(tximport)
library(dplyr)

### Import Data
kallisto_path = "PATH-TO/reprocessing/kallisto/MAYO_TCX/" # Set path
## Import metadata
# https://www.synapse.org/#!Synapse:syn3817650 - MayoRNAseq_RNAseq_TCX_covariates.csv
MayoTCX_metadata <- read.csv("../refs/MayoRNAseq_RNAseq_TCX_covariates.csv", stringsAsFactors=F, header = T)
myDesign <- MayoTCX_metadata %>% select(SampleID, Diagnosis, Braak, RIN, Sex) %>% filter(Diagnosis %in% c("Control", "AD")) %>% arrange(SampleID)
rownames(myDesign) = myDesign$SampleID
myDesign$Diagnosis <- as.factor(myDesign$Diagnosis)
myDesign = myDesign %>% na.omit()
myDesign$condition <- as.factor(plyr::mapvalues(myDesign$Braak, c(0,1,2,2.5,3,4.5,5,5.5,6),c(rep("braak_low",4),rep("braak_mid",2),rep("braak_high",3))))

## Import kallisto files
files = paste0(list.files(kallisto_path, full.names=T), "/abundance.tsv")
files = grep(paste0(paste0("/",myDesign$SampleID), collapse="|"), files, value=T)
tx2gene <- read.table("../refs/tr2g.All.txt.gz", header = T, stringsAsFactors = F)
colnames(tx2gene) <- c("TXNAME", "GENEID", "GENE_NAME")
kallistoQuant <- tximport(files, type = "kallisto", tx2gene = tx2gene[,c(1,2)])
kallistoQuant.tx <- tximport(files, type = "kallisto", txOut = T)

# Create gene-level DESeq2 object
dds <- DESeqDataSetFromTximport(kallistoQuant, myDesign, ~ Sex + RIN + condition)
dds$condition <- relevel(dds$condition, ref = "braak_low") # Reference
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
dres1 <- results(dds, contrast=c("condition","braak_mid","braak_low")) %>% as.data.frame() %>% tibble::rownames_to_column("GENEID") %>% select(GENEID,gene.log2FC=log2FoldChange,gene.padj=padj) %>% mutate(condition="low_mid")
dres2 = results(dds, contrast=c("condition","braak_high","braak_low")) %>% as.data.frame() %>% tibble::rownames_to_column("GENEID") %>% select(GENEID,gene.log2FC=log2FoldChange,gene.padj=padj) %>% mutate(condition="low_high")
dres3 = results(dds, contrast=c("condition","braak_high","braak_mid")) %>% as.data.frame() %>% tibble::rownames_to_column("GENEID") %>% select(GENEID,gene.log2FC=log2FoldChange,gene.padj=padj) %>% mutate(condition="mid_high")

dres <- Reduce(rbind, list(dres1,dres2,dres3))

# Create transcript-level DESeq2 object
#dds.tx <- DESeqDataSetFromTximport(kallistoQuant.tx, myDesign, ~ Sex + RIN + condition)
#dds.tx$condition <- relevel(dds.tx$condition, ref = "braak_low") # Reference
#
#keep.tx <- rowSums(counts(dds.tx)) >= 3
#dds.tx <- dds.tx[keep.tx,]
#dds.tx <- DESeq(dds.tx)
#
#dres.tx1 <- results(dds.tx, contrast=c("condition","braak_mid","braak_low")) %>% as.data.frame() %>% tibble::rownames_to_column("TXNAME") %>% select(TXNAME,transcript.log2FC=log2FoldChange,transcript.padj=padj) %>% mutate(condition="low_mid")
#dres.tx2 = results(dds.tx, contrast=c("condition","braak_high","braak_low")) %>% as.data.frame() %>% tibble::rownames_to_column("TXNAME") %>% select(TXNAME,transcript.log2FC=log2FoldChange,transcript.padj=padj) %>% mutate(condition="low_high")
#dres.tx3 = results(dds.tx, contrast=c("condition","braak_high","braak_mid")) %>% as.data.frame() %>% tibble::rownames_to_column("TXNAME") %>% select(TXNAME,transcript.log2FC=log2FoldChange,transcript.padj=padj) %>% mutate(condition="mid_high")
#
#dres.tx <- Reduce(rbind, list(dres.tx1,dres.tx2,dres.tx3))
#
dres_total <- merge(dres, tx2gene[,2:3]) %>% distinct()
#dres_total <- merge(dres.tx.gene, dres, by = c("GENEID","condition"), all = T)
#dres_total = dres_total[,c(3,1,6,2,7,8,4,5)]
#dres_total <- dres_total %>% arrange(GENE_NAME, condition)

saveRDS(dres_total, file = "braakMayo.TCX.DEseq2.rds")

