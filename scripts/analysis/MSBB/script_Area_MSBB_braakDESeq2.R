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
metadata <- metadata %>% select(SampleID = sampleIdentifier, Diagnosis = NP.1, condition = bbscore, Area = BrodmannArea, CDR, RIN, Sex) %>%
  filter(Diagnosis %in% c(1,2)) %>% arrange(SampleID) %>% distinct() %>% na.omit()

# Change values to nomeclature
metadata$Diagnosis <- plyr::mapvalues(metadata$Diagnosis, c(1,2,3,4), c("Control","AD","prAD","poAD"))
metadata$condition2 <- as.factor(plyr::mapvalues(metadata$condition, c(0,1,2,3,4,5,6),c(rep("braak_low",3),rep("braak_mid",2),rep("braak_high",2)))) # Transform braak scores in 3 different levels

metadata$SampleID <- as.character(metadata$SampleID)

# Select Areas
areas <- levels(metadata$Area)

dds.tabs <- lapply(areas, function(area){

myDesign = metadata %>% filter(Area %in% area) %>% distinct()

# Files to choose
samples <- grep(files, pattern= paste(paste0("/",myDesign$SampleID,"/"), collapse="|"), value = T)
valid <- stringr::str_split_fixed(samples, "/", 5)[,4]

myDesign <- myDesign %>% filter(SampleID %in% valid)
rownames(myDesign) <- myDesign$SampleID
myDesign$condition <- as.factor(myDesign$condition)

# Import kallisto counts gene-leve through tximport
kallistoQuant <- tximport(samples, type = "kallisto", tx2gene = tx2gene[,c(1,2)])
kallistoQuant$braak = myDesign %>% pull(condition)
kallistoQuant$CDR = myDesign %>% pull(CDR)
kallistoQuant$Diagnosis = myDesign %>% pull(Diagnosis)
kallistoQuant$SampleID = myDesign %>% pull(SampleID)  

saveRDS(kallistoQuant, file = paste0("MSBB.",area,".kallisto.tximport.rds"))

myDesign = myDesign %>% select(SampleID, condition2)

# Gene-level DESeq2 analysis
dds <- DESeqDataSetFromTximport(kallistoQuant, myDesign, , ~ Sex + RIN + condition2)
dds$condition <- relevel(dds$condition, ref = "braak_low") # Reference

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

dres1 <- results(dds, contrast=c("condition","braak_mid","braak_low")) %>% as.data.frame() %>% tibble::rownames_to_column("GENEID") %>% select(GENEID,gene.log2FC=log2FoldChange,gene.padj=padj) %>% mutate(condition="low_mid")
dres2 = results(dds, contrast=c("condition","braak_high","braak_low")) %>% as.data.frame() %>% tibble::rownames_to_column("GENEID") %>% select(GENEID,gene.log2FC=log2FoldChange,gene.padj=padj) %>% mutate(condition="low_high")
dres3 = results(dds, contrast=c("condition","braak_high","braak_mid")) %>% as.data.frame() %>% tibble::rownames_to_column("GENEID") %>% select(GENEID,gene.log2FC=log2FoldChange,gene.padj=padj) %>% mutate(condition="mid_high")

dres <- Reduce(rbind, list(dres1,dres2,dres3))

# Rename tx2gene
tx2 = tx2gene[,2:3]
colnames(tx2) = c("GENEID", "GENE_NAME")
dres_total <- merge(dres, tx2, by = "GENEID")
dres_total$Area <- area

print(paste("End of",area,"analysis"))

return(dres_total)

})

merged.dds <- Reduce(rbind, dds.tabs)

saveRDS(merged.dds, file = "braakMSBB.DEseq2.rds")

