#######################
# Necessary libraries #
#######################

library(IsoformSwitchAnalyzeR)
library(dplyr)

### Import Data
kallisto_path = "PATH-TO/reprocessing/kallisto/MSBB/" # Set path
## Import kallisto files
kallistoQuant <- importIsoformExpression(kallisto_path)

## Import metadata
# https://www.synapse.org/#!Synapse:syn6100548 - MSBB_RNAseq_covariates.csv
MSBB_metadata <- read.csv("../refs/MSBB_RNAseq_covariates.csv", header = T)
# https://www.synapse.org/#!Synapse:syn6101474 - MSBB_individual_metadata.csv (deprecated)
MSBB_clinical <- read.csv("../data/MSBB_clinical.csv.gz", header = T)
metadata <- merge(MSBB_metadata, MSBB_clinical, by="individualIdentifier")
metadata <- metadata[metadata$sampleIdentifier %in% colnames(kallistoQuant$abundance),c(3,22,5)]
metadata <- metadata %>% distinct()
# Change values to nomeclature
metadata$NP.1 <- plyr::mapvalues(metadata$NP.1, c(1,2,3,4), c("Control","AD","prAD","poAD"))

## Create design matrix
myDesign <- metadata
colnames(myDesign) = c("sampleID","condition","Area")
rownames(myDesign) = myDesign$sampleID
myDesign <- myDesign[colnames(kallistoQuant$counts)[-1],]
myDesign <- myDesign[myDesign$sampleID %in% colnames(kallistoQuant$abundance),]

### Create switchAnalyzeRlist
# Reference genome: wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# Annotation: wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf.gz 

aSwitchList <- importRdata(
    isoformCountMatrix   = kallistoQuant$counts,
    isoformRepExpression = kallistoQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = "PATH-TO/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf.gz",
    isoformNtFasta       = "PATH-TO/Homo_sapiens.GRCh38.cdna.all.release-94.fa.gz",
    showProgress = T
)
aSwitchList

# Save aSwitchList
saveRDS(aSwitchList, "MSBB.aSwitchList.rds")

# aSwitchList <- readRDS("MSBB.aSwitchList.rds")

aSwitchListFiltered <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 5,
    isoformExpressionCutoff = 1,
    removeSingleIsoformGenes = TRUE
)

aSwitchListAnalyzed <- isoformSwitchTestDRIMSeq(
    switchAnalyzeRlist = aSwitchListFiltered,
    reduceFurtherToGenesWithConsequencePotential=FALSE,
    dIFcutoff=0.05,
    reduceToSwitchingGenes=TRUE
)

saveRDS(aSwitchListAnalyzed, "MSBB.aSwitchListAnalyzed.rds")

# Extract Sequences
extractSequence(aSwitchListAnalyzed, pathToOutput="MSBB_fasta/", outputPrefix="MSBB")




