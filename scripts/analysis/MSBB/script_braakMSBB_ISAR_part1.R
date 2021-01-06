#######################
# Necessary libraries #
#######################

library(IsoformSwitchAnalyzeR)
library(dplyr)

### Import Data
kallisto_path = "PATH-TO/reprocessing/kallisto/MSBB/" # Set path
# Import data files
files = paste0(list.files(kallisto_path, full.names=T),"/abundance.tsv")

# Reference
tx2gene <- read.table("../refs/tr2g.All.txt", header = T, stringsAsFactors = F)


## Import metadata
# https://www.synapse.org/#!Synapse:syn6100548 - MSBB_RNAseq_covariates.csv
MSBB_metadata <- read.csv("../refs/MSBB_RNAseq_covariates.csv", header = T)
# https://www.synapse.org/#!Synapse:syn6101474 - MSBB_individual_metadata.csv (deprecated)
MSBB_clinical <- read.csv("../refs/MSBB_clinical.csv.gz", header = T)
metadata <- merge(MSBB_metadata, MSBB_clinical, by="individualIdentifier")
metadata <- metadata %>% select(SampleID = sampleIdentifier, Diagnosis = NP.1, condition = bbscore, Area = BrodmannArea) %>% filter(Diagnosis %in% c(1,2)) %>% arrange(SampleID) %>% distinct() %>% na.omit()

# Change values to nomeclature
metadata$Diagnosis <- plyr::mapvalues(metadata$Diagnosis, c(1,2,3,4), c("Control","AD","prAD","poAD"))
metadata$condition <- as.factor(plyr::mapvalues(metadata$condition, c(0,1,2,3,4,5,6),c(rep("braak_low",3),rep("braak_mid",2),rep("braak_high",2)))) # Transform braak scores in 3 different levels

metadata$SampleID <- as.character(metadata$SampleID)

# Select Areas
areas <- levels(metadata$Area)

# PART1

ISAR.tabs <- lapply(areas, function(area){

myDesign = metadata %>% filter(Area %in% area) %>% select(SampleID, condition) %>% distinct()

# Files to choose
samples <- grep(files, pattern= paste(paste0("/",myDesign$SampleID,"/"), collapse="|"), value = T)
valid <- stringr::str_split_fixed(samples, "/", 5)[,4]

kallistoQuant <- importIsoformExpression(sampleVector=samples)

myDesign <- myDesign %>% filter(SampleID %in% valid)
rownames(myDesign) <- myDesign$SampleID
myDesign$condition <- as.factor(myDesign$condition)

colnames(kallistoQuant$counts)[-1] = myDesign$SampleID
colnames(kallistoQuant$abundance)[-1] = myDesign$SampleID
colnames(myDesign) = c("sampleID", "condition")

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
saveRDS(aSwitchList, paste0("braakMSBB.",area,".aSwitchList.rds"))

})

ISAR.tabs <- lapply(areas, function(area){

aSwitchList <- readRDS(paste0("braakMSBB.",area,".aSwitchList.rds"))

aSwitchListFiltered <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 10,
    isoformExpressionCutoff = 3,
    removeSingleIsoformGenes = TRUE
)

aSwitchListAnalyzed <- isoformSwitchTestDRIMSeq(
    switchAnalyzeRlist = aSwitchListFiltered,
    reduceFurtherToGenesWithConsequencePotential=FALSE,
    dIFcutoff=0.05,
    reduceToSwitchingGenes=TRUE
)

saveRDS(aSwitchListAnalyzed, file=paste0("braakMSBB.",area,".aSwitchListAnalyzed.rds"))

# Extract Sequences
extractSequence(aSwitchListAnalyzed, pathToOutput="MSBB_fasta/", outputPrefix=paste0("braakMSBB.",area))

})


