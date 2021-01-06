#######################
# Necessary libraries #
#######################

library(IsoformSwitchAnalyzeR)
library(dplyr)
### Import Data

kallisto_path = "PATH-TO/reprocessing/kallisto/ROSMAP/" # Set path
## Import kallisto files
kallistoQuant <- importIsoformExpression(kallisto_path)

## Import metadata
# https://www.synapse.org/#!Synapse:syn3382527 - ROSMAP_IDkey.csv
ROSMap_metadata <- read.csv("../refs/ROSMAP_IDkey.csv", stringsAsFactors=T, header = T)
# https://www.synapse.org/#!Synapse:syn3191087 - ROSMAP_Clinical.csv
ROSMap_clinical <- read.csv("../refs/ROSMAP_Clinical.csv", header = T)

All_metadata <- merge(ROSMap_metadata, ROSMap_clinical, by = "projid")

rename = stringr::str_split_fixed(colnames(kallistoQuant$counts)[-1], "_",3)[,c(1,2)]
rename = paste(rename[,1],rename[,2],sep="_")

# Rename columns
colnames(kallistoQuant$counts)[-1] <- rename
colnames(kallistoQuant$abundance)[-1] <- rename

All_metadata <- All_metadata %>% select(projid, rnaseq_id, braaksc, ceradsc, cogdx) %>% filter(rnaseq_id %in% colnames(kallistoQuant$counts)) %>% filter(cogdx %in% c(1,4,5)) %>% distinct()

keep <- c(T,colnames(kallistoQuant$counts)[-1] %in% All_metadata$rnaseq_id)

# Make sure only annotated data is avaiable
kallistoQuant$counts <- kallistoQuant$counts[,keep]
kallistoQuant$abundance <- kallistoQuant$abundance[,keep]

## Create design matrix
myDesign <- All_metadata[,c(2,5)]
colnames(myDesign) = c("sampleID","condition")
rownames(myDesign) = myDesign$sampleID
#myDesign <- myDesign[colnames(kallistoQuant$counts)[-1],]
#myDesign <- myDesign[myDesign$sampleID %in% colnames(kallistoQuant$abundance),]

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
saveRDS(aSwitchList, "ROSMap.aSwitchList.rds")

# aSwitchList <- readRDS("Mayo.CBE.aSwitchList.rds")

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

saveRDS(aSwitchListAnalyzed, "ROSMap.aSwitchListAnalyzed.rds")

# Extract Sequences
extractSequence(aSwitchListAnalyzed, pathToOutput="ROSMap_fasta/", outputPrefix="ROSMap", dIFcutoff=0.05)

ROSMap <- readRDS("ROSMap.aSwitchListAnalyzed.rds")

ROSMap <- extractSequence(ROSMap, writeToFile = F, dIFcutoff = 0.05)

ROSMap <- isoformSwitchAnalysisPart2(ROSMap,
                                   pathToCPC2resultFile = "ROSMap_fasta/ROSMap_cpc2.txt",
                                   pathToNetSurfP2resultFile = "ROSMap_fasta/ROSMap_AA.csv",
                                   pathToSignalPresultFile = "ROSMap_fasta/ROSMap_summary.signalp5",
                                   pathToPFAMresultFile = "ROSMap_fasta/ROSMap_AA_pfam.txt",
                                   dIFcutoff = 0.05,
                                   removeNoncodinORFs = T)
saveRDS(ROSMap, "ROSMap.fullAnalysis.rds")




