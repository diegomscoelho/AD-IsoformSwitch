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
tx2gene <- read.table("~/diegocoelho/refs/tr2g.All.txt", header = T, stringsAsFactors = F)

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
areas <- levels(metadata$Area)[-c(1:2)]

# PART 2

ISAR.tabs <- lapply(areas, function(area){

BM <- readRDS(paste0("braakMSBB.",area,".aSwitchListAnalyzed.rds"))

BM <- extractSequence(BM, writeToFile = F)

BM <- isoformSwitchAnalysisPart2(BM,
                                 pathToCPC2resultFile = paste0("MSBB_fasta/braakMSBB.",area,"_cpc2.txt"),
                                 pathToNetSurfP2resultFile = paste0("MSBB_fasta/braakMSBB.",area,"_netsurfp2.csv"),
                                 pathToSignalPresultFile = paste0("MSBB_fasta/braakMSBB.",area,"_summary.signalp5"),
                                 pathToPFAMresultFile = paste0("MSBB_fasta/braakMSBB.",area,"_pfam.txt"),
                                 dIFcutoff = 0.05,
                                 removeNoncodinORFs = T)

saveRDS(BM, file = paste0("braakMSBB.",area,".fullAnalysis.rds"))

})


