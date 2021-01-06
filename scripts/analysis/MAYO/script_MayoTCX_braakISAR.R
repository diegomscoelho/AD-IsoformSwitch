#######################
# Necessary libraries #
#######################

library(IsoformSwitchAnalyzeR)
library(dplyr)

## Import metadata
# https://www.synapse.org/#!Synapse:syn3817650 - MayoRNAseq_RNAseq_TCX_covariates.csv
MayoTCX_metadata <- read.csv("../data/MayoRNAseq_RNAseq_TCX_covariates.csv", stringsAsFactors=F, header = T)
myDesign <- MayoTCX_metadata %>% select(SampleID, Diagnosis, Braak) %>% filter(Diagnosis %in% c("Control", "AD")) %>% arrange(SampleID)
rownames(myDesign) = myDesign$SampleID
myDesign$Diagnosis <- as.factor(myDesign$Diagnosis)
myDesign = myDesign %>% na.omit()
myDesign$condition <- as.factor(plyr::mapvalues(myDesign$Braak, c(0,1,2,2.5,3,4.5,5,5.5,6),c(rep("braak_low",4),rep("braak_mid",2),rep("braak_high",3))))
myDesign = myDesign %>% select(SampleID, condition)
colnames(myDesign)[1] = "sampleID" # Rename as library asks for

kallisto_path = "PATH-TO/reprocessing/kallisto/MAYO_TCX/" # Set path
### Import Data
files = list.files(path = kallisto_path, full.names = T, recursive = T) %>%
  grep(pattern = ".tsv", value = T) %>% grep(pattern = paste(paste0("/",myDesign$sampleID,"/"), collapse = "|"), value = T)
## Import kallisto files
kallistoQuant <- importIsoformExpression(sampleVector = files)
colnames(kallistoQuant$counts)[-1] = myDesign$sampleID
colnames(kallistoQuant$abundance)[-1] = myDesign$sampleID

### Create switchAnalyzeRlist
# Reference genome: wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# Annotation: wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf.gz 

aSwitchList <- importRdata(
   isoformCountMatrix   = kallistoQuant$counts,
   isoformRepExpression = kallistoQuant$abundance,
   designMatrix         = myDesign,
   isoformExonAnnoation = "PATH-TO/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf.gz",
   isoformNtFasta       = "PATH-TO/diegocoelho/refs/Homo_sapiens.GRCh38.cdna.all.release-94.fa.gz",
   showProgress = T
)
aSwitchList

# Save aSwitchList
saveRDS(aSwitchList, "braakMayo.TCX.aSwitchList.rds")

# aSwitchList <- readRDS("braakMayo.TCX.aSwitchList.rds")

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

saveRDS(aSwitchListAnalyzed, "braakMayo.TCX.aSwitchListAnalyzed.rds")

extractSequence(preFilter(aSwitchListAnalyzed, IFcutoff=0.05, reduceToSwitchingGenes=T), pathToOutput="TCX_fasta/", outputPrefix="Mayo.TCX")

# TCX <- readRDS("braakMayo.TCX.aSwitchListAnalyzed.rds")

TCX <- extractSequence(TCX, writeToFile = F, dIFcutoff = 0.05)

TCX <- isoformSwitchAnalysisPart2(TCX,
                                  pathToCPC2resultFile = "TCX_fasta/braakMayo.TCX_cpc2.txt",
                                  pathToNetSurfP2resultFile = "TCX_fasta/braakMayo.TCX_AA.csv",
                                  pathToSignalPresultFile = "TCX_fasta/braakMayo.TCX_summary.signalp5",
                                  pathToPFAMresultFile = "TCX_fasta/braakMayo.TCX_pfam.txt",
                                  dIFcutoff = 0.05,
                                  removeNoncodinORFs = T)
saveRDS(TCX, "braakMayo.TCX.fullAnalysis.rds")



