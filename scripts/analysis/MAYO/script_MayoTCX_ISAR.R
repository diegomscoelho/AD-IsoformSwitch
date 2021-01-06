#######################
# Necessary libraries #
#######################

library(IsoformSwitchAnalyzeR)

### Import Data
kallisto_path = "PATH-TO/reprocessing/kallisto/MAYO_TCX/" # Set path
## Import kallisto files
kallistoQuant <- importIsoformExpression(kallisto_path)

## Import metadata
# https://www.synapse.org/#!Synapse:syn3817650 - MayoRNAseq_RNAseq_TCX_covariates.csv
MayoTCX_metadata <- read.csv("../refs/MayoRNAseq_RNAseq_TCX_covariates.csv", stringsAsFactors=T, header = T)

## Create design matrix
myDesign <- MayoTCX_metadata[,c(1,5,6)]
colnames(myDesign) = c("sampleID","condition","Gender")
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
saveRDS(aSwitchList, "Mayo.TCX.aSwitchList.rds")

# aSwitchList <- readRDS("Mayo.TCX.aSwitchList.rds")

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

saveRDS(aSwitchListAnalyzed, "Mayo.TCX.aSwitchListAnalyzed.rds")

extractSequence(preFilter(aSwitchListAnalyzed, IFcutoff=0.05, reduceToSwitchingGenes=T), pathToOutput="TCX_fasta/", outputPrefix="Mayo.TCX")

# TCX <- readRDS("Mayo.TCX.aSwitchListAnalyzed.rds")

TCX <- extractSequence(TCX, writeToFile = F, dIFcutoff = 0.05)

TCX <- isoformSwitchAnalysisPart2(TCX,
                                  pathToCPC2resultFile = "TCX_fasta/Mayo.TCX_cpc2.txt",
                                  pathToNetSurfP2resultFile = "TCX_fasta/Mayo.TCX_AA.csv",
                                  pathToSignalPresultFile = "TCX_fasta/Mayo.TCX_summary.signalp5",
                                  pathToPFAMresultFile = "TCX_fasta/Mayo.TCX_pfam.txt",
                                  dIFcutoff = 0.05,
                                  removeNoncodinORFs = T)
saveRDS(TCX, "Mayo.TCX.fullAnalysis.rds")



