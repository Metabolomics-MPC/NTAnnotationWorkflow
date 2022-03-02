################################################
# Skript loading Slaw Input                    # 
# perform library search on fused MGF file     #
# Store Annotations in Signal Intensity matrix #
################################################
### Set variables => Change for your system  ###
mbank_dir <- "~/Nextcloud/Cloud/9_HGMU/Share/20210203_SF5_MassBankRecords/neg" #directory with mbank records to annotate
datamatrix <- "~/Nextcloud/Cloud/9_HGMU/Share/Slaw_test/Neg/xcms/datamatrix_opt.csv" #slaw file with the MS1 data matrix
fused_mgf <- "~/Nextcloud/Cloud/9_HGMU/Share/Slaw_test/Neg/xcms/fused_mgf_opt.mgf" #slaw fused mgf output
output <- "~/Nextcloud/Cloud/9_HGMU/Share/Slaw_test/Neg/xcms/datamatrix_opt_annotated.csv" #file with the annotated features
annotated_only <- TRUE # Shall only annotated features be reported?
dp_tresh <- 0.4 # Threshold for annotation
mz_ppm <- 5 #mass error in ppm for spectra matching
int_tresh <- 5 #intensity threshold for matching
################################################
### Load or install packages ###################
################################################
{
if (!require("MetaboAnnotation")) devtools::install_github("rformassspectrometry/MetaboAnnotation")
library(MetaboAnnotation)

if (!require("MsBackendMassbank")) devtools::install_github("rformassspectrometry/MsBackendMassbank")
library(MsBackendMassbank)

if (!require("MsBackendMgf")) devtools::install_github("rformassspectrometry/MsBackendMgf")
library(MsBackendMassbank)

if (!require("Spectra")) devtools::install_github("rformassspectrometry/Spectra")
library(Spectra)

}
################################################
### Workflow ###################################
################################################
#Load MassBank Records
fls <- list.files(path=mbank_dir, pattern=".txt$", full.names = T)
library <- Spectra(fls, source = MsBackendMassbank(), backend = MsBackendDataFrame())

#Load Fused MGF file
query <- Spectra(fused_mgf, source = MsBackendMgf(), backend = MsBackendDataFrame())

# remove fragments below x% of base peak intensity
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * (int_tresh / 100 )
}
query <- filterIntensity(query, intensity = low_int)

# normalize intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}

query <- addProcessing(query, norm_int)
library <- addProcessing(library, norm_int)

#Perform library search with the MetaboAnnotation package
prm <- MatchForwardReverseParam(ppm = mz_ppm, #requirePrecursor = TRUE,
                                THRESHFUN = function(x) which(x >= dp_tresh))
mtch <- matchSpectra(query, library, param = prm)
mtch_sub <- mtch[whichQuery(mtch)]

mtches_df <- spectraData(mtch_sub, columns = c( "score", "target_name"))
as.data.frame(mtches_df)

#Load MS1 and store annotations revealed from Library search

