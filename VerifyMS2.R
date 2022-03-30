# ==============================================================================
# Metabolomics, Lipidomics and Exposomics LC-MS Annotation Workflow 
# Add MS2 Annotations to Summarized Experiment
#
# Authors:
# - Carolin Huber, UFZ
# - Michael Witting, HMGU
#
# This data analysis workflow performs an Evaluation of all MS2 MatchedSpectra objects 
# and adds them to the SE
# ==============================================================================
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# load libraries
library(MetaboAnnotation)
library(QFeatures)
library(Spectra)
library(yaml)

# define Add2SE function
addMatches2SE <- function(mtches, se, name){
    # construct rowData from MatchedObject
    rr <- matchedData(mtches)
    aa <- matrix(data=rep(0, nrow(rr)*ncol(assay(se[[1]])) ), nrow=nrow(rr), ncol=ncol(assay(se[[1]])))
    # get Feature row in SE
    idx <- spectraData(mtches)$FEATUREID
    # construct Signal Intensity Assay
    for(y in 1:length(idx)){
        if(!is.na(idx[y])){
            # get the lines with FeatureID from the assay
            rowline <- which(rowData(se[[1]])$id == idx[y])
            aa[y,] <- assay(se[[1]])[rowline,]
        }
    }
    # Add new rownames without duplicated names identical between rowData and Assay
    rownames(aa) <- as.character(c(1:nrow(aa)))
    rownames(rr) <- as.character(c(1:nrow(rr)))
    se <- addAssay(se, SummarizedExperiment(rowData=rr, assay=aa), name=name)
    return(se)
}

# check if args is supplied, else run demote data
if(is.na(args[1])) {
    message("Running demo data!")
    settings_yaml <- "test_input/settings.yaml"
} else {
    settings_yaml <- args[1]
}

# read in setttings
settings <- read_yaml(settings_yaml)

# get the QF to use from settings file
se_path <- list.files(path=paste0(settings$output_dir, "/QFeatures_MS1/"), pattern="MS1annotated.rds$", full.names = T)
se_path_pos <- se_path[grep("pos", se_path)]
se_path_neg <- se_path[grep("neg", se_path)]

if(!length(se_path_pos)){
    message("No QFeature Object found for positive mode")
}

if(!length(se_path_neg)){
    message("No QFeature Object found for negative mode")
}

# Get the Matched Objects to use from settings file
mo_paths_inhouse <- list.files(path=paste0(settings$output_dir, "/Annotation_MS2_inhouse/"), pattern=".rds$", full.names = T)
mo_paths_inhouse_pos <- mo_paths_inhouse[grep("pos", mo_paths_inhouse)]
mo_paths_inhouse_neg <- mo_paths_inhouse[grep("neg", mo_paths_inhouse)]

mo_paths_extern <- list.files(path=paste0(settings$output_dir, "/Annotation_MS2_external/"), pattern=".rds$", full.names = T)
mo_paths_extern_pos <- mo_paths_extern[grep("pos", mo_paths_extern)]
mo_paths_extern_neg <- mo_paths_extern[grep("neg", mo_paths_extern)]

####################################################
# positive mode ####################################
####################################################
# Open SE
message("Start MS2 data positive mode")
se <- readRDS(file = se_path_pos)

# List of MatchedSpectra objects inhouse
for(i in mo_paths_inhouse_pos){
    mo <- readRDS(i)
    mo_q <- mo[whichQuery(mo),]
    mo_v <- MetaboAnnotation::validateMatchedSpectra(mo_q)
    mo_v <- mo_v[whichQuery(mo_v)]
    #Add validated Matches to SE
    se <- addMatches2SE(mo_v, se, basename(i))
}

# List of MatchedSpectra objects inhouse
for(i in mo_paths_extern_pos){
    mo <- readRDS(i)
    mo_q <- mo[whichQuery(mo),]
    mo_v <- MetaboAnnotation::validateMatchedSpectra(mo_q)
    mo_v <- mo_v[whichQuery(mo_v)]
    #Add validated Matches to SE
    se <- addMatches2SE(mo_v, se, basename(i))
}

# Save annotated SE
saveRDS(se, paste0(dirname(se_path_pos),"/",str_replace(basename(se_path_pos), "MS1", "MS2")))
se <- NULL
####################################################
# negative mode ####################################
####################################################
# Open SE
message("Start MS2 data negative mode")
se <- readRDS(file = se_path_neg)

# List of MatchedSpectra objects inhouse
for(i in mo_paths_inhouse_neg){
    mo <- readRDS(i)
    mo_q <- mo[whichQuery(mo),]
    mo_v <- MetaboAnnotation::validateMatchedSpectra(mo_q)
    mo_v <- mo_v[whichQuery(mo_v)]
    #Add validated Matches to SE
    se <- addMatches2SE(mo_v, se, basename(i))
}

# List of MatchedSpectra objects inhouse
for(i in mo_paths_extern_neg){
    mo <- readRDS(i)
    mo_q <- mo[whichQuery(mo),]
    mo_v <- MetaboAnnotation::validateMatchedSpectra(mo_q)
    mo_v <- mo_v[whichQuery(mo_v)]
    #Add validated Matches to SE
    se <- addMatches2SE(mo_v, se, basename(i))
}

# Save annotated SE
saveRDS(se, paste0(dirname(se_path_pos),"/",str_replace(basename(se_path_neg), "MS1", "MS2")))

