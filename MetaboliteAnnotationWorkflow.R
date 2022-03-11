# ==============================================================================
# Metabolomics, Lipidomics and Exposomics LC-MS Annotation Workflow 
# perfoming on two ion modes
#
# Authors:
# - Carolin Huber, UFZ
# - Michael Witting, HMGU
#
# This data analysis workflow perform annotation of untargeted LC-MS data on the
# MS1 and MS2 level using different libraries and matching functions
# ==============================================================================
# get project directory to work on
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#settings_yaml <- args[1]
settings_yaml <- "test_input/settings.yaml"

# ==============================================================================
# 0. Setup 
# ==============================================================================
# source required functions
source("R/00_Setup.R")

# Read in settings of yaml file ------------------------------------------------
settings <- read_yaml(settings_yaml)

# setup output directory with all subfolder ------------------------------------
if(!dir.exists(settings$output_dir)) dir.create(settings$output_dir)

if(!dir.exists(paste0(settings$output_dir, "/Annotation_MS1_external"))) {
  dir.create(paste0(settings$output_dir, "/Annotation_MS1_external"))
}
if(!dir.exists(paste0(settings$output_dir, "/Annotation_MS1_inhouse"))) {
  dir.create(paste0(settings$output_dir, "/Annotation_MS1_inhouse"))
}

# Store yaml file in output directory
write_yaml(settings, paste0(settings$output_dir, "/input_setting.yaml"))

# ==============================================================================
# 1. Read MS1 data
# ==============================================================================
# source required functions ----------------------------------------------------
source("R/01_MS1Import.R")

# read positive and negative mode MS1 data -------------------------------------
ms1_pos_se <- import_ms1(settings$MS1_data_pos,
                         samplegroup = TRUE,
                         studydesign_file = settings$studydesign_pos,
                         prefix = "POS")

ms1_neg_se <- import_ms1(settings$MS1_data_neg,
                         samplegroup = TRUE,
                         studydesign_file = settings$studydesign_neg,
                         prefix = "NEG")

# ==============================================================================
# 2. Read MS2 data
# ==============================================================================
# source required functions ----------------------------------------------------
source("R/02_MS2Import.R")

# read positive and negative mode MS2 data -------------------------------------
ms2_pos_spectra <- import_ms2(settings$MS2_data_pos)
ms2_neg_spectra <- import_ms2(settings$MS2_data_neg)

# add MS1 ID to spectra --------------------------------------------------------
if(!is.na(ms1_pos_se) && !is.na(ms2_pos_spectra)) {
  
}

if(!is.na(ms1_pos_se) && !is.na(ms2_pos_spectra)) {
  
}

# ==============================================================================
# 3. Annotate MS1 data
# ==============================================================================
# source required functions ----------------------------------------------------
source("R/03_MS1Annotation.R")

# perform MS1 annotation for positive mode data --------------------------------
if(!is.na(ms1_pos_se)) {
  
  # perform annotation with in-house libraries
  if(length(list.files(settings$MS1_lib_inhouse))) {
    
    perform_ms1_annotation(ms1_pos_se,
                           settings$MS1_lib_inhouse,
                           adducts = settings$adducts_pos,
                           tolerance = settings$tolerance,
                           ppm = settings$ppm,
                           toleranceRt = settings$toleranceRt,
                           outputdir = settings$output_dir,
                           ionmode = "pos",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }

  # perform annotation with external libraries
  if(length(list.files(settings$MS1_lib_ext))) {
    
    perform_ms1_annotation(ms1_pos_se,
                           settings$MS1_lib_ext,
                           adducts = settings$adducts_pos,
                           tolerance = settings$tolerance,
                           ppm = settings$ppm,
                           toleranceRt = NA,
                           outputdir = settings$output_dir,
                           ionmode = "pos",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }
}

# perform MS1 annotation for positive mode data --------------------------------
if(!is.na(ms1_pos_se)) {
  
  # perform annotation with in-house libraries
  if(length(list.files(settings$MS1_lib_inhouse))) {
    
    perform_ms1_annotation(ms1_neg_se,
                           settings$MS1_lib_inhouse,
                           adducts = settings$adducts_neg,
                           tolerance = settings$tolerance,
                           ppm = settings$ppm,
                           toleranceRt = settings$toleranceRt,
                           outputdir = settings$output_dir,
                           ionmode = "neg",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }
  
  # perform annotation with external libraries
  if(length(list.files(settings$MS1_lib_ext))) {
    
    perform_ms1_annotation(ms1_neg_se,
                           settings$MS1_lib_ext,
                           adducts = settings$adducts_neg,
                           tolerance = settings$tolerance,
                           ppm = settings$ppm,
                           toleranceRt = NA,
                           outputdir = settings$output_dir,
                           ionmode = "neg",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }
}

# # ==============================================================================
# # 4. Annotate MS2 data
# # ==============================================================================
# # source required functions ----------------------------------------------------
# source("R/04_MS2Annotation.R")
# 
# # perform MS2 annotation  for positive mode data -------------------------------
# if(!is.na(ms2_pos_spectra)) {
#   
# }
# 
# # perform MS2 annotation  for positive mode data -------------------------------
# if(!is.na(ms2_neg_spectra)) {
#   
# }
# 
# se <- MS2Annotation(se, query, output_dir, settings)
# 
# #Store annotated SummarizedExperiment in output_dir
# saveRDS(se, file=paste0(output_dir, "/SummarizedExperiment_annotated.rds"))
# # ==============================================================================
# # 5. Generate Output Report
# # ==============================================================================
# #source("R/05_Report.R")
# #ReportMetaboAnnotation(se, output_dir, settings)
# 
# message("Workflow sucessfully finished :)")