# ==============================================================================
# Metabolomics, Lipidomics and Exposomics LC-MS Annotation Workflow 
# performing on two ion modes
#
# Authors:
# - Carolin Huber, UFZ
# - Michael Witting, HMGU
#
# This data analysis workflow perform annotation of untargeted LC-MS data on the
# MS1 and MS2 level using different libraries and matching functions
#
# This is the development version to be in sync with the devel version of SLAW
#
# ==============================================================================
# get project directory to work on
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(warn=-1)

# check if correct command line args are supplied
if(!length(args)) {
  
  message("Running demo data!")
  input <- "new_test_input"
  output <- "new_test_output"
  libraries <- "new_test_library"
  settings_yaml <- paste0(input, "/settings.yaml")
  
} else {
  
  # check if arguments have correct length
  if(!length(args) == 3) {
    stop("Exactly three arguments are required: Input, output and library folder!")
  }
  
  # check if input folder exists
  if(!dir.exists(args[1])) {
    stop(paste0("Input folder ", args[1], " does not exist!"))
  }
  
  # check if settings file is present in input
  if(!file.exists(paste0(args[1], "/settings.yaml"))) {
    stop("Missing settings.yaml in input folder!")
  }
  
  # check for output folder and create if not present
  if(!dir.exists(args[2])) {
    dir.create(args[2])
  }
  
  # check for library folder
  if(!dir.exists(args[3])) {
    stop(paste0("Library folder ", args[3], " does not exist!"))
  }
  
  input <- args[1]
  output <- args[2]
  libraries <- args[3]
  settings_yaml <- paste0(input, "/settings.yaml")
  
}

# ==============================================================================
# 0. Setup 
# ==============================================================================
# source required functions ----------------------------------------------------
source("R/00_Setup.R")

# Read in settings of yaml file ------------------------------------------------
settings <- read_yaml(settings_yaml)

# overwrite data in settings yaml with manually determined values
settings$output_dir <- output
settings$MS1_lib_ext <- paste0(libraries, "/MS1_external")
settings$MS1_lib_inhouse <- paste0(libraries, "/MS1_inhouse")
settings$MS2_lib_pos <- paste0(libraries, "/MS2_inhouse_pos")
settings$MS2_lib_neg <- paste0(libraries, "/MS2_inhouse_neg")
settings$MS2_lib_pos_ext <- paste0(libraries, "/MS2_external_pos")
settings$MS2_lib_neg_ext <- paste0(libraries, "/MS2_external_neg")

# check if folders for annotations exists
if(!dir.exists(settings$MS1_lib_ext)) {
  dir.create(settings$MS1_lib_ext)
}

if(!dir.exists(settings$MS1_lib_inhouse)) {
  dir.create(settings$MS1_lib_inhouse)
}

if(!dir.exists(settings$MS2_lib_pos)) {
  dir.create(settings$MS2_lib_pos)
}

if(!dir.exists(settings$MS2_lib_neg)) {
  dir.create(settings$MS2_lib_neg)
}

if(!dir.exists(settings$MS2_lib_pos_ext)) {
  dir.create(settings$MS2_lib_pos_ext)
}

if(!dir.exists(settings$MS2_lib_neg_ext)) {
  dir.create(settings$MS2_lib_neg_ext)
}

# check for defined format and read accordingly
if(settings$format == "old") {
  
  cat(blue("==================================================================\n"))
  cat(blue("Old SLAW output \n"))
  cat(blue("==================================================================\n"))
  
  # check for positive mode data -----------------------------------------------
  # standard input files
  settings$MS1_data_pos <- list.files(paste0(input, "/output_slaw_pos/datamatrices"),
                                      pattern = "annotated_peaktable_[a-z0-9]*_reduced.csv$",
                                      full.names = TRUE)
  
  settings$MS2_data_pos <- list.files(paste0(input, "/output_slaw_pos/fused_mgf"),
                                      pattern = "fused_mgf_[a-z0-9]*.mgf$",
                                      full.names = TRUE)
  
  # check for study design
  settings$studydesign_pos <- paste0(input, "/output_slaw_pos/studydesign.csv")
  
  # check for full data matrix for isotope pattern reconstruction
  settings$MS1_data_pos_full <- list.files(paste0(input, "/output_slaw_pos/datamatrices"),
                                           pattern = "annotated_peaktable_[a-z0-9]*_full.csv$",
                                           full.names = TRUE)
  
  # check for negative mode data -----------------------------------------------
  # standard input files
  settings$MS1_data_neg <- list.files(paste0(input, "/output_slaw_neg/datamatrices"),
                                      pattern = "annotated_peaktable_[a-z0-9]*_reduced.csv$",
                                      full.names = TRUE)
  
  settings$MS2_data_neg <- list.files(paste0(input, "/output_slaw_neg/fused_mgf"),
                                      pattern = "fused_mgf_[a-z0-9]*.mgf$",
                                      full.names = TRUE)
  # check for study design
  settings$studydesign_neg <- paste0(input, "/output_slaw_neg/studydesign.csv")
  
  # check for full data matrix for isotope pattern reconstruction
  settings$MS1_data_neg_full <- list.files(paste0(input, "/output_slaw_pos/datamatrices"),
                                           pattern = "annotated_peaktable_[a-z0-9]*_full.csv$",
                                           full.names = TRUE)
  
} else if(settings$format == "new") {
  
  cat(blue("==================================================================\n"))
  cat(blue("New SLAW output (devel) \n"))
  cat(blue("==================================================================\n"))
  
  # check for positive mode data -----------------------------------------------
  # standard input files
  settings$MS1_data_pos <- list.files(paste0(input, "/output_slaw_pos/"),
                                      pattern = "data_reduced_[a-z0-9]*.csv$",
                                      full.names = TRUE)
  
  settings$MS2_data_pos <- list.files(paste0(input, "/output_slaw_pos/"),
                                      pattern = "spectra_[a-z0-9]*.mgf$",
                                      full.names = TRUE)
  
  # check for study design
  settings$studydesign_pos <- paste0(input, "/output_slaw_pos/studydesign.csv")
  
  # check for full data matrix for isotope pattern reconstruction
  settings$MS1_data_pos_full <- list.files(paste0(input, "/output_slaw_pos/"),
                                           pattern = "data_full_[a-z0-9]*.csv$",
                                           full.names = TRUE)
  
  # check for negative mode data -----------------------------------------------
  # standard input files
  settings$MS1_data_neg <- list.files(paste0(input, "/output_slaw_neg/"),
                                      pattern = "data_reduced_[a-z0-9]*.csv$",
                                      full.names = TRUE)
  
  settings$MS2_data_neg <- list.files(paste0(input, "/output_slaw_neg/"),
                                      pattern = "spectra_[a-z0-9]*.mgf$",
                                      full.names = TRUE)
  # check for study design
  settings$studydesign_neg <- paste0(input, "/output_slaw_neg/studydesign.csv")
  
  # check for full data matrix for isotope pattern reconstruction
  settings$MS1_data_neg_full <- list.files(paste0(input, "/output_slaw_pos/"),
                                           pattern = "data_full_[a-z0-9]*.csv$",
                                           full.names = TRUE)
}


# validate settings ------------------------------------------------------------
#settings <- validateSettings(settings)

# setup output directory with all subfolder ------------------------------------
if(!dir.exists(settings$output_dir)) dir.create(settings$output_dir)

if(!dir.exists(paste0(settings$output_dir, "/QFeatures_MS1"))) {
  dir.create(paste0(settings$output_dir, "/QFeatures_MS1"))
}

if(!dir.exists(paste0(settings$output_dir, "/Spectra"))) {
  dir.create(paste0(settings$output_dir, "/Spectra"))
}

if(!dir.exists(paste0(settings$output_dir, "/Annotation_MS1_external"))) {
  dir.create(paste0(settings$output_dir, "/Annotation_MS1_external"))
}

if(!dir.exists(paste0(settings$output_dir, "/Annotation_MS1_inhouse"))) {
  dir.create(paste0(settings$output_dir, "/Annotation_MS1_inhouse"))
}

if(!dir.exists(paste0(settings$output_dir, "/Annotation_MS1_ionMode/"))) {
  dir.create(paste0(settings$output_dir, "/Annotation_MS1_ionMode/"))
}

if(!dir.exists(paste0(settings$output_dir, "/Annotation_MS2_external"))) {
  dir.create(paste0(settings$output_dir, "/Annotation_MS2_external"))
}

if(!dir.exists(paste0(settings$output_dir, "/Annotation_MS2_inhouse"))) {
  dir.create(paste0(settings$output_dir, "/Annotation_MS2_inhouse"))
}

if(!dir.exists(paste0(settings$output_dir, "/Sirius"))) {
  dir.create(paste0(settings$output_dir, "/Sirius"))
}

if(!dir.exists(paste0(settings$output_dir, "/FBMN"))) {
  dir.create(paste0(settings$output_dir, "/FBMN"))
}

# setup parallel backend -------------------------------------------------------
if(is.na(settings$cores) | settings$cores == 1) {
  BPParam <- SerialParam()
} else {
  if(.Platform$OS.type == "windows") {
    BPParam <- SnowParam(workers = settings$cores,
                         progressbar = TRUE)
  } else {
    BPParam <- MulticoreParam(workers = settings$cores,
                              progressbar = TRUE)
  }
}

# Store yaml file in output directory
write_yaml(settings, paste0(settings$output_dir, "/input_settings.yaml"))

# settings for export of spectra -----------------------------------------------
# custom mapping for mgf import
custom_mapping_mgf <- c(rtime = "RTINSECONDS",
                        acquisitionNum = "SCANS",
                        precursorMz = "PEPMASS",
                        precursorIntensity = "PEPMASSINT",
                        precursorCharge = "CHARGE",
                        msLevel = "MSLEVEL")

# ==============================================================================
# 1. Read MS1 data
# ==============================================================================
cat(blue("==================================================================\n"))
cat(blue("Read MS1 data...\n"))
cat(blue("==================================================================\n"))
# source required functions ----------------------------------------------------
source("R/01_MS1Import.R")

# read positive and negative mode MS1 data -------------------------------------
if(length(settings$MS1_data_pos)) {
  ms1_pos_se <- import_ms1_data(settings$MS1_data_pos,
                                samplegroup = settings$samplegroup,
                                studydesign_file = settings$studydesign_pos,
                                prefix = "pos",                        
                                outputdir = settings$output_dir,
                                format = settings$format,
                                saveRds = settings$save_rds,
                                saveTsv = settings$save_tsv)
} else {
  ms1_pos_se <- NA
}

if(length(settings$MS1_data_neg)) {
  ms1_neg_se <- import_ms1_data(settings$MS1_data_neg,
                                samplegroup = settings$samplegroup,
                                studydesign_file = settings$studydesign_neg,
                                prefix = "neg",
                                outputdir = settings$output_dir,
                                format = settings$format,
                                saveRds = settings$save_rds,
                                saveTsv = settings$save_tsv)
} else {
  ms1_neg_se <- NA
}

# MS1 spectra dependent on the old or new format -------------------------------
if(settings$format == "old") {
  
  # reconstruct positive and negative mode MS1 spectra (isotope pattern) -------
  if(!is.na(ms1_pos_se)) {
    ms1_pos_spectra <- reconstruct_ms1_spectra(ms1_pos_se,
                                               settings$MS1_data_pos_full,
                                               BPPARAM = BPParam)
  } else {
    ms1_pos_spectra <- NA
  }
  
  if(!is.na(ms1_neg_se)) {
    ms1_neg_spectra <- reconstruct_ms1_spectra(ms1_neg_se,
                                               settings$MS1_data_neg_full,
                                               BPPARAM = BPParam)
  } else {
    ms1_neg_spectra <- NA
  }
} else if(settings$format == "new") {
  
  # reconstruct positive and negative mode MS1 spectra (isotope pattern) -------
  if(check_ms1_spectra(settings$MS2_data_pos)) {
    ms1_pos_spectra <- import_ms1_spectra(settings$MS2_data_pos)
    ms1_pos_spectra <- addFeatureIDMS1(ms1_pos_spectra,
                                       ms1_pos_se,
                                       format = settings$format)
  } else {
    if(!is.na(ms1_pos_se)) {
      ms1_pos_spectra <- reconstruct_ms1_spectra(ms1_pos_se,
                                                 settings$MS1_data_pos_full)
    } else {
      ms1_pos_spectra <- NA
    }
  }
  
  if(check_ms1_spectra(settings$MS2_data_neg)) {
    ms1_neg_spectra <- import_ms1_spectra(settings$MS2_data_neg)
    ms1_neg_spectra <- addFeatureIDMS1(ms1_neg_spectra,
                                       ms1_neg_se,
                                       format = settings$format)
  } else {
    if(!is.na(ms1_neg_se)) {
      ms1_neg_spectra <- reconstruct_ms1_spectra(ms1_neg_se,
                                                 settings$MS1_data_neg_full)
    } else {
      ms1_neg_spectra <- NA
    }
  }
}

# export MS1 spectra -----------------------------------------------------------
if(!is.na(ms1_pos_spectra)) {
  export(ms1_pos_spectra,
         MsBackendMgf(),
         file = paste0(settings$output_dir,
                       "/Spectra/pos_MS1.mgf"),
         mapping = custom_mapping_mgf)
}

if(!is.na(ms1_neg_spectra)) {
  export(ms1_neg_spectra,
         MsBackendMgf(),
         file = paste0(settings$output_dir,
                       "/Spectra/neg_MS1.mgf"),
         mapping = custom_mapping_mgf)
}

# ==============================================================================
# 2. Read MS2 data
# ==============================================================================
cat(blue("==================================================================\n"))
cat(blue("Read MS2 data...\n"))
cat(blue("==================================================================\n"))
# source required functions ----------------------------------------------------
source("R/02_MS2Import.R")

# read positive and negative mode MS2 spectra ----------------------------------
if(length(settings$MS2_data_pos)) {
  ms2_pos_spectra <- import_ms2_spectra(settings$MS2_data_pos)
} else {
  ms2_pos_spectra <- NA
}

if(length(settings$MS2_data_neg)) {
  ms2_neg_spectra <- import_ms2_spectra(settings$MS2_data_neg)
} else {
  ms2_neg_spectra <- NA
}

# add MS1 ID to spectra --------------------------------------------------------
if(!is.na(ms1_pos_se) && !is.na(ms2_pos_spectra)) {
    ms2_pos_spectra <- addFeatureIDMS2(ms2_pos_spectra,
                                    ms1_pos_se,
                                    format = settings$format)
}

if(!is.na(ms1_neg_se) && !is.na(ms2_neg_spectra)) {
    ms2_neg_spectra <- addFeatureIDMS2(ms2_pos_spectra,
                                    ms1_neg_se,
                                    format = settings$format)
}

# export MS2 spectra -----------------------------------------------------------
if(!is.na(ms2_pos_spectra)) {
  export(ms2_pos_spectra,
         MsBackendMgf(),
         file = paste0(settings$output_dir,
                       "/Spectra/pos_MS2.mgf"),
         mapping = custom_mapping_mgf)
}

if(!is.na(ms2_neg_spectra)) {
  export(ms2_neg_spectra,
         MsBackendMgf(),
         file = paste0(settings$output_dir,
                       "/Spectra/neg_MS2.mgf"),
         mapping = custom_mapping_mgf)
}

# ==============================================================================
# 3. Annotate MS1 data
# ==============================================================================
cat(blue("==================================================================\n"))
cat(blue("Annotate MS1 data...\n"))
cat(blue("==================================================================\n"))
# source required functions ----------------------------------------------------
source("R/03_MS1Annotation.R")

# perform MS1 annotation for positive mode data --------------------------------
if(!is.na(ms1_pos_se)) {
  
  # perform annotation with in-house libraries
  if(!is.na(settings$MS1_lib_inhouse) && length(list.files(settings$MS1_lib_inhouse))) {
     
    perform_ms1_annotation(ms1_pos_se,
                           settings$MS1_lib_inhouse,
                           adducts = settings$adducts_pos,
                           tolerance = settings$tolerance_MS1,
                           ppm = settings$ppm_MS1,
                           toleranceRt = settings$toleranceRt_MS1,
                           outputdir = settings$output_dir,
                           ionmode = "pos",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }

  # perform annotation with external libraries
  if(!is.na(settings$MS1_lib_inhouse) && length(list.files(settings$MS1_lib_ext))) {
    
    perform_ms1_annotation(ms1_pos_se,
                           settings$MS1_lib_ext,
                           adducts = settings$adducts_pos,
                           tolerance = settings$tolerance_MS1,
                           ppm = settings$ppm_MS1,
                           toleranceRt = NA,
                           outputdir = settings$output_dir,
                           ionmode = "pos",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }
}

# perform MS1 annotation for negative mode data --------------------------------
if(!is.na(ms1_neg_se)) {
  
  # perform annotation with in-house libraries
  if(!is.na(settings$MS1_lib_inhouse) && length(list.files(settings$MS1_lib_inhouse))) {
    
    perform_ms1_annotation(ms1_neg_se,
                           settings$MS1_lib_inhouse,
                           adducts = settings$adducts_neg,
                           tolerance = settings$tolerance_MS1,
                           ppm = settings$ppm_MS1,
                           toleranceRt = settings$toleranceRt_MS1,
                           outputdir = settings$output_dir,
                           ionmode = "neg",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }
  
  # perform annotation with external libraries
  if(!is.na(settings$MS1_lib_inhouse) && length(list.files(settings$MS1_lib_ext))) {
    
    perform_ms1_annotation(ms1_neg_se,
                           settings$MS1_lib_ext,
                           adducts = settings$adducts_neg,
                           tolerance = settings$tolerance_MS1,
                           ppm = settings$ppm_MS1,
                           toleranceRt = NA,
                           outputdir = settings$output_dir,
                           ionmode = "neg",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
    
  }
}

# ==============================================================================
# 4. Annotate MS2 data
# ==============================================================================
cat(blue("==================================================================\n"))
cat(blue("Annotate MS2 data...\n"))
cat(blue("==================================================================\n"))
# source required functions ----------------------------------------------------
source("R/04_MS2Annotation.R")

#perform MS2 annotation for positive mode --------------------------------------
if(!is.na(ms2_pos_spectra)) {
  
  # perform annotation with in-house libraries
  if(!is.na(settings$MS2_lib_pos) && length(list.files(settings$MS2_lib_pos))) {
    
    print("Performing matching against in-house libraries")
    
    perform_ms2_annotation(ms2_pos_spectra,
                           settings$MS2_lib_pos,
                           tolerance = settings$tolerance_MS2,
                           ppm = settings$ppm_MS2,
                           toleranceRt = settings$toleranceRt_MS2,
                           dpTresh = settings$dp_tresh,
                           relIntTresh = settings$int_tresh,
                           outputdir = settings$output_dir,
                           ionmode = "pos",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv,
                           BPPARAM = SerialParam())
    
  }
  
  # perform annotation with external libraries
  if(!is.na(settings$MS2_lib_pos_ext) && length(list.files(settings$MS2_lib_pos_ext))) {
    
    print("Performing matching against external libraries")
    
    perform_ms2_annotation(ms2_pos_spectra,
                           settings$MS2_lib_pos_ext,
                           tolerance = settings$tolerance_MS2,
                           ppm = settings$ppm_MS2,
                           toleranceRt = NA,
                           dpTresh = settings$dp_tresh,
                           relIntTresh = settings$int_tresh,
                           outputdir = settings$output_dir,
                           ionmode = "pos",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv,
                           BPPARAM = SerialParam())
      
  }
}

#perform MS2 annotation for negative mode --------------------------------------
if(!is.na(ms2_neg_spectra)) {
  
  # perform annotation with in-house libraries
  if(!is.na(settings$MS2_lib_neg) && length(list.files(settings$MS2_lib_neg))) {
    
    print("Performing matching against in-house libraries")
    
    perform_ms2_annotation(ms2_neg_spectra,
                           settings$MS2_lib_neg,
                           tolerance = settings$tolerance_MS2,
                           ppm = settings$ppm_MS2,
                           toleranceRt = settings$toleranceRt_MS2,
                           dpTresh = settings$dp_tresh,
                           relIntTresh = settings$int_tresh,
                           outputdir = settings$output_dir,
                           ionmode = "neg",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv,
                           BPPARAM = SerialParam())
    
  }
  
  # perform annotation with external libraries
  if(!is.na(settings$MS2_lib_neg_ext) && length(list.files(settings$MS2_lib_neg_ext))) {
    
    print("Performing matching against external libraries")
    
    perform_ms2_annotation(ms2_neg_spectra,
                           settings$MS2_lib_neg_ext,
                           tolerance = settings$tolerance_MS2,
                           ppm = settings$ppm_MS2,
                           toleranceRt = NA,
                           dpTresh = settings$dp_tresh,
                           relIntTresh = settings$int_tresh,
                           outputdir = settings$output_dir,
                           ionmode = "neg",
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv,
                           BPPARAM = SerialParam())
    
  }
}


# ==============================================================================
# 5. Perform positive negative matching
# ==============================================================================
cat(blue("==================================================================\n"))
cat(blue("Perform Ionmode matching...\n"))
cat(blue("==================================================================\n"))
# source required functions ----------------------------------------------------
source("R/05_MS1IonModeMatching.R")

# perform MS1 annotation for positive mode data --------------------------------
if(!is.na(ms1_pos_se) && !is.na(ms1_neg_se) && settings$ion_mode_match) {
  
  perform_ionMode_matching(ms1_pos_se,
                           ms1_neg_se,
                           adducts_pos = settings$adducts_pos,
                           adducts_neg = settings$adducts_neg,
                           tolerance = settings$tolerance_MS1,
                           ppm = settings$ppm_MS1,
                           toleranceRt = settings$toleranceRt_MS1,
                           outputdir = settings$output_dir,
                           saveRds = settings$save_rds,
                           saveTsv = settings$save_tsv)
  
}

# ==============================================================================
# 6. Export Sirius files
# ==============================================================================
cat(blue("==================================================================\n"))
cat(blue("Sirius data export...\n"))
cat(blue("==================================================================\n"))
# source required functions ----------------------------------------------------
source("R/06_SiriusExport.R")

# export for Sirius positive mode data -----------------------------------------
if(!is.na(ms2_pos_spectra) && !is.na(ms1_pos_spectra)) {
  
  exportSirius(ms1_pos_se,
               ms1_pos_spectra,
               ms2_pos_spectra,
               ionmode = "pos",
               outputdir = settings$output_dir) 
  
} else if(!is.na(ms2_pos_spectra)) {
  
  exportSirius(ms1_pos_se,
               ms1_spectra = NA,
               ms2_pos_spectra,
               ionmode = "pos",
               outputdir = settings$output_dir) 
  
}

# export for Sirius positive mode data -----------------------------------------
if(!is.na(ms2_neg_spectra) && !is.na(ms1_neg_spectra)) {
  
  exportSirius(ms1_neg_se,
               ms1_neg_spectra,
               ms2_neg_spectra,
               ionmode = "neg",
               outputdir = settings$output_dir) 
  
} else if(!is.na(ms2_neg_spectra)) {
  
  exportSirius(ms1_neg_se,
               ms1_spectra = NA,
               ms2_neg_spectra,
               ionmode = "neg",
               outputdir = settings$output_dir) 
  
}

# ==============================================================================
# 7. Export FBMN files
# ==============================================================================
cat(blue("==================================================================\n"))
cat(blue("FBMN data export...\n"))
cat(blue("==================================================================\n"))
# source required functions ----------------------------------------------------
source("R/07_GnpsFbmn.R")

# export for FBMN positive mode data -------------------------------------------
if(!is.na(ms1_pos_se) && !is.na(ms2_pos_spectra)) {
  
  createFbmnInput(ms1_pos_se,
                  ms2_pos_spectra,
                  ionmode = "pos",
                  outputdir = settings$output_dir,
                  format = settings$format)
  
}

# export for FBMN negative mode data -------------------------------------------
if(!is.na(ms1_neg_se) && !is.na(ms2_neg_spectra)) {
  
  createFbmnInput(ms1_neg_se,
                  ms2_neg_spectra,
                  ionmode = "neg",
                  outputdir = settings$output_dir,
                  format = settings$format)
  
}

# ==============================================================================
# End of Workflow
# ==============================================================================
message("Workflow sucessfully finished")