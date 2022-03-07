# ==============================================================================
# Metabolomics, Lipidomics and Exposomics LC-MS Annotation Workflow 
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
output_dir <- "./testoutput/" #args[1]
settings_yaml <- "./testdata/settings_MetaboAnnotation.yaml" #args[2]

# ==============================================================================
# 0. Setup 
# ==============================================================================
source("R/00_Setup.R")
settings <- prepare_setup(output_dir, settings_yaml)
# ==============================================================================
# 1. Read MS1 data
# ==============================================================================
source("R/01_MS1Import.R")
se <- MS1_export(output_dir, settings)

# ==============================================================================
# 2. Read MS2 data
# ==============================================================================
source("R/02_MS2Import.R")
query <- MS2_export(settings)

# ==============================================================================
# 3. Annotate MS1 data
# ==============================================================================
source("R/03_MS1Annotation.R")
#se <- MS1Annotation(se, output_dir, settings)

# ==============================================================================
# 4. Annotate MS2 data
# ==============================================================================
source("R/04_MS2Annotation.R")
se <- MS2Annotation(se, query, output_dir, settings)

# ==============================================================================
# 5. Generate Output Report
# ==============================================================================
source("R/05_Report.R")
#ReportMetaboAnnotation(se, output_dir, settings)