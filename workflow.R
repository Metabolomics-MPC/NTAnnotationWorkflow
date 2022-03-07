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
project_dir <- args[1]
settings_yaml <- args[2]

# ==============================================================================
# 0. Setup 
# ==============================================================================
source("R/00_Setup.R")
prepare_setup(project_dir, settings_yaml)
# ==============================================================================
# 1. Read MS1 data
# ==============================================================================
source("R/01_MS1Import.R")
MS1_export(project_dir, settings_yaml)

# ==============================================================================
# 2. Read MS2 data
# ==============================================================================
source("R/02_MS2Import.R")
MS2_export(project_dir, settings_yaml)

