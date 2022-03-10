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
settings_yaml <- args[1]

# ==============================================================================
# 0. Setup 
# ==============================================================================
source("R/00_Setup.R")
settings <- suppressPackageStartupMessages(prepare_setup(output_dir, settings_yaml))

# ==============================================================================
# 1. Read MS1 data
# ==============================================================================
source("R/01_MS1Import.R")
se_pos <- MS1_export(settings$pos)
se_neg <- MS1_export(settings$neg)

# ==============================================================================
# 2. Read MS2 data
# ==============================================================================
source("R/02_MS2Import.R")
query_pos <- MS2_export(settings$pos)
query_neg <- MS2_export(settings$neg)

# ==============================================================================
# 3. Annotate MS1 data
# ==============================================================================
source("R/03_MS1Annotation.R")
#se_pos <- MS1Annotation(se_pos, settings$pos)
#se_neg <- MS1Annotation(se_neg, settings$neg)

# ==============================================================================
# 4. Annotate MS2 data
# ==============================================================================
source("R/04_MS2Annotation.R")
se_pos <- MS2Annotation(se_pos, query_pos, settings$pos)
se_neg <- MS2Annotation(se_neg, query_neg, settings$neg)

#Store annotated SummarizedExperiment in output_dir
saveRDS(se_pos, file=paste0(settings$pos$output_dir, "/SummarizedExperiment_annotated.rds"))
saveRDS(se_neg, file=paste0(settings$neg$output_dir, "/SummarizedExperiment_annotated.rds"))
# ==============================================================================
# 5. Generate Output Report
# ==============================================================================
#source("R/05_Report.R")
#ReportMetaboAnnotation(se_pos, settings)
#ReportMetaboAnnotation(se_neg, settings)

message("Workflow sucessfully finished :)")