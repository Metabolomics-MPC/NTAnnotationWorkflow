#' Function for setting up project folder and reading in settings file
#'
#' @param settings_yaml Path to a yaml file which should be used and copied to project directory
#' 
#' @return settings
prepare_setup <- function(settings_yaml){
    
    # Load all required packages
    {
        if (!require("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        
        if (!require("MetaboAnnotation")) devtools::install_github("rformassspectrometry/MetaboAnnotation")
        library(MetaboAnnotation)
        
        if (!require("MsBackendMassbank")) BiocManager::install("MsBackendMassbank")
        library(MsBackendMassbank)
        
        if (!require("MsBackendMgf")) BiocManager::install("MsBackendMgf")
        library(MsBackendMgf)
        
        if (!require("MsBackendMsp")) BiocManager::install("MsBackendMsp")
        library(MsBackendMsp)
        
        if (!require("Spectra")) BiocManager::install("Spectra")
        library(Spectra)
        
        if (!require("QFeatures")) BiocManager::install("QFeatures")
        library(QFeatures)
        
        if (!require("purrr")) install.packages("purrr")
        library(purrr)
        
        if (!require("dplyr")) install.packages("dplyr")
        library(dplyr)
        
        if (!require("ggplot2")) install.packages("ggplot2")
        library(ggplot2)
        
        if (!require("M3C")) BiocManager::install("M3C")
        library(M3C)
        
        if (!require("yaml")) install.packages("yaml")
        library(yaml)
    }
    
    # Load all helperfunctions
    source("R/helperfunctions.R")
    source("R/librarysearch_se.R")
    source("R/slaw2summarizedExperiment.R")
    
    # Read in settings of yaml file
    settings <- read_yaml(settings_yaml)
    
    # check if file contains single ion mode or two ion modes
    if(!is.null(settings$pos)){
        if(!dir.exists(settings$pos$output_dir)) dir.create(settings$pos$output_dir)
        write_yaml(settings$pos, paste0(settings$pos$output_dir, "/settings_MetaboAnnotation.yaml"))
        if(!dir.exists(settings$neg$output_dir)) dir.create(settings$neg$output_dir)
        write_yaml(settings$neg, paste0(settings$neg$output_dir, "/settings_MetaboAnnotation.yaml"))
    }else{
    # Generate output directory in project folder
    if(!dir.exists(settings$output_dir)) dir.create(settings$output_dir)
    
    # Store yaml file in output directory
    write_yaml(settings, paste0(settings$output_dir, "/settings_MetaboAnnotation.yaml"))
    }    
    return(settings)
}