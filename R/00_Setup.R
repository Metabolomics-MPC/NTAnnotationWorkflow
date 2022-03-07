#' Function for setting up project folder and reading in settings file
#'
#' @param project_dir Path to project directory
#' @param settings_yaml Path to a yaml file which should be used and copied to project directory
#' 
#' @return settings
prepare_setup <- function(output_dir, settings_yaml){
    
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
    
    # Read in settings of yaml file
    settings <- read_yaml(settings_yaml)
    
    # Generate output directory in project folder
    if(!dir.exists(output_dir)) dir.create(output_dir)
    
    # Store yaml file in output directory
    write_yaml(settings, paste0(output_dir, "/settings_MetaboAnnotation.yaml"))
        
    return(settings)
}