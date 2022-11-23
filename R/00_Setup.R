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
  
  if (!require("MsBackendMsp")) devtools::install_github("rformassspectrometry/MsBackendMsp")
  library(MsBackendMsp)
  
  if (!require("Spectra")) devtools::install_github("rformassspectrometry/Spectra")
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
  
  if(!require("stringr")) install.packages("stringr")
  library(stringr)
  
  if(!require("tidyr")) install.packages("tidyr")
  library(tidyr)
  
  if(!require("rex")) install.packages("rex")
  library(rex)
  
  if(!require("crayon")) install.packages("crayon")
  library(crayon)

}

# Load all helperfunctions
source("R/helperfunctions.R")
#source("R/slaw2summarizedExperiment.R")
