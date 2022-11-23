# Load all required packages
{
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!require("MetaboAnnotation")) devtools::install_github("rformassspectrometry/MetaboAnnotation")
  suppressMessages(library(MetaboAnnotation))
  
  if (!require("MsBackendMassbank")) BiocManager::install("MsBackendMassbank")
  suppressMessages(library(MsBackendMassbank))
  
  if (!require("MsBackendMgf")) BiocManager::install("MsBackendMgf")
  suppressMessages(library(MsBackendMgf))
  
  if (!require("MsBackendMsp")) devtools::install_github("rformassspectrometry/MsBackendMsp")
  suppressMessages(library(MsBackendMsp))
  
  if (!require("Spectra")) devtools::install_github("rformassspectrometry/Spectra")
  suppressMessages(library(Spectra))
  
  if (!require("QFeatures")) BiocManager::install("QFeatures")
  suppressMessages(library(QFeatures))
  
  if (!require("purrr")) install.packages("purrr")
  suppressMessages(library(purrr))
  
  if (!require("dplyr")) install.packages("dplyr")
  suppressMessages(library(dplyr))
  
  if (!require("ggplot2")) install.packages("ggplot2")
  suppressMessages(library(ggplot2))
  
  if (!require("M3C")) BiocManager::install("M3C")
  suppressMessages(library(M3C))
  
  if (!require("yaml")) install.packages("yaml")
  suppressMessages(library(yaml))
  
  if(!require("stringr")) install.packages("stringr")
  suppressMessages(library(stringr))
  
  if(!require("tidyr")) install.packages("tidyr")
  suppressMessages(library(tidyr))
  
  if(!require("rex")) install.packages("rex")
  suppressMessages(library(rex))
  
  if(!require("crayon")) install.packages("crayon")
  suppressMessages(library(crayon))
}

# Load all helperfunctions
source("R/helperfunctions.R")
#source("R/slaw2summarizedExperiment.R")
