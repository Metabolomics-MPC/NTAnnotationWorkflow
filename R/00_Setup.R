# Load all required packages
{
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!require("MetaboAnnotation")) devtools::install_github("rformassspectrometry/MetaboAnnotation")
  suppressPackageStartupMessages(library(MetaboAnnotation, quietly = TRUE))
  
  if (!require("MsBackendMassbank")) BiocManager::install("MsBackendMassbank")
  suppressPackageStartupMessages(library(MsBackendMassbank, quietly = TRUE))
  
  if (!require("MsBackendMgf")) BiocManager::install("MsBackendMgf")
  suppressPackageStartupMessages(library(MsBackendMgf, quietly = TRUE))
  
  if (!require("MsBackendMsp")) devtools::install_github("rformassspectrometry/MsBackendMsp")
  suppressPackageStartupMessages(library(MsBackendMsp, quietly = TRUE))
  
  if (!require("Spectra")) devtools::install_github("rformassspectrometry/Spectra")
  suppressPackageStartupMessages(library(Spectra, quietly = TRUE))
  
  if (!require("QFeatures")) BiocManager::install("QFeatures")
  suppressPackageStartupMessages(library(QFeatures, quietly = TRUE))
  
  if (!require("purrr")) install.packages("purrr")
  suppressPackageStartupMessages(library(purrr, quietly = TRUE))
  
  if (!require("dplyr")) install.packages("dplyr")
  suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
  
  if (!require("ggplot2")) install.packages("ggplot2")
  suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
  
  if (!require("M3C")) BiocManager::install("M3C")
  suppressPackageStartupMessages(library(M3C, quietly = TRUE))
  
  if (!require("yaml")) install.packages("yaml")
  suppressPackageStartupMessages(library(yaml, quietly = TRUE))
  
  if(!require("stringr")) install.packages("stringr")
  suppressPackageStartupMessages(library(stringr, quietly = TRUE))
  
  if(!require("tidyr")) install.packages("tidyr")
  suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
  
  if(!require("rex")) install.packages("rex")
  suppressPackageStartupMessages(library(rex, quietly = TRUE))
  
  if(!require("crayon")) install.packages("crayon")
  suppressPackageStartupMessages(library(crayon, quietly = TRUE))
  
  if(!require("tibble")) install.packages("tibble")
  suppressPackageStartupMessages(library(crayon, quietly = TRUE))
}

# Load all helperfunctions
source("R/helperfunctions.R")
#source("R/slaw2summarizedExperiment.R")
