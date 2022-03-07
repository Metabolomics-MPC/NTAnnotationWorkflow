################################################
# Skript loading Slaw Input                    # 
# perform library search on fused MGF file     #
# Store Annotations in Signal Intensity matrix #
################################################
### Set variables => Change for your option  ###
mbank_dir <- "../testdata/library/" #Directory with mbank records to annotate
datamatrix <- "../testdata/output_slaw/datamatrix_opt.csv" #Slaw file with the MS1 data matrix
fused_mgf <- "../testdata/output_slaw/fused_mgf_opt.mgf" #Slaw fused mgf output
output <- "../testoutput/" #File with the annotated features

plot_boxplot <- TRUE #Shall boxplot of Signal Intensities be produced for the matches?
plot_PCA <- TRUE #Shall t-SNE/PCA plot be generated for the produced matches?
color_scheme <- c("blank", "QC", "TR") #Patterns in sample names used for coloring.
################################################
### Load or install packages ###################
################################################
{
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    if (!require("MetaboAnnotation")) devtools::install_github("rformassspectrometry/MetaboAnnotation")
    library(MetaboAnnotation)
    
    if (!require("MsBackendMassbank")) BiocManager::install("MsBackendMassbank")
    library(MsBackendMassbank)
    
    if (!require("MsBackendMgf")) BiocManager::install("MsBackendMgf")
    library(MsBackendMassbank)
    
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
}
################################################
source("helperfunctions.R")
source("librarysearch.R")
source("librarysearch_se.R")
source("slaw2summarizedExperiment.R")
################################################
### Workflow ###################################
################################################
mtch <- librarysearch(datamatrix, fused_mgf, library_dir, output, plot_headtail=TRUE)
data_f <- read.csv(paste0(output,"/datamatrix_annotated.csv"))

#Truncate for further calculations
data_f <- data_f[which(!data_f$target_splash==""),]

#Generate signal intensity matrix and groups
if(plot_boxplot | plot_boxplot){

    int_matrix <- data_f[,grep("^intensity", colnames(data_f))] 
    group <- vector(length = ncol(int_matrix))
    for(i in 1:length(color_scheme)){
        group[grep(color_scheme[i], colnames(int_matrix))] <- color_scheme[i]
    }
}

#Boxplot of signal intensities
if(plot_boxplot){
    for(i in 1:nrow(data_f)){
        dl <- tibble(Intensity=as.vector(t(as.matrix(int_matrix[i,]))), Group=as.factor(group))
        p <- ggplot(dl, aes(Group, Intensity)) +
            geom_boxplot() +
            geom_jitter(aes(colour=Group)) +
            labs(title=paste0(data_f$target_name[i])) +
            theme_classic()
        print(p)
        ggsave(paste0(output, "/", i, "_boxplot.png"))
    }
}

#Generate t-SNE and PCA plot for annotations
if(plot_PCA){
    data_log <- log2(int_matrix)
    is.na(data_log)<-sapply(data_log, is.infinite)
    data_log[is.na(data_log)]<-0
    png(paste0(output, "/tsne.png"))
    tsne(data_log, perplex=3, dotsize= 3,labels= group, seed=2000)
    ggsave(paste0(output, "/tsne.png"))
    pca(data_log, labels = group)
    ggsave(paste0(output, "/pca.png"))
}
