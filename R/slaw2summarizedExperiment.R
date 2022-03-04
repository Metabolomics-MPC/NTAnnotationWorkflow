#' Function for storing a  Slaw Output (see \code{\link[Slaw on github]{https://github.com/zamboni-lab/SLAW}}) into a Summarized Experiment class
#'
#' @param datmatrix Path to slaw csv file with the MS1 data matrix and annotations
#' @param group A vector defining different groups of samples
#' 
#' @return SummarizedExperiment Object
slaw2summarizedExperiment <- function(datamatrix, group){
    #Read in the data
    data <- read.delim(datamatrix)
    #Define columns for assay
    int_begin <- grep("^intensity", colnames(data))[1]
    #Start QFeature object of class SummarizedExperiment
    x <- readQFeatures(data, ecol= int_begin:ncol(data), name="intensity", colData=group)
    #Add group information of samples
    x@colData <- DataFrame(name= sub(".csv","", sub("intensity_","",colnames(data)[int_begin:ncol(data)])), group=group)
       
    return(x)  
}

#' Function for adding a fused mgf to a SummarizedExperiment class
#'
#' @param experiment SummarizedExperiment, containing "ms_id" in slaw format
#' @param fused_mgf Path to a fused mgf file
#' 
#' @return SummarizedExperiment Object
fusedMGF2experiment(experiment, fused_mgf){
    
    #TODO
    
    return(experiment)
}