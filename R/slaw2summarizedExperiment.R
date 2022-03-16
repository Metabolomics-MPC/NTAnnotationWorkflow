#' Function for storing a  Slaw Output (see \code{\link[Slaw on github]{https://github.com/zamboni-lab/SLAW}}) into a Summarized Experiment class
#'
#' @param datmatrix Path to slaw csv file with the MS1 data matrix and annotations
#' @param group A vector defining different groups of samples
#' 
#' @return SummarizedExperiment Object
slaw2summarizedExperiment <- function(datamatrix, group=NULL){
    #Read in the data
    data <- read.delim(datamatrix)
    #Define columns for assay
    int_begin <- grep("^intensity", colnames(data))[1]
    #Start QFeature object of class SummarizedExperiment
    x <- readQFeatures(data, ecol= int_begin:ncol(data), name="slaw", colData=group)
    #Add group information of samples
    if(!is.null(group)){
        x@colData <- DataFrame(name= sub(".csv","", sub("intensity_","",colnames(data)[int_begin:ncol(data)])), group=group)
    }else{
        x@colData <- DataFrame(name= sub(".csv","", sub("intensity_","",colnames(data)[int_begin:ncol(data)])))
    }   
       
    return(x)  
}

#' Function for adding Feature Ids from summarized Experiment to a Spectra object
#'
#' @param se SummarizedExperiment
#' @param sps Spectra object
#' 
#' @return Spectra object with new Metadata FeatureID
addFeatureID <- function(sps, se){
    # get id and ms id from summarizedExperiment
    d_idx <- data.frame(id=rowData(se)$slaw$id, ms2=as.numeric(unlist(lapply(strsplit(rowData(se)$slaw$ms2_id, "_"), "[", 1))))
    d_idx <- d_idx[which(!is.na(d_idx$ms2)),]
    # update to length of spectra object
    s_idx <- data.frame(ms2=seq(1,length(sps)))
    idxes <- right_join(d_idx, s_idx, by="ms2")
    # add values to spectra object  
    sps$FEATUREID <- idxes$id
    
    return(sps)
}