#' Function for generating SummarizedExperiment from data in project directory 
#'
#' @param project_dir Path to project directory
#' 
#' @returns A QFeature SummarizedExperiment
MS1_export <- function(output_dir, settings){

    #get name of data used from settingsfile and read in the data
    data <- read.delim(settings$MS1_data)
    #Define columns for assay
    int_begin <- grep("^intensity", colnames(data))[1]
    #get group information
    if(settings$samplegroup){
        group <- vector(length = (ncol(data)-int_begin))
        for(i in 1:length(settings$pattern)){
            group[grep(settings$pattern[i], colnames(data[int_begin:ncol(data)]))] <- settings$pattern[i]
        }
    }else{group=NULL}
    
    #Start QFeature object of class SummarizedExperiment
    x <- readQFeatures(data, ecol= int_begin:ncol(data), name="slaw", colData=group)
    #Add group information of samples
    if(!is.null(group)){
        x@colData <- DataFrame(name= sub(".csv","", sub("intensity_","",colnames(data)[int_begin:ncol(data)])), group=group)
    }else{
        x@colData <- DataFrame(name= sub(".csv","", sub("intensity_","",colnames(data)[int_begin:ncol(data)])))
    }   
    
    #Store Summarized experiment in project_dir
    saveRDS(x, file=paste0(output_dir, "/SummarizedExperiment.rds"))
    
    return(x)  
}