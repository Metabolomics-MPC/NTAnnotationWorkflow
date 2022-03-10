#' Function for generating SummarizedExperiment from data in project directory 
#'
#' @param settings Settings parameter list
#' 
#' @returns A QFeature SummarizedExperiment
MS1_export <- function(settings){
    message("Load MS1 data")
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
    }else{group=rep("sample", ncol(data)-int_begin)}
    
    #Start QFeature object of class SummarizedExperiment
    se <- readQFeatures(data, ecol= int_begin:ncol(data), name="slaw", colData=group)
    #Add group information of samples
    if(!is.null(group)){
        se@colData <- DataFrame(name= sub(".csv","", sub("intensity_","",colnames(data)[int_begin:ncol(data)])), group=group)
    }else{
        se@colData <- DataFrame(name= sub(".csv","", sub("intensity_","",colnames(data)[int_begin:ncol(data)])))
    }   
    
    #Extract ms2_id info and extract ms2 scanIndex in fused_mgf 
    ms2_id <- rowData(se)[[1]]$ms2_id
    ms2_id[which(ms2_id=="")] <- 0
    rowData(se)[[1]]$scanIndex <- as.numeric(unlist(map(strsplit(ms2_id, "_"),  1)))
    
    #Store Summarized experiment in project_dir
    saveRDS(se, file=paste0(settings$output_dir, "/SummarizedExperiment.rds"))
    
    return(se)  
}