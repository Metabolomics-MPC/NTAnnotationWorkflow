#' Function for reading in MS2 Data 
#'
#' @param project_dir Path to project directory
#' 
#' @returns A Spectra object containing the MS2 data 
MS2_export <- function(settings){
    
    #Get filename from settings
    fused_mgf <- settings$fused_mgf
    
    #Load Fused MGF file
    message("Load MS2 data")
    query <- Spectra(fused_mgf, source = MsBackendMgf(), backend = MsBackendDataFrame())
    
    #Generate IDX variable
    query@backend@spectraData$scanIndex <- 1:length(query)
    
    
    return(query)
}