#' Function for reading in MS2 Data 
#'
#' @param project_dir Path to project directory
#' 
#' @returns A Spectra object containing the MS2 data 
import_ms2_spectra <- function(ms2_file){
  
  #Load Fused MGF file
  message("Load MS2 data...")
  
  if(file.exists(ms2_file)) {
    
    if(grepl(".mgf$", ms2_file)) {
      
      ms2_spectra <- Spectra(ms2_file,
                             source = MsBackendMgf(),
                             backend = MsBackendDataFrame())
      
    } else if(grepl(".msp$", ms2_file)) {
      
      ms2_spectra <- Spectra(ms2_file,
                             source = MsBackendMsp(),
                             backend = MsBackendDataFrame())
      
    }
    
    message("... complete")
    return(ms2_spectra[which(ms2_spectra$msLevel == 2L)])
    
  } else {
    
    message("...data not found!")
    return(NA)
  
  }
}
