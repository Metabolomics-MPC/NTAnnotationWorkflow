#' Function for generating SummarizedExperiment from data in project directory 
#'
#' @param output_dir Filepath MS1 data in csv or tsv format
#' 
#' @returns A QFeature SummarizedExperiment
import_ms1 <- function(ms1_file,
                       samplegroup = FALSE,
                       studydesign_file = NA,
                       prefix = "",
                       outputdir = "",
                       saveRds = TRUE,
                       saveTsv = FALSE){
  
  message("Load MS1 data...")
  
  if(file.exists(ms1_file)) {
    
    #get name of data used from settings file and read in the data
    data <- read.delim(ms1_file)
    
    # generate ID for MS1 feature
    strlength <- nchar(nrow(data))
    id <- paste0("FT", prefix, str_pad(1:nrow(data), strlength, side = "left", pad = "0"))
    data <- cbind(id, data)
    
    #Define columns for assay
    int_begin <- grep("^intensity", colnames(data))[1]
    
    #Start QFeature object of class SummarizedExperiment
    se <- readQFeatures(data,
                        ecol = int_begin:ncol(data),
                        name ="slaw")
    
    #get group information
    if(samplegroup & !is.null(studydesign_file)){
      
      studydesign <- read.delim(studydesign_file, row.names = 1)
      se@colData <- DataFrame(studydesign)
      
    } else {
      
      studydesgin <- rep("sample", ncol(data)-int_begin)
      
    }
    
    if(saveRds) {
      
      saveRDS(se,
              paste0(outputdir,
                     "/QFeatures_MS1/",
                     prefix,
                     "_",
                     str_replace(basename(ms1_file), ".tsv$|.csv$", ""),
                     "_qfeatures.rds"))
      
    }

    message("... complete")
    return(se)
    
  } else {
    
    message("...data not found!")
    return(NA)
    
  }
}
