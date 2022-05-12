#' Function for generating SummarizedExperiment from data in project directory 
#'
#' @param output_dir Filepath MS1 data in csv or tsv format
#' 
#' @returns A QFeature SummarizedExperiment
import_ms1_data <- function(ms1_file,
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

#' Function for reading in MS2 Data 
#'
#' @param project_dir Path to project directory
#' 
#' @returns A Spectra object containing the MS2 data 
import_ms1_spectra <- function(se,
                               fulldata){
  
  # reconstruct from MS1 data
  message("Constructing MS1 spectra from peaks...")
  
  se_rowdata <- rowData(se)[[1]]
  
  if(file.exists(fulldata)) {
    
    # create empty Spectra object for storing isotope pettern
    ms1_spectra <- Spectra()
    
    # read data with all peaks and perform reconstruction
    peaks <- read.delim(fulldata)
    pb = txtProgressBar(min = 0, max = nrow(se_rowdata), initial = 0) 
    
    for(i in 1:nrow(se_rowdata)) {
      
      # selected mz and RT of feature for MS1 reconstruction
      selected_mz <- se_rowdata$mz[i]
      selected_rt <- se_rowdata$rt[i]
      selected_id <- se_rowdata$id[i]
      
      # get clique, group and adduct of main peak
      selected_clique <- peaks$clique[which(peaks$mz == selected_mz & peaks$rt == selected_rt)]
      selected_group <- peaks$group[which(peaks$mz == selected_mz & peaks$rt == selected_rt)]
      selected_adduct <- unique(str_replace_all(peaks$annotation[which(peaks$mz == selected_mz & peaks$rt == selected_rt)], " \\+\\d+", ""))
      
      # select corresponding peaks
      selected_peaks <- peaks[which(peaks$clique == selected_clique &
                                      peaks$group == selected_group &
                                      grepl(rex(selected_adduct), peaks$annotation)),]
      
      # get peaks for one adduct
      mz <- selected_peaks$mz
      intensity <- rowSums(selected_peaks[colnames(peaks)[grepl("^intensity.*", colnames(peaks))]])
      
      # create Spectra object
      iso_spd <- DataFrame(msLevel = 1L,
                           FEATUREID = selected_id)
      iso_spd$mz <- list(mz)
      iso_spd$intensity <- list(intensity)
      iso_sps <- Spectra(iso_spd)
      
      plotSpectra(iso_sps)
      
      # add to spectra object
      ms1_spectra <- c(ms1_spectra, iso_sps)
      
      setTxtProgressBar(pb,i)
    }
    
    close(pb)
    message("... complete")
    return(ms1_spectra)
    
  } else {
    
    message("...data not found!")
    return(NA)
    
  }
}

