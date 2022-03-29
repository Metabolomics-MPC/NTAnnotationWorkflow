#' Function for performing spectral library search with MS2 data
#'
#' @param se QFeature SummarizedExperiment
#' @param query Spectra object containing MS2
#' @param output_dir Filepath to output directory
#' @param settings Settings parameter list
#' 
perform_ms2_annotation <- function(spectra,
                                   libpath,
                                   tolerance = 0.000,
                                   ppm = 5,
                                   toleranceRt = NA,
                                   dpTresh = 0.6,
                                   relIntTresh = 1,
                                   outputdir,
                                   ionmode = "",
                                   saveRds = TRUE,
                                   saveTsv = FALSE,
                                   BPPARAM = SerialParam()) {

  message("MS2 Annotation in ", ionmode)
  
  # build param object based on RT selection
  if(is.na(toleranceRt)) {
    
    param <- MatchForwardReverseParam(tolerance = tolerance,
                                      ppm = ppm,
                                      toleranceRt = Inf,
                                      requirePrecursor = TRUE,
                                      THRESHFUN = function(x) which(x >= dpTresh))
    
  } else {
    
    param <- MatchForwardReverseParam(tolerance = tolerance,
                                      ppm = ppm,
                                      toleranceRt = toleranceRt,
                                      requirePrecursor = TRUE,
                                      THRESHFUN = function(x) which(x >= dpTresh))
    
  }
  
  # modify query spectra
  spectra <- addProcessing(spectra, norm_int)
  spectra <- filterIntensity(spectra, intensity = c(relIntTresh, Inf))

  # perform matching for each library in libpath
  ms2_libraries <- list.files(libpath, full.names = TRUE)
  
  for(ms2_library in ms2_libraries) {
    
    print(ms2_library)
    
    # read library data and modify spectra
    if(grepl(".mb$", ms2_library)) {
      
       ms2_lib_data <- Spectra(ms2_library,
                              source = MsBackendMassbank(),
                              backend = MsBackendDataFrame())
      
    } else if(grepl(".msp$", ms2_library)) {
      
       ms2_lib_data <- Spectra(ms2_library,
                              source = MsBackendMsp(),
                              backend = MsBackendDataFrame())
      
    } else if(grepl(".rds$", ms2_library)){
        
        ms2_lib_data <- readRDS(ms2_library)
        
    }
    
    # modify library spectra
    ms2_lib_data <- addProcessing(ms2_lib_data, norm_int)
    
    # perform annotation
    spectra_match <- matchSpectra(spectra,
                                  ms2_lib_data,
                                  param = param,
                                  BPPARAM = BPPARAM)
    
    # print number of matches
    print(spectra_match)
    
    # save results in a rds file
    if(saveRds && is.na(toleranceRt)) {
      
      saveRDS(spectra_match,
              paste0(outputdir,
                     "/Annotation_MS2_external/",
                     ionmode,
                     "_",
                     str_replace(basename(ms2_library), ".msp$|.mb$", ""),
                     "_ms2annotation.rds"))
      
    } else if(saveRds && !is.na(toleranceRt)) {
      
      saveRDS(spectra_match,
              paste0(outputdir,
                     "/Annotation_MS2_inhouse/",
                     ionmode,
                     "_",
                     str_replace(basename(ms2_library), ".msp$|.mb$", ""),
                     "_ms2annotation.rds"))
      
    } 
    
    # save results in a tsv file
    if(saveTsv && is.na(toleranceRt)) {
      
      write.table(matchedData(spectra_match),
                  paste0(outputdir,
                         "/Annotation_MS2_external/",
                         ionmode,
                         "_",
                         str_replace(basename(ms2_library), ".msp$|.mb$", ""),
                         "_ms2annotation.tsv"),
                  sep = "\t", row.names = FALSE)
      
    } else if(saveTsv && !is.na(toleranceRt)) {
      
      write.table(matchedData(spectra_match),
                  paste0(outputdir,
                         "/Annotation_MS2_inhouse/",
                         ionmode,
                         "_",
                         str_replace(basename(ms2_library), ".msp$|.mb$", ""),
                         "_ms2annotation.tsv"),
                  sep = "\t", row.names = FALSE)
      
    }
  }
}
