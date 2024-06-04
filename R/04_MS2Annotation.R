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
                                   rindex = FALSE,
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
    
    # change to rindex if indexing is used
    if(rindex) {
      
      param <- MatchForwardReverseParam(tolerance = tolerance,
                                        ppm = ppm,
                                        toleranceRt = toleranceRt,
                                        requirePrecursor = TRUE,
                                        THRESHFUN = function(x) which(x >= dpTresh),
                                        rtColname = c("rindex", "rindex"))
      
    } else {
      
      param <- MatchForwardReverseParam(tolerance = tolerance,
                                        ppm = ppm,
                                        toleranceRt = toleranceRt,
                                        requirePrecursor = TRUE,
                                        THRESHFUN = function(x) which(x >= dpTresh),
                                        rtColname = c("rt", "rt"))
      
    }
  }
  
  # modify query spectra
  spectra <- addProcessing(spectra, norm_int)
  spectra <- filterIntensity(spectra, intensity = c(relIntTresh, Inf))
  spectra <- applyProcessing(spectra)

  # perform matching for each library in libpath
  ms2_libraries <- list.files(libpath, full.names = TRUE)
  
  for(ms2_library in ms2_libraries) {
    
    try({
      
      print(ms2_library)
      
      # read library data and modify spectra
      if(grepl(".mb$", ms2_library)) {
        ms2_lib_data <- Spectra(ms2_library,
                                source = MsBackendMassbank(),
                                backend = MsBackendDataFrame())
      } else if(grepl(".mblib$", ms2_library)) {
        ms2_lib_data <- Spectra(ms2_library,
                                source = MsBackendMassbank(),
                                backend = MsBackendDataFrame())
      } else if(grepl(".msp$", ms2_library)) {
        ms2_lib_data <- Spectra(ms2_library,
                                source = MsBackendMsp(),
                                backend = MsBackendDataFrame())
      } else if(grepl(".rds$", ms2_library)){
        ms2_lib_data <- readRDS(ms2_library)
        if(!class(ms2_lib_data) == "Spectra") next
      }
      
      # modify library spectra
      ms2_lib_data <- addProcessing(ms2_lib_data, norm_int)
      ms2_lib_data <- applyProcessing(ms2_lib_data)

      # perform annotation
      spectra_match <- matchSpectra(setBackend(spectra, MsBackendMemory()),
                                    setBackend(ms2_lib_data, MsBackendMemory()),
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
        if(nrow(matchedData(spectra_match)) > 0) {
          
          # remove potential list columns for export in text file
          matched_orig <- matchedData(spectra_match)
          matched_dropped <- matched_orig[,which(!sapply(matched_orig, class) == "list")]
          matched_dropped <- as.tibble(matched_dropped)
          
          # add potentially missing columns
          cols <- tibble("msLevel" = NA_integer_,
                          "rtime" = NA_real_,
                          "acquisitionNum" = NA_integer_,
                          "scanIndex" = NA_integer_,
                          "precursorMz" = NA_real_,
                          "precursorCharge" = NA_integer_,
                          "collisionEnergy" = NA_real_,
                          "title" = NA_character_,
                          "ms2_id" = NA_character_,
                          "slaw_id" = NA_integer_,
                          "FEATUREID" = NA_integer_,
                          "target_msLevel" = NA_integer_,
                          "target_rtime" = NA_real_,
                          "target_acquisitionNum" = NA_integer_,
                          "target_scanIndex" = NA_integer_,
                          "target_precursorMz" = NA_real_,
                          "target_precursorCharge" = NA_integer_,
                          "target_collisionEnergy" = NA_real_,
                          "target_title" = NA_character_,
                          "target_accession" = NA_character_,
                          "target_name" = NA_character_,
                          "target_exactmass" = NA_real_,
                          "target_formula" = NA_character_,
                          "target_smiles" = NA_character_,
                          "target_inchi" = NA_character_,
                          "target_inchikey" = NA_character_,
                          "target_adduct" = NA_character_,
                          "score" = NA_real_,
                          "reverse_score" = NA_real_,
                          "presence_ratio" = NA_real_,
                          "matched_peaks_count" = NA_integer_)

          matched_dropped <- add_column(matched_dropped,
                                        cols[!names(cols) %in% names(matched_dropped)])
          
          # sort columns
          matched_dropped <- select(as.data.frame(matched_dropped),
                                    "msLevel",
                                    "rtime",
                                    "acquisitionNum",
                                    "scanIndex",
                                    "precursorMz",
                                    "precursorCharge",
                                    "collisionEnergy",
                                    "title",
                                    "ms2_id",
                                    "slaw_id",
                                    "FEATUREID",
                                    "target_msLevel",
                                    "target_rtime",
                                    "target_acquisitionNum",
                                    "target_scanIndex",
                                    "target_precursorMz",
                                    "target_precursorCharge",
                                    "target_collisionEnergy",
                                    "target_title",
                                    "target_accession",
                                    "target_name",
                                    "target_exactmass",
                                    "target_formula",
                                    "target_smiles",
                                    "target_inchi",
                                    "target_inchikey",
                                    "target_adduct",
                                    "score",
                                    "reverse_score",
                                    "presence_ratio",
                                    "matched_peaks_count")
          
          # write table
          write.table(matched_dropped,
                      paste0(outputdir,
                             "/Annotation_MS2_external/",
                             ionmode,
                             "_",
                             str_replace(basename(ms2_library), ".msp$|.mb$", ""),
                             "_ms2annotation.tsv"),
                      sep = "\t", row.names = FALSE)
        }
      } else if(saveTsv && !is.na(toleranceRt)) {
        if(nrow(matchedData(spectra_match)) > 0) {
          
          
          # remove potential list columns for export in text file
          matched_orig <- matchedData(spectra_match)
          matched_dropped <- matched_orig[,which(!sapply(matched_orig, class) == "list")]
          matched_dropped <- as.tibble(matched_dropped)
          
          # add potentially missing columns
          cols <- tibble("msLevel" = NA_integer_,
                         "rtime" = NA_real_,
                         "acquisitionNum" = NA_integer_,
                         "scanIndex" = NA_integer_,
                         "precursorMz" = NA_real_,
                         "precursorCharge" = NA_integer_,
                         "collisionEnergy" = NA_real_,
                         "title" = NA_character_,
                         "ms2_id" = NA_character_,
                         "slaw_id" = NA_integer_,
                         "FEATUREID" = NA_integer_,
                         "target_msLevel" = NA_integer_,
                         "target_rtime" = NA_real_,
                         "target_acquisitionNum" = NA_integer_,
                         "target_scanIndex" = NA_integer_,
                         "target_precursorMz" = NA_real_,
                         "target_precursorCharge" = NA_integer_,
                         "target_collisionEnergy" = NA_real_,
                         "target_title" = NA_character_,
                         "target_accession" = NA_character_,
                         "target_name" = NA_character_,
                         "target_exactmass" = NA_real_,
                         "target_formula" = NA_character_,
                         "target_smiles" = NA_character_,
                         "target_inchi" = NA_character_,
                         "target_inchikey" = NA_character_,
                         "target_adduct" = NA_character_,
                         "score" = NA_real_,
                         "reverse_score" = NA_real_,
                         "presence_ratio" = NA_real_,
                         "matched_peaks_count" = NA_integer_)

          matched_dropped <- add_column(matched_dropped,
                                        cols[!names(cols) %in% names(matched_dropped)])

          # sort columns
          matched_dropped <- select(as.data.frame(matched_dropped),
                                    "msLevel",
                                    "rtime",
                                    "acquisitionNum",
                                    "scanIndex",
                                    "precursorMz",
                                    "precursorCharge",
                                    "collisionEnergy",
                                    "title",
                                    "ms2_id",
                                    "slaw_id",
                                    "FEATUREID",
                                    "target_msLevel",
                                    "target_rtime",
                                    "target_acquisitionNum",
                                    "target_scanIndex",
                                    "target_precursorMz",
                                    "target_precursorCharge",
                                    "target_collisionEnergy",
                                    "target_title",
                                    "target_accession",
                                    "target_name",
                                    "target_exactmass",
                                    "target_formula",
                                    "target_smiles",
                                    "target_inchi",
                                    "target_inchikey",
                                    "target_adduct",
                                    "score",
                                    "reverse_score",
                                    "presence_ratio",
                                    "matched_peaks_count")
          
          # write table
          write.table(matched_dropped,
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
    ) # end of try
  }
}
