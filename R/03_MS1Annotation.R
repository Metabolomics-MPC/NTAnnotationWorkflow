#' Function for annotating MS1 data
#'
#' @param se A SummarizedExperiment
#' @param output Where to store annotated Output
#' @param settings File containg settings
#' 
perform_ms1_annotation <- function(se,
                                   libpath,
                                   adducts,
                                   tolerance = 0.000,
                                   ppm = 5,
                                   toleranceRt = NA,
                                   rindex = FALSE,
                                   outputdir = NA,
                                   ionmode = "",
                                   saveRds = TRUE,
                                   saveTsv = FALSE) {
    
  message("MS1 Annotation in ", ionmode)
  
  # perform matching for each compound library in libpath
  ms1_libraries <- list.files(libpath,
                              full.names = TRUE)
  
  for(ms1_library in ms1_libraries) {
    
    cat(red(paste0(ms1_library, "\n")))
    #print(se)
      
    # read library data and perform some sanity checks
    ms1_lib_data <- read.delim(ms1_library)
    
    # multiply to seconds
    if("rt" %in% colnames(ms1_lib_data)) {
      ms1_lib_data$rt <- ms1_lib_data$rt * 60
    } else {
      ms1_lib_data$rt <- NA_real_
    }
    
    
    # build param object and perform matching based on input settings
    if(is.na(toleranceRt)) {
      
      # ========================================================================
      # annotation only on m/z
      # ========================================================================
      # sanity check on required columns
      if(!all(c("id", "name", "formula", "exactmass") %in% colnames(ms1_lib_data))) {
        
        message(paste0("Missing one or all required columns id, name, formula, exactmass in library ", basename(ms1_library)))
        next
        
      }
      
      # check if m/z is defined or not
      if("mz" %in% colnames(ms1_lib_data)) {
        
        cat(red("Matching based on precalculated m/z\n"))
        
        # build param object
        param <- MzParam(tolerance = tolerance,
                         ppm = ppm)
        
        # perform matching
        se_match <- matchValues(rowData(se)[[1]],
                                ms1_lib_data,
                                param = param,
                                mzColname = c("mz", "mz"))

      } else {
        
        cat(red("Matching based on m/z from exact mass and adducts\n"))
        
        # build param object
        param <- Mass2MzParam(adducts = adducts,
                              tolerance = tolerance,
                              ppm = ppm)
        
        # perform matching
        se_match <- matchValues(rowData(se)[[1]],
                                ms1_lib_data,
                                param = param,
                                mzColname = "mz",
                                massColname = "exactmass")
        
      }
      
    } else {
      
      # change to rindex if indexing is used
      if(rindex) {
        
        #=======================================================================
        # annotation on m/z and rindex
        # ======================================================================
        # sanity check on required columns
        if(!all(c("id", "name", "formula", "exactmass", "rindex") %in% colnames(ms1_lib_data))) {
          
          message(paste0("Missing one or all required columns id, name, formula, exactmass, rindex in library ", basename(ms1_library)))
          next
          
        }
        
        # check if m/z is defined or not
        if("mz" %in% colnames(ms1_lib_data)) {
          
          cat(red("Matching based on precalculated m/z and RI\n"))
          
          # build param object
          param <- MzRtParam(tolerance = tolerance,
                             ppm = ppm,
                             toleranceRt = toleranceRt)
          
          # perform matching
          se_match <- matchValues(rowData(se)[[1]],
                                  ms1_lib_data,
                                  param = param,
                                  mzColname = c("mz", "mz"),
                                  rtColname = c("rindex", "rindex"))
          
        } else {
          
          cat(red("Matching based on m/z from exact mass and adducts and RI\n"))
          
          # build param object
          param <- Mass2MzRtParam(adducts = adducts,
                                  tolerance = tolerance,
                                  ppm = ppm,
                                  toleranceRt = toleranceRt)
          
          # perform matching
          se_match <- matchValues(rowData(se)[[1]],
                                  ms1_lib_data,
                                  param = param,
                                  mzColname = "mz",
                                  massColname = "exactmass",
                                  rtColname = c("rindex", "rindex"))
          
        }
        
      } else {
        
        #=======================================================================
        # annotation on m/z and rtime
        # ======================================================================
        # sanity check on required columns
        if(!all(c("id", "name", "formula", "exactmass", "rt") %in% colnames(ms1_lib_data))) {
          
          message(paste0("Missing one or all required columns id, name, formula, exactmass, rt in library ", basename(ms1_library)))
          next
          
        }
        
        # check if m/z is defined or not
        if("mz" %in% colnames(ms1_lib_data)) {
          
          cat(red("Matching based on precalculated m/z and RT\n"))
          
          # build param object
          param <- MzRtParam(tolerance = tolerance,
                             ppm = ppm,
                             toleranceRt = toleranceRt)
          
          # perform matching
          se_match <- matchValues(rowData(se)[[1]],
                                  ms1_lib_data,
                                  param = param,
                                  mzColname = c("mz", "mz"),
                                  rtColname = c("rt", "rt"))
          
        } else {
          
          cat(red("Matching based on m/z from exact mass and adducts and RT\n"))
          
          # build param object
          param <- Mass2MzRtParam(adducts = adducts,
                                  tolerance = tolerance,
                                  ppm = ppm,
                                  toleranceRt = toleranceRt)
          
          # perform matching
          se_match <- matchValues(rowData(se)[[1]],
                                  ms1_lib_data,
                                  param = param,
                                  mzColname = "mz",
                                  massColname = "exactmass",
                                  rtColname = c("rt", "rt"))
          
        }
      }
    }

    # print number of matches
    print(se_match)
    
    # save results in a rds file
    if(saveRds && class(param) == "Mass2MzParam") {
      saveRDS(se_match,
              paste0(outputdir,
                     "/Annotation_MS1_external/",
                     ionmode,
                     "_",
                     str_replace(basename(ms1_library), ".tsv$", ""),
                     "_ms1annotation.rds"))
    } else if(saveRds && class(param) == "Mass2MzRtParam") {
      saveRDS(se_match,
              paste0(outputdir,
                     "/Annotation_MS1_inhouse/",
                     ionmode,
                     "_",
                     str_replace(basename(ms1_library), ".tsv$", ""),
                     "_ms1annotation.rds"))
    }
    
    # save results in a tsv file
    if(saveTsv && class(param) == "Mass2MzParam") {
      write.table(matchedData(se_match),
                  paste0(outputdir,
                         "/Annotation_MS1_external/",
                         ionmode,
                         "_",
                         str_replace(basename(ms1_library), ".tsv$", ""),
                         "_ms1annotation.tsv"),
                  sep = "\t", row.names = FALSE)
    } else if(saveTsv && class(param) == "Mass2MzRtParam") {
      write.table(matchedData(se_match),
                  paste0(outputdir,
                         "/Annotation_MS1_inhouse/",
                         ionmode,
                         "_",
                         str_replace(basename(ms1_library), ".tsv$", ""),
                         "_ms1annotation.tsv"),
                  sep = "\t", row.names = FALSE)
    }
  }
}
