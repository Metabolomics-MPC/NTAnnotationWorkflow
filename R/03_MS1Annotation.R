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
                                   outputdir = NA,
                                   ionmode = "",
                                   saveRds = TRUE,
                                   saveTsv = FALSE) {
    
  message("MS1 Annotation in ", ionmode)
  # build param object based on RT selection
  if(is.na(toleranceRt)) {
    param <- Mass2MzParam(adducts = adducts,
                          tolerance = tolerance,
                          ppm = ppm)
  } else {
    param <- Mass2MzRtParam(adducts = adducts,
                            tolerance = tolerance,
                            ppm = ppm,
                            toleranceRt = toleranceRt)
  }
  
  # perform matching for each compound library in libpath
  ms1_libraries <- list.files(libpath, full.names = TRUE)
  
  for(ms1_library in ms1_libraries) {
    
    print(ms1_library)
      
    # read library data and perform some sanity checks
    ms1_lib_data <- read.delim(ms1_libraries)
    
    # check if all required columns are in place
    if(class(param) == "Mass2MzParam" && !all(c("id", "name", "formula", "exactmass") %in% colnames(ms1_lib_data))) {
      message(paste0("Missing one or all required columns id, name, formula, exactmass in library ", basename(ms1_library)))
      next
    } else if(class(param) == "Mass2MzRtParam" && !all(c("id", "name", "formula", "exactmass", "rt") %in% colnames(ms1_lib_data))) {
      message(paste0("Missing one or all required columns id, name, formula, exact_mass, rt in library ", basename(ms1_library)))
      next
    }
    
    # perform annotation
    se_match <- matchMz(rowData(se)[[1]],
                        ms1_lib_data,
                        param = param)
    
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
    
    # Add Annotations to SE
    # construct rowData from MatchedObject
    rr <- matchedData(se_match)[whichQuery(se_match),]
    # construct Signal Intensity Assay
    aa <- assay(se[[1]])[whichQuery(se_match),]
    # Add new rownames without duplicated names identical between rowData and Assay
    rownames(aa) <- as.character(c(1:nrow(aa)))
    rownames(rr) <- as.character(c(1:nrow(rr)))
    se <- addAssay(se, SummarizedExperiment(rowData=rr, assay=aa), name=paste0("MS1_",basename(ms1_library)))
  }
  
  # save features with annotation
  if(saveRds) {
    saveRDS(se,
            paste0(settings$output_dir,
                   "/QFeatures_MS1/",
                   ionmode,
                   "_",
                   str_replace(basename(settings$MS1_data_neg), ".tsv$|.csv$", ""),
                   "_qf_MS1annotated.rds"))
  }
}
