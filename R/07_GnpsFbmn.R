createFbmnInput <- function(ms1_data,
                            ms2_spectra,
                            outputdir = NA,
                            ionmode = "") {
  
  # sanity checks
  # does the MS2 have correct columns
  # if not perform add ID first
  if(!"FEATUREID" %in% spectraVariables(ms2_spectra)) {
    stop("No column FEATUREID found in spectra.")
  }
  
  # check if data contains NAs and remove
  if(any(is.na(ms2_spectra$FEATUREID))) {
    ms2_spectra <- ms2_spectra[which(!is.na(ms2_spectra$FEATUREID))]
  }
  
  # reformat MS1 data to resemble xcms3 FBMN input
  row_anno <- as.data.frame(rowData(ms1_data)[[1]])
  ms_data <- as.data.frame(assay(ms1_data))
  
  # create new dataframe with feature table in XCMS3 format
  feature_table <- data.frame(Row.names = gsub(ionmode, "", row_anno$id),
                              mzmed = row_anno$mz,
                              mzmin = row_anno$min_mz,
                              mzmax = row_anno$max_mz,
                              rtmed = row_anno$rt,
                              rtmin = row_anno$min_rt,
                              rtmax = row_anno$max_rt,
                              npeaks = ncol(ms_data),
                              sample = ncol(ms_data))
  
  # rename features according to XCMS and add intensities
  feature_table <- cbind.data.frame(feature_table, ms_data)
  
  # exchange scan indices
  ms2_spectra$scanIndex <- ms2_spectra$FEATUREID %>% 
    str_extract("\\d+") %>% 
    as.integer()
  
  ms2_spectra$acquisitionNum <- ms2_spectra$FEATUREID %>% 
    str_extract("\\d+") %>% 
    as.integer()
  
  ms2_spectra$precScanNum <- ms2_spectra$FEATUREID %>% 
    str_extract("\\d+") %>% 
    as.integer()
  
  # create add additional PEAKDID
  ms2_spectra$PEAKID <- ms2_spectra$FEATUREID
  
  # write for export
  write.table(feature_table,
              paste0(outputdir,
                     "/FBMN/",
                     ionmode,
                     "_fbmn_ms1.tsv"),
              sep = "\t", row.names = FALSE)
  
  export(ms2_spectra,
         backend = MsBackendMgf(),
         file = paste0(outputdir,
                       "/FBMN/",
                       ionmode,
                       "_fbmn_ms2.mgf"))
  
}
