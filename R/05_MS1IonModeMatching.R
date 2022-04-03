perform_ionMode_matching <- function(ms1_pos_se,
                                     ms1_neg_se,
                                     adducts_pos,
                                     adducts_neg,
                                     tolerance = 0.000,
                                     ppm = 5,
                                     toleranceRt = NA,
                                     outputdir = NA,
                                     saveRds = TRUE,
                                     saveTsv = FALSE) {

  if(is.na(toleranceRt)) {
    param <- Mz2MassParam(queryAdducts = adducts_pos,
                          targetAdducts = adducts_neg,
                          tolerance = tolerance,
                          ppm = ppm)
  } else {
    param <- Mz2MassRtParam(queryAdducts = adducts_pos,
                            targetAdducts = adducts_neg,
                            tolerance = tolerance,
                            ppm = ppm,
                            toleranceRt = toleranceRt)
  }

  # perform ion mode matching
  se_match <- matchMz(rowData(ms1_pos_se)[[1]],
                      rowData(ms1_neg_se)[[1]],
                      param = param)
  
  # print number of matches
  print(se_match)
  
  # save results in a rds file
  if(saveRds) {
    saveRDS(se_match,
            paste0(outputdir,
                   "/Annotation_MS1_ionMode/",
                   "ms1ionmodematching.rds"))
  }
  
  # save results in a tsv file
  if(saveTsv) {
    write.table(matchedData(se_match),
                paste0(outputdir,
                       "/Annotation_MS1_ionMode/",
                       "ms1ionmodematching.tsv"),
                sep = "\t", row.names = FALSE)
  }
}
