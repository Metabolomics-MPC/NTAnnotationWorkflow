exportSirius <- function(ms1_data,
                         ms1_spectra = NA,
                         ms2_spectra,
                         ionmode = "",
                         outputdir = NA) {
  
  # sanity checks on ion mode
  if(ionmode == "pos") {
    adduct <- "[M+?]+"
  } else if(ionmode == "neg") {
    adduct <- "[M+?]-"
  } else {
    stop("Unknown ion mode!")
  }
  
  # get all Feature IDs from MS1 data
  ftids <- rowData(ms1_data)[[1]]$id
  precursor_list <- rowData(ms1_data)[[1]]$mz
  rt_list <- rowData(ms1_data)[[1]]$rt
  
  # create file path
  file_path <- paste0(outputdir,
                      "/Sirius/",
                      ionmode,
                      "_sirius.ms")
  
  # helper function for writing file
  con <- file(description = file_path, open = "wt")
  
  .cat <- function(..., file = con, sep = "", append = FALSE) {
    cat(..., file = file, sep = sep, append = append)
  }
  
  if(!is.na(ms1_spectra)) {
    
    for(ftid in ftids) {
      
      # isolate spectra
      ms1_spectrum <- ms1_spectra[which(ms1_spectra$FEATUREID == ftid)]
      ms2_spectrum <- ms2_spectra[which(ms2_spectra$FEATUREID == ftid)]
      
      # get information if MS1 spectrum is present
      if(length(ms1_spectrum) == 1) {
        # if isotope pattern is present 
        .cat(paste0(">compound ", ftid, "\n"))
        .cat(paste0(">ionization ", adduct, "\n"))
        .cat(paste0(">parentmass ", min(ms1_spectrum$mz), "\n"))
        .cat(paste0(">rt ", rt_list[which(ftids == ftid)] * 60, "\n"))
        .cat("\n>ms1\n")
        .cat(paste(unlist(ms1_spectrum$mz),
                   unlist(ms1_spectrum$intensity),
                   collapse = "\n"))
        .cat("\n")
      } else {
        # if isotope pattern is not present 
        .cat(paste0(">compound ", ftid, "\n"))
        .cat(paste0(">ionization ", adduct, "\n"))
        .cat(paste0(">parentmass ", precursor_list[which(ftids == ftid)], "\n"))
        .cat(paste0(">rt ", rt_list[which(ftids == ftid)] * 60, "\n"))
        .cat("\n>ms1\n")
        .cat(precursor_list[which(ftids == ftid)], " 100")
        .cat("\n")
      }
      
      # get information if MS2 spectrum is present
      if(length(ms2_spectrum) == 1) {
        # if MS2 spectrum is present
        .cat("\n>ms2\n")
        .cat(paste(unlist(ms2_spectrum$mz),
                   unlist(ms2_spectrum$intensity),
                   collapse = "\n"))
        .cat("\n")
        .cat("\n")
      } else {
        # .cat("\n\n>ms2\n")
        # .cat(precursor_list[which(ftids == ftid)], " 100")
        # .cat("\n")
        # .cat("\n")
      }
    }
  } else {
    
    for(ftid in ftids) {
      
      # isolate spectra
      ms2_spectrum <- ms2_spectra[which(ms2_spectra$FEATUREID == ftid)]

      # if isotope pattern is not present 
      .cat(paste0(">compound ", ftid, "\n"))
      .cat(paste0(">ionization ", adduct, "\n"))
      .cat(paste0(">parentmass ", precursor_list[which(ftids == ftid)], "\n"))
      .cat(paste0(">rt ", rt_list[which(ftids == ftid)] * 60, "\n"))
      .cat("\n>ms1\n")
      .cat(precursor_list[which(ftids == ftid)], " 100")
      .cat("\n")
      
      # get information if MS2 spectrum is present
      if(length(ms2_spectrum) == 1) {
        # if MS2 spectrum is present
        .cat("\n>ms2\n")
        .cat(paste(unlist(ms2_spectrum$mz),
                   unlist(ms2_spectrum$intensity),
                   collapse = "\n"))
        .cat("\n")
        .cat("\n")

      } else {
        .cat("\n>ms2\n")
        .cat(precursor_list[which(ftids == ftid)], " 100")
        .cat("\n")
        .cat("\n")
      }
    }
    
  }
  
  close(con)
  
}
