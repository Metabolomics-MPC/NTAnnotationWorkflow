# ==============================================================================
# general functions
# ==============================================================================
# Validation of settings
validateSettings <- function(x) {
  
  # check core setttings
  if(!is.numeric(x$cores)) {
    x$cores <- 1
  }
  
  # check input files ----------------------------------------------------------
  # MS1 data
  if(is.null(x$MS1_data_pos) || is.na(x$MS1_data_pos)) {
    x$MS1_data_pos <- NA_character_
  } else {
    if(!file.exists(x$MS1_data_pos)) {
      stop("MS1 data positive mode input doesn't exist!")
    }
  }

  if(is.null(x$MS1_data_neg) || is.na(x$MS1_data_neg)) {
    x$MS1_data_neg <- NA_character_
  } else {
    if(!file.exists(x$MS1_data_neg)) {
      stop("MS1 data negative mode input doesn't exist!")
    }
  }
  
  # MS1 data
  if(is.null(x$MS2_data_pos) || is.na(x$MS2_data_pos)) {
    x$MS2_data_pos <- NA_character_
  } else {
    if(!file.exists(x$MS2_data_pos)) {
      stop("MS1 data positive mode input doesn't exist!")
    }
  }
  
  if(is.null(x$MS2_data_neg) || is.na(x$MS2_data_neg)) {
    x$MS1_data_neg <- NA_character_
  } else {
    if(!file.exists(x$MS1_data_neg)) {
      stop("MS1 data negative mode input doesn't exist!")
    }
  }
  
  # sample grouping ------------------------------------------------------------
  if(is.null(x$samplegroup) || is.na(x$samplegroup)) {
    x$samplegroup <- FALSE
  } else {
    x$samplegroup <- as.logical(x$samplegroup)
  }
  
  if(is.null(x$studydesign_pos) || is.na(x$studydesign_pos)) {
    x$studydesign_pos <- NA_character_
  } else {
    if(!file.exists(x$studydesign_pos)) {
      stop("File specified in studydesign_pos doesn't exist!")
    }
  }
  
  if(is.null(x$studydesign_neg) || is.na(x$studydesign_neg)) {
    x$studydesign_neg <- NA_character_
  } else {
    if(!file.exists(x$studydesign_neg)) {
      stop("File specified in studydesign_neg doesn't exist!")
    }
  }
  
  # check settings -------------------------------------------------------------
  # check MS1 settings
  if(is.null(x$tolerance_MS1) || is.na(x$tolerance_MS1)) {
    stop("tolerance_MS1 missing!")
  } else {
    x$tolerance_MS1 <- as.numeric(x$tolerance_MS1)
  }
  
  if(is.null(x$ppm_MS1) || is.na(x$ppm_MS1)) {
    stop("ppm_MS1 missing!")
  } else {
    x$ppm_MS1 <- as.numeric(x$ppm_MS1)
  }
  
  if(!all(x$adducts_neg %in% MetaboCoreUtils::adductNames(polarity = "negative"))) {
    stop("wrong negative mode adducts!")
  }
  
  if(!all(x$adducts_pos %in% MetaboCoreUtils::adductNames(polarity = "positive"))) {
    stop("wrong negative mode adducts!")
  }
  
  # check MS2 settings
  if(is.null(x$tolerance_MS2) || is.na(x$tolerance_MS2)) {
    stop("tolerance_MS2 missing!")
  } else {
    x$tolerance_MS2 <- as.numeric(x$tolerance_MS2)
  }
  
  if(is.null(x$ppm_MS2) || is.na(x$ppm_MS2)) {
    stop("ppm_MS2 missing!")
  } else {
    x$ppm_MS2 <- as.numeric(x$ppm_MS2)
  }

  if(is.null(x$dp_tresh) || is.na(x$dp_tresh)) {
    stop("dp_tresh missing!")
  } else {
    x$dp_tresh <- as.numeric(x$dp_tresh)
  }
  
  if(is.null(x$int_tresh) || is.na(x$int_tresh)) {
    stop("dp_tresh missing!")
  } else {
    x$int_tresh <- as.numeric(x$int_tresh)
  }
  
  # check Chromatographic settings
  if(is.null(x$toleranceRt_MS1) || is.na(x$toleranceRt_MS1)) {
    stop("toleranceRt_MS1 missing!")
  } else {
    x$toleranceRt_MS1 <- as.numeric(x$toleranceRt_MS1)
  }
  
  if(is.null(x$toleranceRt_MS2) || is.na(x$toleranceRt_MS2)) {
    stop("toleranceRt_MS2 missing!")
  } else {
    x$toleranceRt_MS2 <- as.numeric(x$toleranceRt_MS2)
  }
  
  # check library settings -----------------------------------------------------
  # MS1 libraries
  if(is.null(x$MS1_lib_ext) || is.na(x$MS1_lib_ext)) {
    x$MS1_lib_ext <- NA_character_
  } else {
    if(!dir.exists(x$MS1_lib_ext)) {
      stop("Directory in MS1_lib_ext doesn't exist")
    }
  }
  
  if(is.null(x$MS1_lib_inhouse) || is.na(x$MS1_lib_inhouse)) {
    x$MS1_lib_inhouse <- NA_character_
  } else {
    if(!dir.exists(x$MS1_lib_inhouse)) {
      stop("Directory in MS1_lib_inhouse doesn't exist")
    }
  }
  
  # MS2 libraries, in-house
  if(is.null(x$MS2_lib_pos) || is.na(x$MS2_lib_pos)) {
    x$MS2_lib_pos <- NA_character_
  } else {
    if(!dir.exists(x$MS2_lib_pos)) {
      stop("Directory in MS2_lib_pos doesn't exist")
    }
  }
  
  if(is.null(x$MS2_lib_neg) || is.na(x$MS2_lib_neg)) {
    x$MS2_lib_neg <- NA_character_
  } else {
    if(!dir.exists(x$MS2_lib_neg)) {
      stop("Directory in MS2_lib_neg doesn't exist")
    }
  }
  
  # MS2 libraries, ext
  if(is.null(x$MS2_lib_pos_ext) || is.na(x$MS2_lib_pos_ext)) {
    x$MS2_lib_pos_ext <- NA_character_
  } else {
    if(!dir.exists(x$MS2_lib_pos_ext)) {
      stop("Directory in MS2_lib_pos doesn't exist")
    }
  }
  
  if(is.null(x$MS2_lib_neg_ext) || is.na(x$MS2_lib_neg_ext)) {
    x$MS2_lib_neg_ext <- NA_character_
  } else {
    if(!dir.exists(x$MS2_lib_neg_ext)) {
      stop("Directory in MS2_lib_neg_ext doesn't exist")
    }
  }
  
  return(x)
  
}

# ==============================================================================
# MS1 related functions
# ==============================================================================
createIsoPatternDb <- function(x) {
  
  isopattern <- Spectra()
  
  for(i in 1:nrow(x)) {
    
  }
  
}

# ==============================================================================
# MS2 related functions
# ==============================================================================
#' Remove fragments below x% of base peak intensity
low_int <- function(x, int_tresh = 1) {
    x > max(x, na.rm = TRUE) * (int_tresh / 100 )
}

#' Normalize intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}

#' Remove precursor ion
removePrecursor <- function(window = 1) {
  
  .between <- function(x, a, b) {
    (x >= a) & (x <= b)
  }  
  
  function(x, precursorMz, ...) {
    x[!.between(x[,1], precursorMz - 0.5*window, precursorMz+0.5*window),,drop=FALSE]
  }
    
}

#' Function for adding Feature Ids from summarized Experiment to a Spectra object
#'
#' @param se SummarizedExperiment
#' @param sps Spectra object
#' 
#' @return Spectra object with new Metadata FeatureID
addFeatureID <- function(sps, se){
  
  # get id and ms id from summarizedExperiment
  d_idx <- data.frame(id = rowData(se)[[1]]$id,
                      ms2_id = rowData(se)[[1]]$ms2_id)
  
  d_idx <- filter(d_idx, ms2_id != "")
  d_idx <- separate_rows(d_idx, ms2_id, sep = "\\|")
  d_idx <- mutate(d_idx, ms2_id = as.integer(str_replace(ms2_id, "_\\(e\\d+\\)", "")))
  
  d_idx <- as.data.frame(d_idx)
  
  sps$FEATUREID <- NA_character_
  
  pb = txtProgressBar(min = 0, max = nrow(d_idx), initial = 0) 
  
  for(i in 1:nrow(d_idx)) {
    sps$FEATUREID[which(sps$number == d_idx$ms2_id[i])] <- d_idx$id[i]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(sps)
}


# ==============================================================================
# GNPS FBMN related functions
# ==============================================================================
createFbmnInput <- function(ms1_data, ms2_spectra) {
  
  # reformat MS1 data to resemble xcms3 FBMN input
  
  # exchange scan indices
  ms2_spectra$scanIndex <- ms2_spectra$CLUSTER_ID %>% 
    str_extract("\\d+") %>% 
    as.integer()
  
  ms2_spectra$acquisitionNum <- ms2_spectra$CLUSTER_ID %>% 
    str_extract("\\d+") %>% 
    as.integer()
  
  ms2_spectra$precScanNum <- ms2_spectra$CLUSTER_ID %>% 
    str_extract("\\d+") %>% 
    as.integer()
  
  # create add additional PEAKDID
  ms2_spectra$PEAKID <- ms2_spectra$FEATUREID
  
  
  # sanity checks
  # does the MS2 have correct columns
  # if not perform add ID first
  if(!"CLUSTER_ID" %in% spectraVariables(spectra)) {
    stop("No column CLUSTER_ID found in spectra. Perform ms2_add_id() first")
  }
  
  # check if data contains NAs and remove
  if(any(is.na(spectra$CLUSTER_ID))) {
    # filter spectra
    spectra <- spectra[which(!is.na(spectra$CLUSTER_ID))]
  }
  
  # create feature table
  # create new DF with feature table in XCMS3 format
  feature_table <- data.frame(Row.names = row.names(row_anno),
                              mzmed = row_anno$`m/z`,
                              mzmin = row_anno$`m/z Min`,
                              mzmax = row_anno$`m/z Max`,
                              rtmed = row_anno$RT * 60,
                              rtmin = row_anno$`RT Min` * 60,
                              rtmax = row_anno$`RT Max` * 60,
                              npeaks = ncol(ms_data),
                              sample = ncol(ms_data))
  
  # rename features according to XCMS and add intensities
  feature_table$Row.names <- gsub("Cluster_", "FT", feature_table$Row.names)
  feature_table <- cbind.data.frame(feature_table, ms_data)
  
  # convert names in MS2 metadata to XCMS3 format
  # create fields required
  spectra$scanIndex <- as.integer(regmatches(spectra$CLUSTER_ID, regexpr("\\d+", spectra$CLUSTER_ID)))
  spectra$FEATURE_ID <- gsub("Cluster_", "FT", spectra$CLUSTER_ID)
  spectra$PEAK_ID <- spectra$FEATURE_ID
  spectra$COMPOUND <- spectra$FEATURE_ID
  
}