# ==============================================================================
# general functions
# ==============================================================================
# Validation of settings
validateSettings <- function(x) {
  
  # check core setttings -------------------------------------------------------
  if(!is.numeric(x$cores)) {
    x$cores <- 1
  }
  
  # check format settings ------------------------------------------------------
  if(!x$format %in% c("old", "new")) {
    stop("format needs to be either 'old' or 'new'")
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
addFeatureIDMS1 <- function(sps,
                            se,
                            format = "old"){
  
  d_idx <- data.frame(id = rowData(se)[[1]]$id,
                      ms1_id = rowData(se)[[1]]$slaw_id)
  
  sps$FEATUREID <- NA_character_

  pb = txtProgressBar(min = 0, max = nrow(d_idx), initial = 0) 
  
  for(i in 1:nrow(d_idx)) {
    sps$FEATUREID[which(sps$SLAW_ID == d_idx$ms1_id[i])] <- d_idx$id[i]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  sps <- sps[which(!is.na(sps$FEATUREID))]
  
  return(sps)
}

#' Function for adding Feature Ids from summarized Experiment to a Spectra object
#'
#' @param se SummarizedExperiment
#' @param sps Spectra object
#' 
#' @return Spectra object with new Metadata FeatureID
addFeatureIDMS2 <- function(sps,
                            se,
                            format = "old"){
  
  # get id and ms id from summarizedExperiment
  if(format == "old") {
    d_idx <- data.frame(id = rowData(se)[[1]]$id,
                        ms2_id = rowData(se)[[1]]$ms2_id)
  } else if(format == "new") {
    d_idx <- data.frame(id = rowData(se)[[1]]$id,
                        ms2_id = rowData(se)[[1]]$mgf_ms2_id)
  }
  
  d_idx <- filter(d_idx, ms2_id != "")
  d_idx <- separate_rows(d_idx, ms2_id, sep = "\\|")
  d_idx <- mutate(d_idx, ms2_id = as.integer(str_replace(ms2_id, "_\\(e\\d+(\\.\\d+)*\\)", "")))
  
  d_idx <- as.data.frame(d_idx)
  
  sps$FEATUREID <- NA_character_
  
  pb = txtProgressBar(min = 0, max = nrow(d_idx), initial = 0) 
  
  for(i in 1:nrow(d_idx)) {
    sps$FEATUREID[which(sps$number == d_idx$ms2_id[i])] <- d_idx$id[i]
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  sps <- sps[which(!is.na(sps$FEATUREID))]
  
  return(sps)
}
