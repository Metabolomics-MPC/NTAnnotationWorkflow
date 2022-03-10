#' Function for performing library search on Slaw Output (see \code{\link[Slaw on github]{https://github.com/zamboni-lab/SLAW}})
#'
#' @param se SummarizedExperiment 
#' @param query Slaw fused mgf output
#' @param library_dir Directory with Library files to use
#' @param library_regex Regex to truncate to specific library
#' @param output Directory where output will be stored (Slaw matrix with with the annotated features)
#' @param annotated_only Shall the output be trucated to only truncated features? Default TRUE
#' @param dp_tresh Threshold for dot product score, default 0.4
#' @param mz_ppm Mass window used for spectra matching in ppm, default 5
#' @param int_tresh Percentage of base peak intensity as threshold that will be removed from Spectra, default 5
#' @param compare_rt Shall annotations be compared on retention times in MassBank records? Use this only if your experminet is performed with the same LC method than your library.
#' @param rt_thresh_min Retention time window, Annotations +- in min of the library will be retained
#' 
#' @return A QFeatures - SummarizedExperiment 
#' The annotated output will be stored in the 
librarysearch_se <- function(se, 
                          query, 
                          library_dir,
                          library_regex,
                          output_dir,
                          output_name,
                          annotated_only=TRUE, 
                          dp_tresh=0.4, 
                          mz_ppm=5,
                          int_tresh=5,
                          compare_rt=FALSE,
                          rt_thresh_min=0.02){
    
    # Load Library, either .mb, .msp or .MSP
    message("Load library data")
    if(nchar(library_regex)>0){
        fls <- list.files(path=library_dir, pattern=library_regex, full.names = T)
    }else{
        fls <- list.files(path=library_dir, pattern=".mb$|.msp$|.MSP$", full.names = T)
    }
    if(length(fls)<1) message("No specral library file found in the applied directory")
    
    # Get file endings for format
    library_format <- unique(unlist(lapply(strsplit(fls, ".", fixed=TRUE), "[[",2)))
    
    if(library_format>0){
        message("Different input formats found for library")
    }
    
    # Load in data
    if(library_format=="mb"){
        message("Library:", paste(fls))
        library <- Spectra(fls, source = MsBackendMassbank(), backend = MsBackendDataFrame())
    }
    
    #Load MSP file
    if(library_format=="msp" | library_format=="MSP"){
        fls <- list.files(path=library_dir, pattern=".msp$|.MSP$", full.names = T)
        if(length(fls)<1) message("No specral library file found in the applied directory")
        message("Library:", paste(fls))
        library <- Spectra(fls, source = MsBackendMsp(), backend = MsBackendDataFrame())
    }
    
    #Modify Spectra objects
    query <- addProcessing(query, norm_int)
    query <- replaceIntensitiesBelow(query, threshold = int_tresh, value = 0)
    
    library <- addProcessing(library, norm_int)
    
    #Perform library search with the MetaboAnnotation package
    prm <- MatchForwardReverseParam(ppm = mz_ppm,
                                    THRESHFUN = function(x) which(x >= dp_tresh))
    message("Performing library search ", output_name, ", this can take time ...")
    mtch <- matchSpectra(query, library, param = prm)
    message("Library search finished")
    
    #Save mtch object
    saveRDS(mtch, paste0(output_dir,"/", output_name, "_librarysearch.rds"))
    
    #Get annotations
    mtches_df <- as.data.frame(spectraData(mtch[whichQuery(mtch)]))
    
    #Convert RT
    mtches_df$target_rtime <- mtches_df$target_rtime/60
    mtches_df$rtime <- mtches_df$rtime/60
    
    #Retention time comparison
    if(compare_rt){
        rt_logical <- which(mtches_df$rtime < (mtches_df$target_rtime + rt_thresh_min) & mtches_df$rtime > (mtches_df$target_rtime - rt_thresh_min))
        mtches_df <- mtches_df[rt_logical,]
    }
    
    #for duplicated Annotations, only keep the one with the highest score (easier?)
    #TODO Change to two possibilities, only keep one or aggregate
    doupleAn <- unique(mtches_df$scanIndex[which(duplicated(mtches_df$scanIndex))])
    if(length(doupleAn)>0){
        remove <- vector()
        for(i in doupleAn){
            x <- mtches_df[which(mtches_df$scanIndex==i),]
            keep <- which.max(x$score)
            remove <- c(remove, as.numeric(rownames(x[-keep,])))
        }
        mtches_df <- mtches_df[-remove, ]
    }
    
    #TODO remove unneeded column from matched object
    
    #Add annotations revealed from library search to summarizedExperiment
    librarySearch <- data.frame(scanIndex=rowData(se)$slaw$scanIndex)
    librarySearch <- dplyr::left_join(librarySearch, mtches_df, by="scanIndex")
    
    #Store library search annotations in SummarizedExperiment
    se@ExperimentList@listData[[length(se@ExperimentList@listData)+1]] <- DataFrame(librarySearch)
    
    #Generate output dataframe from SummarizedExperiment
    data_f <- data.frame(librarySearch,data.frame(rowData(se)$slaw), assay(se))
       
    #Truncate for output
    if(annotated_only){
        data_fi <- data_f[which(!data_f$target_precursorMz==""),]
    }else{
        data_fi <- data_f
    }
    
    #Save output
    write.csv(data_fi, file=paste0(output_dir, "/", output_name, "_datamatrix_annotated.csv"))
    
    return(se)
}