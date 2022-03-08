#' Function for performing library search on Slaw Output (see \code{\link[Slaw on github]{https://github.com/zamboni-lab/SLAW}})
#'
#' @param se SummarizedExperiment 
#' @param query Slaw fused mgf output
#' @param library_dir Directory with MassBank records to annotate
#' @param output Directory where output will be stored (Slaw matrix with with the annotated features)
#' @param annotated_only Shall the output be trucated to only truncated features? Default TRUE
#' @param dp_tresh Threshold for dot product score, default 0.4
#' @param mz_ppm Mass window used for spectra matching in ppm, default 5
#' @param int_tresh Percentage of base peak intensity as threshold that will be removed from Spectra, default 5
#' @param compare_rt Shall annotations be compared on retention times in MassBank records? Use this only if your experminet is performed with the same LC method than your library.
#' @param rt_thresh_min Retention time window, Annotations +- in min of the library will be retained
#' @param plot_headtail Shall be head-tail plots be produced for the spectra matches? Default FALSE
#' 
#' @return A QFeatures - SummarizedExperiment 
#' The annotated output will be stored in the 
librarysearch_se <- function(se, 
                          query, 
                          library_dir,
                          library_format,
                          output_dir,
                          output_name,
                          annotated_only=TRUE, 
                          dp_tresh=0.4, 
                          mz_ppm=5,
                          int_tresh=5,
                          compare_rt=FALSE,
                          rt_thresh_min=0.02,
                          plot_headtail=FALSE){
    
    #Load MassBank Records
    if(library_format=="mbank"){
        fls <- list.files(path=library_dir, pattern=".txt$", full.names = T)
        library <- Spectra(fls, source = MsBackendMassbank(), backend = MsBackendDataFrame())
    }
    
    #Load MSP file
    if(library_format=="msp"){
        fls <- list.files(path=library_dir, pattern=".msp$|.MSP$", full.names = T)
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
    
    #Get annotations
    mtches_df <- as.data.frame(spectraData(mtch[whichQuery(mtch)]))
    
    #Convert RT
    mtches_df$target_rtime <- mtches_df$target_rtime/60
    mtches_df$rtime <- mtches_df$rtime/60
    
    #Retention time comparison
    if(compare_rt){
        rt_logical <- mtches_df$rtime < mtches_df$target_rtime + rt_thresh_min & mtches_df$rtime > mtches_df$target_rtime - rt_thresh_min
        mtches_df <- mtches_df[-which(!rt_logical),]
    }
    
    #Add annotations revealed from library search to summarizedExperiment
    librarySearch <- data.frame(scanIndex=rowData(se)[[1]]$scanIndex)
    librarySearch <- DataFrame(full_join(librarySearch, mtches_df, by="scanIndex"))
    
    #Store library search annotations in SummarizedExperiment
    rowData(se)[[length(rowData(se))+1]] <- librarySearch
    
    #Generate output dataframe from SummarizedExperiment
    data_f <- data.frame(cbind(librarySearch,rowData(se)[[1]], assay(se)))
       
    #Truncate for output
    if(annotated_only){
        data_fi <- data_f[which(!data_f$target_splash==""),]
    }else{
        data_fi <- data_f
    }
    
    #Save output
    write.csv(data_fi, file=paste0(output_dir, "/", output_name, "_datamatrix_annotated.csv"))
    
    #Save mtch object
    saveRDS(mtch, paste0(output_dir,"/", output_name, "_librarysearch.rds"))
    
    #Plot head-tail plots of annotations
    if(plot_headtail){
        mtch_sub <- mtch[whichQuery(mtch)]
        for(i in 1:nrow(mtch_sub@matches)){
            png(paste0(output_dir, "/", output_name, "_", i, ".png"))
            plotSpectraMirror(mtch_sub[i], main=paste0(round(mtch_sub@matches$score[i],2), " for ", data_f$target_name[i]))
            dev.off()
        }
    }
    
    return(se)
}