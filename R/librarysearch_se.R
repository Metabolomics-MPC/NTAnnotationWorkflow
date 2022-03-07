#' Function for performing library search on Slaw Output (see \code{\link[Slaw on github]{https://github.com/zamboni-lab/SLAW}})
#'
#' @param datamatrix datmatrix Slaw file with the MS1 data matrix 
#' @param fused_mgf Slaw fused mgf output
#' @param library_dir Directory with MassBank records to annotate
#' @param output Directory where output will be stored (Slaw matrix with with the annotated features)
#' @param annotated_only Shall the output be trucated to only trunctuated features? Default TRUE
#' @param dp_tresh Threshold for dot product score, default 0.4
#' @param mz_ppm Mass window used for spectra matching in ppm, default 5
#' @param int_tresh Percentage of base peak intensity as threshold that will be removed from Spectra, default 5
#' @param compare_rt Shall annotations be compared on retention times in MassBank records? Use this only if your experminet is performed with the same LC method than your library.
#' @param rt_thresh_min Retention time window, Annotations +- in min of the library will be retained
#' @param plot_headtail Shall be head-tail plots be produced for the spectra matches? Default FALSE
#' 
#' @return A QFeatures - SummarizedExperiment 
#' The annotated output will be stored in the 
librarysearch_se <- function(datamatrix, 
                          fused_mgf, 
                          library_dir, 
                          output, 
                          annotated_only=TRUE, 
                          dp_tresh=0.4, 
                          mz_ppm=5,
                          int_tresh=5,
                          compare_rt=FALSE,
                          rt_thresh_min=0.02,
                          plot_headtail=FALSE){
    
    #Load MassBank Records
    fls <- list.files(path=mbank_dir, pattern=".txt$", full.names = T)
    library <- Spectra(fls, source = MsBackendMassbank(), backend = MsBackendDataFrame())
    
    #Load Fused MGF file
    query <- Spectra(fused_mgf, source = MsBackendMgf(), backend = MsBackendDataFrame())
    
    #Generate IDX variable
    query@backend@spectraData$scanIndex <- 1:length(query)
    
    #Modify Spectra objects
    query <- addProcessing(query, norm_int)
    query <- replaceIntensitiesBelow(query, threshold = int_tresh, value = 0)
    
    library <- addProcessing(library, norm_int)
    
    #Perform library search with the MetaboAnnotation package
    prm <- MatchForwardReverseParam(ppm = mz_ppm,
                                    THRESHFUN = function(x) which(x >= dp_tresh))
    mtch <- matchSpectra(query, library, param = prm)
    
    #Get annotations
    mtches_df <- as.data.frame(spectraData(mtch[whichQuery(mtch)]))
    
    #Convert RT
    mtches_df$target_rtime <- mtches_df$target_rtime /60
    mtches_df$rtime <- mtches_df$rtime /60
    
    #Retention time comparison
    if(compare_rt){
        rt_logical <- mtches_df$rtime < mtches_df$target_rtime + rt_thresh_min & mtches_df$rtime > mtches_df$target_rtime - rt_thresh_min
        mtches_df <- mtches_df[-which(!rt_logical),]
    }
    
    #Load MS1 file
    se <- slaw2summarizedExperiment(datamatrix)
    #Extract ms2_id info and extract ms2 scanIndex in fused_mgf 
    ms2_id <- rowData(se)[[1]]$ms2_id
    ms2_id[which(ms2_id=="")] <- 0
    scanIndex <- as.numeric(unlist(map(strsplit(ms2_id, "_"),  1)))
    
    #Add annotations revealed from library search to summarizedExperiment
    librarySearch <- data.frame(scanIndex)
    librarySearch <- DataFrame(full_join( librarySearch, mtches_df, by="scanIndex"))
    
    #Generate output dataframe from summarizedExperiment
    data_f <- cbind(rowData(se)@listData[[1]],  rowData(se)@listData[[2]],assay(se))
       
    #Truncate for output
    if(annotated_only){
        data_fi <- data_f[which(!data_f$target_splash==""),]
    }else{
        data_fi <- data_f
    }
    
    #Save output
    write.csv(data_fi, file=paste0(output,"datamatrix_annotated.csv"))
    
    #Save mtch object
    saveRDS(mtch, paste0(output,"librarysearch.rds"))
    
    #Plot head-tail plots of annotations
    if(plot_headtail){
        mtch_sub <- mtch[whichQuery(mtch)]
        for(i in 1:nrow(mtch_sub@matches)){
            png(paste0(output, "/", i, ".png"))
            plotSpectraMirror(mtch_sub[i], main=paste0(round(mtch_sub@matches$score[i],2), " for ", data_f$target_name[i]))
            dev.off()
        }
    }
    
    return(se)
}