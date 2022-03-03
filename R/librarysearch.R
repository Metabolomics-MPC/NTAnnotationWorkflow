#' Function for performing library search on Slaw Output (see \code{\link[Slaw on github]{https://github.com/zamboni-lab/SLAW}})
#'
#' @param datmatrix Slaw file with the MS1 data matrix
#' @param fused_mgf Slaw fused mgf output
#' @param library_dir Directory with MassBank records to annotate
#' @param output Name and filepath for the output generated (Slaw matrix with with the annotated features)
#' @param annotated_only Shall the output be trucated to only trunctuated features? Default TRUE
#' @param dp_tresh Threshold for dot product score, default 0.4
#' @param mz_ppm Mass window used for spectra matching in ppm, default 5
#' @param int_tresh Percentage of base peak intensity as threshold that will be removed from Spectra, default 5
#' @param compare_rt Shall annotations be compared on retention times in MassBank records? Use this only if your experminet is performed with the same LC method than your library.
#' @param rt_thresh_min Retention time window, Annotations +- in min of the library will be retained
#' @param plot_headtail Shall be head-tail plots be produced for the spectra matches? Default FALSE
#' @param png_dir Directory where plots will be stored
#' 
#' @return the match object 
#' The annotated output will be stored in the 

librarysearch <- function(datamatrix, 
                          fused_mgf, 
                          library_dir, 
                          output, 
                          annotated_only=TRUE, 
                          dp_tresh=0.4, 
                          mz_ppm=5,
                          int_tresh=5,
                          compare_rt=FALSE,
                          rt_thresh_min=0.02,
                          plot_headtail=FALSE,
                          png_dir="."){
    
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
    
    prm <- MatchForwardReverseParam(ppm = mz_ppm, requirePrecursor = TRUE,
                                    THRESHFUN = function(x) which(x >= dp_tresh))
    
    #Perform library search with the MetaboAnnotation package
    mtch <- matchSpectra(query, library, param = prm)
    
    #Get annotations
    mtches_df <- as.data.frame(spectraData(mtch[whichQuery(mtch)], columns = c( "scanIndex", "score", "target_name", "target_precursorMz","precursorMz", "target_rtime", "target_splash")))
    mtches_df$scanIndex <- as.numeric(mtches_df$scanIndex)
    mtches_df$target_rtime <- mtches_df$target_rtime /60
    
    #Load MS1 file
    data <- read.delim(datamatrix)
    ms2_id <- data$ms2_id
    ms2_id[which(ms2_id=="")] <- 0
    data$scanIndex <- as.numeric(unlist(map(strsplit(ms2_id, "_"),  1)))
    
    #Add annotations revealed from library search
    data_f <- full_join(mtches_df, data, by="scanIndex")
    
    #Retention time comparison
    if(compare_rt){
        rt_logical <- data_f$rt < data_f$target_rtime + rt_thresh_min & data_f$rt > data_f$target_rtime - rt_thresh_min
        data_f <- data_f[-which(!rt_logical),]
    }
    
    #Truncate for output
    if(annotated_only){
        data_fi <- data_f[which(!data_f$target_splash==""),]
    }
    
    #Save output
    write.csv(data_fi, file=output)
    
    #Plot head-tail plots of annotations
    if(plot_headtail){
        mtch_sub <- mtch[whichQuery(mtch)]
        for(i in 1:length(mtch_sub@matches)){
            png(paste0(png_dir, "/", i, ".png"))
            plotSpectraMirror(mtch_sub[i], main=paste0(round(mtch_sub@matches$score[i],2), " for ", library[mtch_sub@matches$target_idx[i]]$name))
            dev.off()
        }
    }
}