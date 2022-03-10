#' Function for performing spectral library search with MS2 data
#'
#' @param se QFeature SummarizedExperiment
#' @param query Spectra object containing MS2
#' @param settings Settings parameter list
#' 
MS2Annotation <- function(se, query, settings){
    
    #Perform external library search without RT
    if(nchar(settings$MS2_lib_ext)>0){
        se <- librarysearch_se(se, query, 
                           library_dir=settings$MS2_lib_ext, 
                           library_regex = settings$MS2_lib_ext_regex,
                           output_dir = settings$output_dir,
                           output_name = "extern",
                           annotated_only=settings$annotated_only, 
                           dp_tresh=settings$dp_tresh, 
                           mz_ppm=settings$mz_ppm,
                           int_tresh=settings$int_tresh,
                           compare_rt=FALSE,
                           rt_thresh_min=settings$rt_thresh_min)
    }
    
    #Perform in-house library search with RT
    if(nchar(settings$MS2_lib_inhouse)>0){
        se <- librarysearch_se(se, query, 
                           library_dir=settings$MS2_lib_inhouse, 
                           library_regex = settings$MS2_lib_inhouse_regex,
                           output_dir = settings$output_dir,
                           output_name = "inhouse",
                           annotated_only=settings$annotated_only, 
                           dp_tresh=settings$dp_tresh, 
                           mz_ppm=settings$mz_ppm,
                           int_tresh=settings$int_tresh,
                           compare_rt=TRUE,
                           rt_thresh_min=settings$rt_thresh_min)
    }
    
    return(se)   
}
