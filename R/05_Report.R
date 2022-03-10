#' Function for generating further evaluation of all annotated features
#'
#' @param se QFeature SummarizedExperiment
#' @param settings Settings parameter list
#' 
ReportMetaboAnnotation <- function(se, settings){
    if(settings$plot_headtail){
        #load in mtch objects inhouse/extern what exists
        
        #perform
    }
    if(settings$plot_boxplot|settings$plot_pca){
        #load in datamatices inhouse/extern what exists
        
        
    }
    if(settings$plot_boxplot){
        #perform boxplot
        
    }
    if(settings$plot_pca){
        #perform pca
        
    }
}

# data_f <- read.csv(paste0(output_dir,"/datamatrix_annotated.csv"))
# 
# 
# #Plot head-tail plots of annotations
# if(plot_headtail){
#     mtch_sub <- mtch[whichQuery(mtch)]
#     for(i in 1:nrow(mtch_sub@matches)){
#         png(paste0(output_dir, "/", output_name, "_", i, ".png"))
#         plotSpectraMirror(mtch_sub[i], main=paste0(round(mtch_sub@matches$score[i],2), " for ", mtch_sub@target$name[i]))
#         dev.off()
#     }
# }
# 
# #Truncate for further calculations
# data_f <- data_f[which(!data_f$target_splash==""),]
# 
# #Generate signal intensity matrix and groups
# if(plot_boxplot | settings$plot_boxplot){
# 
#     int_matrix <- data_f[,grep("^intensity", colnames(data_f))] 
#     group <- vector(length = ncol(int_matrix))
#     for(i in 1:length(color_scheme)){
#         group[grep(color_scheme[i], colnames(int_matrix))] <- color_scheme[i]
#     }
# }
# 
# #Boxplot of signal intensities
# if(plot_boxplot){
#     for(i in 1:nrow(data_f)){
#         dl <- tibble(Intensity=as.vector(t(as.matrix(int_matrix[i,]))), Group=as.factor(group))
#         p <- ggplot(dl, aes(Group, Intensity)) +
#             geom_boxplot() +
#             geom_jitter(aes(colour=Group)) +
#             labs(title=paste0(data_f$target_name[i])) +
#             theme_classic()
#         print(p)
#         ggsave(paste0(output, "/", i, "_boxplot.png"))
#     }
# }
# 
# #Generate t-SNE and PCA plot for annotations
# if(plot_PCA){
#     data_log <- log2(int_matrix)
#     is.na(data_log)<-sapply(data_log, is.infinite)
#     data_log[is.na(data_log)]<-0
#     png(paste0(output, "/tsne.png"))
#     tsne(data_log, perplex=3, dotsize= 3,labels= group, seed=2000)
#     ggsave(paste0(output, "/tsne.png"))
#     pca(data_log, labels = group)
#     ggsave(paste0(output, "/pca.png"))
