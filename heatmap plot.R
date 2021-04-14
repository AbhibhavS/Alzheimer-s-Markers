
library("gplots")
#--------------------------Heat Map---------------------------#
setwd("C:/Protein/AD Gene Expression/Abhibhav/Hippo/results")
#install.packages("gplots")
#
par(mfrow=c(2,2))

dir<-c("C:/Protein/AD Gene Expression/Abhibhav/PFC/results",
       "C:/Protein/AD Gene Expression/Abhibhav/Temproral/results",
       "C:/Protein/AD Gene Expression/Abhibhav/Hippo/results",
       "C:/Protein/AD Gene Expression/Abhibhav/EC/results"
)
brain<-c("Prefrontal cortex", "Middle Temporal Gyrus", 
         "Hippocampus",  "Entorhinal Cortex")

my_palette <- colorRampPalette(c("blue", "white"))(n = 250)

#for(i in 1:4){
for(i in 1:4){
  
  
  setwd(dir[i])
  
  markers<-read.csv("Markers.csv", header = T)
  marker_rf<-markers$RF_marker
  marker_lasso<-markers$lasso_marker
  
  marker_rf<-unique(unlist(c(marker_rf)))[!is.na(unique(unlist(c(marker_rf))))]
  marker_lasso<-unique(unlist(c(marker_lasso)))[!is.na(unique(unlist(c(marker_lasso))))]
  marker_comb<- c(marker_lasso, marker_rf)
  
  flag_rf<-as.numeric(unlist(lapply(strsplit(marker_rf, "x"), function(x){x[2]})))
  flag_lasso<-as.numeric(unlist(lapply(strsplit(marker_lasso, "x"), function(x){x[2]})))
  flag_comb<-as.numeric(unlist(lapply(strsplit(marker_comb, "x"), function(x){x[2]})))
  
  
  data<-read.csv("GSE_GPL_Filtered_with_label_.csv", header = T, row.names = 1)
  
  data_rf<-data.matrix(data[,flag_rf])
  data_lasso<-data.matrix(data[,flag_lasso])
  data_comb<-data.matrix(data[,flag_comb])
  
  meta_data<-read.csv("combine_marker_info.csv", header = T)
  
  colnames(data_comb)<- meta_data$Gene_symbol
  colnames(data_lasso)<- meta_data$Gene_symbol[1:length(marker_lasso)]
  colnames(data_rf)<- meta_data$Gene_symbol[(length(marker_lasso)+1):(length(marker_lasso)+length(marker_rf))]
  
  heat_rf<-cor(data_rf)
  heat_lasso<-cor(data_lasso)
  heat_comb<-cor(data_comb[, !duplicated(t(data_comb))])
  
  cols <- rep('black', nrow(heat_comb))
  dupli<-meta_data$Gene_symbol[which(duplicated(meta_data$Gene_symbol))]
  cols[(length(marker_lasso)+1):(length(marker_lasso)+length(marker_rf))] <- 'red'
  cols[which(colnames(heat_comb) %in% dupli)]<-"orange"
  
  
  
 
  #heatmap.2(heat_rf, trace="none", main="EC", margins = c(10,12))
  #heatmap.2(heat_lasso, trace="none", main="EC", margins = c(10,12))
  heatmap.2(heat_comb, trace="none", 
            main=brain[i], margins = c(8,8), 
            density.info=c("none"),
            Rowv = F, Colv = F,
            key = F,
            cexRow = 0.6, cexCol = 0.9, 
            key.title = "Correlated Gene Density",
            key.xlab = "Correlation combine",
            key.ylab = "Density",
            colCol = cols,
            colRow = cols,
            col= my_palette)
  
}




