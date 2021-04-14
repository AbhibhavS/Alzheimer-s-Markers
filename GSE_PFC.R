setwd("K:/Abhibhav/PFC")
set.seed(123)

library(GEOquery)
library(e1071)
library(randomForest)
library(varSelRF)
library(caret)
library(qpcR)

GPL4372<- read.csv("GPL4372.csv", header =T)

gse1=getGEO(filename="GSE33000_series_matrix.txt.gz",
            GSEMatrix=TRUE,getGPL = TRUE)
GSE33000<- as.matrix(gse1)
GSE33000<-as.data.frame(GSE33000)

colnames(GSE33000)<- c(read.csv("sample info_GSE33000.csv",
                               header =T))[[1]]

gse2=getGEO(filename="GSE44770_series_matrix.txt.gz",
            GSEMatrix=TRUE,getGPL = TRUE)
GSE44770<- as.matrix(gse2)
GSE44770<-as.data.frame(GSE44770)
colnames(GSE44770)<- c(read.csv("sample info_GSE44770.csv",
                                header =T))[[1]]


GSE_GPL4372<- cbind(GSE33000, GSE44770)
GSE_EC<-GSE_GPL4372[,-c(grep("Huntington", colnames(GSE_GPL4372)))]
a<-which(is.na(GSE_EC), arr.ind=TRUE)
GSE_EC<-GSE_EC[-sort(unique(a[,1])),]

##write.csv(GSE_EC, "GSE_PFC.csv")

zscore<-function(X){(X- mean(X))/sqrt(var(X))} #standardization function
GSE_EC<-apply(GSE_EC, MARGIN = 2, zscore) #standardized expression data
GSE_EC<-as.data.frame(GSE_EC)


#dir.create("results")
setwd(paste(getwd(), "/results", sep=""))

#--------------------------------------------------------------------------------------------------------#

gene_name<-as.data.frame(cbind(c(1:nrow(GSE_EC)), rownames(GSE_EC))) #Probe Id
colnames(gene_name)<- c("index", "gene_name")


gene_id<- as.data.frame(cbind(GPL4372[,1:2], GPL4372$GB_ACC ,GPL4372$RosettaGeneModelID)) #meta informations
colnames(gene_id)<- c("gene_name", "entrez_Gene_ID", "ACC_id", "RosettaGeneModelID")

vec <- merge(x = gene_name, y = gene_id, by="gene_name", x.all = TRUE) #prob id mapped with gene ID and other meta informations

info<-vec
#no_info <- which(vec$entrez_Gene_ID == "" & vec$RosettaGeneModelID == "") #avoiding no informative probes
#info <- vec[-no_info,] #informative probes

uni<-which(!(duplicated(info[,2:3]))) #all index of unique probes

#---------------------------------------------------------------------------------------------------------#


for(i in uni){
  
  if(is.na(info$entrez_Gene_ID[i])){ next }
  else if(length(which(info$entrez_Gene_ID==info$entrez_Gene_ID[i]))==1){ next }
  else{
    flag<-which( info$entrez_Gene_ID==info$entrez_Gene_ID[i] ) #identifying duplicate if exist
    x<-apply(GSE_EC[flag,], MARGIN = 1, IQR) #max quantiles are selected
    info[flag[which(x!=max(x))],]<-NA # rest duplicate entries are set to NA in the meta data
    GSE_EC[flag[which(x!=max(x))],]<-NA # rest duplicate entries are set to NA in the gene expression data
  }
}


#GSE_EC[no_info,]<-NA #uniformative probes are set to NA in gene expression data
meta_data<-vec[which(!is.na(GSE_EC[,1])),] # filtering out the NA data from the meta data
GSE<-GSE_EC[which(!is.na(GSE_EC[,1])), ] # filtering out the NA data from the expression data
GSE<-cbind(meta_data, GSE) #combined meta data and expression data
##write.csv(GSE, "GSE_EC_normalized_processed.csv")

#-----------------------------------------------------------------------------------------------

#downsampling
GSE_data1<- as.data.frame(GSE[, 6:ncol(GSE)]) #only choosing the expression data #ncol(GSE_data)
GSE_control1<- c(grep("control", colnames(GSE_data1)),
                grep("Control", colnames(GSE_data1))) #index of control samples
GSE_affected1<- c(1:length(colnames(GSE_data1)))[-GSE_control1]
GSE_affected<-sample(size = length(GSE_control1), GSE_affected1, replace = F) #downsampling
GSE_data<-GSE_data1[,c(GSE_affected, GSE_control1)]

#-------------------------------------------------------------------------------------

GSE_control<- c(grep("control", colnames(GSE_data))) #index of control samples

welch<-function(X){
  a<-t.test(X[GSE_control], X[-GSE_control], var.equal = F)#welch t tailed test
  return(a$p.value)
}                                                                 

wilcox<-function(X){
  a<-wilcox.test(X[GSE_control], X[-GSE_control], var.equal = F)#wilcoxn run tailed test
  return(a$p.value)
}

pvalue_t<-apply(GSE_data, MARGIN = 1, welch) #pvalue generated for expression data through t test
pvalue_wlx<-apply(GSE_data, MARGIN = 1, wilcox) #pvalue generated for expression data through wlx test

my_GSE<- cbind(GSE[,1:5], GSE_data, pvalue_t, pvalue_wlx) #meta data + expression data + p values
my_GSE<- as.data.frame(my_GSE[order(my_GSE$pvalue_t, decreasing = F),]) #ordering the data by t pvalue
my_GSE$pvalue_t<- as.numeric(my_GSE$pvalue_t)
my_GSE$pvalue_wlx<- as.numeric(my_GSE$pvalue_wlx)
my_GSE$adjust<-p.adjust(my_GSE$pvalue_t, method = "bonferroni") #FDR test

Filter_GSE<-my_GSE[which(my_GSE$pvalue_t <= 0.05 
                         & my_GSE$pvalue_wlx <= 0.05),]

Filter_GSE<- Filter_GSE[ which(my_GSE$adjust <= 0.01),] #choosing the appropriate gene expression data

##write.csv(Filter_GSE, "filter_GSE.csv", row.names = F)


#--------------------------------------------------------------------------------------------------------#


data<-Filter_GSE[,6:(ncol(Filter_GSE)-3)] #choosing just the expression data
info<-Filter_GSE[,1:5] #choosing the meta data
#write.csv(info, "meta_data_ml_data.csv", row.names = T)

data<-t(data) #transposed to make sample X gene format

label<-rep(NA, length(row.names(data))) #label vector
label[GSE_control]<- 0 #normal is denoted a 0
label[is.na(label)]<- 1 #affected is denoted a 1

ml_data<-as.data.frame(cbind(data, label)) #machine leaning compatible data
colnames(ml_data)<- c( paste("x",c(1:(ncol(ml_data)-1)), sep=""), "label") #legal naming
shuffle<-sample(1:nrow(ml_data))
#write.table(shuffle, "seed.txt", row.names = F)
#shuffle<-c(read.table("C:/Protein/AD Gene Expression/Abhibhav/PFC/seed.txt", header = T))[[1]] #beacuse my stupid ass has forgot to put seed

ml_data<-ml_data[shuffle,] #shuffling randomly
ml_data$label<-as.factor(ml_data$label) #label factorized

#write.csv(ml_data, "GSE_GPL_Filtered_with_label_.csv", row.names=T)
##write.csv(info, "meta_data_ml_data.csv", row.names = F)
#ml_data<-read.csv("GSE_GPL_Filtered_with_label_.csv", header = T, row.names = 1)
#ml_data$label<-as.factor(ml_data$label)

#_______________________________FEATURE SELECTION and Machine leatning____________________________________

metric<-function(tab){
  TP<-tab[2,2]; FN<-tab[1,2]; TN<-tab[1,1]; FP<-tab[2,1]
  ACC = (sum(diag(tab)))/sum(tab)
  Sen =TP/(TP+FN)
  Spe =TN/(TN+FP) 
  Pre = TP/(TP+FP)
  MCC = (TP*TN - FP*FN)/sqrt((TP+FN)*(TP+FP)*(TN+FN)*(TN+FP))
  return(c(ACC, Sen, Spe, Pre, MCC))
}
row_name<-c(paste(sort(rep(1:5,3)), "CV", rep(c("SVM","RF", "ELasticN"),5)), c("Av.SVM", "Av.RF", "AV.ElasticN"))
Acc<-rep(NA,18); Sen<-rep(NA,18); Spe<-rep(NA,18); Pre<-rep(NA,18); Mcc<-rep(NA,18)
metric_table<-as.data.frame(cbind(Acc, Sen, Spe, Pre, Mcc))
row.names(metric_table)<-row_name


part<- sort(rep(1:5, nrow(ml_data)/5)) #5 CV indexes

#----------------------------------------------------------------------------------------------------------------

######## Random Forest Method ###########
vars<-NULL
for(i in 1:1){
  a<-varSelRF(ml_data[,-ncol(ml_data)], ml_data[,ncol(ml_data)], whole.range = F, verbose = T) #random forest feature selection method
  vars_new<-sort(a$selected.vars) # best set of variables
  vars<-qpcR:::cbind.na(vars, vars_new)
}
vars[which(is.na(vars))]<-0
RF_marker<-rownames(as.matrix(table(vars)[-1]))[which(as.matrix(table(vars)[-1])>=1)]
##write.table(table(vars)[-1], "GSE44770_markers_RF.txt", row.names = F)
imp_data_RF<-ml_data[,c(RF_marker, "label")] #new feature set with label


  

#Classification with RF found marker
imp_data<-imp_data_RF


for(i in 1:5){
  
  train<- imp_data[which(part!=i),]
  test<- imp_data[which(part==i),]
  
  model1<-svm(label~., data= train, kernel="radial", degree=7)
  model2<-randomForest(label~., data= train)
  model3 <- train(label~.,
                  train, 
                  method = "glmnet",
                  tuneGrid = expand.grid(alpha= seq(0,1,length=10),
                                         lambda= seq(0.0001, 1, length = 5)))
  
  pred1<-predict(model1, test)
  pred2<-predict(model2, test)
  pred3 <- predict(model3, test)
  tab1<-table(Predicted= pred1, Actual_svm= test$label)
  tab2<-table(Predicted= pred2, Actual_rf= test$label)
  tab3<-table(Predicted= pred3, Actual_en= test$label)
  
  metric_table[(3*(i-1)+1),]<-metric(tab1)
  metric_table[(3*(i-1)+2),]<-metric(tab2)
  metric_table[(3*(i-1)+3),]<-metric(tab3)
  
  print(tab1)
  print(tab2)
  print(tab3)
  
}
metric_table[16,]<-colMeans(metric_table[seq(1,15,3),])
metric_table[17,]<-colMeans(metric_table[(seq(1,15,3)+1),])
metric_table[18,]<-colMeans(metric_table[(seq(1,15,3)+2),])
t1<-metric_table
##write.csv(t1, "RFmethod_metric.csv")

#---------------------------------------------------------------------------------------------------------
Cstack_info()
############ LASSO Method####################
options(expressions = 5e5)

ml<-apply(ml_data,1,as.numeric)
ml<-t(as.matrix(ml))
colnames(ml)<-colnames(ml_data)

custom <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 5,
                       verboseIter = T)

lasso <- train(x=ml[,-ncol(ml)],
               y=as.factor(ml[,ncol(ml)]),
               data = ml,
               method = "glmnet",
               tuneGrid = expand.grid(alpha=1,
                                      lambda= seq(10^(-4), 1, length = 5)),
               trControl = custom)

imp<-varImp(lasso, scale=T)
imp<-as.matrix(imp$importance)
imp<-imp[order(imp[,1], decreasing = T),]
lasso_marker<-names(which(imp>=10))

#imp<-sort(imp[,1], decreasing = T)
#lasso_marker<-names(imp[1:20])

##write.table(lasso_marker, "GSE44770_markers_lasso.txt", row.names = F)
imp_data_lasso<- ml_data[,c(lasso_marker, "label")]

#Classification with LASSO found marker
imp_data<-imp_data_lasso
#imp_data<-ml_data

for(i in 1:length(unique(part))){
  
  train<- imp_data[which(part!=i),]
  test<- imp_data[which(part==i),]
  
  model1<-svm(label~., data= train, kernel="radial", degree=7)
  model2<-randomForest(label~., data= train)
  model3 <- train(label~.,
                  train, 
                  method = "glmnet",
                  tuneGrid = expand.grid(alpha= seq(0,1,length=10),
                                         lambda= seq(0.0001, 1, length = 5)))
  
  pred1<-predict(model1, test)
  pred2<-predict(model2, test)
  pred3 <- predict(model3, test)
  tab1<-table(Predicted= pred1, Actual_svm= test$label)
  tab2<-table(Predicted= pred2, Actual_rf= test$label)
  tab3<-table(Predicted= pred3, Actual_en= test$label)
  
  metric_table[(3*(i-1)+1),]<-metric(tab1)
  metric_table[(3*(i-1)+2),]<-metric(tab2)
  metric_table[(3*(i-1)+3),]<-metric(tab3)
  
  print(tab1)
  print(tab2)
  print(tab3)
  
}
metric_table[16,]<-colMeans(metric_table[seq(1,15,3),])
metric_table[17,]<-colMeans(metric_table[(seq(1,15,3)+1),])
metric_table[18,]<-colMeans(metric_table[(seq(1,15,3)+2),])
t2<-metric_table
##write.csv(t2, "Lasso_metric.csv")

#---------------------------------------------------------------------------------------------------------


Marker<-qpcR:::cbind.na(RF_marker, lasso_marker)
##write.csv(Marker, "Markers.csv", row.names = F)
#markers<-read.csv("Markers.csv", header = T)
#markers<-unique(unlist(c(markers)))[!is.na(unique(unlist(c(markers))))]
markers<-unique(c(RF_marker, lasso_marker))


imp_data_comb<- ml_data[,c(markers, "label")]
#Classification with combined  marker
imp_data<-imp_data_comb

for(i in 1:5){
  
  train<- imp_data[which(part!=i),]
  test<- imp_data[which(part==i),]
  
  model1<-svm(label~., data= train, kernel="radial", degree=7)
  model2<-randomForest(label~., data= train)
  model3 <- train(label~.,
                  train, 
                  method = "glmnet",
                  tuneGrid = expand.grid(alpha= seq(0,1,length=10),
                                         lambda= seq(0.0001, 1, length = 5)))
  
  pred1<-predict(model1, test)
  pred2<-predict(model2, test)
  pred3 <- predict(model3, test)
  tab1<-table(Predicted= pred1, Actual_svm= test$label)
  tab2<-table(Predicted= pred2, Actual_rf= test$label)
  tab3<-table(Predicted= pred3, Actual_en= test$label)
  
  metric_table[(3*(i-1)+1),]<-metric(tab1)
  metric_table[(3*(i-1)+2),]<-metric(tab2)
  metric_table[(3*(i-1)+3),]<-metric(tab3)
  
  print(tab1)
  print(tab2)
  print(tab3)
  
}
metric_table[16,]<-colMeans(metric_table[seq(1,15,3),])
metric_table[17,]<-colMeans(metric_table[(seq(1,15,3)+1),])
metric_table[18,]<-colMeans(metric_table[(seq(1,15,3)+2),])
t3<-metric_table
##write.csv(t3, "comb_metric.csv")



lasso_marker_info<-info[c(as.numeric(unlist(lapply(strsplit(lasso_marker, "x"), function(x){x[2]})))), ]
Rf_marker_info<-info[c(as.numeric(unlist(lapply(strsplit(RF_marker, "x"), function(x){x[2]})))), ]
#write.csv(lasso_marker_info, "lasso_marker_info.csv", row.names = F)
#write.csv(Rf_marker_info, "rf_marker_info.csv", row.names = F)
#write.csv(rbind(lasso_marker_info, Rf_marker_info), "combine_marker_info.csv", row.names = F)


library(e1071)
?svm()
