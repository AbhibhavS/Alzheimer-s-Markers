setwd("C:/Protein/AD Gene Expression/Abhibhav/EC")
set.seed(123)

library(GEOquery)
library(e1071)
library(randomForest)
library(varSelRF)
library(caret)
library(qpcR)

GPL570<- read.csv("GPL570.csv", header =T)

gse1=getGEO(filename="GSE5281_series_matrix.txt.gz",
            GSEMatrix=TRUE,getGPL = TRUE)
GSE5281<- as.matrix(gse1)
GSE5281<-as.data.frame(GSE5281)
colnames(GSE5281)<- c(read.csv("sample info_GSE5281.csv",
                               header =T))[[1]]

gse2=getGEO(filename="GSE48350_series_matrix.txt.gz",
            GSEMatrix=TRUE,getGPL = TRUE)
GSE48350<- as.matrix(gse2)
GSE48350<-as.data.frame(GSE48350)
colnames(GSE48350)<- c(read.csv("sample info_GSE48350.csv",
                                header =T))[[1]]

gse3=getGEO(filename="GSE4757_series_matrix.txt.gz",
            GSEMatrix=TRUE,getGPL = TRUE)
GSE4757<- as.matrix(gse3)
GSE4757<-as.data.frame(GSE4757)
colnames(GSE4757)<- c(read.csv("sample info_GSE4757.csv",
                                header =T))[[1]]

GSE_GPL570<- cbind(GSE5281, GSE48350, GSE4757)
colnames(GSE_GPL570)

GSE_EC<-GSE_GPL570[,c(grep("Entorhinal", colnames(GSE_GPL570)),
                       grep("entorhinal", colnames(GSE_GPL570)),
                       grep("EC control", colnames(GSE_GPL570)),
                       grep("EC_", colnames(GSE_GPL570)))]
colnames(GSE_EC)
##write.csv(GSE_GPL570, "GSE_GPL570.csv")

zscore<-function(X){(X- mean(X))/sqrt(var(X))} #standardization function
GSE_EC<-apply(GSE_EC, MARGIN = 2, zscore) #standardized expression data
GSE_EC<-as.data.frame(GSE_EC)


#dir.create("results")
setwd(paste(getwd(), "/results", sep=""))

#--------------------------------------------------------------------------------------------------------#

gene_name<-as.data.frame(cbind(c(1:nrow(GSE_EC)), rownames(GSE_EC))) #Probe Id
colnames(gene_name)<- c("index", "gene_name")
gene_id<- as.data.frame(cbind(GPL570[,1:2], GPL570$Gene.Symbol ,GPL570$ENTREZ_GENE_ID, GPL570$RefSeq.Transcript.ID)) #meta informations
colnames(gene_id)<- c("gene_name", "ACC_id", "Gene_symbol", "entrez_Gene_ID", "ref_Gene_ID")

vec <- merge(x = gene_name, y = gene_id, by="gene_name", all = TRUE) #prob id mapped with gene ID and other meta informations


no_info <- which(vec$entrez_Gene_ID == "" & vec$Gene_symbol == "" & vec$ref_Gene_ID == "") #avoiding no informative probes
info <- vec[-no_info,] #informative probes

uni<-which(!(duplicated(info[,4:6]))) #all index of unique probes

#---------------------------------------------------------------------------------------------------------#


for(i in uni){
  
  if(is.na(info$entrez_Gene_ID[i])){ next }
  else if(length(which(info$Gene_symbol==info$Gene_symbol[i] & info$entrez_Gene_ID==info$entrez_Gene_ID[i]))==1){ next }
  else{
    flag<-which(info$Gene_symbol==info$Gene_symbol[i]  & info$entrez_Gene_ID==info$entrez_Gene_ID[i] ) #identifying duplicate if exist
    x<-apply(GSE_EC[flag,], MARGIN = 1, IQR) #max quantiles are selected
    info[flag[which(x!=max(x))],]<-NA # rest duplicate entries are set to NA in the meta data
    GSE_EC[flag[which(x!=max(x))],]<-NA # rest duplicate entries are set to NA in the gene expression data
  }
}
GSE_EC[no_info,]<-NA #uniformative probes are set to NA in gene expression data
meta_data<-vec[which(!is.na(GSE_EC[,1])),] # filtering out the NA data from the meta data
GSE<-GSE_EC[which(!is.na(GSE_EC[,1])), ] # filtering out the NA data from the expression data
GSE<-cbind(meta_data, GSE) #combined meta data and expression data
##write.csv(GSE, "GSE_EC_normalized_processed.csv")

#-----------------------------------------------------------------------------------------------

GSE_data<- as.data.frame(GSE[, 7:ncol(GSE)]) #only choosing the expression data #ncol(GSE_data)


GSE_control<- c(grep("control", colnames(GSE_data)),
                grep("Control", colnames(GSE_data))) #index of control samples
GSE_control
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

my_GSE<- cbind(GSE[,1:6], GSE_data, pvalue_t, pvalue_wlx) #meta data + expression data + p values
my_GSE<- as.data.frame(my_GSE[order(my_GSE$pvalue_t, decreasing = F),]) #ordering the data by t pvalue
my_GSE$pvalue_t<- as.numeric(my_GSE$pvalue_t)
my_GSE$pvalue_wlx<- as.numeric(my_GSE$pvalue_wlx)
my_GSE$adjust<-p.adjust(my_GSE$pvalue_t, method = "bonferroni") #FDR test

Filter_GSE<-my_GSE[which(my_GSE$pvalue_t <= 0.05 
                         & my_GSE$pvalue_wlx <= 0.05
                         & my_GSE$adjust <= 0.01),] #choosing the appropriate gene expression data

##write.csv(Filter_GSE, "GSE_EC_ML.csv", row.names = F)


#--------------------------------------------------------------------------------------------------------#


data<-Filter_GSE[,7:(ncol(Filter_GSE)-3)] #choosing just the expression data
info<-Filter_GSE[,1:6] #choosing the meta data
#write.csv(t(info), "meta_data_ml_data.csv", row.names = T)


data<-t(data) #transposed to make sample X gene format

label<-rep(NA, length(row.names(data))) #label vector
label[GSE_control]<- 0 #normal is denoted a 0
label[is.na(label)]<- 1 #affected is denoted a 1

ml_data<-as.data.frame(cbind(data, label)) #machine leaning compatible data
colnames(ml_data)<- c( paste("x",c(1:(ncol(ml_data)-1)), sep=""), "label") #legal naming
#shuffle<-sample(1:nrow(ml_data))
#write.table(shuffle, "seed.txt", row.names = F)
shuffle<-c(read.table("C:/Protein/AD Gene Expression/Abhibhav/EC/seed.txt", header = T))[[1]] #beacuse my stupid ass has forgot to put seed

ml_data<-ml_data[shuffle,] #shuffling randomly
ml_data$label<-as.factor(ml_data$label) #label factorized

#write.csv(ml_data, "GSE_GPL_Filtered_with_label_.csv", row.names=T)
##write.csv(info, "meta_data_ml_data.csv", row.names = F)

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
for(i in 1:5){
  a<-varSelRF(ml_data[,-ncol(ml_data)], ml_data[,ncol(ml_data)], whole.range = F, verbose = T) #random forest feature selection method
  vars_new<-sort(a$selected.vars) # best set of variables
  vars<-qpcR:::cbind.na(vars, vars_new)
}
vars[which(is.na(vars))]<-0
RF_marker<-rownames(as.matrix(table(vars)[-1]))[which(as.matrix(table(vars)[-1])>=2)]
##write.table(table(vars)[-1], "GSE48350_markers_RF.txt", row.names = F)
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

############ LASSO Method####################

custom <- trainControl(method = "repeatedcv",
                       number = 10,
                       repeats = 5,
                       verboseIter = T)

lasso <- train(label~.,
               ml_data, 
               method = "glmnet",
               tuneGrid = expand.grid(alpha=1,
                                      lambda= seq(10^(-10), 1, length = 5)),
               trControl = custom)

imp<-varImp(lasso, scale=T)
imp<-as.matrix(imp$importance)
lasso_marker<-names(imp[which(imp>=10),])
#imp<-sort(imp[,1], decreasing = T)
#lasso_marker<-names(imp[1:20])

##write.table(lasso_marker, "GSE48350_markers_lasso.txt", row.names = F)
imp_data_lasso<- ml_data[,c(lasso_marker, "label")]

#Classification with LASSO found marker
#imp_data<-imp_data_lasso
imp_data<-imp_data_lasso

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
write.csv(Marker, "Markers.csv", row.names = F)
markers<-read.csv("Markers.csv", header = T)
markers<-unique(unlist(c(markers)))[!is.na(unique(unlist(c(markers))))]

imp_data_comb<- ml_data[,c(markers, "label")]
#Classification with combined  marker
#imp_data<-imp_data_comb
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

t1









