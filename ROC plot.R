

#----------------------------ROC-----------------------------#
library(ROCR)
library(e1071)
library(randomForest)
library(varSelRF)
library(caret)
library(pROC)


par(mfrow=c(2,2))

dir<-c("C:/Protein/AD Gene Expression/Abhibhav/PFC/results",
       "C:/Protein/AD Gene Expression/Abhibhav/Temproral/results",
       "C:/Protein/AD Gene Expression/Abhibhav/Hippo/results",
       "C:/Protein/AD Gene Expression/Abhibhav/EC/results"
       )
brain<-c("Prefrontal cortex", "Middle Temporal Gyrus", 
         "Hippocampus",  "Entorhinal Cortex")


par(mfrow=c(2,2))
for(j in 1:4){
  
  setwd(dir[j])
  
  data<-read.csv("GSE_GPL_Filtered_with_label_.csv", header = T, row.names = 1)
  data$label<-as.factor(data$label)
  
  markers<-read.csv("Markers.csv", header = T)
  markers<-unique(unlist(c(markers)))[!is.na(unique(unlist(c(markers))))]
  
  data_ml<-data[,c(markers,"label")]
  part<- sort(rep(c(1:5), nrow(data_ml)/5))
  
  pred1<-list(NULL)
  pred2<-list(NULL)
  pred3<-list(NULL)
  
  for(i in 1:5){
    
    train<- data_ml[which(part!=i),]
    test<- data_ml[which(part==i),]
    
    model1<-svm(label~., data= train, kernel="radial", degree=7, probability=T)
    model2<-randomForest(label~., data= train)
    model3 <- train(label~.,
                    train, 
                    method = "glmnet",
                    tuneGrid = expand.grid(alpha= seq(0,1,length=10),
                                           lambda= seq(0.0001, 1, length = 5)))
    
    pred<-predict(model1, test, probability = T, type="response")
    pred<-attr(pred, "probabilities")
    pred1[[i]]<- pred[,1]
    
    pred<-predict(model2, test, type = "prob")
    pred2[[i]]<-pred[,1]
    
    pred<-predict(model3, test, type="prob")
    pred3[[i]]<-pred[,1]
    #roc((as.numeric(test$label)-1), 
        #pred1[,1], plot = T)
    #pred1[[i]]<-predict(model1, test, type="prob")
    #pred2[[i]]<-predict(model2, test)
    #pred3[[i]] <- predict(model3, test)
  }
  
  svm_pred<-unlist(pred1)
  lab1<-data_ml$label[1:length(svm_pred)]
  roc_svm<-roc(lab1,svm_pred, percent = T)
  plot(roc_svm, 
       grid=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       main= brain[j],
       xaxt="n", yaxt="n",
       cex.lab=1.3,
       font.lab=9,
       print.auc.y=40,
       print.auc.x=24, 
       print.auc=F,
       legacy.axes=T,
       percent = T,
       col="red",
       lwd=3)
  
  rf_pred<-unlist(pred2)
  lab2<-data_ml$label[1:length(rf_pred)]
  roc_rf<-roc(lab2, rf_pred,percent=T)
  plot(roc_rf,
       add=T,
       xlab="False Positive Rate",
       ylab="True Positive Rate",
       percent = T,
       legacy.axes=T, 
       print.auc.y=30,
       print.auc.x=24, 
       print.auc=F,
       col="darkgoldenrod1", 
       lwd=3)
 
  #par(new=T)
  en_pred<-unlist(pred3)
  lab3<-data_ml$label[1:length(en_pred)]
  roc_en<-roc(lab3, en_pred, percent=T)
  plot(
    roc_en,
    add=T,
    xlab="False Positive Rate",
    ylab="True Positive Rate",
    #lab.cex=5,
    percent = T,
    print.auc.y=20,
    print.auc.x=24, 
    print.auc=F,
    legacy.axes=T, 
    col="deepskyblue4", 
    lwd=4)
  
  text<-c(paste("SVM","(AUC=", as.character(round(roc_svm$auc,2)),"%)"),
          paste("RandomForest","(AUC=", as.character(round(roc_rf$auc,2)),"%)"),
          paste("ElasticNet","(AUC=", as.character(round(roc_en$auc,2)),"%)"))
  
  legend("bottomright",
         legend = text,
         lty=1,
         lwd=4,
         cex=0.8,
         col = c("red","darkgoldenrod1", "deepskyblue4"),
         bty="n",
         text.font = 9)
         
  axis(1, at=seq(100,0,-10), labels =seq(0,1,0.1), tick=F, font=9, cex.axis=1.2)
  axis(2, at=seq(0,100,10), labels =seq(0,1,0.1), tick=F, font=9, cex.axis=1.2)
}

