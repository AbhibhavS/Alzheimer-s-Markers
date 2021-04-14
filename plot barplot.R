


dir<-c("C:/Protein/AD Gene Expression/Abhibhav/PFC/results",
       "C:/Protein/AD Gene Expression/Abhibhav/Temproral/results",
       "C:/Protein/AD Gene Expression/Abhibhav/Hippo/results",
       "C:/Protein/AD Gene Expression/Abhibhav/EC/results"
)
brain<-c("Prefrontal cortex", "Middle Temporal Gyrus", 
         "Hippocampus",  "Entorhinal Cortex")

par(mfrow=c(2,2))

for(i in 1:4){
  setwd(dir[i])
  lasso<-read.csv("Lasso_metric.csv", row.names = 1, header = T)
  varsel<-read.csv("RFmethod_metric.csv", row.names = 1, header = T)
  comb<-read.csv("comb_metric.csv", row.names = 1, header = T)
  
  lasso_svm<-lasso[seq(1,18,3),]
  lasso_RF<-lasso[c(seq(1,18,3)+1),]
  lasso_eln<-lasso[c(seq(1,18,3)+2),]
  
  varsel_svm<-varsel[seq(1,18,3),]
  varsel_RF<-varsel[c(seq(1,18,3)+1),]
  varsel_eln<-varsel[c(seq(1,18,3)+2),]
  
  comb_svm<-comb[seq(1,18,3),]
  comb_RF<-comb[c(seq(1,18,3)+1),]
  comb_eln<-comb[c(seq(1,18,3)+2),]
  
  
  data<-c(varsel_svm[6,1], varsel_RF[6,1], varsel_eln[6,1],
          lasso_svm[6,1], lasso_RF[6,1], lasso_eln[6,1],
          comb_svm[6,1], comb_RF[6,1], comb_eln[6,1])
  data<-matrix(data, 3,3,byrow=F)
  
  colnames(data)<-c("VarSelRF", "Lasso", "VarSelRF + Lasso")
  
  
  barplot(data, beside = T,
          xaxt="n", yaxt="n",
          ylim=c(0.75,1.02), xlim=c(2,20),
          space = c(0.01,3),
          border = c("firebrick1","darkgoldenrod1", "deepskyblue1"),
          col = c("red","darkgoldenrod1", "deepskyblue3"),
          xpd=F)
  
  abline(h=seq(0.75,1.02,0.01), col="grey", lwd=2)
  par(new=T)
  barplot(data, beside = T,
          xaxt="n", yaxt="n",
          ylim=c(0.75,1.02), xlim=c(2,20),
          space = c(0.01,3),
          border = c("firebrick1",
                     "darkgoldenrod1", "deepskyblue4"),
          col = c("red","darkgoldenrod1", "deepskyblue4"),
          xpd=F,
          main = brain[i])
  axis(2, at = seq(1,0.75,-0.05), las=1, font = 9, tick=F,
       labels = round(seq(1,0.75,-0.05),2), col = "black", cex.axis=1.4)
  axis(1, at = c(5,10.5,17.1), las=1, font = 2,
       labels = c("VarSelRF", "LASSO", 
                  "VarSelRF+LASSO"), 
       col = "black", 
       tick=F, cex.axis=1.3)
  
  text<-c("SVM","RandomForest", "ElasticNet")
  legend(#"topleft",
         x = c(1, 16), y = c(1, 1.035),
         legend =text, 
         fill = c("red","darkgoldenrod1", "deepskyblue4"),
         text.width = 1,
         cex=1,
         text.font=9,
         bty="n",
         box.lwd=2,
         ncol = 1)
  
  box(lwd=2)
}










