#' CESur1
#'
#' CESur1 is used for predicting the survival rate of lung adenocarcinoma patient
#' @param x a data.frame containing GeneID and gene expression of a single patient
#'
#' @import GSVA
#'
CESur1<-function(x){
  rownames(x)<-x[,1]
  x<-merge(x,C2,by="GeneID")
  gene_list= split(C1$ENTREZID,C1$Celltype)
  rownames(x)<-x[,1]
  x<-gsva(as.matrix(x[,-1]), gene_list,method='ssgsea',
          kcdf='Gaussian',abs.ranking=TRUE)
  y1<-(1-0.128)^exp(0.30016229*x[,1])
  y3<-(1-0.455)^exp(0.30016229*x[,1])
  y5<-(1-0.876)^exp(0.30016229*x[,1])
  p1<-c(y1,y3,y5)
  p2<-c("1","3","5")
  S<-data.frame(p2,p1)
  colnames(S)<-c("Year","Survival probability")
  S<-as.data.frame(S)
  S
}
