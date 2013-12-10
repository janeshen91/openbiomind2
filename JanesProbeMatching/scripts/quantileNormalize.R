 setwd("/home/janes/Documents/of pigs and men/controls")
 load("Nov26stuff.RData")

source("http://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
GSE28871M<-as.matrix(GSE28871C)
 GSE28871cq<-normalizeBetweenArrays(GSE28871M,method="quantile")
#looks the same so maybe already quantile normalized...check?
 GSE12194<-apply(GSE12149,2,function(x){as.numeric(x)})
GSE12194M<-as.matrix(GSE12194)
 qnGSE12194<-normalizeBetweenArrays(GSE12194M,method="quantile")
#
quantilize<-function(dataset){
datasetM<-as.matrix(dataset)
qn<-normalizeBetweenArrays(datasetM,method="quantile")
}
qnGSE47460GPL14550<-quantilize(GSE47460GPL14550nUM)
qnGSE47460GPL6480<-quantilize(GSE47460GPL6480)

save(qnGSE47460GPL14550,qnGSE47460GPL6480,qnGSE12194,qnGSE28871,file="quantileNormalizedLung.RData")
save.image("lung.RData")



GSM1097985 	
GSM1098028 	
