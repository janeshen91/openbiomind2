#normalizing the alveolar macrophage data
setwd("/home/janes/Documents/of pigs and men/controls")
GSE45145<-read.csv(file="GSE45145_series_matrix.txt",header=T,sep="\t")
GSE45145am<-GSE45145[,c(98:141)]
GSE45145am<-GSE45145am[-c(1:12),]
firstGSE45145am<-GSE45145[-c(1:12),1]
rownames(GSE45145am)<-firstGSE45145am

for(i in 1:44){
j=i+97
colnames(GSE45145am)[i]<-as.character(GSE45145[1,j])
}


GSE45145am<-GSE45145am[-1,]
#binarizedMedian
binarizeMedianFunc<-function(dataset){
#assumes the dataset is already of controls/samples we want only
#make numeric
datasetNum<-apply(dataset,2,function(x){as.numeric(x)})
medianDatasetNum<-median(datasetNum,na.rm=T)
binarizedDataset<-apply(datasetNum,2,function(x){ifelse(x>medianDatasetNum,1,0)})
rownames(binarizedDataset)<-rownames(dataset)
return(binarizedDataset)
}
 
binarizedMedianGSE45145<-binarizeMedianFunc(GSE45145am)

GSE13896<-read.csv(file="GSE13896_series_matrix.txt",header=T,sep="\t")
GSE13896C<-GSE13896[,c(2:12,25:37)]
rownames(GSE13896C)<-GSE13896[,1]

binarizedMedianGSE13896<-binarizeMedianFunc(GSE13896C)

quantilize<-function(dataset){
datasetM<-as.matrix(dataset)
qn<-normalizeBetweenArrays(datasetM,method="quantile")
}
library(limma)
GSE45145am2<-apply(GSE45145am,2,function(x){as.numeric(x)})
qn45145<-quantilize(GSE45145am2)
rownames(qn45145)<-rownames(GSE45145am)

qn13896<-quantilize(GSE13896C)

