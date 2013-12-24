binarizeMedianFunc<-function(dataset){
#assumes the dataset is already of controls/samples we want only
#make numeric
datasetNum<-apply(dataset,2,function(x){as.numeric(x)})
medianDatasetNum<-median(datasetNum,na.rm=T)
binarizedDataset<-apply(datasetNum,2,function(x){ifelse(x>medianDatasetNum,1,0)})
rownames(binarizedDataset)<-rownames(dataset)
return(binarizedDataset)
}
