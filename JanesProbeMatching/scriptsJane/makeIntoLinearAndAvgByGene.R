makeIntoLinearAndAvgByGene<-function(qnDataset,platform){
  max<-max(qnDataset,na.rm=T)
  min<-min(qnDataset,na.rm=T)
  qnDatasetL<-apply(qnDataset,2,function(x){(x-min)/(max-min)})
  qnDataset2<-cbind(rownames(qnDataset),qnDatasetL)
  colnames(qnDataset2)[1]<-"probeID"
  #assuming the platform is 2 columns, probes and names
  namesOfprobe<-colnames(platform)[1]
  qnDataset3<-merge(platform,qnDataset2,by.y="probeID",by.x=namesOfprobe)
  #now we do the tapply
  lengthm<-length(qnDataset3)
  qnDataset4<-apply(qnDataset3[,c(2:lengthm)],2,function(x){as.numeric(x)})
  qnDataset5<-apply(qnDataset4,2,function(x){
   tapply(x,qnDataset4[,1],function(y){mean(y,na.rm=T)}) 
    
  } )
return(qnDataset5)   
}
