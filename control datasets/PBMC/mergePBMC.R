setwd("/home/janes/Documents//of pigs and men/PBMC/")
load("/home/janes/Documents/of pigs and men/controls/PBMC.RData")
#because I want to, I'm going to try doing everythign in data.table
library(data.table)
GSM433321<-read.csv(file="/home/janes/Documents/of pigs and men/controls/gse17320/GSM433321_13685770-5622-CONCY3-LPSCY5-900790.gpr",as.is=T,sep=" ")
setwd("/home/janes/Documents/openbiomind2//openbiomind2/getgeo/")
GSMfiles<-list.files(path="/home/janes/Documents/of pigs and men/controls/gse17320out/",pattern="GSM")
#if even number then 5, if odd number(in the list) then 3 (up to 14)
#15 is 5, 16 is 3 and from 16 onwards even number is 3
#this changes the column number we want to get the data from
# ie 46 for 635(cy5) and 47 for 532 (cy3)
#reading files in and then extracting the right column
#assuming all the rows are in the same order
setwd("/home/janes/Documents/of pigs and men/controls/gse17320out/")
GSMfileValues<-lapply(GSMfiles,function(x){read.csv(file=x)})
GSMfileValues2<-matrix(nrow=nrow(GSMfileValues[[1]]),ncol=length(GSMfiles))
GSMfileValues3<-matrix(nrow=1,ncol=length(GSMfiles))
for (i in 1:length(GSMfiles)){
  if(grep("CONCY3",GSMfiles[1])){
    GSMfileValues2[,i]<-GSMfileValues[[i]]$X.2
  }
  else  if (grep("CONCY5",GSMfiles[1])){
    GSMfileValues2[,i]<-GSMfileValues[[i]]$X.1
  }
  if(i>1){
    GSMfileValues3[i]<-identical(GSMfileValues[[i]]$X,GSMfileValues[[i-1]]$X)}
}
colnames(GSMfileValues2)<-GSMfiles
GSMfileValues2<-GSMfileValues2[,-1]
#ok, so adding in the rest of the pbmc data now
load("/home/janes/Documents/of pigs and men/controls/PBMC.RData")
#adding in the 
rownames(GSMfileValues2)<-GSMfileValues[[1]]$X
GSE17320<-GSMfileValues2
gpl7151<-read.csv(file="/home/janes/Documents//of pigs and men/platforms/GPL7151.csv")
library(limma)
logGSE17320<-apply(GSE17320,2,function(x){log2(x)})
> summary(GSE17320)
GSM433321_13685770-5622-CONCY3-LPSCY5-900790.gpr~
  Min.   : -487.0                                  
1st Qu.:  -15.0                                  
Median :  139.0                                  
Mean   :  875.4                                  
3rd Qu.:  543.2                                  
Max.   :65233.0                                  
#IE note, about 1/4 of all values are below zero, but since it's the value is the foreground
#subtracting the background, it seems reasonable to want it to be above zero

quantilize<-function(dataset){
  datasetM<-as.matrix(dataset)
  qn<-normalizeBetweenArrays(datasetM,method="quantile")
}
qnGSE17320<-quantilize(logGSE17320)

#i guess just read in the GPL3533 AND GPL1
GPL3533<-read.csv(file="/home/janes/Documents//of pigs and men/platforms/GPL3533.txt",header=T,skip=27,sep="\t")
GPL7151<-read.csv(file="/home/janes/Documents//of pigs and men/platforms/GPL7151.csv")

#read.table(text=sub(paste0("[!#]", ".*"), "", readLines("test.tab")), header=TRUE)
geneIdMatch<-match(GPL7151$Gene_IDs,GPL3533$Gene.ID)



#mergeGPL7151WithGPL3533<-merge(GPL7151[,c(1,2)],GPL3533[,c(1,4)],by.x="Gene_IDs",by.y="Gene.ID")
#ok, we're going to merge, but by using the column name of thing we want
#ie we want the ensemble number of the probes in GPL7151
gpl7151WithMatch<-subset(GPL7151,GPL7151[,"Gene_IDs"]==geneIdMatch)
head(geneIdMatch)
library(data.table)
GPL7151DT<-data.table(GPL7151)
GPL3533DT<-data.table(GPL3533[,c(1,4)])
setkey(GPL7151DT,cols="Gene_IDs")
setkey(GPL3533DT,cols="Gene.ID")
new<-GPL7151DT[GPL3533DT,nomatch=0]
# > new<-GPL7151DT[GPL3533DT,nomatch=0]
# Error in vecseq(f__, len__, if (allow.cartesian) NULL else as.integer(max(nrow(x),  : 
# Join results in 247163550 rows; more than 141301 = max(nrow(x),nrow(i)). Check for duplicate key values in i, each of which join to the same group in x over and over again. If that's ok, try including `j` and dropping `by` (by-without-by) so that j runs for each group to avoid the large allocation. If you are sure you wish to proceed, rerun with allow.cartesian=TRUE. Otherwise, please search for this error message in the FAQ, Wiki, Stack Overflow and datatable-help for advice.
# #i think it's all the empty space that's causing trouble
X<-(GPL7151[,3]=="")
GPL7151ii<-GPL7151[!X,]
y<-(GPL3533[,4]=="")
GPL3533ii<-GPL3533[!y,]
GPL7151DT<-data.table(GPL7151ii)
GPL3533DT<-data.table(GPL3533ii[,c(1,4)])
setkey(GPL7151DT,cols="Gene_IDs")
setkey(GPL3533DT,cols="Gene.ID")
gpl7151togpl3533<-GPL7151DT[GPL3533DT,nomatch=0]
rm(new)
#now we can join stuff! like the 3 pig datasets!!
yes<-rownames(qnGSE23503)%in%gpl7151togpl3533$ID
sum(yes)
[1] 5161#same for GSE30874
#WE need to first do the making it between 0 and 1 bit
#ie take it, generate the min and max, and then put it between 0 and 1
is.numeric(qnGSE23503)
min2<-min(qnGSE23503,na.rm=T)
max2<-max(qnGSE23503,na.rm=T)
qnGSE23503ii<-apply(qnGSE23503,2,function(x){(x-min2)/(max2-min2)})

min2<-min(qnGSE30874,na.rm=T)
max2<-max(qnGSE30874,na.rm=T)
qnGSE30874ii<-apply(qnGSE30874,2,function(x){(x-min2)/(max2-min2)})

min2<-min(qnGSE17320,na.rm=T)
max2<-max(qnGSE17320,na.rm=T)
qnGSE17320ii<-apply(qnGSE17320,2,function(x){(x-min2)/(max2-min2)})
rownames(qnGSE17320ii)<-GSMfileValues[[1]]$X
qnGSE23503DT<-data.table(qnGSE23503ii,keep.rownames=T)
qnGSE30874DT<-data.table(qnGSE30874ii,keep.rownames=T)
qnGSE17320DT<-data.table(qnGSE17320ii,keep.rownames=T)
head(qnGSE23503DT)
setkey(qnGSE23503DT,rn)
setkey(qnGSE30874DT,rn)
setkey(qnGSE17320DT,rn)
bothgpl3533<-qnGSE23503DT[qnGSE30874DT]
bothgpl3533<-bothgpl3533[-1,]
#show
# uniqgpl7151togpl3533<-unique(gpl7151togpl3533)
# > dim(uniqgpl7151togpl3533)
# [1] 3430    4
# > dim(gpl7151togpl3533)
# [1] 7911    4
###data.table's unique by default compares only keys
# uniqueGPL7151<-unique(gpl7151togpl3533,by=NULL)
# > dim(uniqueGPL7151)
# [1] 7911    4
# > dim(gpl7151togpl3533)
# [1] 7911    4

gpl7151togpl3533II<-gpl7151togpl3533[,c(3:4),with=FALSE]
mapping<-qnGSE17320DT[gpl7151togpl3533II]
head(qnGSE17320DT)
dim(qnGSE17320DT)
#WHAT ARE THE differences between mapping and qnGSE17320dt?
head(qnGSE17320DT$rn)
x<-qnGSE17320DT$rn %in% qnGSE17320DT[,GPL7151]
# I think this happens to have more levels than the number of elements becuase
#of the way it was constructed: 
# l<-levels(gpl7151togpl3533[,GPL7151])
# > length(l)
# [1] 15479
# > length(gpl7151togpl3533[,GPL7151])
# [1] 7911

# x<- qnGSE17320DT[,rn]%in% mapping[,rn]
# > x<-mapping[,rn] %in% qnGSE17320DT[,rn]
# > head(x)
# [1] TRUE TRUE TRUE TRUE TRUE TRUE
# > sum(x)
# [1] 8453
# > x<- qnGSE17320DT[,rn]%in% mapping[,rn]
# > sum(x)
# [1] 5217
# > length(qnGSE17320DT[,rn])
# [1] 19200
# > length(unique(qnGSE17320DT[,rn]))
# [1] 17116
# > length(unique(gpl7151togpl3533II[,GPL7151]))
# [1] 5072
to<-merge(qnGSE17320DT,gpl7151togpl3533II)
key(qnGSE17320DT)
colnames(qnGSE17320DT)[1]<-"GPL7151"
setkey(gpl7151togpl3533II,GPL7151)
#i think what happens is that we're doing a left join? but merge is supposed
#to do an inner join by default and it might be the same?
setkey(to,"ID")
allpigPBMC<-to[bothgpl3533]
x<-bothgpl3533[,ID] %in% to[,ID]
x<-match(bothgpl3533[,ID],to[,ID])
allpigPBMC<-to[!is.na(x),]
#now pig and human
matchph<-read.csv(file="/home/janes/Documents/openbiomind2/openbiomind2/JanesProbeMatching/Matches/GPL7151AndGPL570.txt")
#make the human between 0 and 1
is.numeric(qnGSE13501)
min2<-min(qnGSE13501,na.rm=T)
max2<-max(qnGSE13501,na.rm=T)
qnGSE13501ii<-apply(qnGSE13501,2,function(x){(x-min2)/(max2-min2)})
head(qnGSE13501ii)
qnGSE13501ii<-cbind(rownames(qnGSE13501ii),qnGSE13501ii)
qnGSE13501III<-data.table(qnGSE13501III,
setnames(qnGSE13501III,"V1","humanProbe")
setkey(qnGSE13501III,humanProbe)
matchph<-data.table(matchph)
setkey(matchph,humanProbe)
qnGSE1301IV<-merge(qnGSE13501III,matchph)
setkey(qnGSE1301IV,pigProbe)     
setnames(allpigPBMC,"GPL7151","pigProbe")
setkey(allpigPBMC,pigProbe)
#remove non matched human rows first
allinQ<-qnGSE1301IV[,pigProbe] %in% allpigPBMC[,pigProbe]
allPBMC<-merge(qnGSE1301IV,allpigPBMC)
#neither hte pig or the human has been averaged by gene :( so redo and 
#average by gene?
#right, waht do we do instead of tapply?
qnGSE23503DT2<-qnGSE23503DT[-1,]
setnames(qnGSE23503DT2,"rn","ID")
#setkey(qnGSE23503DT2,key="rn")
setkey(GPL3533DT,key="ID")
qnGSE23503DT2i<-merge(qnGSE23503DT2,GPL3533DT)
colnames(meanQNGSE23503)[1]<-"Gene.ID"
qnGSE23503DT3<-qnGSE23503DT2i[,lapply(.SD,mean),by=Gene.ID,.SDcol=c(2:15)]
qnGSE23503DT4<-qnGSE23503DT3[,]
qnGSE23503DT<-data.table(qnGSE23503ii,keep.rownames=T)
qnGSE30874DT<-data.table(qnGSE30874ii,keep.rownames=T)
qnGSE17320DT<-data.table(qnGSE17320ii,keep.rowname
meanQNGSE23503<-apply(qnGSE23503DT3[,c(2:15),with=FALSE],1,mean)

setkey(allpigPBMC,ID)
setkey(GPL3533DT,ID)
setnames(GPL3533DT,"pigProbe","ID")
allpigPBMC2<-merge(allpigPBMC,GPL3533DT)
qnGSE23503DT3<-qnGSE23503DT2i[,lapply(.SD,mean),by=Gene.ID,.SDcol=c(2:15)]
dim(allpigPBMC2)
>[1] 8453   53 #first 2 columns are probe ids and so is the last column
allpigPBMC2I<-allpigPBMC2[,lapply(.SD,mean),by=Gene.ID,.SDcol=c(3:52)]
write.table(allpigPBMC2,file="allpigByGene.txt",quote=F)
#the human by gene
head(qnGSE13501)
library(hgu133plus2.db)
columns(hgu133plus2.db)
GSE13501Gene<-select(hgu133plus2.db,keys=qnGSE13501ii[,1],columns="ENTREZID")
GSE13501Genedt<-data.table(GSE13501Gene)
setkey(GSE13501Genedt,PROBEID)
library(data.table)
qnGSE13501IIdt<-data.table(qnGSE13501ii)
setkey(qnGSE13501IIdt,PROBEID)
qnGSE13501iidt<-GSE13501Genedt[qnGSE13501IIdt]
GSE13501Gene2<-merge(qnGSE13501ii,GSE13501Gene,by='PROBEID')
qnGSE13501WITHGENE1<-qnGSE13501iidt[c(1:20),c(2:4),with=FALSE]
qnGSE13501byGene<-qnGSE13501WITHGENE1[,lapply(.SD,mean), .SDcols=c(2:3), by=ENTREZID]
qnGSE13501byGene<-qnGSE13501WITHGENE[,lapply(.SD,function(x){mean(x,na.rm=T)}),by=ENTREZID,.SDcols=c(3:60)]
#removing the NAs from...altho why are there NAs? we shouldn't have merged but should've used
#a left join on GSE13501Gene instead?