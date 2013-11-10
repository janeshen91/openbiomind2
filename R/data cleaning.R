####  code for setting up data sets for pig lung transplant project  ####

## download and set up selected data sets

series <- c("GSE1133",  # human reference
"GSE13134", "GSE28871", "GSE12194", "GSE2335", "GSE2339",   # pig lung # (missing annotation data)
"GSE47460",  # human lung  # (missing annotation data)
"GSE45145",  # pig alveolar macrophage  # (missing annotation data)
"GSE13896",  # human alveolar macrophage
"GSE17320", "GSE23503", "GSE30874",  # pig peripheral blood monocytes # (first missing annotation data)
"GSE13501",  # human PBMC
"GSE8021", "GSE9102",  # human transplant lung tissue
"GSE2018"  # human transplant bronchiolar lavage samples
)
gsenames <- c("href", "plung1", "plung2", "plung3", "plung4", "plung5", "hlung", "palmac", "halmac", "ppbmc1", "ppbmc2", "ppbmc3", "hpbmc", "hlungtx1", "hlungtx2", "hlavtx")

# # check GEOmetadb for platforms used in selected series.  necessary for determining if annotation files exist for automatic construction of annotated expression set objects.
piggpl <- c("GPL1881", "GPL6173", "GPL10162", "GPL16569", "GPL7151", "GPL3533") 

# geoConvert(series, out_type = "gpl", sqlite = 'C:/Users/user/Documents/GitHub/pigRpeople/GEOmetadb.sqlite') 
# $gpl
#    from_acc   to_acc
# 1   GSE1133  GPL1073
# 2   GSE1133  GPL1074
# 3   GSE1133    GPL96
# 4  GSE12194  GPL1881   #  causes error in getGEO(... AnnotGPL = TRUE)
# 5  GSE13134  GPL6173   #  causes error in getGEO(... AnnotGPL = TRUE)
# 6  GSE13501   GPL570
# 7  GSE13896   GPL570
# 8  GSE17320  GPL7151   #  causes error in getGEO(... AnnotGPL = TRUE)
# 9   GSE2018    GPL96
# 10  GSE2335  GPL1881   #  causes error in getGEO(... AnnotGPL = TRUE)
# 11  GSE2339  GPL1881   #  causes error in getGEO(... AnnotGPL = TRUE)
# 12 GSE23503  GPL3533
# 13 GSE28871 GPL10162   #  causes error in getGEO(... AnnotGPL = TRUE)
# 14 GSE30874  GPL3533
# 15 GSE45145 GPL16569   #  causes error in getGEO(... AnnotGPL = TRUE)
# 16 GSE47460 GPL14550   #  causes error in getGEO(... AnnotGPL = TRUE)
# 17 GSE47460  GPL6480
# 18  GSE8021  GPL5356
# 19  GSE9102  GPL5909   #  causes error in getGEO(... AnnotGPL = TRUE)

ann_file <- c(T,F,F,F,F,F,F,F,T,F,T,T,T,T,F,T)  # eliminates errors when downloading vector of GEO objects with GEOquery

library("GEOquery")

#  generate list of series matrix expression sets and store compressed downloaded files in directory "getgeo" in pwd.
gsematrix <- list(NULL)
length(gsematrix) <- 16
names(gsematrix) <- gsenames
for(i in 1:16){gsematrix[[i]] <- getGEO(GEO = series[i], destdir = 'getgeo', AnnotGPL = ann_file[i])}

## apply median norm to object

med.normalize <- function(mat) {
  out <- mat
  for (i in seq(dim(mat)[2])) { 
    vect <- mat[,i]
    med <- median(vect, na.rm = TRUE)
    out[,i] <- as.numeric(vect >= med)
  }
  return(out)
}

