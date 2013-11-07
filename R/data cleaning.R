####  code for setting up data sets for pig lung transplant project  ####

## set up GEOmetadb GEO metadata database for GEO data exploration  --  useful for screening GEO data sets
##  skip down to the next double hashmarks if you just want the code to download the data for the expression sets...
library("GEOmetadb")
library("stringr")
sqldb <- 'C:/Users/user/Documents/GitHub/pigRpeople/GEOmetadb.sqlite' # download from 11/1/13
con <- dbConnect(SQLite(),sqldb) 

# set up auxiliary data & access functions

# needed for access functions
colDesc <- columnDescriptions(sql = 'C:/Users/user/Documents/GitHub/pigRpeople/GEOmetadb.sqlite')
gseDesc <- colDesc[colDesc[[1]] == 'gse', 2:3]  # for reference - variable descriptions for series
gsmDesc <- colDesc[colDesc[[1]] == 'gsm', 2:3]  # for reference - variable descriptions for samples
gplDesc <- colDesc[colDesc[[1]] == 'gpl', 2:3]  # for reference - variable descriptions for platforms

# get vector of values from vector of variable indeces and geo accession number.
# TODO: add other types (gpl, gds, others?)

GEOVar <- function(geo, vars) {
  type <- tolower(str_sub(geo, end = 3))
  if(is.na(type)) return("missing GEO accession number")
  sql <- paste("SELECT ",
    paste(colDesc[colDesc$TableName == type,2][vars], collapse = ", "),
    " FROM ", type,
    " WHERE ", type, " = '", geo, "'",
    sep = "")
  df <- dbGetQuery(con, sql)
  nval <- as.character(df)
  names(nval) <- colnames(df)
  return(nval)
}

# get dataframe of variables from vector of geo accession #s.  assumes records are all same type.
# TODO: check codes are same type and return codes & variable names as dimnames

getVars <- function(geo, vars) {
  type <- tolower(str_sub(geo, end = 3))
  if(is.na(type)) return("missing GEO accession number")
  sql <- paste("SELECT ",
    paste(colDesc[colDesc$TableName == type,2][vars], collapse = ", "),
    " FROM ", type,
    " WHERE ", type, " IN ('", paste(geo, collapse = "', '"), "')",
    sep = "")
  df <- dbGetQuery(con, sql)
  row.names(df) <- geo
  return(df)
}

selVar <- function(char, var) {
  val <- matrix(nrow = length(char), ncol = length(var))
#    dimnames = list(char, ))
  for(i in seq_along(char)) val[i,] <- GEOVar(char[i], var)
  return(val)
} 

# get vector of sample codes associated with a series code

gse2gsm <- function(gse) {
  out <- list()
  for(i in seq_along(gse)) {
    samples <- geoConvert(gse[i], sqlite = sqldb)$gsm[,2]
    out <- c(out, list(samples))
    print(length(samples))
  }
  return(out)
}

# get vector of dataset codes associated with a series code

gse2gds <- function(gse) {
  for(i in seq_along(gse)){
    sql <- paste("SELECT gds.gse, gds.gds, gds.title, gds.sample_count",
    " FROM gds WHERE gds.gse = '",
    gse[i], "'", sep = "")
    print(dbGetQuery(con, sql))
  }
}

# parse string of variable/value pairs to value vector with named as variable

string2char <- function(str, sp = ";\t") {
  lns <- str_count(str, sp) +1
  vars <- character(lns)
  for(i in seq_along(vars)) vars[i] <- word(str, i, sep = sp)
  var <- matrix(unlist(str_split(vars, ": ")), nrow = 2)
  vars <- var[2,]
  names(vars) <- var[1,]
  return(vars)
}

# get variable from variable string by indexes from geo record name

getSubVar <- function(str, var, sub, sp = ";\t") string2char(GEOVar(str, var), sp)[sub]

# vectorization of above function (assumes all records have same structure)

extSubVar <- function(str, var, sub, sp = ";\t") {
  val <- character(length(str))
  for(i in seq_along(val)) val[i] <- getSubVar(str[i], var, sub, sp)
  return(val)
}

# separate string names by substring into string matrix and optionally return one column

split_char <- function(chr, sp, rw = 0) {
  vars <- str_split_fixed(chr, sp, n = str_count(chr[1], sp) + 1)
  if(rw > 0 && rw <= length(vars[1,])) return(vars[,rw]) else return(vars) 
}

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
ann_file <- c(T,F,F,F,F,F,F,F,T,F,T,T,T,T,F,T)

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
# 
# dbDisconnect(con) # run when finished searching GEOmetadatadb

library("GEOquery")

#  generate list of series matrix expression sets and store compressed downloaded files in directory "getgeo" in pwd.
gsematrix <- list(NULL)
length(gsematrix) <- 16
names(gsematrix) <- gsenames
for(i in 1:16){gsematrix[[i]] <- getGEO(GEO = series[i], destdir = 'getgeo', AnnotGPL = ann_file[i])}

gsmlists <- gse2gsm(series)
names(gsmlists) <- gsenames
