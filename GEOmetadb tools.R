## set up GEOmetadb GEO metadata database for GEO data exploration  --  useful for screening GEO data sets

library("GEOmetadb")
library("stringr")
sqldb <- 'C:/Users/user/Documents/GitHub/pigRpeople/GEOmetadb.sqlite' # download from 11/1/13
con <- dbConnect(SQLite(),sqldb) 

## set up auxiliary data & access functions

# needed for access functions

colDesc <- columnDescriptions(sql = 'C:/Users/user/Documents/GitHub/pigRpeople/GEOmetadb.sqlite')
gseDesc <- colDesc[colDesc[[1]] == 'gse', 2:3]  # for reference - variable descriptions for series
gsmDesc <- colDesc[colDesc[[1]] == 'gsm', 2:3]  # for reference - variable descriptions for samples
gplDesc <- colDesc[colDesc[[1]] == 'gpl', 2:3]  # for reference - variable descriptions for platforms

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

# make a list of of lists of samples named for GEO series

gsmlists <- gse2gsm(series)
names(gsmlists) <- gsenames

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
