####################################################
# map probes from hlungtx1 to hlungtx2
library("org.Hs.eg.db")

#add manufacturer annotation
H19kdump <- read.delim("C:/Users/user/Desktop/biomind/chimpiggie/H19kdump.tab/H19kdump.csv", colClasses = "character")
H19k <- H19kdump[, c(3,6,7,11)]
mode(H19k$Accession) <- "character"

# list of unigene named vectors of entregenes
x <- org.Hs.egUNIGENE2EG
mapped_genes <- mappedkeys(x)
ug <- as.list(x[mapped_genes])

# list of entregene named vectors of refseqs
y <- org.Hs.egREFSEQ
mapped_genes <- mappedkeys(y)
rs <- as.list(y[mapped_genes])

# merge on entregene id
hlungtx.map <- merge(toTable(x), toTable(y))

# make mapping table
# list unigene labels in hlungtx2 because genebank acc #s only map for 90 probes ?!?
tx2unigene <- str_split_fixed(fData(hlungtx2)[[3]], fixed("_"), 2)
hlungtx <- cbind(fData(hlungtx2)[[2]], as.data.frame(tx2unigene, stringsAsFactors = FALSE))
names(hlungtx) <- c("gb_acc", "unigene", "symbol")
hlungtx$gb_acc <- as.character(hlungtx$gb_acc)
hlungtx <- merge(hlungtx, H19k, by = NULL, by.x = "gb_acc", by.y = "Accession", all.x = TRUE)
hlungtx <- hlungtx[hlungtx$gb_acc != "",]

# # get first refseq mapped to unigene id
# hlungtx$refseq.map <- character(16677)
# for(i in seq_along(hlungtx$unigene)) {
#   if(hlungtx$unigene[i] == "") next
#   if(sum(hlungtx.map$unigene_id == hlungtx$unigene[i]) == 0) next
#   hlungtx$refseq.map[i] <- 
#     hlungtx.map[match(hlungtx.map[hlungtx.map$unigene_id == hlungtx$unigene[i], 3], colnames(hlungtx1.moses), no = 10),3]
# }


# map unigene in hlungtx2 to entregene numbers
for(i in seq_along(hlungtx$unigene)){     # replace with map_1 ?
  if(length(ug[[hlungtx$unigene[i]]]) == 0) hlungtx$eg[i] <- "chk" 
  else hlungtx$eg[i] <- str_c(ug[[hlungtx$unigene[i]]], collapse = ", ")
}
# list all the refseq numbers matching probes with single unigene name
hlungtx$refseq.all <- map_n(rs, hlungtx$eg, dimnames(hlungtx1.moses)[[2]][-1])
# pick first one
for(i in seq_along(hlungtx$unigene)) hlungtx$refseq[i] <- str2char(hlungtx$refseq.all[i])

# fill in manufacturer's missing refseq with unigene -> refseq mapping
hlungtx$tx2name <- hlungtx$RefSeq.ID
hlungtx$tx2name[is.na(hlungtx$tx2name)] <- ""
hlungtx$tx2name[hlungtx$tx2name == ""] <- hlungtx$refseq[hlungtx$tx2name == ""]

# list the entregene numbers for probes with multiple matches skipped above
dup.eg <- str_split(hlungtx$eg[str_detect(hlungtx$eg, fixed(", "))], fixed(", "))
ex.eg <- unlist(lapply(str_split(hlungtx$eg[str_detect(hlungtx$eg, fixed(", "))], fixed(", ")), tail, n = -1))

dup.refseq <- map_n(rs, ex.eg, dimnames(hlungtx1.moses)[[2]][-1])
dup.refseq <- dup.refseq[dup.refseq != "chk"]
dup.refseq <- unlist(str_split(dup.refseq, fixed(", ")))

# > summary(hlungtx$tx2name %in% colnames(hlungtx1.moses))
#    Mode   FALSE    TRUE    NA's 
# logical    8791    7886       0 
# > summary(colnames(hlungtx1.moses) %in% hlungtx$tx2name)
#    Mode   FALSE    TRUE    NA's 
# logical   12649    4925       0 
# > summary(dup.refseq %in% hlungtx$tx2name)
#    Mode   FALSE    TRUE    NA's 
# logical     219      17       0 
