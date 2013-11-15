# human lung transplant data set number 2 (not annotated)

################################################################
# setup & description

hlungtx2 <- gsematrix$hlungtx2
hlungtx2 <- hlungtx2[[1]]

# $hlungtx2
# $hlungtx2$GSE9102_series_matrix.txt.gz
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 18746 features, 26 samples 
#   element names: exprs 
# protocolData: none
# phenoData
#   sampleNames: GSM230361 GSM230362 ... GSM230386 (26 total)
#   varLabels: title geo_accession ... data_row_count (43 total)
#   varMetadata: labelDescription
# featureData
#   featureNames: 1 2 ... 18746 (18746 total)
#   fvarLabels: ID GB_ACC Name SPOT_ID
#   fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: GPL5909 

# > pData(hlungtx2)[1,] # GSM230361 
#                                                  title geo_accession                status submission_date last_update_date type
# Control (Good outcome) -matched to 127_PGD_1     Public on Sep 15 2008     Sep 19 2007      Aug 14 2011  RNA
#           channel_count                                                       source_name_ch1 organism_ch1 characteristics_ch1
#             2 Total RNA from human donor lung tissues labeled with Cyanine-5 (red). Homo sapiens          Human lung
#              characteristics_ch1.1 molecule_ch1
# End cold ischemic period    total RNA
#                                                                                                                                                                                                                                                                    extract_protocol_ch1
# Total RNA (10 Î¼g) was prepared from lung tissue sample using Trizol Reagent (Invitrogen, Carlsbad, CA) and further purified using an RNeasy kit (Qiagen, Valencia, CA) with a DNase treatment (RNase-free DNase Set, Qiagen) according to the manufacturerâ€™s instructions.
#           label_ch1
#       Cy5
#                                                                                                                                                                                                                                                                                                                                                           label_protocol_ch1
# Total RNA was added to 8 Âµl of 5Ã— first strand reaction buffer, 0.75 Âµl of 200 ÂµM AncT mRNA primer (5â€™-T20VN-3â€™) and 4 Âµl of 0.1 M DTT. Three Âµl of 20 mM ddNTP (dATP, dGTP, and dTTP mix), 1 Âµl of 2 mM dCTP, 1 ng of control RNA and 1 Âµl of 1 mM Cyanine 3 or Cyanine 5 dCTP were added to each RNA sample to a total volume of 40 Âµl in the dark.
#           taxid_ch1                                                             source_name_ch2 organism_ch2 characteristics_ch2
#      9606 Human Universal Reference RNA (Stratagene), labeled with Cyanine-3 (green). Homo sapiens       reference RNA
#           molecule_ch2
#    total RNA
#                                                                                                                                                                                                                                                                    extract_protocol_ch2
# Total RNA (10 Î¼g) was prepared from lung tissue sample using Trizol Reagent (Invitrogen, Carlsbad, CA) and further purified using an RNeasy kit (Qiagen, Valencia, CA) with a DNase treatment (RNase-free DNase Set, Qiagen) according to the manufacturerâ€™s instructions.
#           label_ch2
#       Cy3
#                                                                                                                                                                                                                                                                                                                                                           label_protocol_ch2
# Total RNA was added to 8 Âµl of 5Ã— first strand reaction buffer, 0.75 Âµl of 200 ÂµM AncT mRNA primer (5â€™-T20VN-3â€™) and 4 Âµl of 0.1 M DTT. Three Âµl of 20 mM ddNTP (dATP, dGTP, and dTTP mix), 1 Âµl of 2 mM dCTP, 1 ng of control RNA and 1 Âµl of 1 mM Cyanine 3 or Cyanine 5 dCTP were added to each RNA sample to a total volume of 40 Âµl in the dark.
#           taxid_ch2
#      9606
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           hyb_protocol
# The reaction mixture was heated to 65Â°C for 5 minutes and then placed at 42Â°C for 5 minutes. Subsequently, reverse transcription was performed at 42Â°C for 2 hours by adding 2 Âµl of reverse transcriptase (Superscript II, Life Technologies, Inc., Burlington, Canada) to each sample. RNA was hydrolyzed by adding 4 Âµl of 50 mM EDTA and 2 Âµl of 10 N NaOH, and then incubated at 65Â°C for 20 minutes. Samples were neutralized with 4 Âµl of 50 M of acetic acid and 50 Âµl of isopropanol. At this point, the two reactions to be hybridized to the same slide (i.e. control and experimental) were combined together. Air dried pellets of labeled cDNA were obtained after removing isopropanol; the samples were then resuspended with 5 Âµl of distilled water. Samples were resuspended in 100 Âµl of DIGeasy hybridization buffer (Roche, Basel, Switzerland) combined with 5 Âµl of 10 mg/ml yeast tRNA (Life Technologies, Inc.) and 5 Âµl of 10 mg/ml salmon sperm (boiled 10 minutes prior to use), denatured at 65Â°C for 2 minutes, then hybridized to the arrays overnight at 37Â°C in humid hybridization chambers. Arrays were washed next day, three times for 15 minutes each at 50Â°C in slide staining boxes containing pre-warmed 1Ã— SSC, 0.1% SDS, and then rinsed two times with 1Ã— SSC. Slides were dried by centrifugation at 1000 rpm for 5 minutes.
#                                                                                                                                                                                                              scan_protocol
# Hybridized arrays were scanned using a GenePix 4000 scanner (Axon Instruments Inc., Union City, CA), and fluorescent images were analyzed with GenePix Pro version 5.0 (Molecular Devices Corp., Sunnywale, CA).
#                    description
# patients without PGD
#                                                                                                                                                                                                        data_processing
# LOWESS normalized, background subtracted data obtained from log2 of processed Red signal/processed Green signal. Microarray Data Analysis System version 2.19 (The Institute for Genomic Research) was used.
#           platform_id   contact_name        contact_email contact_phone  contact_fax        contact_laboratory contact_department
#     GPL5909 Masaki,,Anraku m.anraku@utoronto.ca  416-581-7503 416-581-7702 Thoracic Surgery Research            Surgery
#                            contact_institute                contact_address contact_city contact_state contact_zip/postal_code
# Toronto General Research Institute 101 College Street, TMDT 2-816      Toronto       Ontario                  M5G1L7
#           contact_country                                                                                 supplementary_file
#          Canada ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM230nnn/GSM230361/GSM230361.gpr.gz
#           data_row_count
#          18746

#####   Differential expression analysis with limma
library(limma)

# make proper column names to match toptable 
fvarLabels(hlungtx2) <- make.names(fvarLabels(hlungtx2))

# group names for all samples
sml <- c("control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","control","case","case","case","case","case","case","case","case","case","case");

################################################################
#   differential gene expression by empirical bayes linear models for 2 channel competitive binding array

# set up the data and  and apply lmFit and eBayes with contrasts

fl <- as.factor(sml)
hlungtx2$description <- fl
design <- model.matrix(~ description + 0, hlungtx2)
colnames(design) <- levels(fl)
fit <- lmFit(hlungtx2, design)
cont.matrix <- makeContrasts(control-case, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
hlungtx2.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=500)#[c(-5, -6, -7, -11, -12, -13)] remove toptable columns?

# export file
write.csv(hlungtx2.tT, file = "results/transplant_samples/hlungtx2_diffex500.csv", row.names = FALSE)

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# group names for all samples in a series
sml <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)

# order samples by group
ex <- exprs(hlungtx2)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("primary graft disease (PGD)","no PGD")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(hlungtx2)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(hlungtx2)))/2),4,2,1))
title <- paste ("GSE9102", '/', annotation(hlungtx2), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, inset = c(.1, .2), fill=palette())

################################################################
#   construct moses dataset

# get probe x sample log2 normalized expression level matrix from expression set
hlungtx2.moses <- exprs(hlungtx2)
rownames(hlungtx2.moses) <- fData(hlungtx2)$GB_ACC

# remove control spots
hlungtx2.moses <- hlungtx2.moses[!(fData(hlungtx2)$GB_ACC == ""),]

#median normalize with med.normalize() from "data cleaning.R"
hlungtx2.moses <- med.normalize(hlungtx2.moses)

# add control binary (cases are transplants resulting in primary graft disfunction)
controls <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)
hlungtx2.moses <- t(rbind(controls, hlungtx2.moses))

# export file
write.csv(hlungtx2.moses, file = "results/transplant_samples/hlungtx2_moses.csv")

####################################################
# map probes from hlungtx1 to hlungtx2

x <- org.Hs.egUNIGENE2EG
mapped_genes <- mappedkeys(x)
ug <- as.list(x[mapped_genes])

x <- org.Hs.egREFSEQ
mapped_genes <- mappedkeys(x)
rs <- as.list(x[mapped_genes])

tx2unigene <- str_split_fixed(fData(hlungtx2)[[3]], fixed("_"), 2)
hlungtx <- cbind(fData(hlungtx2)[[2]], as.data.frame(tx2unigene, stringsAsFactors = FALSE),)
hlungtx <- hlungtx[hlungtx$gb_acc != "",]
names(hlungtx) <- c("gb_acc", "unigene", "symbol")
for(i in seq_along(hlungtx$unigene)){
  if(length(ug[[hlungtx$unigene[i]]]) == 0) hlungtx$eg[i] <- "chk" 
  else hlungtx$eg[i] <- str_c(ug[[hlungtx$unigene[i]]], collapse = ", ")
}
for(i in seq_along(hlungtx$unigene)){
  if(length(rs[[hlungtx$eg[i]]]) == 0) hlungtx$refseq[i] <- "chk"
  else hlungtx$refseq[i] <- str_c(rs[[hlungtx$eg[i]]][rs[[hlungtx$eg[i]]] %in% dimnames(hlungtx1.moses)[[2]]], collapse = ", ")
}
hlungtx$refseq <- map_n(rs, hlungtx$eg, dimnames(hlungtx1.moses)[[2]][-1])
