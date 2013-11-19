# human lung transplant bronciolar lavage data set

################################################################
# setup & description

hlavtx <- gsematrix$hlavtx
hlavtx <- hlavtx[[1]]

# $hlavtx
# $hlavtx$GSE2018_series_matrix.txt.gz
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 22283 features, 34 samples 
#   element names: exprs 
# protocolData: none
# phenoData
#   sampleNames: GSM36423 GSM36424 ... GSM36456 (34 total)
#   varLabels: title geo_accession ... data_row_count (27 total)
#   varMetadata: labelDescription
# featureData
#   featureNames: 1007_s_at 1053_at ... AFFX-TrpnX-M_at (22283 total)
#   fvarLabels: ID Gene title ... GO:Component ID (21 total)
#   fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: GPL96 

# > pData(hlavtx)[1:2,]
#                 title geo_accession                status submission_date last_update_date type channel_count
# GSM36423  BAL_37_A1B1      GSM36423 Public on Nov 30 2004     Nov 26 2004      Jun 26 2007  RNA             1
# GSM36424 BAL_15a_A0B1      GSM36424 Public on Nov 30 2004     Nov 26 2004      Jun 26 2007  RNA             1
#                       source_name_ch1 organism_ch1 molecule_ch1 taxid_ch1
# GSM36423 Bronchoalveolar lavage cells Homo sapiens    total RNA      9606
# GSM36424 Bronchoalveolar lavage cells Homo sapiens    total RNA      9606
#                                                                                                           description
# GSM36423      Bronchoalveolar lavage cells from a patient with an ISHLT biopsy score of A=1 and B=1 (Acute rejection)
# GSM36424 Bronchoalveolar lavage cells from a patient with an ISHLT biopsy score of A=0 and B=1 (Non-acute rejection).
#                                                        description.1 platform_id contact_name    contact_email contact_phone
# GSM36423                                                                   GPL96  Jeff,,Lande land0038@umn.edu  612-625-5628
# GSM36424 One of two samples from the same patient, different visits.       GPL96  Jeff,,Lande land0038@umn.edu  612-625-5628
#          contact_department       contact_institute          contact_address contact_city contact_state
# GSM36423           Medicine University of Minnesota 420 Washington Avenue SE  Minneapolis            MN
# GSM36424           Medicine University of Minnesota 420 Washington Avenue SE  Minneapolis            MN
#          contact_zip/postal_code contact_country
# GSM36423                   55455             USA
# GSM36424                   55455             USA
#                                                                                       supplementary_file
# GSM36423 ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM36nnn/GSM36423/GSM36423.CEL.gz
# GSM36424 ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM36nnn/GSM36424/GSM36424.CEL.gz
#                                                                                     supplementary_file.1 data_row_count
# GSM36423 ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM36nnn/GSM36423/GSM36423.EXP.gz          22283
# GSM36424 ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM36nnn/GSM36424/GSM36424.EXP.gz          22283

#####   Differential expression analysis with limma
library(limma)

# make proper column names to match toptable (no white space)
fvarLabels(hlavtx) <- make.names(fvarLabels(hlavtx))

# group names for all samples
sml <- c("case","control","control","control","control","control","control","control","control","case","case","control","control","case","control","case","control","control","control","control","control","control","control","control","control","control","control","control","control","case","case","control","control","control")

################################################################
#   differential gene expression by empirical bayes linear models for 1 channel array

# set up the data and apply lmFit and eBayes with contrasts
fl <- as.factor(sml)
hlavtx$description <- fl  # see above for original description
design <- model.matrix(~ description + 0, hlavtx)
colnames(design) <- levels(fl)
fit <- lmFit(hlavtx, design)
cont.matrix <- makeContrasts(control-case, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
hlavtx.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=500)[c(-5, -6, -7, -11, -12, -13)]

# export file
write.csv(hlavtx.tT, file = "results/transplant_samples/hlavtx_diffex500.csv", row.names = FALSE)

################################################################
#   Boxplot for selected GEO samples (note:  the .jpg & .svg files in ~/openbiomind2/results/transplant_samples were exported from RStudio environment)

# order samples by group
ex <- exprs(hlavtx)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("acute rejection","no acute rejection")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(hlavtx)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(hlavtx)))/2),4,2,1))
title <- paste ("GSE2018"," (transplant bronchiolar lavage samples)", " log2 transformed expression levels", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, inset = c(.1, .2), fill=palette())

################################################################
#   construct moses dataset

# get probe x sample log2 normalized expression level matrix from expression set
hlavtx.moses <- exprs(hlavtx)

#median normalize with med.normalize() from "data cleaning.R"
hlavtx.moses <- med.normalize(hlavtx.moses)

# add control binary (cases are transplants resulting in primary graft disfunction)
controls <- c(0,1,1,1,1,1,1,1,1,0,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1)
hlavtx.moses <- t(rbind(controls, hlavtx.moses))

# remove control spots
gpl96.controls <- as.character(fData(hlavtx)$ID[fData(hlavtx)$Platform_SPOTID == "--Control"])
hlavtx.moses <- hlavtx.moses[,!(colnames(hlavtx.moses) %in% gpl96.controls)]

# export file
write.csv(hlavtx.moses, file = "results/transplant_samples/hlavtx_moses.csv")

