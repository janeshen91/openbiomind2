# human lung transplant data set number 1

################################################################
# setup & description

hlungtx1 <- gsematrix$hlungtx1
hlungtx1 <- hlungtx1[[1]]

# > hlungtx1
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 17635 features, 50 samples 
#   element names: exprs 
# protocolData: none
# phenoData
#   sampleNames: GSM198616 GSM198617 ... GSM198665 (50 total)
#   varLabels: title geo_accession ... data_row_count (32 total)
#   varMetadata: labelDescription
# featureData
#   featureNames: AFFX-BioB-3_at AFFX-BioB-5_at ... XR_001531_at (17635 total)
#   fvarLabels: ID Gene title ... GO:Component ID (21 total)
#   fvarMetadata: Column Description labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: GPL5356 

# > pData(hlungtx1)[1,]    # geo sample GSM198616
#                    title                     status submission_date last_update_date type channel_count
#  6416_Lung developed PGD      Public on Jun 06 2007     Jun 06 2007      Jun 06 2007  RNA             1

#                                                                                   source_name_ch1          organism_ch1
#  lung biopsies from the anterior right middle lobe or lingula immediately prior to cold-flushing.          Homo sapiens
#           characteristics_ch1 
#  P/F ratio at T0:180    total RNA
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             extract_protocol_ch1
#  Single isolates of donor lung samples were homogenised in the presence of RNAzolB and finally dissolved in RNase-free H2O.  25 micro-grams of total RNA was treated with DNase using the Qiagen RNase-free DNase kit and samples were further purified using RNeasy spin columns (Qiagen, Valencia, CA).  Total RNA treated with DNase was dissolved in RNase-free H2O to a final concentration of 0.2 g/l.  RNA quality was assessed by 1% agarose gel electrophoresis in the presence of ethidium bromide.  Samples that did not reveal intact and approximately equal 18S and 28S ribosomal bands were excluded from further study.
#  label_ch1                                                                                          label_protocol_ch1
#     biotin Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 25 micrograms total RNA.

#  hyb_protocol
#  Following fragmentation, 10 micrograms of cRNA were hybridized for 16 hr at 45C on GeneChip Human Genome Array (HGU133Av2). GeneChips were washed and stained in the Affymetrix Fluidics Station 450.
#                                                                                           scan_protocol
#  Affymetrix GeneArray scanner 3000. Image analysis was performed with the Affymetrix GeneChip software.

#   description
#  Samples were immediately snap-frozen in liquid nitrogen and then stored in a -70Â° Celsius freezer until used for analysis.  An area of lung tissue approximately 1 x 1 cm was isolated and excised using 2 staple lines from a 30 mm EndoGIA stapler (US Surgical, Norwalk, CT).

#                                                data_processing platform_id contact_name      contact_email contact_phone
#  The data were normalised using gcRMA.  Log2-transformed data.     GPL5356  Monika,,Ray mray@cse.wustl.edu    3149358788
# contact_fax                  contact_institute                contact_address contact_city contact_state
#   3149356160 Washington University in St. Louis 1 Brookings Deive P O Box 1045    St. Louis            MO

#  contact_zip/postal_code contact_country
#                    63112             USA
#                                                                supplementary_file data_row_count
#  ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM198nnn//.CEL.gz          17635

#####   Differential expression analysis with limma
library(limma)

# make proper column names to match toptable (no white space)
fvarLabels(hlungtx1) <- make.names(fvarLabels(hlungtx1))

# group names for all samples
sml <- c("case","case","control","case","case","case","control","control","case","case","control","control","control","control","control","control","control","control","control","control","control","case","case","control","case","control","control","control","control","control","case","case","case","control","case","control","case","control","control","control","control","control","control","control","control","case","control","control","control","control")

# set up the data and proceed with analysis
fl <- as.factor(sml)
hlungtx1$description <- fl  # see above for original description
design <- model.matrix(~ description + 0, hlungtx1)
colnames(design) <- levels(fl)
fit <- lmFit(hlungtx1, design)
cont.matrix <- makeContrasts(control-case, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
hlungtx1.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=500)[c(-5, -6, -7, -11, -12, -13)]

################################################################
#   Boxplot for selected GEO samples

# order samples by group
ex <- exprs(hlungtx1)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("primary graft disease (PGD)","no PGD")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(hlungtx1)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(hlungtx1)))/2),4,2,1))
title <- paste ("GSE8021","(transplanted human lung tissue)", " log2 transformed expression levels", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, inset = c(.1, .2), fill=palette())

################################################################
#   construct moses dataset

# get probe x sample log2 normalized expression level matrix from expression set
hlungtx1.moses <- as.data.frame(exprs(hlungtx1))

# median normalize with med.normalize() from "data cleaning.R" and transpose to moses form
hlungtx1.moses <- t(lapply(hlungtx1.moses, med.normalize))

# add control binary (cases are transplants resulting in primary graft disfunction)
controls <- c(0,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,1,1,0,0,0,1,0,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1)
hlungtx1.moses <- as.data.frame(cbind(controls, hlungtx1.moses))

# remove control spots
gpl5356.controls <- as.character(fData(hlungtx1)$ID[fData(hlungtx1)$Platform_SPOTID == "--CONTROL"])
hlungtx1.moses <- hlungtx1.moses[!(names(hlungtx1.moses) %in% gpl5356.controls)]

