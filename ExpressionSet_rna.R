dat<-read.table("/Users/carlasanchezalmirall/Downloads/LaCHMI_RNASEQ.808046.csv", header=TRUE, sep=";")
pheno<-read.table("/Users/carlasanchezalmirall/01_v20.05.08_clinical_chmi_dat_rna.csv", header=TRUE, sep=";")
dat1<-dat[,-1]
assayData <- dat1
rownames(assayData) <- dat$feature
pData <- pheno
rownames(pData) <- colnames(assayData)
phenoData <- new('AnnotatedDataFrame', data = pData)
rnalachmi<-new("ExpressionSet", exprs=as.matrix(assayData), phenoData = phenoData)

save(rnalachmi, file = 'rnalachmi.Rdata')
