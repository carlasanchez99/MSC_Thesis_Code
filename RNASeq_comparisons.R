#RNA-SEQ DATA

#Load dataset and explore with PCA

load("/Users/carlasanchezalmirall/Desktop/TFM/rnalachmi copia.RData")

#Compute variance to select top 100 genes
V <- apply(exprs(rnalachmi), 1, var)
selectedGenes <- names(V[order(V, decreasing = T)][1:100])
M <- t(exprs(rnalachmi)[selectedGenes,])
# transform to log2 scale 
M <- log2(M + 1)
#PCA 
pcaResults <- prcomp(M)
GROUP<-pData(rnalachmi)$gender
pcplot<-data.frame(GROUP, pcaResults$x[,1:2])
qplot(x=PC1, y=PC2, data=pcplot, colour=GROUP) 


#1. COMPARISON BETWEEN AFRICANS AND EUROPEANS

###########
##1.1 DEG##
###########

baseline<-rnalachmi[,rnalachmi$t_point=="C-1"]
colDatab<-pData(baseline)
countsb<-dat[,rownames(colDatab)]
countDatab <- as.matrix(countsb)

mod1 <- model.matrix( ~ immune_status+age+gender+hb_status, data=colDatab)
mod0 <- model.matrix( ~ age+gender+hb_status, data=colDatab)   

sv <- svaseq(countDatab, mod1, mod0)
temp <- DataFrame(colDatab, sv$sv)
designFormula<- as.formula( ~ immune_status+age+gender+hb_status+V1+V2+V3+V4)

#create a DESeq dataset object from the count matrix and the colData 
dds <- DESeqDataSetFromMatrix(countData = countDatab, 
                              colData = temp, 
                              design = designFormula)

dds <- DESeq(dds)
ans.sv <- results(dds, contrast = c("immune_status", 'naive', 'semi_immune'))
genes<-ans.sv[ans.sv$padj<0.05  & abs(ans.sv$log2FoldChange) > log2(2) & !is.na(ans.sv$padj),] 


###1.2 CAMERA+VOOM

design1<- model.matrix(~0+immune_status+gender+age+hb_status+V1+V2+V3+V4, data=colDatab)
v <- voom(countDatab, design=design1, plot=TRUE)
geneIDsgood <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(v), columns = c("SYMBOL","GENEID"))
idx <- ids2indices(btm$genesets,id=geneIDsgood$SYMBOL)
contr.matrix <- makeContrasts(
  NaivevsSemi = immune_statusnaive-immune_statussemi_immune, 
  levels = colnames(design1))
contr.matrix
cam.BasalvsLP <- camera(countDatab,idx,design1,contrast=contr.matrix[,1])
sigeneset<-cam.BasalvsLP[cam.BasalvsLP$FDR<0.25,]
sigeneset<-sigeneset[sigeneset$NGenes>4,]
numgeneset<-as.numeric(rownames(sigeneset))
sigs<-cbind(sigeneset,btm$geneset.names[numgeneset])
colnames(sigs)<-c(colnames(sigeneset), "name")
sigs<-sigs%>%dplyr::filter(!str_detect(sigs$name, "TBA"))


#D5 COMPARISON BETWEEN AFRICANS AND EUROPEANS (adjusting for baseline)


d5<-rnalachmi[,rnalachmi$t_point=="C-1"|rnalachmi$t_point=="D5"]
colDatad5<-pData(d5)
counts5<-dat[,rownames(colData5)]
countData5 <- as.matrix(counts5)

d5.2<-rnalachmi[,rnalachmi$t_point=="D5"]
pData(d5.2)$t_point<-droplevels(pData(d5.2)$t_point)
colData5.2<-pData(d5.2)
countsd5.2<-dat[,rownames(colData5.2)]
countDatad5.2<- as.matrix(countsd5.2)
mod1 <- model.matrix( ~ immune_status+age+gender+hb_status, data=colData5.2)
mod0 <- model.matrix( ~ age+gender+hb_status, data=colData5.2)
sv <- svaseq(countDatad5.2, mod1, mod0)

temp5 <- DataFrame(colData5.2, sv$sv)
designFormula5<- as.formula( ~ immune_status+t_point+immune_status:t_point+age+gender+hb_status+V1+V2+V3+V4)

dds5 <- DESeqDataSetFromMatrix(countData = countData5,
                                colData =temp5,
                                design = designFormula5)

ddsTC5 <- DESeq(dds5, test="LRT", reduced = ~ immune_status+t_point+age+gender+hb_status+V1+V2+V3+V4)
resTC5 <- results(ddsTC5)

genes<-resTC5[resTC5$padj<0.05  & abs(resTC5$log2FoldChange) > log2(2) & !is.na(resTC5$padj),] 

#D11 COMPARISON BETWEEN AFRICANS AND EUROPEANS (adjusting for baseline)


d11<-rnalachmi[,rnalachmi$t_point=="C-1"|rnalachmi$t_point=="D11"]
colDatad11<-pData(d11)
counts11<-dat[,rownames(colData11)]
countData11 <- as.matrix(counts11)

d11.2<-rnalachmi[,rnalachmi$t_point=="D11"]
pData(d11.2)$t_point<-droplevels(pData(d11.2)$t_point)
colData11.2<-pData(d11.2)
countsd11.2<-dat[,rownames(colData11.2)]
countDatad11.2<- as.matrix(countsd11.2)
mod1 <- model.matrix( ~ immune_status+age+gender+hb_status, data=colData11.2)
mod0 <- model.matrix( ~ age+gender+hb_status, data=colData11.2)
sv <- svaseq(countDatad11.2, mod1, mod0)

temp11 <- DataFrame(colData11.2, sv$sv)
designFormula11<- as.formula( ~ immune_status+t_point+immune_status:t_point+age+gender+hb_status+V1+V2+V3+V4)

dds11 <- DESeqDataSetFromMatrix(countData = countData11,
                              colData =temp11,
                              design = designFormula11)

ddsTC11 <- DESeq(dds11, test="LRT", reduced = ~ immune_status+t_point+age+gender+hb_status+V1+V2+V3+V4)
resTC11 <- results(ddsTC11)

genes<-resTC11[resTC11$padj<0.05  & abs(resTC11$log2FoldChange) > log2(2) & !is.na(resTC11$padj),] 



#2. AMONG AFRICANS COMPARING TBS+/TBS-  

baseline<-rnalachmi[,rnalachmi$t_point=="C_1"]
baseline<-baseline[, baseline$immune_status=="semi_immune"]
pData(baseline)$t_point<-droplevels(pData(baseline)$t_point)
colDatab<-pData(baseline)
countsb<-dat[,rownames(colDatab)]
countDatab <- as.matrix(countsb)

mod1 <- model.matrix( ~ tbs_pos+gender+age+hb_status, data=colDatab)
mod0 <- model.matrix( ~ gender+age+hb_status, data=colDatab)   

sv <- svaseq(countDatab, mod1, mod0)
temp <- DataFrame(colDatab, sv$sv)
designFormula<- as.formula( ~ tbs_pos+gender+age+hb_status+V1+V2+V3)

dds <- DESeqDataSetFromMatrix(countData = countDatab, 
                              colData = temp, 
                              design = designFormula)

dds <- DESeq(dds)
ans.sv <- results(dds, contrast = c("tbs_pos", '0', '1'))
genes<-ans.sv[ans.sv$pvalue<0.05  &!is.na(ans.sv$pvalue),] 



###1.2 CAMERA+VOOM

design1<- model.matrix(~0+tbs_pos+gender+age+hb_status+V1+V2+V3, data=colDatab)
v <- voom(countDatab, design=design1, plot=TRUE)
geneIDsgood <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(v), columns = c("SYMBOL","GENEID"))
idx <- ids2indices(btm$genesets,id=geneIDsgood$SYMBOL)
contr.matrix <- makeContrasts(
  NegvsPos = tbs_pos0-tbs_pos1, 
  levels = colnames(design1))
contr.matrix
cam.BasalvsLP <- camera(countDatab,idx,design1,contrast=contr.matrix[,1])
sigeneset<-cam.BasalvsLP[cam.BasalvsLP$FDR<0.25,]
sigeneset<-sigeneset[sigeneset$NGenes>4,]
numgeneset<-as.numeric(rownames(sigeneset))
sigs<-cbind(sigeneset,btm$geneset.names[numgeneset])
colnames(sigs)<-c(colnames(sigeneset), "name")
sigs<-sigs%>%dplyr::filter(!str_detect(sigs$name, "TBA"))



#D11 COMPARISON BETWEEN AFRICANS TBS+ AND TBS- (adjusting for baseline)


d11<-rnalachmi[,rnalachmi$t_point=="C-1"|rnalachmi$t_point=="D11"]
d11<-d11[, d11$immune_status=="semi_immune"]
colDatad11<-pData(d11)
counts11<-dat[,rownames(colData11)]
countData11 <- as.matrix(counts11)

d11.2<-rnalachmi[,rnalachmi$t_point=="D11"]
pData(d11.2)$t_point<-droplevels(pData(d11.2)$t_point)
colData11.2<-pData(d11.2)
countsd11.2<-dat[,rownames(colData11.2)]
countDatad11.2<- as.matrix(countsd11.2)
mod1 <- model.matrix( ~ tbs_pos+age+gender+hb_status, data=colData11.2)
mod0 <- model.matrix( ~ age+gender+hb_status, data=colData11.2)
sv <- svaseq(countDatad11.2, mod1, mod0)

temp11 <- DataFrame(colData11, sv$sv)
designFormula11<- as.formula( ~ tbs_pos+t_point+tbs_pos:t_point+age+gender+hb_status+V1+V2+V3+V4)

dds11 <- DESeqDataSetFromMatrix(countData = countData11,
                                colData =temp11,
                                design = designFormula11)

ddsTC11 <- DESeq(dds11, test="LRT", reduced = ~ tbs_pos+t_point+age+gender+hb_status+V1+V2+V3+V4)
resTC11 <- results(ddsTC11)

genes<-resTC11[resTC11$pvalue<0.011  &!is.na(resTC11$pvalue),]



###1.2 CAMERA+VOOM

colDatad11$tt<-paste0(colDatad11$t_point, sep="",colDatad11$tbs_pos)
design1<- model.matrix(~0+tt+gender+age+hb_status+V1+V2+V3+V4, data=temp11)
colnames(design1)<-c("tbsnegc1", "tbsposc1", "tbsnegd11", "tbsposd11", "genderM", "age", "hb_status","V1",  "V2", "V3", "V4")
v <- voom(countDatab, design=design1, plot=TRUE)
geneIDsgood <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(v), columns = c("SYMBOL","GENEID"))
idx <- ids2indices(btm$genesets,id=geneIDsgood$SYMBOL)
contr.matrix <- makeContrasts(
  NegvsPos = (tbsnegd11-tbsnegc1)- (tbsposd11-tbsposc1), 
  levels = colnames(design1))
contr.matrix
cam.BasalvsLP <- camera(countDatab11,idx,design1,contrast=contr.matrix[,1])
sigeneset<-cam.BasalvsLP[cam.BasalvsLP$FDR<0.25,]
sigeneset<-sigeneset[sigeneset$NGenes>4,]
numgeneset<-as.numeric(rownames(sigeneset))
sigs<-cbind(sigeneset,btm$geneset.names[numgeneset])
colnames(sigs)<-c(colnames(sigeneset), "name")
sigs<-sigs%>%dplyr::filter(!str_detect(sigs$name, "TBA"))


#1.3 GSEA with BTMs
foldch <-genes$log2FoldChange
my.symbols <- geneIDs$SYMBOL
my.symbols<-gsub(" ", "", my.symbols)
names(foldch) <- my.symbols

foldch<- sort(foldch, decreasing = TRUE)
head(foldch)

msigC_1D5 <- GSEA(foldch, TERM2GENE=pathw, verbose=FALSE, pvalueCutoff = 0.1, minGSSize = 2)


#1.4 Enrichment GO database
neg<-genes[genes$log2FoldChange>0,]
pos<-genes[genes$log2FoldChange<0,]

negg <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(neg), columns = c("SYMBOL","GENEID"))
posg <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(pos), columns = c("SYMBOL","GENEID"))

hs<-org.Hs.eg.db
syneg<-AnnotationDbi::select(hs, 
                          keys =  negg$SYMBOL,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

sypos<-AnnotationDbi::select(hs, 
                             keys =  posg$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


ans.go <- enrichGO(gene = sy$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)



