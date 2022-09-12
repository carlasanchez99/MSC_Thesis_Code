#Microarrays

#1. Load the data

load("/Users/carlasanchezalmirall/Desktop/lachmi_ExpressionSet.Rdata")

#Boxplot to ensure normalization was correctly performed 

boxplot(exprs(lachmi))

#Delete positive and negative probes

lachmimicro<-lachmi
lachmimicro<-lachmimicro[!fData(lachmimicro)$`mRna.-.Description` == 'neg_control ']
lachmimicro<-lachmimicro[!fData(lachmimicro)$`mRna.-.Description` == 'pos_control ']

dgr<-lachmimicro
pData(dgr)$tbs<-as.factor(pData(dgr)$tbs)
fData(dgr)$Gene.Symbol<-gsub("---", NA, fData(dgr)$Gene.Symbol)
nan<-is.na(fData(dgr)$Gene.Symbol)
dgr<-dgr[!nan]
fData(dgr)$Gene.Symbol<-gsub(" ", "", fData(dgr)$Gene.Symbol)

#Perform baseline analysis
dgr<-dgr[,dgr$t_point=="C-1"]

mod1 <- model.matrix( ~ tbs+donor, data=pData(dgr))
mod0 <- model.matrix( ~ donor, data=pData(dgr))
sv <- svaseq(exprs(dgr), mod1, mod0)

V1<-sv$sv[,1]
V2<-sv$sv[,2]
design1 <- model.matrix(~0+pData(dgr)$tbs+pData(dgr)$donor+V1+V2)
colnames(design1) <- c("tbsneg", "tbspos", "donor","V1","V2")
fit1 <- lmFit(exprs(dgr),design1)
cm1 <- makeContrasts(cont=tbsneg-tbspos, levels=design1)
# contrasts fit and Bayesian adjustment
fit2.1 <- contrasts.fit(fit1, cm1)
fit2.1 <- eBayes(fit2.1)
top.table1<-topTable(fit2.1, coef=1, number = Inf, adjust.method = "fdr")

genesmicro<-top.table1[top.table1$P.Val<0.05,] 
micros<-dgr[rownames(genesmicro)]
negmic<-genesmicro[genesmicro$logFC>0,]
negmic<-fData(dgr)%>% dplyr::filter(Probe.Set.ID %in% rownames(negmic))
posmic<-genesmicro[genesmicro$logFC<0,]
posmic<-fData(dgr)%>% dplyr::filter(Probe.Set.ID %in% rownames(posmic))

hs<-org.Hs.eg.db
syneg<-AnnotationDbi::select(hs, 
                             keys =  negmic$Gene.Symbol,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")

sypos<-AnnotationDbi::select(hs, 
                             keys =  posg$SYMBOL,
                             columns = c("ENTREZID", "SYMBOL"),
                             keytype = "SYMBOL")


ans.gopos <- enrichGO(gene = syneg$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)


ans.goneg <- enrichGO(gene = syneg$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)


