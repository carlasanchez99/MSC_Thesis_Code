#DIABLO BASELINE AFRICANS vs EUROPEANS (RNA-SEQ, CYTOS, ANTIBODIES)

genesdiab<-read.csv("baselinenaivesemiimmunecbo2.csv", header = TRUE, row.names = 1)

#ab
abbaseline<-pcabaseline %>% dplyr::filter(fullname %in% sigbaseline)

vesctbaseline<-c()
for(i in unique(abbaseline$fullname)){
  abanti<-abbaseline[abbaseline$fullname==i,]
  vesctbaseline<-cbind(vesctbaseline, abanti$log10_mfi)
  vesctbaseline<-vesctbaseline
}
colnames(vesctbaseline)<-c(unique(abbaseline$fullname))
rownames( vesctbaseline)<-unique(abbaseline$original_id)

#cytos
cytossig<-c("egf", "eotaxin")
cytobas<-dat[dat$t_point=="C_1",]
mfisd7<-c()
for(i in unique(cytossig)){
  cytobas1<-cytobas[cytobas$analyte==i]
  mfisd7=cbind(mfisd7, cytobas1$log10_mfi)
  mfisd7<-mfisd7
}

colnames(mfisd7)<-unique(cytossig)
rownames(mfisd7)<-unique(cytobas$original_id)

ens<-strsplit(colnames(genesdiab), "[_]")
ori<-sapply(ens, "[", 1)
colnames(genesdiab) <-ori

r<-rownames(vesctbaseline)
r<-gsub("_", ".", r)
rownames(vesctbaseline)<-r
genesdiab<-genesdiab[,r]

rownames(mfisd7)<-gsub("_", ".", rownames(mfisd7))

#DIABLO
X=list(RNASeq=t(genesdiab),  Antibodies=vesctbaseline, Cytokines=mfisd7)
cond<-abanti$immune_status
Y=cond

list.keepX <- list(RNASeq = c(309, 309), Antibodies = c(47, 47), Cytokines=c(2,2))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo,  legend=TRUE)

plotIndiv(MyResult.diablo2, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2))

circ<-circosPlot(MyResult.diablo, cutoff=0.8,  line = TRUE,var.adj = 100)
pdf("cimnaivesemi.pdf",  width = 20, height = 8)
cimDiablo(MyResult.diablo, legend.position = "right",size.legend =3, scale="col",margins=c(2, 38), col.names = FALSE, row.names = FALSE)
dev.off() 

abcor<-c()
for (i in colnames(circ)){
  if (sum(circ[i,]>0.8) >0 | sum(circ[i,]<(-0.8))>0){
    abcor<-c(i, abcor)
    abcor<-abcor
  }
}

sigscor<-vesctbaseline[,abcor[3:20]]
genescor<-genesdiab[abcor[20:68],]
X=list(RNASeq=t(genescor),  antibodies=sigscor, cytokines=mfisd7)
Y=cond

list.keepX <- list(RNASeq = c(48, 48), antibodies = c(18, 18), cytokines=c(2,2))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
png('circosPlotbaseline.png')
circ<-circosPlot(MyResult.diablo,  cutoff=0.8,  line = TRUE,var.adj = 0.9, size.variables=0.7, size.labels=2)
dev.off()

#Eotaxin correlations


syeot<-AnnotationDbi::select(hs, 
                          keys =rownames(circ[circ[,"eotaxin"]>0.8,]),
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

ans.goeot <- enrichGO(gene = syeot$ENTREZID, ont = "BP",
                   OrgDb ="org.Hs.eg.db", 
                   readable = TRUE,
                   pvalueCutoff = 0.05)


#BASELINE DEG MICROARRAYS AFRICANS TBS+/TBS- AND ANTIBODIES (FDR<0.05)

microdif<-read.csv("microc1tbs+tbs-.csv")
donorslach<-unique(lachmi$original_id)
colnames(microdif)<-c("probe", donorslach)
rownames(microdif)<-microdif[,1]
microdif<-microdif[,-1]


#ab: we have previously analyzed and obtained that the differentially expressed are: "IgG4 AMA-1[FVO]","IgG4 AMA-1[3D7]", "IgG3 CelTOS", "IgG3 AMA-1[FVO]" "IgG3 AMA-1[3D7]"

abaseline<-ab[ab$t_point=="C_1",]
abaseline<-abaseline[abaseline$immune_status=="semi_immune",]
abbaseline<-pcabaseline %>% dplyr::filter(fullname %in% c("IgG4 AMA-1[FVO]","IgG4 AMA-1[3D7]", "IgG3 CelTOS", "IgG3 AMA-1[FVO]", "IgG3 AMA-1[3D7]"))

vesctbaseline<-c()
for(i in unique(abbaseline$fullname)){
  abanti<-abbaseline[abbaseline$fullname==i,]
  vesctbaseline<-cbind(vesctbaseline, abanti$log10_mfi)
  vesctbaseline<-vesctbaseline
}

colnames(vesctbaseline)<-c(unique(abbaseline$fullname))
rownames( vesctbaseline)<-unique(abbaseline$original_id)

rownames(vesctbaseline)<-gsub("_", "-", rownames(vesctbaseline))
vesctbaseline<-vesctbaseline[colnames(microdif),]


#DIABLO

X=list(mRNA=t(microdif),  ab=vesctbaseline)
abantidon<-abanti%>%dplyr::filter(gsub("_", "-", abanti$original_id)%in%donorslach)
Y<-abantidon$tbs_pos
list.keepX <- list(mRNA = c(1819, 1819), ab = c(5, 5))
MyResult.diablo2 <- block.plsda(X, Y)
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
circ<-circosPlot(MyResult.diablo, cutoff=0.8,  line = TRUE,var.adj = 100)


abcor<-c()
for (i in colnames(circ)){
  if (sum(circ[i,]>0.8) >0 | sum(circ[i,]<(-0.8))>0){
    abcor<-c(i, abcor)
    abcor<-abcor
  }
}

sigscor<-vesctbaseline[,abcor[1:5]]
genescor<-genesdiab[abcor[6:700],]
X=list(RNASeq=t(genescor),  antibodies=sigscor, cytokines=mfisd7)
Y=cond

list.keepX <- list(RNASeq = c(695, 695), antibodies = c(5, 5))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
circ<-circosPlot(MyResult.diablo, cutoff=0.8,  line = TRUE,var.adj = 100)

gencyr<-fData(lachmi)%>%filter(fData(lachmi)$Probe.Set.ID %in% rownames(circ[circ[,"IgG4 AMA-1[3D7]"]>0.8,]))
hs<-org.Hs.eg.db
sy<-AnnotationDbi::select(hs,
                          keys =gencyr$Gene.Symbol,
                          columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")

ans.go <- enrichGO(gene = sy$ENTREZID, ont = "BP",OrgDb ="org.Hs.eg.db",readable = TRUE,pvalueCutoff = 0.05)
#Pick the genes that are involved in immune pathways

g<-c("BST2 ", "CD80 ", "EIF2AK2 ", "EXOSC4 ", "MX2 ", "OAS1 ","OAS3 ", "OAS2 ", "CCR8 ","IL12B ", "TRIM25 ", "UBE2J1 ")


ls<-fData(lachmi)%>%filter(Gene.Symbol %in% g)
newd<-c()
for (i in ls$Probe.Set.ID){
  i<-as.character(i)
  newd<-cbind(newd, t(microdif[i,]))
}

newd<-cbind(newd, vesctbaseline[,4])
colnames(newd)<-c(ls$Gene.Symbol, "log10_mfi")
datgat<-newd %>%
  as_tibble() %>%
  gather(key = "variable", value = "exprs",
         -log10_mfi)
ggplot(datgat, aes(x = exprs, y = log10_mfi))  +
  facet_wrap(~variable,  scales = "free")+
  geom_point( color="#69b3a2")+geom_smooth(method=lm , color="red")+
  stat_cor(method = "pearson", p.digits=2,  label.y = 5)
library(ggplot2)
