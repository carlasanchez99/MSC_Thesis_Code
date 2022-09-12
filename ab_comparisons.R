#Antibodies

#Obtaining antibodies dataset and filtering 
ab <- chmi.02_phen.argument_ab_data(dat_imp)
ab<-chmi.update_labels.dat_imp_ab_lb(ab, fold_change = FALSE)
ab <- ab %>%
  as_tibble()
ab<-ab[ab$dataset=="L1",]
ab$fullname<-paste(ab$t_igg, ab$antigen_lb, sep=" ")

#Using raw data for the PCA
abpca<-read.csv("~/02_v20.11.16_chmi_ab_dat copia.csv", sep=";")
rownames(abpca)<-paste(abpca$original_id, abpca$t_point, abpca$t_igg, sep="_")
pcadata<-abpca[7:27]
pcadata<-select(pcadata, -pf_rh4, -pf_rh5, -pf_rh1)
pcadata<-na.omit(pcadata)
pcaResults <- prcomp(pcadata)

abpcafilt<-abpca[rownames(pcadata),]
GROUP<-abpcafilt$tbs_pos
pcplot<-data.frame(GROUP, pcaResults$x[,1:2])
qplot(x=PC1, y=PC2, data=pcplot, colour=GROUP)

fviz_pca_biplot(pcaResults, 
                col.ind = abpcafilt$t_igg,palette = "lancet", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                geom.ind = "point", pointshape = 20, title="", legend.title = "", pointsize = 2) 


#COMPARING AFRICANS AND EUROPEANS

#Chose the t_point we want to compare("C_1", "D7", "D13", "D19", "D28")
time<-
abaseline<-ab[ab$t_point==time,]
ab_names<-unique(abaseline$fullname)
sigbaseline<-c()
for(i in ab_names){
  datnew<-abaseline[abaseline$fullname==i,]
  pshap<-shapiro.test(datnew$log10_mfi)$p.value
  if (pshap <0.05){
    wilc<-wilcox.test(datnew$log10_mfi ~ datnew$immune_status)$p.value
    wilc<-p.adjust(wilc)
    if (wilc< 0.05){
      sigbaseline<-c(unique(datnew$fullname),sigbaseline)
      sigbaseline<-sigbaseline
    }
  }
  else{
    variance<-var.test(datnew$log10_mfi ~ datnew$immune_status)$p.value 
    if (variance< 0.05){
      ttest<-t.test(datnew$log10_mfi ~ datnew$immune_status, var.equal=FALSE)$p.value 
      ttest<-p.adjust(ttest)
    }
    else{
      ttest<-t.test(datnew$log10_mfi ~ datnew$immune_status, var.equal=TRUE)$p.value
      ttest<-p.adjust(ttest)}
    if (ttest< 0.05){
      sigbaseline<-c(unique(datnew$fullname),sigbaseline)
      sigbaseline<-sigbaseline}
  }}


#AMONG AFRICANS COMPARING TBS+/TBS- 

#Chose the t_point we want to compare("C_1", "D7", "D13", "D19", "D28")

time<-
abaseline<-ab[ab$t_point==time,]
abaseline<-abaseline[abaseline$immune_status=="semi_immune",]
ab_names<-unique(abaseline$fullname)
abafricans<-c()
for(i in ab_names){
  datnew<-abaseline[abaseline$fullname==i,]
  pshap<-shapiro.test(datnew$log10_mfi)$p.value
  if (pshap <0.05){
    wilc<-wilcox.test(datnew$log10_mfi ~ datnew$tbs_pos)$p.value
    wilc<-p.adjust(wilc)
    if (wilc< 0.05){
      abafricans<-c(unique(datnew$fullname),abafricans)
      abafricans<-abafricans
    }
  }
  else{
    variance<-var.test(datnew$log10_mfi ~ datnew$tbs_pos)$p.value 
    if (variance< 0.05){
      ttest<-t.test(datnew$log10_mfi ~ datnew$tbs_pos, var.equal=FALSE)$p.value 
      ttest<-p.adjust(ttest)
    }
    else{
      ttest<-t.test(datnew$log10_mfi ~ datnew$tbs_pos, var.equal=TRUE)$p.value
      ttest<-p.adjust(ttest)}
    if (ttest< 0.05){
      abafricans<-c(unique(datnew$fullname),abafricans)
      abafricans<-abafricans}
  }}

