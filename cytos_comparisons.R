#CYTOKINES

#Obtaining cytokines dataset and filtering 
cyto <- chmi.02_phen.argument_cyto_data(dat_imp)
cyto<-cyto[cyto$cytoaset=="L1",]
cyto_names<-unique(cyto$analyte)

#PCA WITH RAW DATA
cytop<-read.csv("/Users/carlasanchezalmirall/datasets/chmiddbb/03_chmi_cyto/01_v16.12.19_cyto_chmi_dat.csv", sep=";")
cytop<-cytop[cytop$dataset=="L1",]
cytokines<-cytop[,10:39]
cytokines<-cytokines %>% mutate_all(as.numeric)  
# transform the counts to log2 scale 
M <- log2(cytokines + 1)
#compute PCA 
pcaResults <- prcomp(M)
pcplot<-data.frame(cytop$t_point, pcaResults$x[,1:2])
GROUP<-cyto$gender
qplot(x=PC1, y=PC2, data=pcplot, colour=GROUP)

#COMPARING AFRICANS AND EUROPEANS

#Chose the t_point we want to compare("C_1", "D7", "D13", "D19", "D28")
#For time point C_1 threshold is 0.1 but at other points it is changed to 0.05.

time<-
cytonew<-cyto[cyto$t_point==time,]
sigcytoc1<-c()
for(i in cyto_names){
  cytonew<-cyto[cyto$analyte==i]
  pshap<-shapiro.test(cytonew$log10_mfi)$p.value
  if (pshap <0.05){
    wilc<-wilcox.test(cytonew$log10_mfi ~ cytonew$immune_status.x)$p.value
    wilc<-p.adjust(wilc)
    if (wilc< 0.1){
      sigcytoc1<-c(i,sigcytoc1)
      sigcytoc1<-sigcytoc1}
  }
  else{
    variance<-var.test(cytonew$log10_mfi ~ cytonew$immune_status.x)$p.value 
    if (variance< 0.05){
      ttest<-t.test(cytonew$log10_mfi ~ cytonew$immune_status.x, var.equal=FALSE)$p.value
      ttest<-p.adjust(ttest)
    }
    else{
      ttest<-t.test(cytonew$log10_mfi ~ cytonew$immune_status.x, var.equal=TRUE)$p.value
      ttest<-p.adjust(ttest)}
    if (ttest< 0.1){
      sigcytoc1<-c(i,sigcytoc1)
      sigcytoc1<-sigcytoc1}
  }}


#AMONG AFRICANS COMPARING TBS+/TBS- 

#Chose the t_point we want to compare("C_1", "D7", "D13", "D19", "D28")


time<-
cyto<-cyto[cyto$immune_status.x=="semi_immune",]
cytonew<-cyto[cyto$t_point==time,]
sigcytoc1<-c()
for(i in cyto_names){
  cytonew<-cyto[cyto$analyte==i]
  pshap<-shapiro.test(cytonew$log10_mfi)$p.value
  if (pshap <0.05){
    wilc<-wilcox.test(cytonew$log10_mfi ~ cytonew$tbs_pos.x)$p.value
    wilc<-p.adjust(wilc)
    if (wilc< 0.05){
      sigcytoc1<-c(i,sigcytoc1)
      sigcytoc1<-sigcytoc1}
  }
  else{
    variance<-var.test(cytonew$log10_mfi ~ cytonew$tbs_pos.x)$p.value 
    if (variance< 0.05){
      ttest<-t.test(cytonew$log10_mfi ~ cytonew$tbs_pos.x, var.equal=FALSE)$p.value
      ttest<-p.adjust(ttest)
    }
    else{
      ttest<-t.test(cytonew$log10_mfi ~ cytonew$tbs_pos.x, var.equal=TRUE)$p.value
      ttest<-p.adjust(ttest)}
    if (ttest< 0.05){
      sigcytoc1<-c(i,sigcytoc1)
      sigcytoc1<-sigcytoc1}
  }}


