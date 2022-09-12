#Hematologcial data

#Obtaining hemo dataset and filtering 

hemo <- chmi.02_phen.argument_hemo_data(dat_imp)
hemo <- hemo %>%
  as_tibble()

#PCA

hemo2<-hemo[,57:62]
hemo2<-cbind(hemo$original_id, hemo2)
hemo2<-na.omit(hemo2)
colnames(hemo2)<-c("original_id", colnames(hemo2)[-1])

# transform the counts to log2 scale 
M <- log2(hemo2[-1] + 1)
#compute PCA 
pcaResults <- prcomp(M)
pcplot<-data.frame(hemo2$original_id , pcaResults$x[,1:2])
Original_id<-hemo2$original_id
qplot(x=PC1, y=PC2, data=pcplot, colour=Original_id)

#COMPARING AFRICANS AND EUROPEANS

#Chose the t_point we want to compare("C_1", "D13", "D28")

time<-
hemon<-hemo[hemo$t_point==time,]
hemo_names<-c("hb", "leucos", "neutros", "lymphos", "platelets", "eos")
for(i in hemo_names){
  dat<-hemon[i]
  dat<-cbind(hemon$immune_status, dat)
  pshap<-shapiro.test(dat[,2])$p.value
  if (pshap <0.05){
    wilc<-wilcox.test(dat[,2]~dat[,1])$p.value
    wilc<-p.adjust(wilc)
    if (wilc< 0.05){
      print(i)}
    else{
      variance<-var.test(dat[,2]~dat[,1])$p.value 
      if (variance< 0.05){
        ttest<-t.test(dat[,2]~dat[,1], var.equal=FALSE)$p.value 
        ttest<-p.adjust(ttest)
      }
      else{
        ttest<-t.test(dat[,2]~dat[,1], var.equal=TRUE)$p.value
        ttest<-p.adjust(ttest)
      }
      if (ttest< 0.05){
        print(i)}
    }}
}

#AMONG AFRICANS COMPARING TBS+/TBS- 

#Chose the t_point we want to compare("C_1", "D13", "D28")

time<-
hemon<-hemo[hemo$t_point==time,]
hemon<-hemo[hemo$immune_status=="semi_immune",]
hemo_names<-c("hb", "leucos", "neutros", "lymphos", "platelets", "eos")
for(i in hemo_names){
  dat<-hemon[i]
  dat<-cbind(hemon$tbs_pos, dat)
  pshap<-shapiro.test(dat[,2])$p.value
  if (pshap <0.05){
    wilc<-wilcox.test(dat[,2]~dat[,1])$p.value
    wilc<-p.adjust(wilc)
    if (wilc< 0.05){
      print(i)}
    else{
      variance<-var.test(dat[,2]~dat[,1])$p.value 
      if (variance< 0.05){
        ttest<-t.test(dat[,2]~dat[,1], var.equal=FALSE)$p.value 
        ttest<-p.adjust(ttest)
      }
      else{
        ttest<-t.test(dat[,2]~dat[,1], var.equal=TRUE)$p.value
        ttest<-p.adjust(ttest)
      }
      if (ttest< 0.05){
        print(i)}
    }}
}

