#CLINICAL VARIABLES


#Obtain clinical dataset and filter 
clin <- dat_imp[[1]]
clin <- clin %>%
  as_tibble()
clinical<-clin[clin$dataset=="L1",]

for(i in unique(clinical$immune_status)){
  datnew<-clinical[clinical$immune_status==i,]
  cat(i, "( n=", nrow(datnew) ,"(", nrow(datnew)/nrow(clinical)*100, "%", ")", ")")
  cat(" \n") 
  cat("Gender:", " ")
  cat("Males:", nrow(datnew[datnew$gender=="M",]), "(",  nrow(datnew[datnew$gender=="M",])/nrow(datnew)*100, "%", ")")
  cat(" \n")  
  cat("Age: ", " \t")
  cat("-Mean: ", mean(datnew$age), "       -Sd: ", sd(datnew$age))
  cat(" \n")  
  cat("HbAA:"," \t", nrow(datnew[datnew$hb_status=="AA",]), "(", nrow(datnew[datnew$hb_status=="AA",])/nrow(datnew)*100, "%", ")")
  cat(" \n") 
  cat("TBS+:"," \t", nrow(datnew[datnew$tbs_pos=="1",]), "(", nrow(datnew[datnew$tbs_pos=="1",])/nrow(datnew)*100, "%", ")")
  cat(" \n") 
  cat("Days TBS+:"," \t", "-Mean:", mean(datnew$ppp_tbs), "-Min:", min(datnew$ppp_tbs), "-Max:", max(datnew$ppp_tbs))
  cat(" \n") 
  cat("PCR+:"," \t", nrow(datnew[datnew$pcr_pos=="1",]), "(", nrow(datnew[datnew$pcr_pos=="1",])/nrow(datnew)*100, "%", ")")
  cat(" \n") 
  cat("Days PCR+:"," \t", "-Mean:", mean(datnew$ppp_mal_rtqpcr), "-Min:", min(datnew$ppp_mal_rtqpcr), "-Max:", max(datnew$ppp_tbs))
  cat(" \n") 
  cat(" \n") 
}


vars<-c("age","ppp_tbs", "ppp_mal_rtqpcr")
clinical<-as.data.frame(clinical)


for(i in vars){
  var<-clinical[,i]
  pshap<- shapiro.test(var)$p.value
  if (pshap <0.05){
    wilc<-wilcox.test(var ~ clinical$immune_status)$p.value
    wilc<-p.adjust(wilc)
    if (wilc< 0.05){
      cat(i, " \t", " wilcox padj:",wilc)}
  }
  else{
    variance<-var.test(var ~ clinical$immune_status)$p.value 
    if (variance< 0.05){
      ttest<-t.test(var ~ clinical$immune_status, var.equal=FALSE)$p.value
      ttest<-p.adjust(ttest)
    }
    else{
      ttest<-t.test(var ~ clinical$immune_status, var.equal=TRUE)$p.value 
      ttest<-p.adjust(ttest)}
    if (ttest< 0.05){
      cat(i, " \t", "t.test padj:",ttest)}
  }
}

fisher.test(table(clinical$gender,clinical$immune_status))
fisher.test(table(clinical$hb_status,clinical$immune_status))
fisher.test(table(clinical$tbs_pos,clinical$immune_status))
fisher.test(table(clinical$pcr_pos,clinical$immune_status))
