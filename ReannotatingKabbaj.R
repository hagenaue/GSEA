#Reannotating Kabbaj & Evans (2004) Table 1 of HR/LR differences:

KabbaJTable1<-read.csv("Kabbaj_Table1.csv", header=T, stringsAsFactors = F)

library(org.Rn.eg.db)

xx <- as.list(org.Rn.egACCNUM2EG)

EntrezGeneID<-unlist(xx, use.names=FALSE)
GenBankAccNum<-rep(names(xx), lengths(xx))
table(EntrezGeneID)
#but each EntrezGeneID maps to many accession numbers

GenBankAccNumToEntrezGeneID<-data.frame(GenBankAccNum, EntrezGeneID, stringsAsFactors=F)

NumberOfEntrezPerAccNum<-lengths(xx)
GenBankAccNum<-names(xx)

AccNumVsNumberOfEntrezID<-data.frame(GenBankAccNum, NumberOfEntrezPerAccNum, stringsAsFactors=F)
table(NumberOfEntrezPerAccNum)

library(plyr)
AccNumToEntrezID_NumOfDups<-join(GenBankAccNumToEntrezGeneID, AccNumVsNumberOfEntrezID, by="GenBankAccNum", type="left")

colnames(KabbaJTable1)[5]<-"GenBankAccNum"
is.data.frame(KabbaJTable1)
str(KabbaJTable1)

KabbaJTable1_ReannotatedProbes<-join(KabbaJTable1, AccNumToEntrezID_NumOfDups, by="GenBankAccNum", type="left")


x <- org.Rn.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

xx[1]
[1] "Slc39a4l"
xx$"100360501"
[1] "Rnh1"

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))
# 1 
# 42306
#A 1:1 mapping. Nice...

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)


KabbaJTable1_ReannotatedProbes2<-join(KabbaJTable1_ReannotatedProbes, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")


write.csv(KabbaJTable1_ReannotatedProbes2, "KabbaJTable1_ReannotatedProbes2.csv")