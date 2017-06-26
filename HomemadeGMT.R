#Reading in Module information .csv files
setwd("C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Co-Expression")

humanCoexpressionModules<-read.csv("HumanCoexpressionModules.csv")
colnames(humanCoexpressionModules)
#[1] "X"                 "Illumina.Probe.ID" "Ensembl.Gene.ID"   "GeneSymbol_Human"  "Gene.Description" 
#[6] "Chromosome"        "Gene.Start..bp."   "Gene.End..bp."     "Module"  

mouseCoexpressionModules<-read.csv("MouseCoexpressionModules.csv")
colnames(mouseCoexpressionModules)
#[1] "X"                   "Search_Key"          "Transcript"          "GeneSymbol_Human"    "Source_Reference_ID"
#[6] "RefSeq_ID"           "Probe_Id"            "cis_minP_hippo"      "cis_minP_striatum"   "chromosome"         
#[11] "chrStart"            "moduleHippo"         "moduleStriatum"      "kTotalHippo"         "kWithinHippo"       
#[16] "kOutHippo"           "kDiffHippo"          "kTotalStriatum"      "kWithinStriatum"     "kOutStriatum"       
#[21] "kDiffStriatum" 


###
#Making list of genes according to module

humanModules<-split(humanCoexpressionModules$GeneSymbol_Human, humanCoexpressionModules$Module)


#splitting according to hippocampal modules
mouseModules<-split(mouseCoexpressionModules$GeneSymbol_Human, mouseCoexpressionModules$moduleHippo)


#Making gmt files for both

#human gmt
str(humanModules)
names(humanModules[1])
humanModules[1]

temp<-vapply(humanModules$M1, paste, collapse = ", ", character(1L))
temp<-vapply(humanModules[10], paste, collapse = ", ", character(1L))

humanModulesGMT<-matrix(0, 25, 3)
for(i in 1:25){
  humanModulesGMT[i,1]<-names(humanModules[i]) #module names
  humanModulesGMT[,2]<-"NA" #dummy row
  humanModulesGMT[i,3]<-vapply(humanModules[i], paste, collapse="\t", character(1L)) #gene list per module
}

colnames(humanModulesGMT)<-c("Modules", "Dummy", "Genes")
write.table(humanModulesGMT, "humanCoexpressionModulesGMT.gmt", col.names = FALSE, row.names = FALSE, sep="\t")


#mouse gmt
length(mouseModules)
#[1] 30

mouseModulesGMT<-matrix(0, 30, 3)
for(i in 1:30){
  mouseModulesGMT[i,1]<-names(mouseModules[i]) #module names
  mouseModulesGMT[,2]<-"NA" #dummy row
  mouseModulesGMT[i,3]<-vapply(mouseModules[i], paste, collapse="\t", character(1L)) #gene list per module
}

colnames(mouseModulesGMT)<-c("Mouse Modules", "Dummy", "Genes")
write.table(mouseModulesGMT, "mouseCoexpressionModulesGMT.gmt", col.names = FALSE, row.names = FALSE, sep="\t")




###################################################################

# Marker Genes

###################################################################

setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Marker Genes")

dendrogramMarkers<-read.csv("HippoSeq_DendrogramMarkerGenes2.csv")
colnames(dendrogramMarkers)
#[1] "X"                    "enriched"             "depleted"             "gene_id"             
#[5] "GeneSymbol_Mouse"     "uniprot_protein_name" "coronalIshAvailable"  "coronalIshConsistent"
#[9] "GeneSymbol_Human"     "EnrichedVSDepleted"
#
dgGranuleMarkers<-read.csv("HippoSeq_DGGranuleCellMarkerGenes.csv")
colnames(dgGranuleMarkers)
#[1] "X"                    "gene_id"              "GeneSymbol_Mouse"     "uniprot_protein_name" "dg_d_fpkm"           
#[6] "dg_v_fpkm"            "ratio"                "coronalIshExists"     "coronalIshCorrect"  
#
dgMossyMarkers<-read.csv("HippoSeq_DGMossyCellMarkerGenes.csv")
colnames(dgMossyMarkers)
#[1] "X"                    "enriched"             "gene_id"              "GeneSymbol_Mouse"     "uniprot_protein_name"
#[6] "novel"                "coronalIshAvailable"  "coronalIshConsistent"
#

###################################
#Extracting enriched vs depleted for dendrogram makers with gene lists for each value

colnames(dendrogramMarkers)
#[1] "X"                    "enriched"             "depleted"            
#[4] "gene_id"              "GeneSymbol_Mouse"     "uniprot_protein_name"
#[7] "coronalIshAvailable"  "coronalIshConsistent" "GeneSymbol_Human"    
#[10] "EnrichedVSDepleted" 


#Splitting genes according to which module they fall in (enriched vs depleted category)
dendrogramEvsD_genes<-split(dendrogramMarkers$GeneSymbol_Human, dendrogramMarkers$EnrichedVSDepleted)
length(dendrogramEvsD_genes)
#[1] 12

#Splitting genes using mouse gene symbols (first letter uppercase rest lowercase)
dendrogramEvsD_mousegenes<-split(dendrogramMarkers$GeneSymbol_Mouse, dendrogramMarkers$EnrichedVSDepleted)
length(dendrogramEvsD_mousegenes)
#[1] 12

#Creating matrix with correct gmt format (first column is module/pathway names, second is a dummy row or whatever you want, and third is the genes)
dendrogramMarkersEvsD<-matrix(0, 12, 3)
for(i in 1:12){
  dendrogramMarkersEvsD[i,1]<-names(dendrogramEvsD_genes[i]) #enriched versus depleted names
  dendrogramMarkersEvsD[i,2]<-vapply(dendrogramEvsD_mousegenes[i], paste, collapse=", ", character(1L)) #mouse gene list per module
  dendrogramMarkersEvsD[i,3]<-vapply(dendrogramEvsD_genes[i], paste, collapse="\t", character(1L)) #gene list per module
}

colnames(dendrogramMarkersEvsD)<-c("EnrichedVSDepleted", "Mouse Symbol", "Genes")
str(dendrogramMarkersEvsD)

write.table(dendrogramMarkersEvsD, "dendrogramMarkersEnrichedvsDepleted.gmt", col.names = FALSE, row.names = FALSE, sep="\t")


################################
#testing homemade gmt with fgsea

library(fgsea)

temp<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/dendrogramMarkersEnrichedvsDepleted.gmt")
names(temp)
#correct .gmt format, read in correctly

temp[[1]] #Maybe?


## Reading in .rnk file of top adult meta genes

setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA")

adult_rnkinput<-read.csv("MetaAdult_RankedList.csv", header=F)

rankedGeneList_upper<-adult_rnkinput[[2]] #remaking ranked gene list
names(rankedGeneList_upper)<-toupper(adult_rnkinput[[1]]) #naming with uppercase gene symbols to match rat pathways
names(rankedGeneList_upper)[1:5]

temp1<-fgsea(temp, rankedGeneList_upper, nperm=2000, minSize = 1, maxSize = 500)
#okay that worked, think we're good to go


############################################################################
#
#              Running fgsea with all homemade gmts
#
#
setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA")

#Mouse Coexpression
homemadeMouseCoexpression<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Co-Expression/mouseCoexpressionModulesGMT.gmt")

adultMetaGenes_MouseCoexpressiongsea<-fgsea(homemadeMouseCoexpression, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)

adultMetaGenes_MouseCoexpressiongsea$leadingEdge<-vapply(adultMetaGenes_MouseCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(adultMetaGenes_MouseCoexpressiongsea, "fgsea_Adult_MouseCoexpression.csv")


#Human Coexpression
homemadeHumanCoexpression<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Co-Expression/humanCoexpressionModulesGMT.gmt")

adultMetaGenes_HumanCoexpressiongsea<-fgsea(homemadeHumanCoexpression, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)

sigadultHuman<-subset(adultMetaGenes_HumanCoexpressiongsea, adultMetaGenes_HumanCoexpressiongsea$pval < 0.05)

adultMetaGenes_HumanCoexpressiongsea$leadingEdge<-vapply(adultMetaGenes_HumanCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(adultMetaGenes_HumanCoexpressiongsea, "fgsea_Adult_HumanCoexpression.csv")


#Dendrogram Markers Encrichment vs Depletion
homemadeDendrogramMarkers<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/dendrogramMarkersEnrichedvsDepleted.gmt")

adultMetaGenes_DendrogramMarkersgsea<-fgsea(homemadeDendrogramMarkers, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)



#######################################################################################
####### Remaking dendrogram marker module adding in markers for granule and mossy cells

#Condensing genes to list for both dg granule and dg mossy files

colnames(dgGranuleMarkers)
#[1] "X"                    "gene_id"              "GeneSymbol_Mouse"     "uniprot_protein_name"
#[5] "dg_d_fpkm"            "dg_v_fpkm"            "ratio"                "coronalIshExists"    
#[9] "coronalIshCorrect"  

colnames(dgMossyMarkers)
#[1] "X"                    "enriched"             "gene_id"              "GeneSymbol_Mouse"    
#[5] "uniprot_protein_name" "novel"                "coronalIshAvailable"  "coronalIshConsistent"

mossy<-split(toupper(dgMossyMarkers$GeneSymbol_Mouse), "Mossy")
(mossy[[1]])[1:5]
#[1] "6330406I15RIK" "AJAP1"         "AP2S1"         "ASS1"          "B230216N24RIK"

granule<-split(toupper(dgGranuleMarkers$GeneSymbol_Mouse), "Granule")
(granule[[1]])[1:5]
#[1] "AKAP5"  "DES"    "TNNT2"  "TAGLN2" "IFITM2"

#######################

#  Remaking dendrogram gmt with granule and mossy cell markers

dendrogramandCellTypeMarkers<-matrix(0, 14, 3)
for(i in 1:12){
  dendrogramandCellTypeMarkers[i,1]<-names(dendrogramEvsD_genes[i]) #enriched versus depleted names
  dendrogramandCellTypeMarkers[13,1]<-names(mossy)
  dendrogramandCellTypeMarkers[14,1]<-names(granule)
  dendrogramandCellTypeMarkers[(1:14),2]<-"NA"  #dummy column
  dendrogramandCellTypeMarkers[i,3]<-vapply(dendrogramEvsD_genes[i], paste, collapse="\t", character(1L)) #gene list per module
  dendrogramandCellTypeMarkers[13,3]<-vapply(mossy, paste, collapse="\t", character(1L))
  dendrogramandCellTypeMarkers[14,3]<-vapply(granule, paste, collapse="\t", character(1L))
  }

colnames(dendrogramandCellTypeMarkers)<-c("EnrichedVSDepleted", "dummy", "Genes")
#
write.table(dendrogramandCellTypeMarkers, "dendrogram_CellTypeMarkers.gmt", col.names = FALSE, row.names = FALSE, sep="\t")

###################################################################################
#
#         Re-rerunning fgsea with updated dendrogram cell marker gmt
#
###################################################################################

#Dendrogram Cell Markers Encrichment vs Depletion
homemadeDendrogram_CellMarkers<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/dendrogram_CellTypeMarkers.gmt")

adultMetaGenes_Dendrogram_CellMarkersgsea<-fgsea(homemadeDendrogram_CellMarkers, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)

#testing a second run to see if it changes
temp<-fgsea(homemadeDendrogram_CellMarkers, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)
#It does, okay that's to be expected

str(adultMetaGenes_Dendrogram_CellMarkersgsea)

#outputting fgsea results for cell marker enrichment

adultMetaGenes_Dendrogram_CellMarkersgsea$leadingEdge<-vapply(adultMetaGenes_Dendrogram_CellMarkersgsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(adultMetaGenes_Dendrogram_CellMarkersgsea, "GSEA_AdultMetaGenes_DendrogramCellMarkers.csv")



########## Making tables of top coexpressin modules/cell enrichment for each fgsea output


topMarkers <- adultMetaGenes_Dendrogram_CellMarkersgsea[head(order(pval), n=15)][order(NES), pathway]
## 

png("MetaAdult_TopCellMarkerEnrichmentvsDepletion.png", width=800)
plotGseaTable(homemadeDendrogram_CellMarkers[topMarkers], rankedGeneList_upper,
              adultMetaGenes_Dendrogram_CellMarkersgsea, gseaParam=0.5)
dev.off()



######## Mouse coexpression
adultMetaGenes_MouseCoexpressiongsea<-fgsea(homemadeMouseCoexpression, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)


topMouseCoexpression <- adultMetaGenes_MouseCoexpressiongsea[head(order(pval), n=15)][order(NES), pathway]
## 

png("MetaAdult_TopMouseCoexpressionModules.png")
plotGseaTable(homemadeMouseCoexpression[topMouseCoexpression], rankedGeneList_upper,
              adultMetaGenes_MouseCoexpressiongsea, gseaParam=0.5)
dev.off()


######## Human coexpression
adultMetaGenes_HumanCoexpressiongsea<-fgsea(homemadeHumanCoexpression, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)


topHumanCoexpression <- adultMetaGenes_HumanCoexpressiongsea[head(order(pval), n=15)][order(NES), pathway]
## 

png("MetaAdult_TopHumanCoexpressionModules.png")
plotGseaTable(homemadeHumanCoexpression[topHumanCoexpression], rankedGeneList_upper,
              adultMetaGenes_HumanCoexpressiongsea, gseaParam=0.5)
dev.off()

######################################################
#   Outputting picture of most sig pathway/module enricment
#  

### region/cell type enrichment

#Looking at most significantly enriched gene sets
topMarkers <- adultMetaGenes_Dendrogram_CellMarkersgsea[head(order(pval), n=15)]


(topMarkers$pathway[1])
#[1] "\"dg_v__VS__dg_d\""


png("dg_v__VS__dg_d_AdultEnrichmentvsDepletion.png")
plotEnrichment(homemadeDendrogram_CellMarkers[["\"dg_v__VS__dg_d\""]], rankedGeneList_upper)
dev.off()


### region/cell type enrichment

#Looking at highly enriched gene sets
topMarkers <- adultMetaGenes_Dendrogram_CellMarkersgsea[head(order(pval), n=15)]


(topMarkers$pathway[1])
#[1] "\"dg_v__VS__dg_d\""


png("dg_v__VS__dg_d_AdultEnrichmentvsDepletion.png")
plotEnrichment(homemadeDendrogram_CellMarkers[["\"dg_v__VS__dg_d\""]], rankedGeneList_upper)
dev.off()


### Top Mouse Co-expression module

#Looking at highly enriched gene sets
topMouse <- adultMetaGenes_MouseCoexpressiongsea[head(order(pval), n=15)]


(topMouse$pathway[1])
#[1] "\"paleturquoise\""


png("paleturquoise_AdultMouseCoexpressionModule.png")
plotEnrichment(homemadeMouseCoexpression[["\"paleturquoise\""]], rankedGeneList_upper)
dev.off()

### Top Human Co-expression module

#Looking at highly enriched gene sets
topHuman <- adultMetaGenes_HumanCoexpressiongsea[head(order(pval), n=15)]


(topHuman$pathway[1])
#[1] "\"M11\""


png("M11_AdultHumanCoexpressionModule.png")
plotEnrichment(homemadeHumanCoexpression[["\"M11\""]], rankedGeneList_upper)
dev.off()




#########################################
########################################
########### Looking at overlap between pathways/modules of various .gmts

str(mossy) 
str(granule)
#both character lists of 1

#Using VennDiagram package

#install.packages("VennDiagram")
library("VennDiagram")

#putting two lists together
temp<-c(mossy, granule)

#using venndiagram function to calculate overlap
temp1<-calculate.overlap(temp)
str(temp1)
#List of 3
#$ a1: chr [1:61] "6330406I15RIK" "AJAP1" "AP2S1" "ASS1" ...
#$ a2: chr [1:128] "AKAP5" "DES" "TNNT2" "TAGLN2" ...
#$ a3: chr(0) 
#
# cool, and there's no overlap which should be correct?

#However the VennDiagram package can't handle very large lists like the ones I have...


############################ observing overlap through intersect() and %in%

#######
# Reading in original rat gmt file used to run gsea

rat_Pathways<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Rattus_norvegicus_GSEA_GO_sets_bp_symbols_highquality_April_2015.gmt")


temp<-intersect(rat_Pathways[[i]], homemadeMouseCoexpression[[i]]) 
#This works for comparing two module/pathways but I want to compare all at once
#This calls for a forloop



#Using %in% to get sum of genes in each
i<-1
rat_Pathways[[i]]%in%homemadeMouseCoexpression[[i]]
rat_Pathways[[i]][rat_Pathways[[i]]%in%homemadeMouseCoexpression[[i]]]
sum(rat_Pathways[[i]]%in%homemadeMouseCoexpression[[i]])
#[1] 1
sum(rat_Pathways[[i]]%in%homemadeMouseCoexpression[[i]])/length(rat_Pathways[[i]])
#[1] 0.09090909

sum(homemadeMouseCoexpression[[i]]%in%rat_Pathways[[i]])/length(homemadeMouseCoexpression[[i]])
#[1] 0.001923077

length(rat_Pathways)
#[1] 2761



############################## Performing intersect for each mouse coexpression module against the entire dataset of rat pathways

#forlooping intersect to compare each mouse module gene list to each rat_pathway gene list


################### Test Run

#making list of vectors with length of rat pathway
temp<-vector("list", 5)
for(i in 1:5){
     temp[[i]]<-c(intersect(homemadeMouseCoexpression[[1]], rat_Pathways[[i]]))
   }  #test run and it worked

#Test # 2, implementing a matrix so that each column represents a Mouse Co-expression module
temp1<-matrix(0, 5, 1)
for(i in 1:5){
  temp[[i]]<-(intersect(homemadeMouseCoexpression[[1]], rat_Pathways[[i]]))
  temp1[i,]<-vapply(temp[i], paste, collapse=",", character(1L))
  }  #All good should be ready to run

length(homemadeMouseCoexpression)
#[1] 30 #how many columns will be needed


#########################################################
#        ()_()
#        =\"/=   Mouse Coexpression/Rat pathway overlap
#          O
#########################################################

temp<-vector("list", length(rat_Pathways))
intersectionMouseCoexp_RatPathway<-matrix(0, length(rat_Pathways), 30)

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeMouseCoexpression)){
  temp[[i]]<-(intersect(homemadeMouseCoexpression[[l]], rat_Pathways[[i]]))
  intersectionMouseCoexp_RatPathway[i,l]<-vapply(temp[i], paste, collapse=",", character(1L))
}  
}

colnames(intersectionMouseCoexp_RatPathway)<-names(homemadeMouseCoexpression)
row.names(intersectionMouseCoexp_RatPathway)<-names(rat_Pathways)

write.csv(intersectionMouseCoexp_RatPathway, "intersectionMouseCoexp_RatGOPathways.csv")

##### Making dataframe of sums of overlapping gene

intersectionMouseCoexp_RatPathway_SUMS<-matrix(0, length(rat_Pathways), length(homemadeMouseCoexpression))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeMouseCoexpression)){
    intersectionMouseCoexp_RatPathway_SUMS[i,l]<-sum(rat_Pathways[[i]]%in%homemadeMouseCoexpression[[l]])
  }  
}

colnames(intersectionMouseCoexp_RatPathway_SUMS)<-names(homemadeMouseCoexpression)
row.names(intersectionMouseCoexp_RatPathway_SUMS)<-names(rat_Pathways)

write.csv(intersectionMouseCoexp_RatPathway_SUMS, "SUMSintersectionMouseCoexp_RatGOPathways.csv")


#################################################################################
#  Now running for Human co-expression modules
#     
#     www
#    ('o')   
#    (/_\)
#     | |
#     

temp<-vector("list", length(rat_Pathways))
intersectionHumanCoexp_RatPathway<-matrix(0, length(rat_Pathways), length(homemadeHumanCoexpression))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeHumanCoexpression)){
    temp[[i]]<-(intersect(homemadeHumanCoexpression[[l]], rat_Pathways[[i]]))
    intersectionHumanCoexp_RatPathway[i,l]<-vapply(temp[i], paste, collapse=",", character(1L))
  }  
}

colnames(intersectionHumanCoexp_RatPathway)<-names(homemadeHumanCoexpression)
row.names(intersectionHumanCoexp_RatPathway)<-names(rat_Pathways)

write.csv(intersectionHumanCoexp_RatPathway, "intersectionHumanCoexp_RatGOPathways.csv")

######### Making dataframe of sums of overlapping gene

intersectionHumanCoexp_RatPathway_SUMS<-matrix(0, length(rat_Pathways), length(homemadeHumanCoexpression))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeHumanCoexpression)){
    intersectionHumanCoexp_RatPathway_SUMS[i,l]<-sum(rat_Pathways[[i]]%in%homemadeHumanCoexpression[[l]])
  }  
}

colnames(intersectionHumanCoexp_RatPathway_SUMS)<-names(homemadeHumanCoexpression)
row.names(intersectionHumanCoexp_RatPathway_SUMS)<-names(rat_Pathways)

write.csv(intersectionHumanCoexp_RatPathway_SUMS, "SUMSintersectionHumanCoexp_RatGOPathways.csv")


#############################################################
#
#####################################        
#                                              
############################ Region/Cell type module e vs d compared to Rat Pathways gmt

temp<-vector("list", length(rat_Pathways))
intersectionHippoSeqModules_RatPathway<-matrix(0, length(rat_Pathways), length(homemadeDendrogram_CellMarkers))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeDendrogram_CellMarkers)){
    temp[[i]]<-(intersect(homemadeDendrogram_CellMarkers[[l]], rat_Pathways[[i]]))
    intersectionHippoSeqModules_RatPathway[i,l]<-vapply(temp[i], paste, collapse=",", character(1L))
  }  
}

colnames(intersectionHippoSeqModules_RatPathway)<-names(homemadeDendrogram_CellMarkers)
row.names(intersectionHippoSeqModules_RatPathway)<-names(rat_Pathways)

write.csv(intersectionHippoSeqModules_RatPathway, "intersectionHippoSeqModules_RatGOPathways.csv")


###############################
#Making dataframe with number of genes present overlapping between cell marker/region enrichment vs rat pathway sum

intersectionHippoSeqModules_RatPathway_SUMS<-matrix(0, length(rat_Pathways), length(homemadeDendrogram_CellMarkers))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeDendrogram_CellMarkers)){
    intersectionHippoSeqModules_RatPathway_SUMS[i,l]<-sum(rat_Pathways[[i]]%in%homemadeDendrogram_CellMarkers[[l]])
  }  
}

colnames(intersectionHippoSeqModules_RatPathway_SUMS)<-names(homemadeDendrogram_CellMarkers)
row.names(intersectionHippoSeqModules_RatPathway_SUMS)<-names(rat_Pathways)
write.csv(intersectionHippoSeqModules_RatPathway_SUMS, "SUMSintersectionHippoSeqModules_RatGOPathways.csv")



###################################################################################

######## Running homemade .gmts with P14 meta genes

rankedGeneList_P14input<-read.csv("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Gene Lists (rnk files)/MetaP14_RankedList.csv", header=F)
rankedGeneList_P14<-rankedGeneList_P14input[[2]]
names(rankedGeneList_P14)<-toupper(rankedGeneList_P14input[[1]])
rankedGeneList_P14[1:5]
#RGD1310311  LOC501089       ABP1     OLR257    OLR1566 
#2.148696   2.069144   2.043404   2.003478   1.878333

#Mouse Coexpression
homemadeMouseCoexpression<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Co-Expression/mouseCoexpressionModulesGMT.gmt")

P14MetaGenes_MouseCoexpressiongsea<-fgsea(homemadeMouseCoexpression, rankedGeneList_P14, nperm=10000, minSize = 1, maxSize = 1000)

P14MetaGenes_MouseCoexpressiongsea$leadingEdge<-vapply(P14MetaGenes_MouseCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(P14MetaGenes_MouseCoexpressiongsea, "fgsea_P14_MouseCoexpression.csv")


#Human Coexpression
homemadeHumanCoexpression<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Co-Expression/humanCoexpressionModulesGMT.gmt")

P14MetaGenes_HumanCoexpressiongsea<-fgsea(homemadeHumanCoexpression, rankedGeneList_P14, nperm=10000, minSize = 1, maxSize = 1000)

sigP14Human<-subset(P14MetaGenes_HumanCoexpressiongsea, P14MetaGenes_HumanCoexpressiongsea$pval < 0.05)
sigP14Human

P14MetaGenes_HumanCoexpressiongsea$leadingEdge<-vapply(P14MetaGenes_HumanCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(P14MetaGenes_HumanCoexpressiongsea, "fgsea_P14_HumanCoexpression.csv")


#Dendrogram Cell Markers/region Encrichment vs Depletion
homemadeDendrogram_CellMarkers<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/dendrogram_CellTypeMarkers.gmt")

P14MetaGenes_Dendrogram_CellMarkersgsea<-fgsea(homemadeDendrogram_CellMarkers, rankedGeneList_P14, nperm=10000, minSize = 1, maxSize = 1000)


#outputting fgsea results for cell marker enrichment

P14MetaGenes_Dendrogram_CellMarkersgsea$leadingEdge<-vapply(P14MetaGenes_Dendrogram_CellMarkersgsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(P14MetaGenes_Dendrogram_CellMarkersgsea, "fgsea_P14MetaGenes_DendrogramCellMarkers.csv")

######

############# Creating Tables for P14 homemade gmt fgsea output

#Subsetting most significant pathways and ordering by NES
topMarkersP14 <- P14MetaGenes_Dendrogram_CellMarkersgsea[head(order(pval), n=15)][order(NES), pathway]
## 

png("MetaP14_TopCellMarkerEnrichmentvsDepletion.png", width=800)
plotGseaTable(homemadeDendrogram_CellMarkers[topMarkersP14], rankedGeneList_P14,
              P14MetaGenes_Dendrogram_CellMarkersgsea, gseaParam=0.5)
dev.off()



######## Mouse coexpression

topMouseCoexpressionP14 <- P14MetaGenes_MouseCoexpressiongsea[head(order(pval), n=15)][order(NES), pathway]
## 

png("MetaP14_TopMouseCoexpressionModules.png")
plotGseaTable(homemadeMouseCoexpression[topMouseCoexpressionP14], rankedGeneList_P14,
              P14MetaGenes_MouseCoexpressiongsea, gseaParam=0.5)
dev.off()


######## Human coexpression
P14MetaGenes_HumanCoexpressiongsea<-fgsea(homemadeHumanCoexpression, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 1000)


topHumanCoexpression <- adultMetaGenes_HumanCoexpressiongsea[head(order(pval), n=15)][order(NES), pathway]
## 

png("MetaAdult_TopHumanCoexpressionModules.png")
plotGseaTable(homemadeHumanCoexpression[topHumanCoexpression], rankedGeneList_upper,
              adultMetaGenes_HumanCoexpressiongsea, gseaParam=0.5)
dev.off()



