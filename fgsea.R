#Installing fgsea
source("https://bioconductor.org/biocLite.R")
biocLite("fgsea")
library("fgsea")


##### Example stuff from fgsea package
#data(exampleRanks)
#data(examplePathways)
#ranks <- sort(exampleRanks, decreasing=TRUE)
#es <- calcGseaStat(ranks, na.omit(match(examplePathways[[1]], names(ranks))))
#
#fgsea(pathways, stats, nperm, minSize = 1, maxSize = Inf, nproc = 0,
#      gseaParam = 1, BPPARAM = NULL)
#
####################################################################


## Reading in .gmt file for Rat GO pathways

setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA")

rat_Pathways<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Rat_AllPathways_June_01_2017_symbol.gmt")
names(rat_Pathways)

#Reading in ranked list

adult_rnkinput<-read.csv("MetaAdult_RankedList.csv", header=F)

#Making rnk input into a named numeric vector

rankedGeneList<-adult_rnkinput[[2]] #subsetting numeric info from rnk input to rankedGeneList
names(rankedGeneList)<-(adult_rnkinput[[1]]) #adding gene names from rnk input to rankedGeneList


#Running test fgsea

test<-fgsea(rat_Pathways, rankedGeneList, nperm=1000, minSize = 1, maxSize = 500) #generic setting for minsize and maxsize for this test, probably better to use more permutations than 1000
#Yes quite

str(test$leadingEdge) #Observing output column that contains all genes enriched in any particular pathway
#chr [1:440] "C3" "RGD1308923" "RGD1308923" "RGD1308923" "C3" "C3" "RGD621098" ...
###
#Because the leadingEdge column contains long strings of characters (gene names) it must be collapsed before being output

test$leadingEdge<-vapply(test$leadingEdge, paste, collapse = ",", character(1L)) #collapsing leadingEdge column so gene symbols are separated by comma and made into string smalle enough to fit in .csv file column
#output gsea data
write.csv(test, "testFGSEA_topmetagene.csv")

##################################################################################################
#    Validating fgsea method through circadian gene list (known enrichment in circadian pathways)
#
# Testing with circadian data
#
##################################################################################################

setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA") #setting directory to GSEA folder

#reading in .gmt file with human gene symbol info
human_Pathways<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/c2.all.v6.0.symbols.gmt")

#observing subset of human pathways
names(human_Pathways)[1:5]
#[1] "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"              
#[2] "KEGG_CITRATE_CYCLE_TCA_CYCLE"                 
#[3] "KEGG_PENTOSE_PHOSPHATE_PATHWAY"               
#[4] "KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS"
#[5] "KEGG_FRUCTOSE_AND_MANNOSE_METABOLISM" 

#Reading in ranked list

circadian_rnkinput<-read.csv("Circadian_FishersStatCrossRegion.csv", header=F)

rankedCircadianGeneList<-circadian_rnkinput[[2]] #subsetting numeric values from .rnk file
names(rankedCircadianGeneList)<-(circadian_rnkinput[[1]]) #assigning gene names to numeric values

#running fgsea on circadian genes using generic settings
circadianTest<-fgsea(human_Pathways, rankedCircadianGeneList, nperm=1000, minSize = 1, maxSize = 500)

#subsetting fgsea output with pvalues less than 0.05
sigCircadian<-subset(circadianTest, circadianTest$pval < 0.05)

circadianTest$leadingEdge<-vapply(circadianTest$leadingEdge, paste, collapse = ",", character(1L)) #collapsing for output

#output circadian fgsea test as .csv
write.csv(circadianTest, "testFGSEA_CircadianGenes.csv")



### Re-running top meta adult genes with different .gmt file

## Reading in .gmt file 

setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA")

rat_Pathways<-gmtPathways("file:///C:/Users/Izzy/Documents/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA/Rattus_norvegicus_GSEA_GO_sets_bp_symbols_highquality_April_2015.gmt")
names(rat_Pathways)
rat_Pathways[[1]]
#[1] "AADAT" "DLD"   "DLST"  "GOT1"  "GOT2"  "IDH1"  "IDH2"  "IDH3A" "IDH3B" "IDH3G"
#[11] "OGDH"           
#They're all caps -_-

#forcing rankedGeneList to be uppercase

rankedGeneList_upper<-adult_rnkinput[[2]] #remaking ranked gene list
names(rankedGeneList_upper)<-toupper(adult_rnkinput[[1]]) #naming with uppercase gene symbols to match rat pathways
names(rankedGeneList_upper)[1:5]
#[1] "ASB15"     "CAR9"      "SLC39A12"  "HMGN5B"    "LOC691921"
#all caps now ~('O')~

#Running fgsea with new pathway
gseaMetaGenes<-fgsea(rat_Pathways, rankedGeneList_upper, nperm=1000, minSize = 1, maxSize = 500)

gseaMetaGenes$leadingEdge<-vapply(gseaMetaGenes$leadingEdge, paste, collapse = ",", character(1L))
write.csv(gseaMetaGenes, "FGSEA_topmetagene.csv")

#subsetting fgsea output with pvals < 0.05
topMetaGSEA<-subset(gseaMetaGenes, temp$pval < 0.05)



###################################################
#Testing new .gmt file with circadian data again (since circadian genes are already uppercase)

#setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA")
#Reading in ranked list

#rankedCircadianGeneList<-circadian_rnkinput[[2]]
#names(rankedCircadianGeneList)<-(circadian_rnkinput[[1]])

#circadianTest2<-fgsea(rat_Pathways, rankedCircadianGeneList, nperm=1000, minSize = 1, maxSize = 500)

#sigCircadian2<-subset(circadianTest2, circadianTest2$pval < 0.05)
#View(sigCircadian2)
#About the same results as before so both .gmt files should be relatively similar

#circadianTest2$leadingEdge<-vapply(circadianTest2$leadingEdge, paste, collapse = ",", character(1L))
#write.csv(circadianTest2, "testFGSEA_CircadianGenes2.csv")

################################################################################
###Running a third time with 10,000 permutations instead of 1,000

circadianTest3<-fgsea(rat_Pathways, rankedCircadianGeneList, nperm=10000, minSize = 1, maxSize = 500)

#Subsetting enriched pathways with pvals < 0.05
sigCircadian3<-subset(circadianTest3, circadianTest3$pval < 0.05)
View(sigCircadian3)
#Now circadian rhythm and circadian regulation are both significant after BH
#Using 10,000 permutations results in much stronger p-values

#circadianTest3$leadingEdge<-vapply(circadianTest3$leadingEdge, paste, collapse = ",", character(1L))
#write.csv(circadianTest3, "testFGSEA_CircadianGenes3.csv")


#Testing plot function with Circadian genes

#Choose pathway of interest to plot against list of ranked genes
plotEnrichment(rat_Pathways[["circadian_rhythm(3)"]], rankedCircadianGeneList)
#Each pathway must be run separately
plotEnrichment(rat_Pathways[["ATP_biosynthetic_process(9)"]], rankedCircadianGeneList)
plotEnrichment(rat_Pathways[["gluconeogenesis(8)"]], rankedCircadianGeneList)

#Output for enrichment plot of circadian rhythm pathway for circadian genes
png("Circadian Rhythm Pathway Enrichment.png")
plotEnrichment(rat_Pathways[["circadian_rhythm(3)"]], rankedCircadianGeneList)
dev.off()

#cool


##Testing table function with circadian genes

#Subsetting the top 15 pathways with greatest significance from circadian gsea output
topPathways <- circadianTest3[head(order(pval), n=15)][order(NES), pathway]
## 

#Plotting a table of top 15 pathways for circadian genes
plotGseaTable(rat_Pathways[topPathways], rankedCircadianGeneList,
              circadianTest3, gseaParam=0.5)
#also cool

#Output .png of table containing info for top 15 pathways in circadian genes
png("CircadianGSEATable.png")
plotGseaTable(rat_Pathways[topPathways], rankedCircadianGeneList,
              circadianTest3, gseaParam=0.5)
dev.off()

#  End of testing with circadian genes
####################################################################
####################################################################

### Need to rerun analyses using greater permutations

###
### Re-running fgsea for top genes from adult meta-analysis with greater permutations
###
setwd("~/Phenotype Project/MetaAnalysis/MetaAnalysis 3.0 (newest)/GSEA")

rat_Pathways[[1]]
#[1] "AADAT" "DLD"   "DLST"  "GOT1"  "GOT2"  "IDH1"  "IDH2"  "IDH3A" "IDH3B" "IDH3G"
#[11] "OGDH"           


names(rankedGeneList_upper)[1:5]
#[1] "ASB15"     "CAR9"      "SLC39A12"  "HMGN5B"    "LOC691921"


#Running fgsea with 10,000 permutations

gseaMetaGenes2<-fgsea(rat_Pathways, rankedGeneList_upper, nperm=10000, minSize = 1, maxSize = 500)


gseaMetaGenes2$leadingEdge<-vapply(gseaMetaGenes2$leadingEdge, paste, collapse = ",", character(1L))
write.csv(gseaMetaGenes2, "FGSEA_topmetagene_10000perm.csv")


topMetaGSEA2<-subset(gseaMetaGenes2, gseaMetaGenes2$pval < 0.05)


#############################################################
### Re-running fgsea with 100,000 permutations

gseaMetaGenes3<-fgsea(rat_Pathways, rankedGeneList_upper, nperm=100000, minSize = 1, maxSize = 500)


gseaMetaGenes3$leadingEdge<-vapply(gseaMetaGenes3$leadingEdge, paste, collapse = ",", character(1L))
write.csv(gseaMetaGenes3, "FGSEA_topmetagene_100,000perm.csv")


#################### Visualizing Adult gsea output

##Creating table of top pathways for Adult Meta genes

topAdultMetaPathways <- gseaMetaGenes3[head(order(pval), n=15)][order(NES), pathway]

#Output table

png("MetaAdult_TopPathways.png", width=900)
plotGseaTable(rat_Pathways[topAdultMetaPathways], rankedGeneList_upper,
              gseaMetaGenes3, gseaParam=0.5)
dev.off()


#Looking at highly enriched gene sets
(topAdultMetaPathways[1:5])
#[1] "negative_regulation_of_programmed_cell_death(5)"
#[2] "negative_regulation_of_apoptotic_process(6)"    
#[3] "positive_regulation_of_cell_proliferation(4)"   
#[4] "tube_development(4)"                            
#[5] "tissue_morphogenesis(4)"


#plotting a few of the top Adult Meta pathways

png(paste(topAdultMetaPathways[1], "Pathway Enrichment for Adult.png"))
plotEnrichment(rat_Pathways[["negative_regulation_of_programmed_cell_death(5)"]], rankedGeneList_upper)
dev.off()

png(paste(topAdultMetaPathways[2], "Pathway Enrichment for Adult.png"))
plotEnrichment(rat_Pathways[[paste(topAdultMetaPathways[2])]], rankedGeneList_upper)
dev.off()




#########################################################
###########           P14 data           ################

### Now reading in P14 Meta ranked genes

metaP14_RankedGenes<-read.csv("MetaP14_RankedList.csv", header=F)

P14_rankedGeneList<-metaP14_RankedGenes[[2]]
names(P14_rankedGeneList)<-toupper(metaP14_RankedGenes[[1]]) #making them upper case

names(P14_rankedGeneList)[1:5]
#[1] "RGD1310311" "LOC501089"  "ABP1"       "OLR257"     "OLR1566" 

rat_Pathways[[1]]
#[1] "AADAT" "DLD"   "DLST"  "GOT1"  "GOT2"  "IDH1"  "IDH2"  "IDH3A" "IDH3B" "IDH3G"
#[11] "OGDH"           


###### Running fgsea on P14 genes with 100,000 perms

gsea_P14MetaGenes<-fgsea(rat_Pathways, P14_rankedGeneList, nperm=100000, minSize = 1, maxSize = 500)


str(gsea_P14MetaGenes$leadingEdge)
gsea_P14MetaGenes$leadingEdge<-vapply(gsea_P14MetaGenes$leadingEdge, paste, collapse = ",", character(1L))
write.csv(gsea_P14MetaGenes, "FGSEA_topP14metagene_100,000perms.csv")


sigMetaP14GSEA<-subset(gsea_P14MetaGenes, gsea_P14MetaGenes$pval < 0.001)

# 21 pathways with pvals < 0.001


#Visualizing gsea output


##Creating table of top pathways for P14

topP14Pathways <- gsea_P14MetaGenes[head(order(pval), n=15)][order(NES), pathway]
## 

png("MetaP14_TopPathways.png", width=1000)
plotGseaTable(rat_Pathways[topP14Pathways], P14_rankedGeneList,
              gsea_P14MetaGenes, gseaParam=0.5)
dev.off()

(topP14Pathways[1])
#[1] "positive_regulation_of_stem_cell_proliferation(5)"



#plotting top P14 pathway

png("Regulation of Stem Cell Proliferation Pathway Enrichment for P14.png")
plotEnrichment(rat_Pathways[["positive_regulation_of_stem_cell_proliferation(5)"]], P14_rankedGeneList)
dev.off()


