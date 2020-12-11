## GO Enrichment Tests for dN/dS
library(topGO)
library(dplyr)

# pairwise dN/dS enrichment 

NDUV <- read.delim("NDUV_all_contigs_terms.txt", header = F)
NDUV <- NDUV %>% 
  dplyr::select(V1, V2, V4, V6, V8, V10, V12, V14, V16) # get columns with GO terms only
  
SKYV <- read.delim("SKYV_all_contigs_terms", header = F)
SKYV <- SKYV %>% 
  dplyr::select(V1, V2, V4, V6, V8, V10, V12, V14, V16, V18, V20, V22, V24, V26, V28, V30, V32, V34, V36, V38) # get columns with GO terms only 

NDUV$V18 <- NA; NDUV$V20 <- NA; NDUV$V22 <- NA; NDUV$V24 <- NA; NDUV$V26 <- NA;NDUV$V28 <- NA;
NDUV$V30 <- NA; NDUV$V32 <- NA; NDUV$V34 <- NA; NDUV$V36 <- NA; NDUV$V38 <- NA
NDUV_SKYV <- rbind(SKYV,NDUV)

purifying_dnds <- read.delim("all_dn_ds_pur_contigs.txt", header = F)

purifying_GO <- NDUV_SKYV %>% 
  filter(V1 %in% purifying_dnds$V1)

write.table(x= purifying_GO, "purifying_to_edit.txt")

pos_outliers <- read.delim("postive_outliers_contigs.txt", header = F)

pos_outliers_GO <- NDUV_SKYV %>% 
  filter(V1 %in% pos_outliers$V1)

write.table(x= pos_outliers_GO, "positive_outliers_to_edit.txt")

neg_outliers <- read.delim("negative_outliers_contigs.txt", header = F)

neg_outliers_GO <- NDUV_SKYV %>% 
  filter(V1 %in% neg_outliers$V1)

write.table(x = neg_outliers_GO, "negative_outliers_to_edit.txt")
  
pos_GO <- read.delim("pos_outliers_GO_terms.txt", header = F)
neg_GO <- read.delim("negative_outliers_GO_terms.txt", header = F)
geneID2GO <- readMappings(file = "purifying_to_edit.txt")

str(geneID2GO)

# define gene universe 
all_genes <- names(geneID2GO)

# read in GO annotations for genes of interest 
genes_of_interest <- as.character(pos_outliers$V1)
# define genes of interest as a subset of the total dataset 
geneList <- factor(as.integer(all_genes %in% genes_of_interest))
names(geneList) <- all_genes

#############################################################################################

## BP POSITIVE dN/dS DIFFERENCE 

pos_GO_data_BP <- new("topGOdata", description = "dN/dS positive outliers", 
                      ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_BP <- sigGenes(pos_GO_data_BP)
str(sig_genes_BP)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(pos_GO_data_BP)
#16

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
BP_result_Fisher <- runTest(pos_GO_data_BP, algorithm = "classic", statistic = "fisher")
BP_result_Fisher
#16 terms have raw p<0.01

BP_list <- usedGO(object = pos_GO_data_BP)

BP_fisher_only <- GenTable(pos_GO_data_BP, classicFisher = BP_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(BP_list))

BP_corrected <- p.adjust(BP_fisher_only$classicFisher, method = "fdr")
BP_corrected <- as.data.frame(BP_corrected)
#no terms are significantly enriched following FDR adjustment

showSigOfNodes(pos_GO_data_BP, score(BP_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(pos_GO_data_BP, BP_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')

## CC POSITIVE dN/dS DIFFERENCE

pos_GO_data_CC <- new("topGOdata", description = "dN/dS positive outliers", 
                      ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_CC <- sigGenes(pos_GO_data_CC)
str(sig_genes_CC)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(pos_GO_data_CC)
#19

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
CC_result_Fisher <- runTest(pos_GO_data_CC, algorithm = "classic", statistic = "fisher", )
CC_result_Fisher
#3 terms have raw p<0.01

CC_list <- usedGO(object = pos_GO_data_CC)

CC_fisher_only <- GenTable(pos_GO_data_CC, classicFisher = CC_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(CC_list))
#write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

CC_corrected <- p.adjust(CC_fisher_only$classicFisher, method = "fdr")
CC_corrected <- as.data.frame(CC_corrected)
#no terms are significantly enriched following FDR adjustment

showSigOfNodes(pos_GO_data_CC, score(CC_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(pos_GO_data_CC, CC_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')

## MF POSITIVE dN/dS DIFFERENCE

pos_GO_data_MF <- new("topGOdata", description = "dN/dS positive outliers", 
                      ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_MF <- sigGenes(pos_GO_data_MF)
str(sig_genes_MF)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(pos_GO_data_MF)
#17

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
MF_result_Fisher <- runTest(pos_GO_data_MF, algorithm = "classic", statistic = "fisher")
MF_result_Fisher
#1 terms have raw p<0.01

MF_list <- usedGO(object = pos_GO_data_MF)

MF_fisher_only <- GenTable(pos_GO_data_MF, classicFisher = MF_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(MF_list))
#write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

MF_corrected <- p.adjust(MF_fisher_only$classicFisher, method = "fdr")
MF_corrected <- as.data.frame(MF_corrected)
#no terms are significantly enriched following FDR adjustment

showSigOfNodes(pos_GO_data_MF, score(MF_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(pos_GO_data_MF, MF_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')

#############################################################################################################

## GENES WITH NEG dN/dS DIFF

# define gene universe 
all_genes <- names(geneID2GO)

# read in GO annotations for genes of interest 
genes_of_interest <- as.character(neg_outliers$V1)
# define genes of interest as a subset of the total dataset 
geneList <- factor(as.integer(all_genes %in% genes_of_interest))
names(geneList) <- all_genes

## BP POSITIVE dN/dS DIFFERENCE 

neg_GO_data_BP <- new("topGOdata", description = "dN/dS positive outliers", 
                      ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_BP <- sigGenes(neg_GO_data_BP)
str(sig_genes_BP)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(neg_GO_data_BP)
#10

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
BP_result_Fisher <- runTest(neg_GO_data_BP, algorithm = "classic", statistic = "fisher", )
BP_result_Fisher
#13 terms have raw p<0.01

BP_list <- usedGO(object = neg_GO_data_BP)

BP_fisher_only <- GenTable(neg_GO_data_BP, classicFisher = BP_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(BP_list))
#write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

BP_corrected <- p.adjust(BP_fisher_only$classicFisher, method = "fdr")
BP_corrected <- as.data.frame(BP_corrected)
#no terms are significantly enriched following FDR adjustment

showSigOfNodes(neg_GO_data_BP, score(BP_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(neg_GO_data_BP, BP_result_Fisher, firstSigNodes = 5, fn.prefix = "tGO", useInfo = 'all')

## CC NEGATIVE dN/dS DIFFERENCE

neg_GO_data_CC <- new("topGOdata", description = "dN/dS positive outliers", 
                      ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_CC <- sigGenes(neg_GO_data_CC)
str(sig_genes_CC)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(neg_GO_data_CC)
#8

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
CC_result_Fisher <- runTest(neg_GO_data_CC, algorithm = "classic", statistic = "fisher", )
CC_result_Fisher
#1 terms have raw p<0.01

CC_list <- usedGO(object = neg_GO_data_CC)

CC_fisher_only <- GenTable(neg_GO_data_CC, classicFisher = CC_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(CC_list))
#write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

CC_corrected <- p.adjust(CC_fisher_only$classicFisher, method = "fdr")
CC_corrected <- as.data.frame(CC_corrected)
#no terms are significantly enriched following FDR adjustment

showSigOfNodes(neg_GO_data_CC, score(CC_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(neg_GO_data_CC, CC_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')


## MF POSITIVE dN/dS DIFFERENCE

neg_GO_data_MF <- new("topGOdata", description = "dN/dS positive outliers", 
                      ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_MF <- sigGenes(pos_GO_data_MF)
str(sig_genes_MF)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(pos_GO_data_MF)
#17

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
MF_result_Fisher <- runTest(neg_GO_data_MF, algorithm = "classic", statistic = "fisher", )
MF_result_Fisher
#3 terms have raw p<0.01

MF_list <- usedGO(object = neg_GO_data_MF)

MF_fisher_only <- GenTable(neg_GO_data_MF, classicFisher = MF_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(MF_list))
#write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

MF_corrected <- p.adjust(MF_fisher_only$classicFisher, method = "fdr")
MF_corrected <- as.data.frame(MF_corrected)
#no terms are significantly enriched following FDR adjustment

showSigOfNodes(pos_GO_data_MF, score(MF_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(pos_GO_data_MF, MF_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')


