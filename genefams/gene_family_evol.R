## Gene family analysis
## Import results from Count under Wagner Parsimony, gain penalty = 1.2

library(dplyr)
library(topGO)

gene_families <- read.table("Orthogroups.txt", heade =F, sep = ':')
gf_count <- read.delim("count_genefamilies_wagner.txt")

gf_count_losses <- gf_count %>% 
  filter(Losses > 0) %>% 
  filter(NDUV.v.appalachiana < 1) %>% 
  filter(NDUV.v.appalachiana < X3)

# get contigs corresponding to lost OGs
OGs_lost_contigs <- gene_families %>% 
  filter(V1 %in% gf_count_losses$Family)
write.table(OGs_lost_contigs, "OGs_lost.txt")

gf_count_con <- gf_count %>% #checked against GUI results 
  filter(Contractions > 0) %>% 
  dplyr::select(Family, NDUV.v.appalachiana, X3) %>% 
  filter(NDUV.v.appalachiana <= (0.5*X3)) %>%
  filter(NDUV.v.appalachiana == 1)
  
# get contigs corresponding to contracted OGs
OGs_contracted_contigs <- gene_families %>% 
  filter(V1 %in% gf_count_con$Family)
write.table(OGs_contracted_contigs, "OGs_contracted.txt")

gf_count_gain <- gf_count %>% 
  filter(Gains > 0) %>%
  filter(X3 == 0) %>% 
  filter(NDUV.v.appalachiana > 0)

# get contigs corresponding to gained OGs
OGs_gained_contigs <- gene_families %>% 
  filter(V1 %in% gf_count_gain$Family)
write.table(OGs_gained_contigs, "OGs_gained.txt")

gf_count_exp <- gf_count %>% #Checked against GUI results 
  filter(Expansions > 0) %>%
  filter(X3 == 1) %>% 
  filter(NDUV.v.appalachiana >= 1.5*X3)

# get contigs corresponding to expanded OGs
OGs_expanded_contigs <- gene_families %>% 
  filter(V1 %in% gf_count_exp$Family)
write.table(OGs_expanded_contigs, "OGs_expanded.txt")

#in bash: grep "scaffold\-NDUV\-[0-9A-Za-z\_\.\-]*" OGs_XX.txt -o > OGs_XX_NDUV_contigs.txt
OGs_gained_NDUV_contigs <- read.delim("OGs_gained_NDUV_contigs.txt", header = F)
#OGs_gained_GOs <- inner_join(all_NDUV_transcriptome, OGs_gained_NDUV_contigs)

geneID2GO <- readMappings(file = "NDUV_all_contigs_terms_mapping.txt")

# define gene universe 
all_genes <- names(geneID2GO)
genes_gained <- as.character(OGs_gained_NDUV_contigs$V1)
# define genes of interest as a subset of the total dataset 
geneList <- factor(as.integer(all_genes %in% genes_gained))
names(geneList) <- all_genes

## BP GAINED OGs

gained_GO_data_BP <- new("topGOdata", description = "gained OGs", 
                      ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_BP <- sigGenes(gained_GO_data_BP)
str(sig_genes_BP)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(gained_GO_data_BP)
#1025

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
BP_result_Fisher <- runTest(gained_GO_data_BP, algorithm = "classic", statistic = "fisher")
BP_result_Fisher
#504 terms have raw p<0.01

BP_list <- usedGO(object = gained_GO_data_BP)

BP_fisher_only <- GenTable(gained_GO_data_BP, classicFisher = BP_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(BP_list))

BP_corrected <- p.adjust(BP_fisher_only$classicFisher, method = "fdr")
# NA induced by coercion, all NAs should be p.adjust(p=1e-30, n = 7933, method= "fdr") = 7.933e-27
BP_corrected <- as.data.frame(BP_corrected)
#BP_corrected_significant <- BP_corrected %>% 
#  filter(BP_corrected < 0.01)
# 226 GO terms significant after FDR correction (p < .01)
BP_fisher_top <- BP_fisher_only[1:226,]
BP_corrected$classicFisher = BP_fisher_only$classicFisher
BP_table <- dplyr::left_join(BP_fisher_top, BP_corrected)
BP_table_uni <- unique(BP_table)
write.csv(BP_table_uni, file = "BP_corrected_gained_OGs.csv") 

showSigOfNodes(gained_GO_data_BP, score(BP_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(gained_GO_data_BP, BP_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO_gained_OGs_BP", useInfo = 'all')

## CC gained OGs

gained_GO_data_CC <- new("topGOdata", description = "gained OGs", 
                      ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_CC <- sigGenes(gained_GO_data_CC)
str(sig_genes_CC)
#Number of significantly enriched genes 
numSigGenes(gained_GO_data_CC)
#1205

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
CC_result_Fisher <- runTest(gained_GO_data_CC, algorithm = "classic", statistic = "fisher", )
CC_result_Fisher
#91 terms have raw p<0.01

CC_list <- usedGO(object = gained_GO_data_CC)

CC_fisher_only <- GenTable(gained_GO_data_CC, classicFisher = CC_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(CC_list))

CC_corrected <- p.adjust(CC_fisher_only$classicFisher, method = "fdr")
CC_corrected <- as.data.frame(CC_corrected)
#67 terms p<0.01 significantly enriched following FDR adjustment
# NA induced by coercion, all NAs should be p.adjust(p=1e-30, n = 1234, method= "fdr") = 1.234e-27

CC_fisher_top <- CC_fisher_only[1:67,]
CC_corrected$classicFisher = CC_fisher_only$classicFisher
CC_table <- dplyr::left_join(CC_fisher_top, CC_corrected)
CC_table_uni <- unique(CC_table)
write.csv(CC_table_uni, file = "CC_corrected_gained_OGs.csv") 

showSigOfNodes(gained_GO_data_CC, score(CC_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(gained_GO_data_CC, CC_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO_gained_CC_OGs", useInfo = 'all')

## MF gained OGs 

gained_GO_data_MF <- new("topGOdata", description = "gained_OGs", 
                      ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_MF <- sigGenes(gained_GO_data_MF)
str(sig_genes_MF)
#Number of significantly enriched genes
numSigGenes(gained_GO_data_MF)
#1057

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
MF_result_Fisher <- runTest(gained_GO_data_MF, algorithm = "classic", statistic = "fisher")
MF_result_Fisher
#90 terms have raw p<0.01

MF_list <- usedGO(object = gained_GO_data_MF)

MF_fisher_only <- GenTable(gained_GO_data_MF, classicFisher = MF_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(MF_list))

MF_corrected <- p.adjust(MF_fisher_only$classicFisher, method = "fdr")
MF_corrected <- as.data.frame(MF_corrected)
# 39 terms p < 0.01, significantly enriched following FDR adjustment
# NAs introduced by coercion, NAs should be p.adjust(p=1e-30, n=2818, method = "fdr") = 2.818e-27

MF_fisher_top <- MF_fisher_only[1:39,]
MF_corrected$classicFisher = MF_fisher_only$classicFisher
MF_table <- dplyr::left_join(MF_fisher_top, MF_corrected)
MF_table_uni <- unique(MF_table)
write.csv(MF_table_uni, file = "MF_corrected_gained_OGs.csv") 

showSigOfNodes(gained_GO_data_MF, score(MF_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(gained_GO_data_MF, MF_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO_gained_OGs_MF", useInfo = 'all')

################ GO ENRICHMENT FOR EXPANDED GENE FAMILIES #######################################
## Change genes of interest to expanded gene families 
OGs_expanded_NDUV_contigs <- read.delim("OGs_expanded_NDUV_contigs.txt", header = F)
OGs_expanded <- as.character(OGs_expanded_NDUV_contigs$V1)
# define genes of interest as a subset of the total dataset 
geneList <- factor(as.integer(all_genes %in% OGs_expanded))
names(geneList) <- all_genes

## BP EXPANDED OGs

expanded_GO_data_BP <- new("topGOdata", description = "expnded OGs", 
                         ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_BP <- sigGenes(expanded_GO_data_BP)
str(sig_genes_BP)
#Number of significantly enriched genes 
numSigGenes(gained_GO_data_BP)
#1045

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
BP_result_Fisher <- runTest(expanded_GO_data_BP, algorithm = "classic", statistic = "fisher")
BP_result_Fisher
#182 terms have raw p<0.01

BP_list <- usedGO(object = expanded_GO_data_BP)

BP_fisher_only <- GenTable(expanded_GO_data_BP, classicFisher = BP_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(BP_list))

BP_corrected <- p.adjust(BP_fisher_only$classicFisher, method = "fdr")
BP_corrected <- as.data.frame(BP_corrected)
#BP_corrected_significant <- BP_corrected %>% 
#  filter(BP_corrected < 0.01)
# 8 GO terms significant after FDR correction (p < .01)
BP_fisher_top <- BP_fisher_only[1:8,]
BP_corrected$classicFisher = BP_fisher_only$classicFisher
BP_table <- dplyr::left_join(BP_fisher_top, BP_corrected)
BP_table_uni <- unique(BP_table)
write.csv(BP_table_uni, file = "BP_corrected_expanded_OGs.csv") 

showSigOfNodes(expanded_GO_data_BP, score(BP_result_Fisher), firstSigNodes = 8, useInfo = 'all')
printGraph(expanded_GO_data_BP, BP_result_Fisher, firstSigNodes = 8, fn.prefix = "tGO_expanded_OGs_BP", useInfo = 'all')

## CC expanded OGs

expanded_GO_data_CC <- new("topGOdata", description = "exapnded OGs", 
                         ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_CC <- sigGenes(expanded_GO_data_CC)
str(sig_genes_CC)
#Number of significantly enriched genes 
numSigGenes(expanded_GO_data_CC)
#1170

#Run Fisher's Exact Test
#Note that this does not correct for multiple comparisons 
CC_result_Fisher <- runTest(expanded_GO_data_CC, algorithm = "classic", statistic = "fisher", )
CC_result_Fisher
#35 terms have raw p<0.01

CC_list <- usedGO(object = expanded_GO_data_CC)

CC_fisher_only <- GenTable(expanded_GO_data_CC, classicFisher = CC_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(CC_list))
#write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

CC_corrected <- p.adjust(CC_fisher_only$classicFisher, method = "fdr")
CC_corrected <- as.data.frame(CC_corrected)
#0 terms p<0.01 significantly enriched following FDR adjustment

## MF expanded OGs 

expanded_GO_data_MF <- new("topGOdata", description = "expanded_OGs", 
                         ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_MF <- sigGenes(expanded_GO_data_MF)
str(sig_genes_MF)
#Number of significantly enriched genes 
numSigGenes(expanded_GO_data_MF)
#1042

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
MF_result_Fisher <- runTest(expanded_GO_data_MF, algorithm = "classic", statistic = "fisher")
MF_result_Fisher
#97 terms have raw p<0.01

MF_list <- usedGO(object = expanded_GO_data_MF)

MF_fisher_only <- GenTable(expanded_GO_data_MF, classicFisher = MF_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(MF_list))

MF_corrected <- p.adjust(MF_fisher_only$classicFisher, method = "fdr")
MF_corrected <- as.data.frame(MF_corrected)
# 3 terms p < 0.01, significantly enriched following FDR adjustment

MF_fisher_top <- MF_fisher_only[1:3,]
MF_corrected$classicFisher = MF_fisher_only$classicFisher
MF_table <- dplyr::left_join(MF_fisher_top, MF_corrected)
MF_table_uni <- unique(MF_table)
write.csv(MF_table_uni, file = "MF_corrected_expanded_OGs.csv") 

showSigOfNodes(expanded_GO_data_MF, score(MF_result_Fisher), firstSigNodes = 3, useInfo = 'all')
printGraph(expanded_GO_data_MF, MF_result_Fisher, firstSigNodes = 3, fn.prefix = "tGO_expanded_OGs_MF", useInfo = 'all')

################ GO ENRICHMENT FOR CONTRACTED GENE FAMILIES #######################################
## Change genes of interest to contracted gene families 
OGs_contracted_NDUV_contigs <- read.delim("OGs_contracted_NDUV_contigs.txt", header = F)
OGs_contracted <- as.character(OGs_contracted_NDUV_contigs$V1)
# define genes of interest as a subset of the total dataset 
geneList <- factor(as.integer(all_genes %in% OGs_contracted))
names(geneList) <- all_genes

## BP contraced OGs

contracted_GO_data_BP <- new("topGOdata", description = "contracted OGs", 
                           ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_BP <- sigGenes(contracted_GO_data_BP)
str(sig_genes_BP)
#Number of significantly enriched genes 
numSigGenes(contracted_GO_data_BP)
#381

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
BP_result_Fisher <- runTest(contracted_GO_data_BP, algorithm = "classic", statistic = "fisher")
BP_result_Fisher
#130 terms have raw p<0.01

BP_list <- usedGO(object = contracted_GO_data_BP)

BP_fisher_only <- GenTable(contracted_GO_data_BP, classicFisher = BP_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(BP_list))

BP_corrected <- p.adjust(BP_fisher_only$classicFisher, method = "fdr")
BP_corrected <- as.data.frame(BP_corrected)
#BP_corrected_significant <- BP_corrected %>% 
#  filter(BP_corrected < 0.01)
# 5 GO terms significant after FDR correction (p < .01)
BP_fisher_top <- BP_fisher_only[1:5,]
BP_corrected$classicFisher = BP_fisher_only$classicFisher
BP_table <- dplyr::left_join(BP_fisher_top, BP_corrected)
BP_table_uni <- unique(BP_table)
write.csv(BP_table_uni, file = "BP_corrected_contracted_OGs.csv") 

showSigOfNodes(expanded_GO_data_BP, score(BP_result_Fisher), firstSigNodes = 5, useInfo = 'all')
printGraph(expanded_GO_data_BP, BP_result_Fisher, firstSigNodes = 5, fn.prefix = "tGO_contracted_OGs_BP", useInfo = 'all')

## CC contracted OGs

contracted_GO_data_CC <- new("topGOdata", description = "contracted OGs", 
                           ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_CC <- sigGenes(contracted_GO_data_CC)
str(sig_genes_CC)
#Number of significantly enriched genes 
numSigGenes(contracted_GO_data_CC)
#422

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
CC_result_Fisher <- runTest(contracted_GO_data_CC, algorithm = "classic", statistic = "fisher", )
CC_result_Fisher
#22 terms have raw p<0.01

CC_list <- usedGO(object = contracted_GO_data_CC)

CC_fisher_only <- GenTable(contracted_GO_data_CC, classicFisher = CC_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(CC_list))
#write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

CC_corrected <- p.adjust(CC_fisher_only$classicFisher, method = "fdr")
CC_corrected <- as.data.frame(CC_corrected)
#0 terms p<0.01 significantly enriched following FDR adjustment

## MF contracted OGs 

contracted_GO_data_MF <- new("topGOdata", description = "contracted OGs", 
                           ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#Review significantly enriched genes  
sig_genes_MF <- sigGenes(contracted_GO_data_MF)
str(sig_genes_MF)
#Number of significantly enriched genes 
numSigGenes(contracted_GO_data_MF)
#378

#Run Fisher's Exact Test
#Note that this does not correct for multiple comparisons 
MF_result_Fisher <- runTest(contracted_GO_data_MF, algorithm = "classic", statistic = "fisher")
MF_result_Fisher
#41 terms have raw p<0.01

MF_list <- usedGO(object = contracted_GO_data_MF)

MF_fisher_only <- GenTable(contracted_GO_data_MF, classicFisher = MF_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(MF_list))

MF_corrected <- p.adjust(MF_fisher_only$classicFisher, method = "fdr")
MF_corrected <- as.data.frame(MF_corrected)
# 0 terms p < 0.01, significantly enriched following FDR adjustment

################ GO ENRICHMENT FOR LOST GENE FAMILIES #######################################
## NOTE: CHANGE BACKGROUND TO INCLUDE LOST GENE FAMILIES (ANNOTATED IN V. LINEATA)
lost_SKYV <- read.table("OGs_lost_SKYV_contigs.txt", header = F)
all_SKYV <- read.delim("SKYV_all_contigs_terms", header = F)
lost_SKYV_contig_GO <- all_SKYV %>% 
  filter(V1 %in% lost_SKYV$V1)
write.table(lost_SKYV_contig_GO, "SKYV_lost_contig_go.txt", sep = '\t')
lost_contig_mapping <- readMappings("lost_contigs_mapping.txt")

# define gene universe 
all_genes <- names(lost_contig_mapping)
genes_lost <- as.character(lost_SKYV$V1)
# define genes of interest as a subset of the total dataset 
geneList <- factor(as.integer(all_genes %in% genes_lost))
names(geneList) <- all_genes

## BP lost OGs

lost_GO_data_BP <- new("topGOdata", description = "Lost OGs", 
                         ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = lost_contig_mapping)

#Review significantly enriched genes  
sig_genes_BP <- sigGenes(lost_GO_data_BP)
str(sig_genes_BP)
#Number of significantly enriched genes 
numSigGenes(lost_GO_data_BP)
#540

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
BP_result_Fisher <- runTest(lost_GO_data_BP, algorithm = "classic", statistic = "fisher")
BP_result_Fisher
#104 terms have raw p<0.01

BP_list <- usedGO(object = lost_GO_data_BP)

BP_fisher_only <- GenTable(lost_GO_data_BP, classicFisher = BP_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(BP_list))

BP_corrected <- p.adjust(BP_fisher_only$classicFisher, method = "fdr")
BP_corrected <- as.data.frame(BP_corrected)
# 7 GO terms significant after FDR correction (p < .01)
BP_fisher_top <- BP_fisher_only[1:7,]
BP_corrected$classicFisher = BP_fisher_only$classicFisher
BP_table <- dplyr::left_join(BP_fisher_top, BP_corrected)
BP_table_uni <- unique(BP_table)
write.csv(BP_table_uni, file = "BP_corrected_lost_OGs.csv") 

showSigOfNodes(lost_GO_data_BP, score(BP_result_Fisher), firstSigNodes = 7, useInfo = 'all')
printGraph(lost_GO_data_BP, BP_result_Fisher, firstSigNodes = 7, fn.prefix = "tGO_lost_OGs_BP", useInfo = 'all')

## CC lost OGs

lost_GO_data_CC <- new("topGOdata", description = "lost OGs", 
                         ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = lost_contig_mapping)

#Review significantly enriched genes  
sig_genes_CC <- sigGenes(lost_GO_data_CC)
str(sig_genes_CC)
#Number of significantly enriched genes 
numSigGenes(lost_GO_data_CC)
#627

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
CC_result_Fisher <- runTest(lost_GO_data_CC, algorithm = "classic", statistic = "fisher")
CC_result_Fisher
#19 terms have raw p<0.01

CC_list <- usedGO(object = lost_GO_data_CC)

CC_fisher_only <- GenTable(lost_GO_data_CC, classicFisher = CC_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(CC_list))

CC_corrected <- p.adjust(CC_fisher_only$classicFisher, method = "fdr")
CC_corrected <- as.data.frame(CC_corrected)
#13 terms p<0.01 significantly enriched following FDR adjustment

CC_fisher_top <- CC_fisher_only[1:13,]
CC_corrected$classicFisher = CC_fisher_only$classicFisher
CC_table <- dplyr::left_join(CC_fisher_top, CC_corrected)
CC_table_uni <- unique(CC_table)
write.csv(CC_table_uni, file = "CC_corrected_lost_OGs.csv") 

showSigOfNodes(lost_GO_data_CC, score(CC_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(lost_GO_data_CC, CC_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO_lost_OGs_CC", useInfo = 'all')

## MF gained OGs 

lost_GO_data_MF <- new("topGOdata", description = "lost OGs", 
                         ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = lost_contig_mapping)

#Review significantly enriched genes  
sig_genes_MF <- sigGenes(lost_GO_data_MF)
str(sig_genes_MF)
#Number of significantly enriched genes
numSigGenes(lost_GO_data_MF)
#516

#Run Fisher's Exact Test 
#Note that this does not correct for multiple comparisons 
MF_result_Fisher <- runTest(lost_GO_data_MF, algorithm = "classic", statistic = "fisher")
MF_result_Fisher
#68 terms have raw p<0.01

MF_list <- usedGO(object = lost_GO_data_MF)

MF_fisher_only <- GenTable(lost_GO_data_MF, classicFisher = MF_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(MF_list))

MF_corrected <- p.adjust(MF_fisher_only$classicFisher, method = "fdr")
MF_corrected <- as.data.frame(MF_corrected)
# 6 terms p < 0.01, significantly enriched following FDR adjustment

MF_fisher_top <- MF_fisher_only[1:6,]
MF_corrected$classicFisher = MF_fisher_only$classicFisher
MF_table <- dplyr::left_join(MF_fisher_top, MF_corrected)
MF_table_uni <- unique(MF_table)
write.csv(MF_table_uni, file = "MF_corrected_lost_OGs.csv") 

showSigOfNodes(lost_GO_data_MF, score(MF_result_Fisher), firstSigNodes = 6, useInfo = 'all')
printGraph(lost_GO_data_MF, MF_result_Fisher, firstSigNodes = 6, fn.prefix = "tGO_lost_OGs_MF", useInfo = 'all')
