# miRNA target Gene Application tools
# Publication : https://www.frontiersin.org/articles/10.3389/fgene.2019.01330/full
# Package Instruction for miRNA conversion name : 
# Package Instruction for target genes: https://www.bioconductor.org/packages/release/bioc/vignettes/miRNAtap/inst/doc/miRNAtap.pdf
# Targets are aggregated from 5 most commonly cited prediction algorithms:
    #1. DIANA (Maragkakis et al., 2011), 
    #2. Miranda (Enright et al., 2003), 
    #3. PicTar (Lall et al., 2006), 
    #4. TargetScan (Friedman et al., 2009), and 
    #5. miRDB (Wong and Wang,2015)
# Package Instruction ggplot : https://r-graph-gallery.com/web-horizontal-barplot-with-labels-the-economist.html


#TIPS : Make sure install & package correctly

#BiocManager::install("miRNAtap")
#BiocManager::install("miRNAtap.db")
#BiocManager::install("topGO")            #ERROR: lazy loading failed for package 'rstatix' SOLUTION : update vctrs to 0.5.2
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("miRBaseConverter")

library(miRNAtap)
library(topGO)        
library(org.Hs.eg.db)
library(miRBaseConverter)
library("AnnotationDbi") #gene Entrez ID to Gene symbol

library(ggplot2)
library(dplyr)
library(writexl)

# The colors for GGPLOT2
BLUE <- "#076fa2"
RED <- "#E3120B"
BLACK <- "#202020"
GREY <- "grey50"

data(miRNATest) #raw datasests
Accessions = miRNATest$Accession #accession data
input=c("MIMAT0004494",
        "MIMAT0000450",
        "MIMAT0002171",
        "MIMAT0000280",
        "MIMAT0019880",
        "MIMAT0002807",
        "MIMAT0000318",
        "MIMAT0003218",
        "MIMAT0000254",
        "MIMAT0000689",
        "MIMAT0000073",
        "MIMAT0000099",
        "MIMAT0015066",
        "MIMAT0014995")

#For converting function
  #result2=miRNA_PrecursorToMature(miRNANames)                        #pma-mir-100a  ---> pma-miR-100a-5p pma-miR-100a-3p
  #result1=miRNA_MatureToPrecursor(miRNANames)                        # hsa-miR-196a-5p ---> hsa-mir-196a-1
  #result1 = miRNA_NameToAccession(miRNANames,version = version)      #hsa-miR-520b-3p ---> MIMAT0002843
  #result2 = miRNA_AccessionToName(result1[,2],targetVersion = "v22") #MIMAT0002843 ---> hsa-miR-520b-3p

result2 = miRNA_AccessionToName(input,targetVersion = "v22") #converting to miRNA name

#find target genes for each miRNA
#QC Test 
  #mir = 'miR-10b'
  #predictions = getPredictedTargets(mir, species= 'hsa',method = 'geom', min_src = 2)

allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL, mapping="org.Hs.eg.db", ID = "entrez") #get DATABASE of the Gene

#TEST USING SINGLE MIRNA, see:  result2
# NOTE : IF re-run data, make sure library activate before repeat these Line
##-----------------------miR-21-3p--------------------------1---##
mir_1 = 'miR-21-3p' #remove "hsa"
predictions_1 = getPredictedTargets(mir_1, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_1 = predictions_1[,'rank_product']
selection = function(x) TRUE
allGO2genes = annFUN.org(whichOnto='BP', feasibleGenes = NULL, mapping="org.Hs.eg.db", ID = "entrez") #get DATABASE of the Gene
GOdata = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_1,
             annot = annFUN.GO2genes, GO2genes = allGO2genes,
             geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks = runTest(GOdata, algorithm = "classic", statistic = "ks")
results.ks
allRes = GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 20) #only 20
sum_mir1=allRes[,c('GO.ID','Term','KS')] 
sum_mir1$KS=as.numeric(sum_mir1$KS) #Convert data structure into numeric from char
sum_mir1= filter(sum_mir1, sum_mir1$KS < 0.05)
# barplot() function is used to
plt <- ggplot(sum_mir1) +
  geom_col(aes(sum_mir1$KS,sum_mir1$Term), fill = "#076fa2", width = 0.5) 
plt + ggtitle("miR-21-3p target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=12))
plt #to show the plot

# Save Prediction miR-21-3p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_1), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_1_save= data.frame(geneSymbols, predictions_1) #ALL PREDICTION
sheets1 <- list("sheet1Name" = predictions_1_save, "sheet2Name" = sum_mir1) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets1, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-21-3p.xlsx")



##-----------------------miR-149-5p--------------------------2---##
mir_2 = 'miR-149-5p' #remove "hsa"
predictions_2 = getPredictedTargets(mir_2, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_2 = predictions_2[,'rank_product']
GOdata_2 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_2,
             annot = annFUN.GO2genes, GO2genes = allGO2genes,
             geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_2 = runTest(GOdata_2, algorithm = "classic", statistic = "ks")
results.ks_2
allRes_2 = GenTable(GOdata_2, KS = results.ks_2, orderBy = "KS", topNodes = 60) # CHECK always how many significant TRY:20,50,60
sum_mir_2=allRes_2[,c('GO.ID','Term','KS')] 
sum_mir_2$KS=as.numeric(sum_mir_2$KS) #Convert data structure into numeric from char
sum_mir_2= filter(sum_mir_2, sum_mir_2$KS < 0.05)
# barplot() function is used to
plt2 <- ggplot(sum_mir_2) +
  geom_col(aes(sum_mir_2$KS,sum_mir_2$Term), fill = "#E3120B", width = 0.5) 
plt2 + ggtitle("miR-149-5p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=8))
plt2 #to show the plot

# Save Prediction miR-149-5p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_2), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_2_save= data.frame(geneSymbols, predictions_2) #ALL PREDICTION
sheets2 <- list("sheet1Name" = predictions_2_save, "sheet2Name" = sum_mir_2) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets2, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-149-5p.xlsx")



##-----------------------hsa-miR-410-3p---------------------------3--##
mir_3 = 'miR-410-3p' #remove "hsa"
predictions_3 = getPredictedTargets(mir_3, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_3 = predictions_3[,'rank_product']
GOdata_3 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_3,
               annot = annFUN.GO2genes, GO2genes = allGO2genes,
               geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_3 = runTest(GOdata_3, algorithm = "classic", statistic = "ks")
results.ks_3
allRes_3 = GenTable(GOdata_3, KS = results.ks_3, orderBy = "KS", topNodes = 10) # in This case not SIG, show 10 
sum_mir_3=allRes_3[,c('GO.ID','Term','KS')] 
sum_mir_3$KS=as.numeric(sum_mir_3$KS) #Convert data structure into numeric from char
#sum_mir_3= filter(sum_mir_3, sum_mir_3$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt3 <- ggplot(sum_mir_3) +
  geom_col(aes(sum_mir_3$KS,sum_mir_3$Term), fill = "#E3120B", width = 0.5) 
plt3 + ggtitle("miR-410-3p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=12))
plt3 #to show the plot

# Save Prediction miR-410-3p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_3), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_3_save= data.frame(geneSymbols, predictions_3) #ALL PREDICTION
sheets3 <- list("sheet1Name" = predictions_3_save, "sheet2Name" = sum_mir_3) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets3, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-410-3p.xlsx")




##-----------------------miR-223-3p---------------------------4--##
mir_4 = 'miR-223-3p' #remove "hsa"
predictions_4 = getPredictedTargets(mir_4, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_4 = predictions_4[,'rank_product']
GOdata_4 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_4,
               annot = annFUN.GO2genes, GO2genes = allGO2genes,
               geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_4 = runTest(GOdata_4, algorithm = "classic", statistic = "ks")
results.ks_4
allRes_4 = GenTable(GOdata_4, KS = results.ks_4, orderBy = "KS", topNodes = 50) # 50 is significant
sum_mir_4=allRes_4[,c('GO.ID','Term','KS')] 
sum_mir_4$KS=as.numeric(sum_mir_4$KS) #Convert data structure into numeric from char
#sum_mir_4= filter(sum_mir_4, sum_mir_4$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt4 <- ggplot(sum_mir_4) +
  geom_col(aes(sum_mir_4$KS,sum_mir_4$Term), fill = "#E3120B", width = 0.5) 
plt4 + ggtitle("miR-223-3p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=8))
plt4 #to show the plot AND Save it manually

# Save Prediction miR-223-3p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_4), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_4_save= data.frame(geneSymbols, predictions_4) #ALL PREDICTION
sheets4 <- list("sheet1Name" = predictions_4_save, "sheet2Name" = sum_mir_4) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets4, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-223-3p.xlsx")




##-----------------------miR-4746-5p---------------------------5--##
mir_5 = 'miR-4746-5p' #remove "hsa"
predictions_5 = getPredictedTargets(mir_5, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_5 = predictions_5[,'rank_product']
GOdata_5 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_5,
               annot = annFUN.GO2genes, GO2genes = allGO2genes,
               geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_5 = runTest(GOdata_5, algorithm = "classic", statistic = "ks")
results.ks_5
allRes_5 = GenTable(GOdata_5, KS = results.ks_5, orderBy = "KS", topNodes = 10) # 50 is significant
sum_mir_5=allRes_5[,c('GO.ID','Term','KS')] 
sum_mir_5$KS=as.numeric(sum_mir_5$KS) #Convert data structure into numeric from char
#sum_mir_5= filter(sum_mir_5, sum_mir_5$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt5 <- ggplot(sum_mir_5) +
  geom_col(aes(sum_mir_5$KS,sum_mir_5$Term), fill = "#202020", width = 0.5) 
plt5 + ggtitle("miR-4746-5p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=8))
plt5 #to show the plot AND Save it manually

# Save Prediction miR-4746-5p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_5), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_5_save= data.frame(geneSymbols, predictions_5) #ALL PREDICTION
sheets5 <- list("sheet1Name" = predictions_5_save, "sheet2Name" = sum_mir_5) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets5, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-4746-5p.xlsx")



##-----------------------miR-491-5p---------------------------6--##
mir_6 = 'miR-491-5p' #remove "hsa"
predictions_6 = getPredictedTargets(mir_6, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_6 = predictions_6[,'rank_product']
GOdata_6 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_6,
               annot = annFUN.GO2genes, GO2genes = allGO2genes,
               geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_6 = runTest(GOdata_6, algorithm = "classic", statistic = "ks")
results.ks_6
allRes_6 = GenTable(GOdata_6, KS = results.ks_6, orderBy = "KS", topNodes = 65) # 50 is significant
sum_mir_6=allRes_6[,c('GO.ID','Term','KS')] 
sum_mir_6$KS=as.numeric(sum_mir_6$KS) #Convert data structure into numeric from char
#sum_mir_6= filter(sum_mir_6, sum_mir_6$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt6 <- ggplot(sum_mir_6) +
  geom_col(aes(sum_mir_6$KS,sum_mir_6$Term), fill = "grey50", width = 0.5) 
plt6 + ggtitle("miR-491-5p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=7))
plt6 #to show the plot AND Save it manually

# Save Prediction miR-491-5p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_6), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_6_save= data.frame(geneSymbols, predictions_6) #ALL PREDICTION
sheets6 <- list("sheet1Name" = predictions_6_save, "sheet2Name" = sum_mir_6) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets6, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-491-5p.xlsx")



##-----------------------miR-200b-3p---------------------------7--##
mir_7 = 'miR-200b-3p' #remove "hsa"
predictions_7 = getPredictedTargets(mir_7, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_7 = predictions_7[,'rank_product']
GOdata_7 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_7,
               annot = annFUN.GO2genes, GO2genes = allGO2genes,
               geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_7 = runTest(GOdata_7, algorithm = "classic", statistic = "ks")
results.ks_7
allRes_7 = GenTable(GOdata_7, KS = results.ks_7, orderBy = "KS", topNodes = 70) # 50 is significant
sum_mir_7=allRes_7[,c('GO.ID','Term','KS')] 
sum_mir_7$KS=as.numeric(sum_mir_7$KS) #Convert data structure into numeric from char
sum_mir_7= filter(sum_mir_7, sum_mir_7$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt7 <- ggplot(sum_mir_7) +
  geom_col(aes(sum_mir_7$KS,sum_mir_7$Term), fill = "grey50", width = 0.5) 
plt7 + ggtitle("miR-200b-3p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=7))
plt7 #to show the plot AND Save it manually

# Save Prediction miR-200b-3p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_7), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_7_save= data.frame(geneSymbols, predictions_7) #ALL PREDICTION
sheets7 <- list("sheet1Name" = predictions_7_save, "sheet2Name" = sum_mir_7) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets7, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-200b-3p.xlsx")





##-----------------------miR-92b-3p---------------------------8--##
mir_8 = 'miR-92b-3p' #remove "hsa"
predictions_8 = getPredictedTargets(mir_8, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_8 = predictions_8[,'rank_product']
GOdata_8 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_8,
               annot = annFUN.GO2genes, GO2genes = allGO2genes,
               geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_8 = runTest(GOdata_8, algorithm = "classic", statistic = "ks")
results.ks_8
allRes_8 = GenTable(GOdata_8, KS = results.ks_8, orderBy = "KS", topNodes = 90) # 86 is significant to many to show
sum_mir_8=allRes_8[,c('GO.ID','Term','KS')] 
sum_mir_8$KS=as.numeric(sum_mir_8$KS) #Convert data structure into numeric from char
sum_mir_8= filter(sum_mir_8, sum_mir_8$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
sum_mir_8_1=sum_mir_8[1:40,]
sum_mir_8_2=sum_mir_8[41:86,]

plt8.1 <- ggplot(sum_mir_8_1) +
  geom_col(aes(sum_mir_8_1$KS,sum_mir_8_1$Term), fill = "#E3120B", width = 0.5) 
plt8.1 + ggtitle("miR-92b-3p Target Gene Part-1") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=7))

plt8.2 <- ggplot(sum_mir_8_2) +
  geom_col(aes(sum_mir_8_2$KS,sum_mir_8_2$Term), fill = "#E3120B", width = 0.5) 
plt8.2 + ggtitle("miR-92b-3p Target Gene Part-2") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=7))

plt8 #to show the plot AND Save it manually

# Save Prediction miR-92b-3p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_8), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_8_save= data.frame(geneSymbols, predictions_8) #ALL PREDICTION
sheets8 <- list("sheet1Name" = predictions_8_save, "sheet2Name" = sum_mir_8) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets8, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-92b-3p.xlsx")



##-----------------------miR-10b-5p---------------------------9--##
mir_9 = 'miR-10b-5p' #remove "hsa"
predictions_9 = getPredictedTargets(mir_9, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_9 = predictions_9[,'rank_product']
GOdata_9 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_9,
               annot = annFUN.GO2genes, GO2genes = allGO2genes,
               geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_9 = runTest(GOdata_9, algorithm = "classic", statistic = "ks")
results.ks_9
allRes_9 = GenTable(GOdata_9, KS = results.ks_9, orderBy = "KS", topNodes = 70) # 50 is significant
sum_mir_9=allRes_9[,c('GO.ID','Term','KS')] 
sum_mir_9$KS=as.numeric(sum_mir_9$KS) #Convert data structure into numeric from char
sum_mir_9= filter(sum_mir_9, sum_mir_9$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt9 <- ggplot(sum_mir_9) +
  geom_col(aes(sum_mir_9$KS,sum_mir_9$Term), fill = "#202020", width = 0.5) 
plt9 + ggtitle("miR-10b-5p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=10))
plt9 #to show the plot AND Save it manually

# Save Prediction miR-10b-5p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_9), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_9_save= data.frame(geneSymbols, predictions_9) #ALL PREDICTION
sheets9 <- list("sheet1Name" = predictions_9_save, "sheet2Name" = sum_mir_9) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets9, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-10b-5p.xlsx")



##-----------------------miR-99b-5p---------------------------10--##
mir_10 = 'miR-99b-5p' #remove "hsa"
predictions_10 = getPredictedTargets(mir_10, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_10 = predictions_10[,'rank_product']
GOdata_10 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_10,
                annot = annFUN.GO2genes, GO2genes = allGO2genes,
                geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_10 = runTest(GOdata_10, algorithm = "classic", statistic = "ks")
results.ks_10
allRes_10 = GenTable(GOdata_10, KS = results.ks_10, orderBy = "KS", topNodes = 5) # only one sig
sum_mir_10=allRes_10[,c('GO.ID','Term','KS')] 
sum_mir_10$KS=as.numeric(sum_mir_10$KS) #Convert data structure into numeric from char
sum_mir_10= filter(sum_mir_10, sum_mir_10$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt10 <- ggplot(sum_mir_10) +
  geom_col(aes(sum_mir_10$KS,sum_mir_10$Term), fill = "grey50", width = 0.5) 
plt10 + ggtitle("miR-99b-5p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=15))
plt10 #to show the plot AND Save it manually

# Save Prediction miR-99b-5p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_10), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_10_save= data.frame(geneSymbols, predictions_10) #ALL PREDICTION
sheets10 <- list("sheet1Name" = predictions_10_save, "sheet2Name" = sum_mir_10) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets10, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-99b-5p.xlsx")




##-----------------------miR-19a-3p---------------------------11--##
mir_11 = 'miR-19a-3p' #remove "hsa"
predictions_11 = getPredictedTargets(mir_11, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_11 = predictions_11[,'rank_product']
GOdata_11 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_11,
                annot = annFUN.GO2genes, GO2genes = allGO2genes,
                geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_11 = runTest(GOdata_11, algorithm = "classic", statistic = "ks")
results.ks_11
allRes_11 = GenTable(GOdata_11, KS = results.ks_11, orderBy = "KS", topNodes = 175) # only one sig
sum_mir_11=allRes_11[,c('GO.ID','Term','KS')] 
sum_mir_11$KS=as.numeric(sum_mir_11$KS) #Convert data structure into numeric from char
sum_mir_11= filter(sum_mir_11, sum_mir_11$KS < 0.05) # SIGnificant is 111 GO
# barplot() function is used to
sum_mir_11_1=sum_mir_11[1:40,]
sum_mir_11_2=sum_mir_11[41:80,]
sum_mir_11_3=sum_mir_11[81:111,]

plt11.1 <- ggplot(sum_mir_11_1) +
  geom_col(aes(sum_mir_11_1$KS,sum_mir_11_1$Term), fill = "#E3120B", width = 0.5) 
plt11.1 + ggtitle("miR-19a-3p Target Gene part 1") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=10))

plt11.2 <- ggplot(sum_mir_11_2) +
  geom_col(aes(sum_mir_11_2$KS,sum_mir_11_2$Term), fill = "#E3120B", width = 0.5) 
plt11.2 + ggtitle("miR-19a-3p Target Gene part 2") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=10))

plt11.3 <- ggplot(sum_mir_11_3) +
  geom_col(aes(sum_mir_11_3$KS,sum_mir_11_3$Term), fill = "#E3120B", width = 0.5) 
plt11.3 + ggtitle("miR-19a-3p Target Gene part 3") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=10))

plt11 #to show the plot AND Save it manually

# Save Prediction miR-19a-3p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_11), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_11_save= data.frame(geneSymbols, predictions_11) #ALL PREDICTION
sheets11 <- list("sheet1Name" = predictions_11_save, "sheet2Name" = sum_mir_11) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets11, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-19a-3p.xlsx")



##-----------------------miR-101-3p---------------------------12--##
mir_12 = 'miR-101-3p' #remove "hsa"
predictions_12 = getPredictedTargets(mir_12, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_12 = predictions_12[,'rank_product']
GOdata_12 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_12,
                annot = annFUN.GO2genes, GO2genes = allGO2genes,
                geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_12 = runTest(GOdata_12, algorithm = "classic", statistic = "ks")
results.ks_12
allRes_12 = GenTable(GOdata_12, KS = results.ks_12, orderBy = "KS", topNodes = 150) # 142 is significant
sum_mir_12=allRes_12[,c('GO.ID','Term','KS')] 
sum_mir_12$KS=as.numeric(sum_mir_12$KS) #Convert data structure into numeric from char
sum_mir_12= filter(sum_mir_12, sum_mir_12$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
sum_mir_12_1=sum_mir_12[1:40,]
sum_mir_12_2=sum_mir_12[41:80,]
sum_mir_12_3=sum_mir_12[81:120,]
sum_mir_12_4=sum_mir_12[121:142,]

plt12.1 <- ggplot(sum_mir_12_1) +
  geom_col(aes(sum_mir_12_1$KS,sum_mir_12_1$Term), fill = "#202020", width = 0.5) 
plt12.1 + ggtitle("miR-101-3p Target Gene part 1") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=10))

plt12.2 <- ggplot(sum_mir_12_2) +
  geom_col(aes(sum_mir_12_2$KS,sum_mir_12_2$Term), fill = "#202020", width = 0.5) 
plt12.2 + ggtitle("miR-101-3p Target Gene part 2") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=8))

plt12.3 <- ggplot(sum_mir_12_3) +
  geom_col(aes(sum_mir_12_3$KS,sum_mir_12_3$Term), fill = "#202020", width = 0.5) 
plt12.3 + ggtitle("miR-101-3p Target Gene part 3") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=8))

plt12.4 <- ggplot(sum_mir_12_4) +
  geom_col(aes(sum_mir_12_4$KS,sum_mir_12_4$Term), fill = "#202020", width = 0.5) 
plt12.4 + ggtitle("miR-101-3p Target Gene part 4") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=10))

plt12 #to show the plot AND Save it manually

# Save Prediction miR-101-3p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_12), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_12_save= data.frame(geneSymbols, predictions_12) #ALL PREDICTION
sheets12 <- list("sheet1Name" = predictions_12_save, "sheet2Name" = sum_mir_12) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets12, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-101-3p.xlsx")



##-----------------------miR-3065-5p---------------------------13--##
mir_13 = 'miR-3065-5p' #remove "hsa"
predictions_13 = getPredictedTargets(mir_13, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_13 = predictions_13[,'rank_product']
GOdata_13 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_13,
                annot = annFUN.GO2genes, GO2genes = allGO2genes,
                geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_13 = runTest(GOdata_13, algorithm = "classic", statistic = "ks")
results.ks_13
allRes_13 = GenTable(GOdata_13, KS = results.ks_13, orderBy = "KS", topNodes = 55) # 
sum_mir_13=allRes_13[,c('GO.ID','Term','KS')] 
sum_mir_13$KS=as.numeric(sum_mir_13$KS) #Convert data structure into numeric from char
sum_mir_13= filter(sum_mir_13, sum_mir_13$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt13 <- ggplot(sum_mir_13) +
  geom_col(aes(sum_mir_13$KS,sum_mir_13$Term), fill = "grey50", width = 0.5) 
plt13 + ggtitle("miR-3065-5p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=10))
plt13 #to show the plot AND Save it manually

# Save Prediction miR-3065-5p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_13), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_13_save= data.frame(geneSymbols, predictions_13) #ALL PREDICTION
sheets13 <- list("sheet1Name" = predictions_13_save, "sheet2Name" = sum_mir_13) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets13, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-3065-5p.xlsx")




##-----------------------miR-3130-5p---------------------------14--##
mir_14 = 'miR-3130-5p' #remove "hsa"
predictions_14 = getPredictedTargets(mir_14, species = 'hsa', method = 'geom', min_src = 2)
rankedGenes_14 = predictions_14[,'rank_product']
GOdata_14 = new('topGOdata', ontology = 'BP', allGenes = rankedGenes_14,
                annot = annFUN.GO2genes, GO2genes = allGO2genes,
                geneSel = selection, nodeSize=10)
# Kolomonogorov Smirnov (K-S) test were used to analyzed rank information
results.ks_14 = runTest(GOdata_14, algorithm = "classic", statistic = "ks")
results.ks_14
allRes_14 = GenTable(GOdata_14, KS = results.ks_14, orderBy = "KS", topNodes = 55) # 
sum_mir_14=allRes_14[,c('GO.ID','Term','KS')] 
sum_mir_14$KS=as.numeric(sum_mir_14$KS) #Convert data structure into numeric from char
sum_mir_14= filter(sum_mir_14, sum_mir_14$KS < 0.05) # NOT SIGNIFICANT
# barplot() function is used to
plt14 <- ggplot(sum_mir_14) +
  geom_col(aes(sum_mir_14$KS,sum_mir_14$Term), fill = "#076fa2", width = 0.5) 
plt14 + ggtitle("miR-3130-5p Target Gene") +    #\n is enter
  xlab("KS P-value") + ylab("Gene Ontology") +theme(axis.text.y = element_text(size=12))
plt14 #to show the plot AND Save it manually

# Save Prediction miR-3130-5p
geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(predictions_14), column="SYMBOL", keytype="ENTREZID", multiVals="first") # Get the Gene Symbol
predictions_14_save= data.frame(geneSymbols, predictions_14) #ALL PREDICTION
sheets14 <- list("sheet1Name" = predictions_14_save, "sheet2Name" = sum_mir_14) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets14, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-3130-5p.xlsx")






# NOTE to replace
#       1. Change the miRNA according to arrangement in result2 [miRNA name, plot, save predict, file name ]
#       2. use Notepad++ replace the number ex: _3 to _4
#       3. editing in plt4, sheets4, 
#       4. check the allRes_(n) to see significant miRNA exists (try big and then small number)

# The colors for GGPLOT2
BLUE <- "#076fa2"
RED <- "#E3120B"
BLACK <- "#202020"
GREY <- "grey50"


# QC-1:  need to check accuracy of each code
#Warning messages:
#  1: Use of `sum_mir_11_1$KS` is discouraged.
#ℹ Use `KS` instead. 
#2: Use of `sum_mir_11_1$Term` is discouraged.
#ℹ Use `Term` instead. 




#check pathway for miRNA targets genes
# source material : https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
rm()

