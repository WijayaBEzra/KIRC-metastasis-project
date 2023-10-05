### Kidney Cancer ####
# WE INVESTIGATE Kidney Cancer USING BIOINFORMATICS MACHINE LEARNING PIPELINE
# REVISION : 0
# YEAR : 2022

library(mlbench)
library(caret)
library(dplyr)
library(stringr)

library(tibble)
library(UBL) #for data imbalance
library(themis) #for data imbalance
library(edgeR)

library( "DESeq2" )
library(ggplot2)

###DATA EXTRACTION###
# Only choose read_count criteria
myFile = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your file manually
#mirna_name=data.frame(myFile[,1])
#tmyFile=t(myFile) #transpose matrix
#selectedRows <- tmyFile[grep("read_count", tmyFile), ] #select only Row naming read_count REMOVE miRNA name
#myFile=t(selectedRows)
#File = data.frame(mirna_name,myFile) #combine mirna name and expression data
#colnames(File)[1] <- "miRNA" #Set first col name Header
#File <- File[-c(1), ]  #Remove first row but Not Header

# only miRNA col
#mirna_name_=data.frame(File[,1])
Gene_name=data.frame(myFile[,1])
colnames(Gene_name)[1] <- "Gene" #Set first col name Header

#Seperate Normal-Primary Cancer
Normal=myFile[,grep(".11", colnames(myFile))] 
Primary=myFile[,grep(".01", colnames(myFile))] #need to grep patient of metastasis name from clinical

#Get Clinical data metastasis
myFile2 = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your clinical file manually
tmyFile2 = t(myFile2)
tmyFile2 = data.frame(tmyFile2) #converting into matrix
tmyFile2 <- tmyFile2[-c(1), ] #removing first row

new_row= data.frame(rownames(tmyFile2),tmyFile2[,10]) #only choose TNM (M) stage
#new_row= data.frame(rownames(tmyFile2),tmyFile2[,11]) #CESC-----only choose TNM (M) stage

cc=dplyr::filter(tmyFile2, grepl('m1', X10)) #Grab row berdasarkan fitur di col x10
#cc=dplyr::filter(tmyFile2, grepl('m1', X11)) #CESC---------Grab row berdasarkan fitur di col x10

tcc=t(cc) #transpose back
getcolname=data.frame(row.names(cc))
a=t(getcolname) #Show m1 Patients
aa=substring(a, 9) #Extract only patient name
da=(paste(aa, collapse="|")) #To make it as "OR"
da=toupper(da) #Prepare input for m1 data in miRNA file

Metastasis=myFile %>% dplyr::select(matches(da)) #calling metastasis data
#NOTE : error in get metastasis datasets for CESC #FINDING : location TNM in clinical data CESC is different
# WARNING : Patients name in PRIMARY still contain metastasis Patients, so we should remove

#set table for Limma package
# note 1: design matrix
# ncol(x) to calculate number of table
design_n=data.frame(colnames(Normal),code=0) ; colnames(design_n)=c("patient_ID","code")
design_p=data.frame(colnames(Primary),code=1) ; colnames(design_p)=c("patient_ID","code")
design_m=data.frame(colnames(Metastasis),code=2); colnames(design_m)=c("patient_ID","code")
design_nvp=rbind(design_n, design_p)
design_pvm=rbind(design_p, design_m)

#combined Data Frame N vs P
df_normalvprimary=data.frame(Gene_name,Normal,Primary)
#df_normalvprimary <- df_normalvprimary[,-c(1)] #remove miRNA col; usually in first col
#colnames(df_normalvprimary)[1] <- "Gene"
df_normalvprimary=data.matrix(df_normalvprimary)
#df_normalvprimary=as.matrix(df_normalvprimary) #note usage: as.matrix https://bioinformatics.stackexchange.com/questions/8808/why-does-deseq2-convert-numeric-columns-to-factor-during-differential-expression
df_normalvprimary[, 'Gene'] <- as.factor(df_normalvprimary[, 'Gene']) #for character error in desq2

#combined Data Frame P vs M
df_primaryvmetas=data.frame(Gene_name,Primary,Metastasis)
#colnames(df_primaryvmetas)[1] <- "miRNA"
df_primaryvmetas=data.matrix(df_primaryvmetas)
#df_primaryvmetas=as.matrix(df_primaryvmetas)
df_primaryvmetas[, 'Gene'] <- as.factor(df_primaryvmetas[, 'Gene']) #try to converting miRNA name into number (arranged)


#### construct DESQ datas set object Normal vs Primary ####
df_normal_vs_Primary <- DESeqDataSetFromMatrix(countData=round(df_normalvprimary), 
                              colData=design_nvp, 
                              design=~code, tidy = TRUE)
#NOTE : Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are not integers
# Try : https://support.bioconductor.org/p/133326/


#Summarized Results
np <- DESeq(df_normal_vs_Primary) #wait for 1 min
res_np <- results(np)
#head(results(res_np, tidy=TRUE))
summary(res_np)
res_np <- res_np[order(res_np$padj),] #order based on p-val
head(res_np)
result_np=data.frame(res_np)
result_np$index <- as.numeric(row.names(result_np)) #arranged row names
result_np=result_np[order(result_np$index),]
#result_np=data.frame(mirna_name_,result_np)
rownames(result_np)=Gene_name[,1] #set gene names
result_np=result_np[,-7]
result_np <- result_np[order(result_np$padj),]
# still bug because gene result not the same with internet : https://lashlock.github.io/compbio/R_presentation.html

## Volcano Plot normal vs primary tumor ##
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res_np, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_np, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_np, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#### construct DESQ datas set object Primary vs Metastasis ####
###############################################################
df_Primary_vs_metastasis <- DESeqDataSetFromMatrix(countData=round(df_primaryvmetas), 
                                               colData=design_pvm, 
                                               design=~code, tidy = TRUE)
#Summarized Results
pm <- DESeq(df_Primary_vs_metastasis) #wait for 1 min
res_pm <- results(pm)
#head(results(res_np, tidy=TRUE))
summary(res_pm)
res_pm <- res_pm[order(res_pm$padj),] #order based on p-val
head(res_pm)
result_pm=data.frame(res_pm)
result_pm$index <- as.numeric(row.names(result_pm)) #arranged row names
result_pm=result_pm[order(result_pm$index),]
#result_np=data.frame(mirna_name_,result_np)
rownames(result_pm)=Gene_name[,1] #set gene names
result_pm=result_pm[,-7]
result_pm <- result_pm[order(result_pm$padj),]
# still bug because gene result not the same with internet : https://lashlock.github.io/compbio/R_presentation.html

## Volcano Plot  primary tumor vs metastasis ##
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res_pm, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res_np, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res_np, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#Filtering Significant Result
padjustsig<- which(result_np$padj < 0.05) #filtering data which padj should less than 0.05
padjustsig_np<- data.frame (result_np[c(padjustsig),])
padjustsig_np1<- data.frame (rownames(padjustsig_np),padjustsig_np)

padjustsig_pm<- which(result_pm$padj < 0.05) #filtering data which padj should less than 0.05
padjustsig_pm<- data.frame (result_pm[c(padjustsig_pm),])
padjustsig_pm1<- data.frame (rownames(padjustsig_pm),padjustsig_pm)

#find common miRNA
common_mirna <- data.frame(intersect(row.names(padjustsig_pm), row.names(padjustsig_np)))  #testing code :  padjustsig_np[grep("1250", rownames(padjustsig_np)), ]
setdiff_np <- data.frame(setdiff(row.names(padjustsig_np),common_mirna[1,])) #find set diff np vs common
setdiff_pm <- data.frame(setdiff(row.names(padjustsig_pm),common_mirna[1,]))  #find set diff pm vs common                    

#$# Save
library(writexl)
write_xlsx(common_mirna,"D:\\Dissertation 4.08.2022\\Pan-Cancer miRNA\\Common_Gene_KIRC.xlsx")
write_xlsx(padjustsig_np1,"D:\\Dissertation 4.08.2022\\Pan-Cancer miRNA\\DEG_Gene_NvP_KIRC.xlsx")
write_xlsx(padjustsig_pm1,"D:\\Dissertation 4.08.2022\\Pan-Cancer miRNA\\DEG_Gene_PvM_KIRC.xlsx")

