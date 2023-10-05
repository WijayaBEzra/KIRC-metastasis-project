### Kidney Cancer Project ####
# WE INVESTIGATE Kidney cancer USING BIOINFORMATICS MACHINE LEARNING PIPELINE
# REVISION : 0
# YEAR : 2022 - 2023
# note for multi mir (updated: 13-01.2023: multiMiR :https://www.bioconductor.org/packages/devel/bioc/vignettes/multiMiR/inst/doc/multiMiR.html)
# Re-Fix RPM data when downloaded as excel dosent represent the million value; Note for RPM & RNA-seq analysis : https://holab-hku.github.io/R-workshop/rna-seq-analysis-in-r.html
# Re-Fis RPM data analysis limma : https://rpubs.com/jrgonzalezISGlobal/transcriptomic_analyses
# Issue-1 : XENA provide log unit for RPM in this testing we converting back to RPM


library(mlbench)
library(caret)
library(dplyr)
library(stringr)

library(tibble)
#library(UBL) #for data imbalance
#library(themis) #for data imbalance
library(edgeR)

library( "DESeq2" )
library(ggplot2)

###DATA EXTRACTION###
# Only choose read_count criteria
myFile = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your file manually
mirna_name=data.frame(myFile[,1])

##### OLD DATA #######################################
#tmyFile=t(myFile) #transpose matrix
#selectedRows <- tmyFile[grep("read_count", tmyFile), ] #select only Row naming read_count REMOVE miRNA name
#myFile=t(selectedRows)
#File = data.frame(mirna_name,myFile) #combine mirna name and expression data
#colnames(File)[1] <- "miRNA" #Set first col name Header
#File <- File[-c(1), ]  #Remove first row but Not Header
#######################################################

####### Check Quality of Missing Data miRNA #######
count_na <- data.frame(rowSums(is.na(myFile)))
count_na <- data.frame(mirna_name,count_na)

QC_count_na <- (rowSums(is.na(myFile))/length (myFile))
QC_count_na_ <- data.frame(mirna_name,QC_count_na) #Check Percentage of missing value per-Row miRNA

#Show barplot of 100 row data
barplot(QC_count_na_[1:100,2],
main = "Missing data in 100 miRNA",
xlab = "Missing Data in %",
ylab = "miRNA",
names.arg = QC_count_na_[1:100,1],
col = "red",
horiz = TRUE)
#Show hetmap plot of missing value of 100 miRNA & 100 patients
myFile_h=myFile[,-1]
myFile_h=data.matrix(myFile_h)
myFile_h[is.na(myFile_h)] = 0# convert NA into zero
pheatmap(myFile_h[1:200,1:200])

#Check missing value higher to certain %
sum_NA=sum(QC_count_na_$QC_count_na > 0.5) #counting missing value >0.5
sum_NA_7=sum(QC_count_na_$QC_count_na > 0.7) #counting missing value >0.7
sum_NA_9=sum(QC_count_na_$QC_count_na > 0.9) #counting missing value >0.9
data_mir=QC_count_na_[QC_count_na_$QC_count_na < 0.5,] # Show you miRNA with less than 50% missing value

#Choose only miRNA data with <50% missing Value
myFile_=myFile[row.names(data_mir),] #since Row names in data_mir is fixed (use fixed index instead variable name)
miRNA_name=data.frame(myFile_[,1]) #mirna name considerable in <50% missing value

##### Divide Data into normal - Primary - Metastasis #####
#Seperate Normal-Primary Cancer
Normal=myFile_[,grep(".11", colnames(myFile_))] #miRNA name is lost in this data
Primary=myFile_[,grep(".01", colnames(myFile_))] #need to grep patient of metastasis name from clinical

#Get Clinical data metastasis
myFile2 = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your clinical file manually
new_row= data.frame(myFile2$sampleID,myFile2$pathologic_M) #only choose col SampleID & TNM (M) stage

cc=dplyr::filter(new_row, grepl('M1', myFile2$pathologic_M)) #Grab row only M1
cc_=str_replace_all(cc$myFile2.sampleID,"-",".") # edited TCGA name in from clinical data which use '-' in the patient's name
#Note : patient from Primay need to subtract

#aa=substring(a, 9) #Extract only patient name
da=(paste(cc_, collapse="|")) #To make it as "OR"
#da=toupper(da) #Prepare input for m1 data in miRNA file

Metastasis=Primary %>% dplyr::select(matches(da)) #calling metastasis data from Primary Data

#Update Primary M0 only
cp=dplyr::filter(new_row, grepl('M0', myFile2$pathologic_M)) #Grab patients row only M0
cp_=str_replace_all(cp$myFile2.sampleID,"-",".")
dap=(paste(cp_, collapse="|")) #To make it as "OR"
Primary_=Primary %>% dplyr::select(matches(dap)) #take Patients name in clinical with M0 cathegory

#Converting log2(RPM+1)into original value RPM = 2^(value)-1
Normal=(2^(Normal))-1
Primary_=(2^(Primary_))-1
Metastasis=(2^(Metastasis))-1
#QC negative values
Normal[Normal<0] = 0
Primary_[Primary_<0]=0
Metastasis[Metastasis<0]=0

##### Imputation #####
library(imputeR)
#Imputation info from https://cran.r-project.org/web/packages/imputeR/imputeR.pdf
# Note : some publication DOI:10.1186/s12859-015-0853-0
#                         https://doi.org/10.1186/s12859-022-04656-4
#                         https://doi.org/10.1038/s41598-021-03438-x
#impdata_N <- impute(Normal, lmFun = "lassoR") #using Lasso so in Normal can be inputed (Filled the data)
#myFile_Normal <- data.frame(impdata_N$imp) #Data Frame Filled with lasso imputation
#Simulation in R : Impute_KRIC.R
#Note-1 : Experiment consist of 8 Algorithms and 10 condition, repeat 5 time 
#Note-2 : output will be comparison between RMSE and computation time
#Note-3 : Result shows plsR is good enough

# After We got the best Performance we apply Such Imputation
impdata_P <- impute(Primary_, lmFun = "plsR") #using method based on Optimize performance
impdata_M <- impute(Metastasis,lmFun= "plsR") #using method based on Optimize performance
impdata_N <- impute(Normal, lmFun = "plsR")
myFile_Normal <- data.frame(impdata_N$imp)
myFile_Primary <- data.frame(impdata_P$imp)
myFile_Metastasis <- data.frame(impdata_M$imp)

#QC
#a=data.frame(myFile_Primary,myFile_Normal)
count_na <- data.frame(rowSums(is.na(myFile_Normal))) #evaluation missing data NA
count_min <- data.frame(rowSums(myFile_Normal<0)) #evaluation missing data NA
#Note: minus value in inputation due to huge missing values


#### DEG USING LIMMA #####
#set table for Limma package
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://rpubs.com/jrgonzalezISGlobal/transcriptomic_analyses
# note 1: design matrix
# ncol(x) to calculate number of table
library(limma)
library(EnhancedVolcano)

design_n=data.frame(colnames(myFile_Normal),code=0) ; colnames(design_n)=c("patient_ID","code")
design_p=data.frame(colnames(myFile_Primary),code=1) ; colnames(design_p)=c("patient_ID","code")
design_m=data.frame(colnames(myFile_Metastasis),code=2); colnames(design_m)=c("patient_ID","code")

# Design Matrix
design_pvn=rbind(design_p,design_n)
design_pvn=data.matrix(design_pvn, rownames.force = NA) #convert DF to numeric matrix
design_mvp=rbind(design_m, design_p)
design_mvp=data.matrix(design_mvp, rownames.force = NA) #convert DF to numeric matrix

#combined Data Frame P vs N #
df_primary_v_normal=data.frame(miRNA_name,myFile_Primary,myFile_Normal) #dont input miRNA yet because affecting dimension
rownames(df_primary_v_normal) <- df_primary_v_normal[,1] #col1 as rownames
df_primary_v_normal = df_primary_v_normal[,-1]

#QC : edit negative value to zero
df_primary_v_normal[df_primary_v_normal < 0] <- 0
df_primary_v_normal=data.matrix(df_primary_v_normal,rownames.force = NA)

#count_na <- data.frame(rowSums(is.na(df_primary_v_normal))) # evaluation code for missing data
#df_primary_v_normal[df_primary_v_normal < 0] <- 0.000001 #Removing negative value 
#see=df_primary_v_normal[df_primary_v_normal < 0] #show you negative value row

y <- voom(df_primary_v_normal, design_pvn, plot = T)
#Error 19.01.2023 : Error in lmFit(y, design, block = block, correlation = correlation, weights = weights) : row dimension of design doesn't match column dimension of data ob
#QC for dimension problem error
dim(df_primary_v_normal)
dim(design_pvn)
#Error:   Negative counts not allowed its assume data is read count with (0 -> inf+) not float; due to log2 scale and ZERo(0)not allowed (https://support.bioconductor.org/p/58673/)
#Option-1 : convert back to count/ or antilog2 for each value (from beginning) data=2^(score)
#Option-2 : process usual formula (re-construct your own)

#QC Variance trends not good, it may need filters on low RPM genes
# for filter code: https://rdrr.io/bioc/edgeR/man/filterByExpr.html
QC_expression = data.frame(rowSums2(df_primary_v_normal))
keep <- filterByExpr(df_primary_v_normal, design_pvn)
keep_ = data.frame(df_primary_v_normal[keep,]) # Found 470 miRNA with good expression
#Second Test for limma-voom
z <- voom(keep_, design_pvn, plot = T) #Result little better somehow

# LogFC analysis(version filtered low RPM genes)
fit <- lmFit(z, design_pvn)
fit <- eBayes(fit)
topTable(fit)

# logFC original
fit_ <- lmFit(y, design_pvn)
fit_ <- eBayes(fit_)
topTable(fit_)

top.table <- topTable(fit, sort.by = "F", n = Inf) #with Filter
top.table2 <- topTable(fit_, sort.by = "F", n = Inf) #without Filter
length(which(top.table$adj.P.Val < 0.05)) #counting DEGs

#QC Call row names
#top.table2[c("MIMAT0000073"),]
#top.table2[c("MIMAT0000099"),]
#top.table2[c("MIMAT0015066"),]
#top.table2[c("MIMAT0014995"),]


#plot Patient sample PvN
plotMDS(z, col = as.numeric(design_pvn))

#Check Gene Function
#BiocManager::install("MEAL")
#library(MEAL)
#need to update 'cli' to min version 3.4.1
#fit.meal <- getProbeResults(ans, coef=2, 
#                            fNames=c("probe", "Gene.symbol", "Chromosome.location"))
#head(fit.meal)

#Use volcano plot to represent
names(toptable) #logFC still using "code" as name

volcanoplot(fit, coef=2, names = fit$F, 
            highlight = 5)

toptable <- topTable(fit, n = Inf)
toptable=add_column(toptable, row.names(toptable), .after = 0) #adding miRNA name in toptable

# Use enchanced volcano plot to show
EnhancedVolcano(toptable,
                lab = toptable$`row.names(toptable)`, #show gene/mir name
                x = 'code',
                y = 'P.Value')

#more filter on volcanoplot
# link : https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
EnhancedVolcano(toptable,
                lab = row.names(toptable),
                x = 'code',
                y = 'P.Value',
                title = 'Primary vs Normal miRNA',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0, #dot font size
                labSize = 3.0) #mirna name font size

#fixed the table names
toptable_=toptable
colnames(toptable_)=c("miRNA_ID", "patient_ID",  "logFC",       "AveExpr",    "F" ,         "P.Value",    "adj.P.Val" )
Result_pvn=data.frame(toptable_) #Table results

#NEED TO FIX REPRESENTATION, need to check from design the toptable not fit (OVERALL XENA IS MESS)
#library(writexl)
#check_result=data.frame(fit)
#write_xlsx(a,"D:\\Dissertation 4.08.2022\\check Results.xlsx")


# NOTE: In Fit data we convert into table
# logFC == Coefficients.code
# AveExpr == Amean
# t == t.code
# P.Value == p.val.code
# B == lods.code
# only adjusted pval is not recorded

#combined Data Frame M vs P #
df_metastasis_v_primary=data.frame(miRNA_name,myFile_Metastasis,myFile_Primary) #dimensional problem is fixed using code below
rownames(df_metastasis_v_primary) <- df_metastasis_v_primary[,1] #col1 as rownames
df_metastasis_v_primary = df_metastasis_v_primary[,-1] #remove col 1

#QC : edit negative value to zero
df_metastasis_v_primary[df_metastasis_v_primary < 0] <- 0
df_metastasis_v_primary=data.matrix(df_metastasis_v_primary,rownames.force = NA)

x <- voom(df_metastasis_v_primary, design_mvp, plot = T)

#QC Variance trends not good, it may need filters on low RPM genes
# for filter code: https://rdrr.io/bioc/edgeR/man/filterByExpr.html 
# library(edgeR)
QC_expression_mvp = data.frame(rowSums2(df_metastasis_v_primary))
keep_qc_mvp <- filterByExpr(df_metastasis_v_primary, design_mvp)
keep_qc_mvp = data.frame(df_metastasis_v_primary[keep,]) # Found 470 miRNA with good expression

x1 <- voom(keep_qc_mvp, design_mvp, plot = T)

# LogFC analysis with Filter quality expression
fit_mvp <- lmFit(x1, design_mvp)
fit_mvp_ <- eBayes(fit_mvp)
topTable(fit_mvp_)

# LogFC analysis Raw [630 miRNA]
fit_mvp2 <- lmFit(x, design_mvp)
fit_mvp2_ <- eBayes(fit_mvp2)
topTable(fit_mvp2_)
toptable_mvp2 <- topTable(fit_mvp2_, sort.by = "F", n = Inf)

#check Result
toptable_mvp2[c("MIMAT0002171"),]
toptable_mvp2[c("MIMAT0019880"),]
toptable_mvp2[c("MIMAT0002807"),]

#QC for results
check_table_2=data.frame(fit_mvp_) #check the table

#TOP table metastasis vs primary
toptable_mvp <- topTable(fit_mvp_, sort.by = "F", n = Inf)
length(which(toptable_mvp$adj.P.Val < 0.05)) #counting DEGs

#plot Patient sample MvP
plotMDS(x1, col = as.numeric(design_mvp))

# volcano plot M v P
volcanoplot(fit_mvp_, coef=2, names = fit_mvp_$rank, 
            highlight = 5)

#toptable <- topTable(fit_mvp_, n = Inf)
toptable_mvp_1=add_column(toptable_mvp, row.names(toptable_mvp), .after = 0) #adding miRNA name in toptable

EnhancedVolcano(toptable_mvp_1,
                lab = toptable_mvp_1$`row.names(toptable_mvp)`, #show gene/mir name
                x = 'code',
                y = 'P.Value')

#more filter on volcanoplot
# link : https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
EnhancedVolcano(toptable_mvp_1,
                lab = row.names(toptable_mvp_1),
                x = 'code',
                y = 'P.Value',
                title = 'Metastasis vs Primary miRNA',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0, #dot font size
                labSize = 3.0) #mirna name font size

#fixed the table names
toptable_mvp_1x=toptable_mvp_1
colnames(toptable_mvp_1x)=c("miRNA_ID", "patient_ID", "logFC",       "AveExpr",    "F" ,         "P.Value",    "adj.P.Val" )
Result_mvp=data.frame(toptable_mvp_1x) #Table results

#Evaluation Note: Question why the results all upregulated??? does using TCGA or remove the imputation will help in this analysis?
# or because we removed missing value since beginning?


#find common miRNA
# Evaluation : all common miRNA is belong to metastasis with also belongs to 470/630 PvN mirna
common_mirna <- data.frame(intersect(row.names(Result_mvp), row.names(Result_pvn)))  #testing code :  padjustsig_np[grep("1250", rownames(padjustsig_np)), ]
setdiff_pn <- data.frame(setdiff(row.names(Result_pvn),common_mirna$intersect.row.names.Result_mvp...row.names.Result_pvn..)) #find set diff np vs common
setdiff_mp <- data.frame(setdiff(row.names(Result_mvp),common_mirna$intersect.row.names.Result_mvp...row.names.Result_pvn..))  #find set diff pm vs common                    

#$# Save
library(writexl)
write_xlsx(common_mirna,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\common_mirna.xlsx")
write_xlsx(Result_mvp,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Result_mvp.xlsx")
write_xlsx(Result_pvn,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Result_pvn.xlsx")

#Save Imputed data Normal-Primary-Metastasis
miRNA_Normal=data.frame(miRNA_name, myFile_Normal)
miRNA_Normal[miRNA_Normal < 0] <- 0 #remove negative value!

miRNA_Primary=data.frame(miRNA_name, myFile_Primary)
miRNA_Primary[miRNA_Primary < 0] <- 0 #remove negative value!

miRNA_Metastasis=data.frame(miRNA_name, myFile_Metastasis)
miRNA_Metastasis[miRNA_Metastasis< 0] <- 0 #remove negative value!

write_xlsx(miRNA_Normal,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Normal.xlsx")
write_xlsx(miRNA_Primary,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Primary.xlsx")
write_xlsx(miRNA_Metastasis,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Metastasis.xlsx")

#EVALUATION NOTES
#Evaluation Note: Question why the results all upregulated??? does using TCGA or remove the imputation will help in this analysis?
# or because we imput missing value AFTER seperated into normal-Primary-metastasis? so it lock to its features
# Evaluation : all common miRNA is belong to metastasis with also belongs to 470/630 PvN mirna

#Request 12/07/2023
# save full miRNA list without filter
pvn=data.frame(row.names(top.table2),top.table2)
mvp=data.frame(row.names(toptable_mvp2),toptable_mvp2)
write_xlsx(top.table2,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_DEG_PvN.xlsx")
write_xlsx(toptable_mvp2,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_DEG_MvP.xlsx")
