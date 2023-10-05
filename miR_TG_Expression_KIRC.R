# miRNA-Gene Expression relation
# Purpose :
# 1. Checking correlation between miRNA and list of target genes from miRTAP analysis
# 2. Side-analysis will be checking how many of those Significant-Gene are member in Pathway analysis
# 3. remapped using Subcluster analysis
# 4. Final best miRNA which covering the function of Metastasis Pathway will be selected for Deep Learning

library(dplyr)
library(stringr)
library(tibble)
library("readxl")
options(scipen = 99) #avoid unit as exponential number as 1e+100

### --------------------------------------- DATA EXTRACTION Gene Expression -------------------------------------------###
# Only choose read_count criteria
myFile = read.delim(file.choose(), header = TRUE) #WARNING: check your tab, choose your file manually to Gene Expression
Gene_name=data.frame(myFile[,1])

### Divide Data into normal - Primary - Metastasis ###
#Seperate Normal-Primary Cancer
Normal_gene=myFile[,grep(".11", colnames(myFile))] #miRNA name is lost in this data
Primary_gene=myFile[,grep(".01", colnames(myFile))] #need to grep patient of metastasis name from clinical

#Get Clinical data metastasis
myFile2 = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your clinical file manually
new_row= data.frame(myFile2$sampleID,myFile2$pathologic_M) #only choose col SampleID & TNM (M) stage
cc=dplyr::filter(new_row, grepl('M1', myFile2$pathologic_M)) #Grab row only M1
cc_=str_replace_all(cc$myFile2.sampleID,"-",".") # edited TCGA name in from clinical data which use '-' in the patient's name
da=(paste(cc_, collapse="|")) #To make it as "OR"
Metastasis_gene=Primary_gene %>% dplyr::select(matches(da)) #calling metastasis data from Primary Data
row.names(Metastasis_gene)=Gene_name[,1]

#Update Primary M0 only
cp=dplyr::filter(new_row, grepl('M0', myFile2$pathologic_M)) #Grab patients row only M0
cp_=str_replace_all(cp$myFile2.sampleID,"-",".")
dap=(paste(cp_, collapse="|")) #To make it as "OR"
Primary_gene=Primary_gene %>% dplyr::select(matches(dap)) #take Patients name in clinical with M0 cathegory
row.names(Primary_gene)=Gene_name[,1]

### --------------------------------------- DATA EXTRACTION miRNA -------------------------------------------###
input_mirna_M=c("MIMAT0004494",
        "MIMAT0000450",
        "MIMAT0002171",
        "MIMAT0000280",
        "MIMAT0019880",
        "MIMAT0002807",
        "MIMAT0000318",
        "MIMAT0003218",
        "MIMAT0000254",
        "MIMAT0000689")
        
input_mirna_P=c("MIMAT0000073",
        "MIMAT0000099",
        "MIMAT0015066",
        "MIMAT0014995")

# Choose excel File imputed miRNA Primary & Metastasis
Metastasis_mirna=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Metastasis.xlsx")
Primary_mirna=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Primary.xlsx")

#since dplyr filter only works in col-wise
Metastasis_mirna1=t(Metastasis_mirna)
Primary_mirna1=t(Primary_mirna)
colnames(Metastasis_mirna1)=Metastasis_mirna1[1,] ; Metastasis_mirna1 <- Metastasis_mirna1[-1, ] # first row into col names
colnames(Primary_mirna1)=Primary_mirna1[1,] ; Primary_mirna1 <- Primary_mirna1[-1, ] # first row into col names
dmp=(paste(input_mirna_P, collapse="|")) #get miRNA Cox primary
dmm=(paste(input_mirna_M, collapse="|")) #get miRNA Cox metastasis
#rechanged into data.frame dplyr required
Metastasis_mirna1=data.frame(Metastasis_mirna1)
Primary_mirna1=data.frame(Primary_mirna1)

Primary_mirna1=Primary_mirna1 %>% dplyr::select(matches(dmp)) #take miRNA expression according to miRNA cox NOTE: produced CHR data type
Metastasis_mirna1=Metastasis_mirna1 %>% dplyr::select(matches(dmm)) #take miRNA expression according to miRNA cox NOTE: produced CHR data type

# The patients list in miRNA same with patients list in Gene Expression
pp=row.names(Primary_mirna1);pp=(paste(pp, collapse="|")) #Get patients name in Primary
Primary_patients=Primary_gene %>% dplyr::select(matches(pp)) #take

mp=row.names(Metastasis_mirna1);mp=(paste(mp, collapse="|"))
Metastasis_patients=Metastasis_gene %>% dplyr::select(matches(mp)) #take
#Note : no problems

### --------------------------------------- Match the Target Genes-------------------------------------------###
#see supplementary Target Genes (Metastasis)
MIMAT0000254_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-10b-5p.xlsx")
MIMAT0000689_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-99b-5p.xlsx")
MIMAT0004494_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-21-3p.xlsx")
MIMAT0000450_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-149-5p.xlsx")
MIMAT0002171_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-410-3p.xlsx")
MIMAT0000280_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-223-3p.xlsx")
MIMAT0019880_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-4746-5p.xlsx")
MIMAT0002807_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-491-5p.xlsx")
MIMAT0000318_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-200b-3p.xlsx")
MIMAT0003218_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-92b-3p.xlsx")
#see supplementary Target Genes (Primary)
MIMAT0000073_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-19a-3p.xlsx")
MIMAT0000099_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-101-3p.xlsx")
MIMAT0015066_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-3065-5p.xlsx")
MIMAT0014995_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miR-3130-5p.xlsx")
#_tg=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\.xlsx")

Primary_patients1=t(Primary_patients) # COL as gene name
Metastasis_patients1=t(Metastasis_patients) #COL as gene name

#get Target Genes NAMES for each miRNA and set up as list to find in Gene Expression datasets
#Metastasis
mir10b=MIMAT0000254_tg[,1];mir10b=t(mir10b); mir10b=(paste(mir10b, collapse="|")) #Need to transpose to create 'list' data frame; then combined as OR
mir99b=MIMAT0000689_tg[,1];mir99b=t(mir99b); mir99b=(paste(mir99b, collapse="|")) 
mir21=MIMAT0004494_tg[,1];mir21=t(mir21); mir21=(paste(mir21, collapse="|")) 
mir149=MIMAT0000450_tg[,1];mir149=t(mir149); mir149=(paste(mir149, collapse="|")) 
mir410=MIMAT0002171_tg[,1];mir410=t(mir410); mir410=(paste(mir410, collapse="|")) 
mir223=MIMAT0000280_tg[,1];mir223=t(mir223); mir223=(paste(mir223, collapse="|")) 
mir4746=MIMAT0019880_tg[,1];mir4746=t(mir4746); mir4746=(paste(mir4746, collapse="|")) 
mir491=MIMAT0002807_tg[,1];mir491=t(mir491); mir491=(paste(mir491, collapse="|")) 
mir200b=MIMAT0000318_tg[,1];mir200b=t(mir200b); mir200b=(paste(mir200b, collapse="|")) 
mir92b=MIMAT0003218_tg[,1];mir92b=t(mir92b); mir92b=(paste(mir92b, collapse="|")) 
#Primary
mir19a=MIMAT0000073_tg[,1];mir19a=t(mir19a); mir19a=(paste(mir19a, collapse="|"))  
mir101=MIMAT0000099_tg[,1];mir101=t(mir101); mir101=(paste(mir101, collapse="|"))
mir3065=MIMAT0015066_tg[,1];mir3065=t(mir3065); mir3065=(paste(mir3065, collapse="|"))
mir3130=MIMAT0014995_tg[,1];mir3130=t(mir3130); mir3130=(paste(mir3130, collapse="|"))

#Filter Data FROM miRNA target genes to Primary Genes or Metastasis Genes patients data
#Filter Data in Primary Patients Genes
mir19a_patients=data.frame(Primary_patients1) %>% dplyr::select(matches(mir19a)) #take from gene expression data selected mirna-target genes
mir101_patients=data.frame(Primary_patients1) %>% dplyr::select(matches(mir101))
mir3065_patients=data.frame(Primary_patients1) %>% dplyr::select(matches(mir3065))
mir3130_patients=data.frame(Primary_patients1) %>% dplyr::select(matches(mir3130))

#Filter Data in Metastasis
mir10b_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir10b))
mir99b_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir99b))
mir21_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir21))
mir149_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir149))
mir410_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir410))
mir223_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir223))
mir4746_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir4746))
mir491_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir491))
mir200b_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir200b))
mir92b_patients=data.frame(Metastasis_patients1) %>% dplyr::select(matches(mir92b))

# QC check value type; Gene Expression = log2(norm_count+1) where miRNA =log2(RPM+1)
# Constrcut matrix for correlation analysis
# QC-1 : check the order of patients data
Primary_mirna1=Primary_mirna1[order(row.names(Primary_mirna1)), ]
mir19a_patients=mir19a_patients[order(row.names(mir19a_patients)), ]
mir101_patients=mir19a_patients[order(row.names(mir101_patients)), ]
mir3065_patients=mir19a_patients[order(row.names(mir3065_patients)), ]
mir3130_patients=mir19a_patients[order(row.names(mir3130_patients)), ]

Metastasis_mirna1=Metastasis_mirna1[order(row.names(Metastasis_mirna1)), ]
mir10b_patients=mir10b_patients[order(row.names(mir10b_patients)),]
mir99b_patients=mir99b_patients[order(row.names(mir99b_patients)),]
mir21_patients=mir21_patients[order(row.names(mir21_patients)),]
mir149_patients=mir149_patients[order(row.names(mir149_patients)),]
mir410_patients=mir410_patients[order(row.names(mir410_patients)),]
mir223_patients=mir223_patients[order(row.names(mir223_patients)),]
mir4746_patients=mir4746_patients[order(row.names(mir4746_patients)),]
mir491_patients=mir491_patients[order(row.names(mir491_patients)),]
mir200b_patients=mir200b_patients[order(row.names(mir200b_patients)),]
mir92b_patients=mir92b_patients[order(row.names(mir92b_patients)),]

#Combined miRNA & Gene Expression
#Primary
mir19a_tg_patients=data.frame(mir19a_patients,Primary_mirna1[,1]);
colnames(mir19a_tg_patients)[1937] ='mir19a_3p'
mir101_tg_patients=data.frame(mir101_patients,Primary_mirna1[,2]);
colnames(mir101_tg_patients)[1937]='mir101_3p'
mir3065_tg_patients=data.frame(mir3065_patients,Primary_mirna1[,3]);
colnames(mir3065_tg_patients)[1937]='mir3065_5p'
mir3130_tg_patients=data.frame(mir3130_patients,Primary_mirna1[,4]);
colnames(mir3130_tg_patients)[1937]='mir3130_5p'
#QC : why TG of these miRNA have same number? is it because some of gene names is similar so it will be repeated?

#Metastasis
mir10b_tg_patients=data.frame(mir10b_patients,Metastasis_mirna1[,9]);
colnames(mir10b_tg_patients)[656] = 'mir10b_3p' #change col names based on last col
mir99b_tg_patients=data.frame(mir99b_patients,Metastasis_mirna1 [,10]);
colnames(mir99b_tg_patients)[100] = 'mir99b_5p'
mir21_tg_patients=data.frame(mir21_patients,Metastasis_mirna1 [,1]);
colnames(mir21_tg_patients)[1474] = 'mir21_3p'
mir149_tg_patients=data.frame(mir149_patients,Metastasis_mirna1[,2]);
colnames(mir149_tg_patients)[1358]= 'mir149_5p'
mir410_tg_patients=data.frame(mir410_patients,Metastasis_mirna1[,3]);
colnames(mir410_tg_patients)[73]= 'mir410_3p'
mir223_tg_patients=data.frame(mir223_patients,Metastasis_mirna1[,4]);
colnames(mir223_tg_patients)[902]= 'mir223_3p'
mir4746_tg_patients=data.frame(mir4746_patients,Metastasis_mirna1[,5]);
colnames(mir4746_tg_patients)[476]= 'mir4746_5p'
mir491_tg_patients=data.frame(mir491_patients,Metastasis_mirna1[,6]);
colnames(mir491_tg_patients)[1624]= 'mir491_5p'
mir200b_tg_patients=data.frame(mir200b_patients,Metastasis_mirna1[,7]);
colnames(mir200b_tg_patients)[837]= 'mir200b_3p'
mir92b_tg_patients=data.frame(mir92b_patients,Metastasis_mirna1[,8]);
colnames(mir92b_tg_patients)[942]= 'mir92b_3p'

### --------------------------------------- correlation miRNA & TG -------------------------------------------###
# This can be reconstruct without combined matrix, since the corr matrix will be to big
library(Hmisc)
#https://www.tutorialspoint.com/how-to-find-the-correlation-matrix-with-p-values-for-an-r-data-frame
# PRIMARY
cor_mir19a_tg=rcorr(as.matrix(mir19a_tg_patients1),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir19a_tg=data.frame(cor_mir19a_tg$r['mir19a_3p',],cor_mir19a_tg$P['mir19a_3p',]); colnames(cor_mir19a_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir19a_tg=data.frame(row.names(cor_mir19a_tg),cor_mir19a_tg)

cor_mir101_tg=rcorr(as.matrix(mir101_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir101_tg=data.frame(cor_mir101_tg$r['mir101_3p',],cor_mir101_tg$P['mir101_3p',]); colnames(cor_mir101_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir101_tg=data.frame(row.names(cor_mir101_tg),cor_mir101_tg)

cor_mir3065_tg=rcorr(as.matrix(mir3065_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir3065_tg=data.frame(cor_mir3065_tg$r['mir3065_5p',],cor_mir3065_tg$P['mir3065_5p',]); colnames(cor_mir3065_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir3065_tg=data.frame(row.names(cor_mir3065_tg),cor_mir3065_tg)

#cor_mir3065_tg=rcorr(as.matrix(mir3065_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
#cor_mir3065_tg=data.frame(cor_mir3065_tg$r['mir3065_5p',],cor_mir3065_tg$P['mir3065_5p',]); colnames(cor_mir3065_tg)= c('rho','p_val') #extract only r-value and p-value

cor_mir3130_tg=rcorr(as.matrix(mir3130_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir3130_tg=data.frame(cor_mir3130_tg$r['mir3130_5p',],cor_mir3130_tg$P['mir3130_5p',]); colnames(cor_mir3130_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir3130_tg=data.frame(row.names(cor_mir3130_tg),cor_mir3130_tg)

# METASTASIS
#Error out of bound, meaning the mirna col name hasn't changed correctly 
cor_mir10b_tg=rcorr(as.matrix(mir10b_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir10b_tg=data.frame(cor_mir10b_tg$r['mir10b_3p',],cor_mir10b_tg$P['mir10b_3p',]); colnames(cor_mir10b_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir10b_tg=data.frame(row.names(cor_mir10b_tg),cor_mir10b_tg)

cor_mir99b_tg=rcorr(as.matrix(mir99b_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir99b_tg=data.frame(cor_mir99b_tg$r['mir99b_5p',],cor_mir99b_tg$P['mir99b_5p',]); colnames(cor_mir99b_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir99b_tg=data.frame(row.names(cor_mir99b_tg),cor_mir99b_tg)

cor_mir21_tg=rcorr(as.matrix(mir21_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir21_tg=data.frame(cor_mir21_tg$r['mir21_3p',],cor_mir21_tg$P['mir21_3p',]); colnames(cor_mir21_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir21_tg=data.frame(row.names(cor_mir21_tg),cor_mir21_tg)

cor_mir149_tg=rcorr(as.matrix(mir149_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir149_tg=data.frame(cor_mir149_tg$r['mir149_5p',],cor_mir149_tg$P['mir149_5p',]); colnames(cor_mir149_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir149_tg=data.frame(row.names(cor_mir149_tg),cor_mir149_tg)

cor_mir410_tg=rcorr(as.matrix(mir410_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir410_tg=data.frame(cor_mir410_tg$r['mir410_3p',],cor_mir410_tg$P['mir410_3p',]); colnames(cor_mir410_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir410_tg=data.frame(row.names(cor_mir410_tg),cor_mir410_tg)

cor_mir223_tg=rcorr(as.matrix(mir223_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir223_tg=data.frame(cor_mir223_tg$r['mir223_3p',],cor_mir223_tg$P['mir223_3p',]); colnames(cor_mir223_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir223_tg=data.frame(row.names(cor_mir223_tg),cor_mir223_tg)

cor_mir4746_tg=rcorr(as.matrix(mir4746_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir4746_tg=data.frame(cor_mir4746_tg$r['mir4746_5p',],cor_mir4746_tg$P['mir4746_5p',]); colnames(cor_mir4746_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir4746_tg=data.frame(row.names(cor_mir4746_tg),cor_mir4746_tg)

cor_mir491_tg=rcorr(as.matrix(mir491_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir491_tg=data.frame(cor_mir491_tg$r['mir491_5p',],cor_mir491_tg$P['mir491_5p',]); colnames(cor_mir491_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir491_tg=data.frame(row.names(cor_mir491_tg),cor_mir491_tg)

cor_mir200b_tg=rcorr(as.matrix(mir200b_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir200b_tg=data.frame(cor_mir200b_tg$r['mir200b_3p',],cor_mir200b_tg$P['mir200b_3p',]); colnames(cor_mir200b_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir200b_tg=data.frame(row.names(cor_mir200b_tg),cor_mir200b_tg)

cor_mir92b_tg=rcorr(as.matrix(mir92b_tg_patients),type="spearman") #wait for 20s if data too big; output will be large matrix with 8 list
cor_mir92b_tg=data.frame(cor_mir92b_tg$r['mir92b_3p',],cor_mir92b_tg$P['mir92b_3p',]); colnames(cor_mir92b_tg)= c('rho','p_val') #extract only r-value and p-value
cor_mir92b_tg=data.frame(row.names(cor_mir92b_tg),cor_mir92b_tg)

# SAVE-1 no cut-off performed
library(writexl)
sheets_p <- list("mir19a" = cor_mir19a_tg, "mir101" = cor_mir101_tg, "mir3065" = cor_mir3065_tg, "mir3130" = cor_mir3130_tg ) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets_p, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Primary_correlation_mir_tg.xlsx")
sheets_m <- list("mir10b" = cor_mir10b_tg, "mir99b" = cor_mir99b_tg, "mir21" = cor_mir21_tg, "mir149" = cor_mir149_tg, 
                 "mir410"=cor_mir410_tg,  "mir223"=cor_mir223_tg,"mir4746"=cor_mir4746_tg,"mir491"=cor_mir491_tg,
                 "mir200b"=cor_mir200b_tg,"mir92b"=cor_mir92b_tg ) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets_m, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Metastasis_correlation_mir_tg.xlsx")

### CUT-OFF
#primary
cor_mir19a_tg_f <- filter(cor_mir19a_tg, p_val< 0.05)
cor_mir101_tg_f <- filter(cor_mir101_tg,  p_val < 0.05)
cor_mir3065_tg_f <- filter(cor_mir3065_tg,  p_val < 0.05)
cor_mir3130_tg_f = filter(cor_mir3130_tg,  p_val < 0.05)
#Metastasis
cor_mir10b_tg_f= filter(cor_mir10b_tg,  p_val < 0.05)
cor_mir99b_tg_f= filter(cor_mir99b_tg,  p_val < 0.05)
cor_mir21_tg_f= filter(cor_mir21_tg,  p_val < 0.05)
cor_mir149_tg_f= filter(cor_mir149_tg,  p_val < 0.05)
cor_mir410_tg_f= filter(cor_mir410_tg,  p_val < 0.05)
cor_mir223_tg_f= filter(cor_mir223_tg,  p_val < 0.05)
cor_mir4746_tg_f= filter(cor_mir4746_tg,  p_val < 0.05)
cor_mir491_tg_f= filter(cor_mir491_tg,  p_val < 0.05)
cor_mir200b_tg_f= filter(cor_mir200b_tg,  p_val < 0.05)
cor_mir92b_tg_f = filter(cor_mir92b_tg,  p_val < 0.05)

# SAVE-2 no cut-off performed
sheets_pf <- list("mir19a" = cor_mir19a_tg_f, "mir101" = cor_mir101_tg_f, "mir3065" = cor_mir3065_tg_f, "mir3130" = cor_mir3130_tg_f ) #assume sheet1 and sheet2 are data frames
write_xlsx(sheets_pf, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Primary_Filtered_correlation_mir_tg.xlsx")
sheets_mf <- list("mir10b" = cor_mir10b_tg_f, "mir99b" = cor_mir99b_tg_f, "mir21" = cor_mir21_tg_f, "mir149" = cor_mir149_tg_f, 
                  "mir410"=cor_mir410_tg_f,  "mir223"=cor_mir223_tg_f,"mir4746"=cor_mir4746_tg_f,"mir491"=cor_mir491_tg_f,
                  "mir200b"=cor_mir200b_tg_f,"mir92b"=cor_mir92b_tg_f ) #assume sheet1 and sheet2 are data frames

write_xlsx(sheets_mf, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Metastasis_Filtered_correlation_mir_tg.xlsx")

#Meaning there might be indication those miRNA affecting those Target Genes lists 

#NOTE 1 miRNA & Gene Patients needs to be aligned
#Note 2 only gene listed in Target genes
#where colnames is patient name; other col will be gene & miRNA feature?