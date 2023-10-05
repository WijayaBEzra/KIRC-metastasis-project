### Kidney Cancer Project ####
# WE INVESTIGATE Kidney cancer USING BIOINFORMATICS MACHINE LEARNING PIPELINE
# REVISION : 2
# Purposes : Imputation Evaluation
# YEAR : 2022
# Addition code in 42 -50

library(caret)
library(dplyr)
library(stringr)

library(tibble)
library(edgeR)

library("DESeq2")
library(ggplot2)

###DATA EXTRACTION###
# Only choose read_count criteria
myFile = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your file manually E:\KIRC_Ez_2022 in WIndown PC
mirna_name=data.frame(myFile[,1])

####### Check Quality of Missing Data miRNA #######
count_na <- data.frame(rowSums(is.na(myFile)))
count_na <- data.frame(mirna_name,count_na)

QC_count_na <- (rowSums(is.na(myFile))/length (myFile))
QC_count_na_ <- data.frame(mirna_name,QC_count_na) #Check Percentage of missing value per-Row miRNA

#Show barplot
barplot(QC_count_na_[1:100,2],
        main = "Missing data in 100 miRNA",
        xlab = "Missing Data in %",
        ylab = "miRNA",
        names.arg = QC_count_na_[1:100,1],
        col = "red",
        horiz = TRUE)

# Check missing value higher to certain %
sum_NA=sum(QC_count_na_$QC_count_na > 0.5) #counting missing value >0.5
sum_NA_7=sum(QC_count_na_$QC_count_na > 0.7) #counting missing value >0.7
sum_NA_9=sum(QC_count_na_$QC_count_na > 0.9) #counting missing value >0.9
data_mir=QC_count_na_[QC_count_na_$QC_count_na < 0.5,] # Show you miRNA with less than 50% missing value
data_mir_clear=QC_count_na_[QC_count_na_$QC_count_na < 0.001,] #show miRNA no need imputation

#Choose only miRNA data with <50% missing Value
myFile_=myFile[row.names(data_mir),] #since Row names in data_mir is fixed (use fixed index instead variable name)
miRNA_name=data.frame(myFile_[,1]) #mirna name considerable in <50% missing value

##### Divide Data into normal - Primary - Metastasis #####
#Seperate Normal-Primary Cancer
Normal=myFile_[,grep(".11", colnames(myFile_))] #miRNA name is lost in this data
Primary=myFile_[,grep(".01", colnames(myFile_))] #Have both Primary & Metastasis

#Get Clinical data metastasis
myFile2 = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your clinical file manually
new_row= data.frame(myFile2$sampleID,myFile2$pathologic_M) #only choose col SampleID & TNM (M) stage

cc=dplyr::filter(new_row, grepl('M1', myFile2$pathologic_M)) #Grab row only M1
cc_=str_replace_all(cc$myFile2.sampleID,"-",".") # edited TCGA name in from clinical data which use '-' in the patient's name
da=(paste(cc_, collapse="|")) #To make it as "OR"
Metastasis=Primary %>% dplyr::select(matches(da)) #calling metastasis data from Primary Data

#Update Primary M0 only
cp=dplyr::filter(new_row, grepl('M0', myFile2$pathologic_M)) #Grab patients row only M0
cp_=str_replace_all(cp$myFile2.sampleID,"-",".")
dap=(paste(cp_, collapse="|")) #To make it as "OR"
Primary_=Primary %>% dplyr::select(matches(dap)) #take Patients name in clinical with M0 cathegory


##### Imputation #####
library(imputeR)
#Imputation info from https://cran.r-project.org/web/packages/imputeR/imputeR.pdf
# Note : some publication DOI:10.1186/s12859-015-0853-0
#                         https://doi.org/10.1186/s12859-022-04656-4
#                         https://doi.org/10.1038/s41598-021-03438-x
impdata_N <- impute(Normal, lmFun = "lassoR") #using Lasso so in Normal can be inputed (Filled the data)
myFile_Normal <- data.frame(impdata_N$imp) #Data Frame Filled with lasso imputation

#Simulation
library("Cubist")
library("gbm")
library("mboost")
library("pls")
library("ridge")
#need to check imputation method one by one

miss_percent = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5) #nature paper + our 50% missing data problem
Imputation_method = list("lassoR",
                         "CubistR", 
                         "glmboostR",
                         "pcrR",  
                         "plsR",
                         "ridgeR", 
                         "stepBackR",
                         "stepBothR",
                         "stepForR")

#### FOR LOOP for assing Missing data ####
#NOTE : For TESTING ONLY
#Testing Performance using imputed Normal Datasets (MANUAL)
missdata=SimIm(myFile_Normal, 0.1) # using Missing data
impute_test=impute(missdata, lmFun = "lassoR")
lasso_eval=SimEval(myFile_Normal, p=0.5, n.sim=1, method = "lassoR")
Rmse(impute_test$imp, missdata, myFile_Normal, norm = TRUE)

#Loop Function to test simulation (AUTOMATIC)
for (q in Imputation_method[3:9]){
  for (x in miss_percent){
    assign(paste("eval",q,x,sep="_"),SimEval(myFile_Normal, p=x, n.sim=5, method = q)) #Looping evaluation for (5 try iteration)
  }
}

#Imputation_method[4] cubit on hold need 'cli' (internet problems)
#UNDER DEVELOPMENT : Loop Function to test simulation (AUTOMATIC with Error calling)
for (q in Imputation_method[3:9]){
  for (x in miss_percent){
    tryCatch({
      assign(paste("eval",q,x,sep="_"),SimEval(myFile_Normal, p=x, n.sim=1, method = q))
    },error=function(e){cat("ERROR:",e)})
  }
}
#NOte: CubistR freeze

Simulation_lassoR_RMSE=data.frame(eval_lassoR_0.05$error,eval_lassoR_0.1$error,eval_lassoR_0.15$error,eval_lassoR_0.2$error,eval_lassoR_0.25$error,
                                  eval_lassoR_0.3$error,eval_lassoR_0.35$error,eval_lassoR_0.4$error,eval_lassoR_0.45$error,eval_lassoR_0.5$error)
Simulation_lassoR_Time=data.frame(eval_lassoR_0.05$time,eval_lassoR_0.1$time,eval_lassoR_0.15$time,eval_lassoR_0.2$time,eval_lassoR_0.25$time,
                                  eval_lassoR_0.3$time,eval_lassoR_0.35$time,eval_lassoR_0.4$time,eval_lassoR_0.45$time,eval_lassoR_0.5$time)


Simulation_glmbosstR_RMSE=data.frame(eval_glmboostR_0.05$error,eval_glmboostR_0.1$error,eval_glmboostR_0.15$error,eval_glmboostR_0.2$error,eval_glmboostR_0.25$error,
                                     eval_glmboostR_0.3$error,eval_glmboostR_0.35$error,eval_glmboostR_0.4$error,eval_glmboostR_0.45$error,eval_glmboostR_0.5$error)
Simulation_glmbosstR_Time=data.frame(eval_glmboostR_0.05$time,eval_glmboostR_0.1$time,eval_glmboostR_0.15$time,eval_glmboostR_0.2$time,eval_glmboostR_0.25$time, eval_glmboostR_0.3$time,
                                     eval_glmboostR_0.35$time,eval_glmboostR_0.4$time,eval_glmboostR_0.45$time,eval_glmboostR_0.5$time)

Simulation_pcrR_RMSE=data.frame(eval_pcrR_0.05$error,eval_pcrR_0.1$error,eval_glmboostR_0.15$error,eval_glmboostR_0.2$error,eval_glmboostR_0.25$error,
                                eval_glmboostR_0.3$error,eval_glmboostR_0.35$error,eval_glmboostR_0.4$error,eval_glmboostR_0.45$error,eval_glmboostR_0.5$error)
Simulation_glmbosstR_Time=data.frame(eval_glmboostR_0.05$time,eval_glmboostR_0.1$time,eval_glmboostR_0.15$time,eval_glmboostR_0.2$time,eval_glmboostR_0.25$time, eval_glmboostR_0.3$time,
                                     eval_glmboostR_0.35$time,eval_glmboostR_0.4$time,eval_glmboostR_0.45$time,eval_glmboostR_0.5$time)

Simulation_pcrR_RMSE=data.frame(eval_pcrR_0.05$error,eval_pcrR_0.1$error,eval_pcrR_0.15$error,eval_pcrR_0.2$error,eval_pcrR_0.25$error,
                                eval_pcrR_0.3$error,eval_pcrR_0.35$error,eval_pcrR_0.4$error,eval_pcrR_0.45$error,eval_pcrR_0.5$error)
Simulation_pcrR_Time=data.frame(eval_pcrR_0.05$time,eval_pcrR_0.1$time,eval_pcrR_0.15$time,eval_pcrR_0.2$time,eval_pcrR_0.25$time, eval_pcrR_0.3$time,
                                eval_pcrR_0.35$time,eval_pcrR_0.4$time,eval_pcrR_0.45$time,eval_pcrR_0.5$time)

Simulation_pcrR_RMSE=data.frame(eval_pcrR_0.05$error,eval_pcrR_0.1$error,eval_pcrR_0.15$error,eval_pcrR_0.2$error,eval_pcrR_0.25$error,
                                eval_pcrR_0.3$error,eval_pcrR_0.35$error,eval_pcrR_0.4$error,eval_pcrR_0.45$error,eval_pcrR_0.5$error)
Simulation_pcrR_Time=data.frame(eval_pcrR_0.05$time,eval_pcrR_0.1$time,eval_pcrR_0.15$time,eval_pcrR_0.2$time,eval_pcrR_0.25$time, eval_pcrR_0.3$time,
                                eval_pcrR_0.35$time,eval_pcrR_0.4$time,eval_pcrR_0.45$time,eval_pcrR_0.5$time)

Simulation_plsR_RMSE=data.frame(eval_plsR_0.05$error,eval_plsR_0.1$error,eval_plsR_0.15$error,eval_plsR_0.2$error,eval_plsR_0.25$error,
                                eval_plsR_0.3$error,eval_plsR_0.35$error,eval_plsR_0.4$error,eval_plsR_0.45$error,eval_plsR_0.5$error)
Simulation_plsR_Time=data.frame(eval_plsR_0.05$time,eval_plsR_0.1$time,eval_plsR_0.15$time,eval_plsR_0.2$time,eval_plsR_0.25$time, eval_plsR_0.3$time,
                                eval_plsR_0.35$time,eval_plsR_0.4$time,eval_plsR_0.45$time,eval_plsR_0.5$time)

Simulation_ridgeR_RMSE=data.frame(eval_ridgeR_0.05$error,eval_ridgeR_0.1$error,eval_ridgeR_0.15$error,eval_ridgeR_0.2$error,eval_ridgeR_0.25$error,
                                  eval_ridgeR_0.3$error,eval_ridgeR_0.35$error,eval_ridgeR_0.4$error,eval_ridgeR_0.45$error,eval_ridgeR_0.5$error)
Simulation_ridgeR_Time=data.frame(eval_ridgeR_0.05$time,eval_ridgeR_0.1$time,eval_ridgeR_0.15$time,eval_ridgeR_0.2$time,eval_ridgeR_0.25$time, eval_ridgeR_0.3$time,
                                  eval_ridgeR_0.35$time,eval_ridgeR_0.4$time,eval_ridgeR_0.45$time,eval_ridgeR_0.5$time)

Simulation_stepBackR_RMSE=data.frame(eval_stepBackR_0.05$error,eval_stepBackR_0.1$error,eval_stepBackR_0.15$error,eval_stepBackR_0.2$error,eval_stepBackR_0.25$error,
                                     eval_stepBackR_0.3$error,eval_stepBackR_0.35$error,eval_stepBackR_0.4$error,eval_stepBackR_0.45$error,eval_stepBackR_0.5$error)
Simulation_stepBackR_Time=data.frame(eval_stepBackR_0.05$time,eval_stepBackR_0.1$time,eval_stepBackR_0.15$time,eval_stepBackR_0.2$time,eval_stepBackR_0.25$time, eval_stepBackR_0.3$time,
                                     eval_stepBackR_0.35$time,eval_stepBackR_0.4$time,eval_stepBackR_0.45$time,eval_stepBackR_0.5$time)

Simulation_stepBothR_RMSE=data.frame(eval_stepBothR_0.05$error,eval_stepBothR_0.1$error,eval_stepBothR_0.15$error,eval_stepBothR_0.2$error,eval_stepBothR_0.25$error,
                                     eval_stepBothR_0.3$error,eval_stepBothR_0.35$error,eval_stepBothR_0.4$error,eval_stepBothR_0.45$error,eval_stepBothR_0.5$error)
Simulation_stepBothR_Time=data.frame(eval_stepBothR_0.05$time,eval_stepBothR_0.1$time,eval_stepBothR_0.15$time,eval_stepBothR_0.2$time,eval_stepBothR_0.25$time, eval_stepBothR_0.3$time,
                                     eval_stepBothR_0.35$time,eval_stepBothR_0.4$time,eval_stepBothR_0.45$time,eval_stepBothR_0.5$time)

Simulation_stepForR_RMSE=data.frame(eval_stepForR_0.05$error,eval_stepForR_0.1$error,eval_stepForR_0.15$error,eval_stepForR_0.2$error,eval_stepForR_0.25$error,
                                    eval_stepForR_0.3$error,eval_stepForR_0.35$error,eval_stepForR_0.4$error,eval_stepForR_0.45$error,eval_stepForR_0.5$error)
Simulation_stepForR_Time=data.frame(eval_stepForR_0.05$time,eval_stepForR_0.1$time,eval_stepForR_0.15$time,eval_stepForR_0.2$time,eval_stepForR_0.25$time, eval_stepForR_0.3$time,
                                    eval_stepForR_0.35$time,eval_stepForR_0.4$time,eval_stepForR_0.45$time,eval_stepForR_0.5$time)

##combined all Result (standarized all data frame col name to joint the row)##

#Simulation_lassoR_RMSE_=Simulation_lassoR_RMSE
#avg=data.frame(colMeans(Simulation_lassoR_RMSE_))
#Simulation_glmbosstR_RMSE_=Simulation_glmbosstR_RMSE
#avg_=data.frame(colMeans(Simulation_glmbosstR_RMSE_))
#test_combine=data.frame(avg,avg_)

#Standarized COl names RMSE
colnames(Simulation_lassoR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_glmbosstR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
#colnames(Simulation_CubistR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_pcrR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_plsR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_ridgeR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_stepBackR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_stepBothR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_stepForR_RMSE)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
# for Time
colnames(Simulation_lassoR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_glmbosstR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
#colnames(Simulation_CubistR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_pcrR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_plsR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_ridgeR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_stepBackR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_stepBothR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')
colnames(Simulation_stepForR_Time)=c('0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5')

#Average RMSE each Experiment
Avr_lassoR_RMSE=data.frame(colMeans(Simulation_lassoR_RMSE))
Avr_glmbosstR_RMSE=data.frame(colMeans(Simulation_glmbosstR_RMSE))
Avr_pcrR_RMSE=data.frame(colMeans(Simulation_pcrR_RMSE))
Avr_plsR_RMSE=data.frame(colMeans(Simulation_plsR_RMSE))
Avr_ridgeR_RMSE=data.frame(colMeans(Simulation_ridgeR_RMSE))
Avr_stepBackR_RMSE=data.frame(colMeans(Simulation_stepBackR_RMSE))
Avr_stepBothR_RMSE=data.frame(colMeans(Simulation_stepBothR_RMSE))
Avr_stepForR_RMSE=data.frame(colMeans(Simulation_stepForR_RMSE))

All_Avr_RMSE=data.frame(Avr_lassoR_RMSE,Avr_glmbosstR_RMSE,Avr_pcrR_RMSE,Avr_plsR_RMSE,Avr_ridgeR_RMSE,
                        Avr_stepBackR_RMSE,Avr_stepBothR_RMSE,Avr_stepForR_RMSE)
col.names(All_Avr_RMSE) <- c('lassoR','glmbosstR','pcrR','plsR','ridgeR','stepBackR','stepBothR','stepForR')

#All_Simulation_RMSE=rbind(Simulation_lassoR_RMSE,Simulation_glmbosstR_RMSE,Simulation_pcrR_RMSE,Simulation_plsR_RMSE,Simulation_ridgeR_RMSE,
#Simulation_stepBackR_RMSE,Simulation_stepBothR_RMSE,Simulation_stepForR_RMSE)
#Simulation_lassoR_Time_=data.frame(t(Simulation_lassoR_Time))

All_Simulation_Time=rbind(Simulation_lassoR_Time,Simulation_glmbosstR_Time,Simulation_pcrR_Time,Simulation_plsR_Time,Simulation_ridgeR_Time,
                          Simulation_stepBackR_Time,Simulation_stepBothR_Time,Simulation_stepForR_Time)
row.names(All_Simulation_Time) <- c('lassoR','glmbosstR','pcrR','plsR','ridgeR','stepBackR','stepBothR','stepForR')

#$# Save [note: this directory using WIN PC i636]
library(writexl)
write_xlsx(Simulation_lassoR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_lassoR_RMSE.xlsx")
write_xlsx(Simulation_glmbosstR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_glmbosstR_RMSE.xlsx")
write_xlsx(Simulation_pcrR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_pcrR_RMSE.xlsx")
write_xlsx(Simulation_plsR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_plsR_RMSE.xlsx")
write_xlsx(Simulation_ridgeR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_ridgeR_RMSE.xlsx")
write_xlsx(Simulation_stepBackR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_stepBackR_RMSE.xlsx")
write_xlsx(Simulation_stepBothR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_stepBothR_RMSE.xlsx")
write_xlsx(Simulation_stepForR_RMSE,"E:\\KIRC_Ez_2022\\Simulation_stepForR_RMSE.xlsx")

write_xlsx(Simulation_lassoR_Time,"E:\\KIRC_Ez_2022\\Simulation_lassoR_Time.xlsx")
write_xlsx(Simulation_glmbosstR_Time,"E:\\KIRC_Ez_2022\\Simulation_glmbosstR_Time.xlsx")
write_xlsx(Simulation_pcrR_Time,"E:\\KIRC_Ez_2022\\Simulation_pcrR_Time.xlsx")
write_xlsx(Simulation_plsR_Time,"E:\\KIRC_Ez_2022\\Simulation_plsR_Time.xlsx")
write_xlsx(Simulation_ridgeR_Time,"E:\\KIRC_Ez_2022\\Simulation_ridgeR_Time.xlsx")
write_xlsx(Simulation_stepBackR_Time,"E:\\KIRC_Ez_2022\\Simulation_stepBackR_Time.xlsx")
write_xlsx(Simulation_stepBothR_Time,"E:\\KIRC_Ez_2022\\Simulation_stepBothR_Time.xlsx")
write_xlsx(Simulation_stepForR_Time,"E:\\KIRC_Ez_2022\\Simulation_stepForR_Time.xlsx")

write_xlsx(All_Simulation_Time,"E:\\KIRC_Ez_2022\\All_Simulation_Time.xlsx")
write_xlsx(All_Avr_RMSE,"E:\\KIRC_Ez_2022\\All_Avr_RMSE.xlsx")
# NOTE : when save to excel we need to rename again

#[22-06-2023] update for reviewer
write_xlsx(QC_count_na_,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Supplementary Information RENAMED\\Review & Revision\\QC_count_na_.xlsx")
write_xlsx(data_mir,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Supplementary Information RENAMED\\Review & Revision\\miRNA_less50.xlsx")
write_xlsx(data_mir_clear,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Supplementary Information RENAMED\\Review & Revision\\miRNA_no missing value.xlsx")
