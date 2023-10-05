## Validation Code ##
# Option
# 1. Using CCLE to verify Expression level of miRNA [metas cell & primary]
#     mirNA : https://depmap.org/portal/download/all/
#              [select data view : miRNA expression]
#     Gene  : CCLE_Expression_Entrez_2012-10-18.res
#
# 2. Using survival analysis independent database (GEO)

# Tools
# 1. Barplot analysis
#.2. Cluster analysis

# Additional
# 1.Drugs targeting miRNA or Genes

library("ggplot2")
library("readxl")
library(dplyr)

### VALIDATION BY KM-PLOT GEO STUDY ###
library(ggfortify) #need to install
library(haven) #For importing STATA files.
library(survival)

url <- "http://web1.sph.emory.edu/dkleinb/allDatasets/surv2datasets/addicts.dta"
addicts <- data.frame(read_dta(url))
#only status need to be as.Factor 
addicts$id <- as.numeric(addicts$id)
addicts$clinic <- as.factor(addicts$clinic)
addicts$status <- as.numeric(addicts$status)
addicts$survt <- as.numeric(addicts$survt)
addicts$prison <- as.numeric(addicts$prison)
addicts$dose <- as.numeric(addicts$dose)

Y = Surv(addicts$survt, addicts$status == 1)
kmfit = survfit(Y ~ addicts$clinic)
model_fit <- survfit(Surv(survt, status) ~ clinic, data = addicts)

model_fit <- survfit(Surv(survt, status) ~ clinic, data = addicts) #survt = time in days; #Status=life or death [1=die] #Clinic = group to compore [can be more than 1]

autoplot(model_fit) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n miRNA-Cox Primary-Normal KIRC with High & Low  \n") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))

# Test 3 condition [KM plot only achieve 2 condition]
status2=sample(1:3,238, replace=TRUE) #create random value 1 : 3
addicts=data.frame(addicts,status2)

model_fit2 <- survfit(Surv(survt, status2) ~ clinic, data = addicts)
autoplot(model_fit2) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n miRNA-Cox Primary-Normal KIRC with High & Low  \n") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))

### USING TCGA miRNA ###
### Input data miRNA with Imputation Set [maybe from Cox?] ###
library(ggfortify) #need to install
library(haven) #For importing STATA files.
library(survival)
library(survminer)

library("readxl")
library(dplyr)
library(stringr)
Metastasis=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Metastasis.xlsx")
Primary=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Primary.xlsx")
Normal=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Normal.xlsx")

All_miRNA_expression = data.frame(Metastasis, Primary, Normal)
All_miRNA_expression_ = t(All_miRNA_expression)
colnames(All_miRNA_expression_)=All_miRNA_expression_ [1,]; All_miRNA_expression_ = All_miRNA_expression_[-c(1),]

#Get only miRNA in Primary Cox
All_miRNA_expression_Primary = data.frame(All_miRNA_expression_[,c("MIMAT0000073", "MIMAT0000099", "MIMAT0015066", "MIMAT0014995")])
All_miRNA_expression_Metastasis = data.frame(All_miRNA_expression_[,c("MIMAT0004494", "MIMAT0000450", "MIMAT0002171", "MIMAT0000280", "MIMAT0019880", 
                                                           "MIMAT0002807", "MIMAT0000318", "MIMAT0003218", "MIMAT0000254", "MIMAT0000689")])

#Get Clinical data
myFile2 = read.delim(file.choose(), header = TRUE) #load file manually
time=data.frame(myFile2$days_to_death,myFile2$days_to_last_followup) #Combined time days_to_death & Days to follow up by choosing the bigger dayss
time[is.na(time)] <- 0 #NA set as zero
time=pmax(time$myFile2.days_to_death, time$myFile2.days_to_last_followup)
time=data.frame(time)
patient_id=str_replace_all(myFile2$sampleID,"-",".")
vital_stat=str_replace_all(myFile2$vital_status,c(LIVING="0", DECEASED = "1"))
vital_stat=data.frame(vital_stat)
Censored_=data.frame(patient_id,time,vital_stat) #final clinical survival time
colnames(Censored_)=c("patient_ID","time","")

row.names(Censored_)=Censored_[,1] 
Censored_ = Censored_[,-c(1)]

#New Method Using Merge [ since the main focus is to create new data frame]
# Merge required each variable to set row.name as patients data; it will cause a lot of NA but removed easily
  #test=merge(as.data.frame(All_miRNA_expression_Metastasis), as.data.frame(Censored_), by='row.names', all=TRUE)
  #test=na.omit(test)

# Combine Clinical data with miRNA Expression Data in METASTASIS
miRNA_Metastasis=merge(as.data.frame(All_miRNA_expression_Metastasis), as.data.frame(Censored_), by='row.names', all=TRUE)
miRNA_Metastasis=na.omit(miRNA_Metastasis)

# Combine Clinical data with miRNA Expression Data in PRIMARY
miRNA_Primary=merge(as.data.frame(All_miRNA_expression_Primary), as.data.frame(Censored_), by='row.names', all=TRUE)
miRNA_Primary=na.omit(miRNA_Primary)
miRNA_Primary$MIMAT0000073 <- as.numeric(miRNA_Primary$MIMAT0000073) # test converting into numeric

#NEXT check either 1 or 0 for censored
## Primary [median]
#median(miRNA_Primary$MIMAT0000073)            # "2.068246e+01"
#mean(as.numeric(miRNA_Primary$MIMAT0000073))  # 15.70791
miRNA_Primary = miRNA_Primary %>% 
  mutate(x_MIMAT0000073 = if_else(MIMAT0000073 > median(miRNA_Primary$MIMAT0000073)  , 2, 1))
miRNA_Primary = miRNA_Primary %>% 
  mutate(x_MIMAT0000099 = if_else(MIMAT0000099 > median(miRNA_Primary$MIMAT0000099)  , 2, 1))
miRNA_Primary = miRNA_Primary %>% 
  mutate(x_MIMAT0015066 = if_else(MIMAT0015066 > median(miRNA_Primary$MIMAT0015066)  , 2, 1))
miRNA_Primary = miRNA_Primary %>% 
  mutate(x_MIMAT0014995 = if_else(MIMAT0014995 > median(miRNA_Primary$MIMAT0014995)  , 2, 1))

## Metastasis
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0004494 = if_else(MIMAT0004494 > median(miRNA_Metastasis$MIMAT0004494)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0000450 = if_else(MIMAT0000450 > median(miRNA_Metastasis$MIMAT0000450)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0002171 = if_else(MIMAT0002171 > median(miRNA_Metastasis$MIMAT0002171)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0000280 = if_else(MIMAT0000280 > median(miRNA_Metastasis$MIMAT0000280)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0019880 = if_else(MIMAT0019880 > median(miRNA_Metastasis$MIMAT0019880)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0002807 = if_else(MIMAT0002807 > median(miRNA_Metastasis$MIMAT0002807)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0000318 = if_else(MIMAT0000318 > median(miRNA_Metastasis$MIMAT0000318)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0003218 = if_else(MIMAT0003218 > median(miRNA_Metastasis$MIMAT0003218)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0000254 = if_else(MIMAT0000254 > median(miRNA_Metastasis$MIMAT0000254)  , 2, 1))
miRNA_Metastasis = miRNA_Metastasis %>% 
  mutate(x_MIMAT0000689 = if_else(MIMAT0000318 > median(miRNA_Metastasis$MIMAT0000318)  , 2, 1))

## Data Structuring before KM-Plot Survival
names(miRNA_Primary)[7] = "status" #replace the name for status censored & uncensored in variable
miRNA_Primary$status <- as.numeric(miRNA_Primary$status)
miRNA_Primary$time <- as.numeric(miRNA_Primary$time)

names(miRNA_Metastasis)[13] = "status" #replace the name for status censored & uncensored in variable
miRNA_Metastasis$status <- as.numeric(miRNA_Metastasis$status)
miRNA_Metastasis$time <- as.numeric(miRNA_Metastasis$time)

### KM-Plot ###
# MIMAT0000073 OR hsa-miR-19a-3p
Y = Surv(miRNA_Primary$time, miRNA_Primary$status == 1)
kmfit_Y = survfit(Y ~ miRNA_Primary$x_MIMAT0000073) #MIMAT0000073
model_fit_Y <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = miRNA_Primary)
autoplot(model_fit_Y) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-19a-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Y) #0.084

# MIMAT0000099 or hsa-miR-101-3p
#Y2 = Surv(miRNA_Primary$time, miRNA_Primary$status == 1)
kmfit_Y2 = survfit(Y ~ miRNA_Primary$x_MIMAT0000099) #MIMAT0000099
model_fit_Y2 <- survfit(Surv(time, status) ~ x_MIMAT0000099, data = miRNA_Primary)
autoplot(model_fit_Y2) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-101-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Y2) #0.23

# MIMAT0015066		OR		hsa-miR-3065-5p
#Y2 = Surv(miRNA_Primary$time, miRNA_Primary$status == 1)
kmfit_Y3 = survfit(Y ~ miRNA_Primary$x_MIMAT0015066) #MIMAT0015066
model_fit_Y3 <- survfit(Surv(time, status) ~ x_MIMAT0015066, data = miRNA_Primary)
autoplot(model_fit_Y3) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-3065-5p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Y3) #0.099

# MIMAT0014995	OR	hsa-miR-3130-5p
#Y2 = Surv(miRNA_Primary$time, miRNA_Primary$status == 1)
kmfit_Y4 = survfit(Y ~ miRNA_Primary$x_MIMAT0014995) #MIMAT0014995
model_fit_Y4 <- survfit(Surv(time, status) ~ x_MIMAT0014995, data = miRNA_Primary)
autoplot(model_fit_Y4) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-3130-5p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Y4) #0.6

#Metastasis miRNA
# MIMAT0004494 OR hsa-miR-21-3p
Z = Surv(miRNA_Metastasis$time, miRNA_Metastasis$status == 1)
kmfit_Z = survfit(Z ~ miRNA_Metastasis$x_MIMAT0004494) #MIMAT0004494
model_fit_Z <- survfit(Surv(time, status) ~ x_MIMAT0004494, data = miRNA_Metastasis)
autoplot(model_fit_Z) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-21-3p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z) #0.86

#	MIMAT0000450 OR hsa-miR-149-5p
kmfit_Z1 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0000450) #MIMAT0000450
model_fit_Z1 <- survfit(Surv(time, status) ~ x_MIMAT0000450, data = miRNA_Metastasis)
autoplot(model_fit_Z1) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-149-5p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z1) #0.7

#	MIMAT0002171 OR hsa-miR-410-3p
kmfit_Z2 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0002171) #MIMAT0002171
model_fit_Z2 <- survfit(Surv(time, status) ~ x_MIMAT0002171, data = miRNA_Metastasis)
autoplot(model_fit_Z2) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-410-3p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z2) #0.87

#	MIMAT0000280 OR hsa-miR-223-3p
kmfit_Z3 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0000280) #MIMAT0000280
model_fit_Z3 <- survfit(Surv(time, status) ~ x_MIMAT0000280, data = miRNA_Metastasis)
autoplot(model_fit_Z3) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-223-3p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z3) #0.11

#	MIMAT0019880 OR hsa-miR-4746-5p
kmfit_Z4 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0019880) #MIMAT0019880
model_fit_Z4 <- survfit(Surv(time, status) ~ x_MIMAT0019880, data = miRNA_Metastasis)
autoplot(model_fit_Z4) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-4746-5p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z4) #0.81

#	MIMAT0002807 OR hsa-miR-491-5p
kmfit_Z5 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0002807) #MIMAT0019880
model_fit_Z5 <- survfit(Surv(time, status) ~ x_MIMAT0002807, data = miRNA_Metastasis)
autoplot(model_fit_Z5) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-491-5p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z5) #0.92

#	MIMAT0000318 OR hsa-miR-200b-3p
kmfit_Z6 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0000318) #MIMAT0000318
model_fit_Z6 <- survfit(Surv(time, status) ~ x_MIMAT0000318, data = miRNA_Metastasis)
autoplot(model_fit_Z6) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-200b-3p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z6) #0.19

#	MIMAT0003218 OR	hsa-miR-92b-3p
kmfit_Z7 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0003218) #MIMAT0003218
model_fit_Z7 <- survfit(Surv(time, status) ~ x_MIMAT0003218, data = miRNA_Metastasis)
autoplot(model_fit_Z7) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-92b-3p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z7) #0.00045

#	MIMAT0000254 OR	hsa-miR-10b-5p
kmfit_Z8 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0000254) #MIMAT0000254
model_fit_Z8 <- survfit(Surv(time, status) ~ x_MIMAT0000254, data = miRNA_Metastasis)
autoplot(model_fit_Z8) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-10b-5p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z8)#0.84

#	MIMAT0000689 OR	hsa-miR-99b-5p
kmfit_Z9 = survfit(Z ~ miRNA_Metastasis$x_MIMAT0000689) #MIMAT0000689
model_fit_Z9 <- survfit(Surv(time, status) ~ x_MIMAT0000689, data = miRNA_Metastasis)
autoplot(model_fit_Z9) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-99b-5p Metastasis KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
surv_pvalue(model_fit_Z9)#0.19


##-------------- Combined miRNA model by selecting patients with criteria High & Low Expression------------------##
# Note : we check addition combination, NOT all possible combination Yet
library(gtools)
list_miR_primary = c("MIMAT0000073", "MIMAT0000099", "MIMAT0015066", "MIMAT0014995")
combination_miR_primary_2key=combinations(n=4,r=2,v=list_miR_primary,set=TRUE,repeats.allowed=F) # n=possible item, r=2 combination key without replacement
combination_miR_primary_3key=combinations(n=4,r=3,v=list_miR_primary,set=TRUE, repeats.allowed=F) #set= duplicate will be removed
  
print(nrow(combination_miR_primary_2key)) #6 combination
print(nrow(combination_miR_primary_3key)) #4 combination

#Primary MIMAT0000073 & MIMAT0000099
        see= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==2 & miRNA_Primary$x_MIMAT0000099 == 2 ,]  #Get the set fo 
        see2= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==1 & miRNA_Primary$x_MIMAT0000099 == 1 ,] #
        see_combine=rbind(see,see2) #combined base on col names
        
        Q = Surv(see_combine$time, see_combine$status == 1)
        kmfit_Q = survfit(Q ~ see_combine$x_MIMAT0000073) # it doesn't matter since the number is same 
        model_fit_Q <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = see_combine)
        autoplot(model_fit_Q) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-19a-3p & hsa-miR-101-3p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        # surv_pvalue(model_fit_Q) #Pval = 0.022

# Primary MIMAT0000073 & MIMAT0014995 [31.03.2023]        
        see_73_14995= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==2 & miRNA_Primary$x_MIMAT0014995 == 2 ,]  #Get the set fo 
        see2_73_14995_2= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==1 & miRNA_Primary$x_MIMAT0014995 == 1 ,] #
        see_combine_73_14995=rbind(see_73_14995,see2_73_14995_2) #combined base on col names
        
        Q_73_14995 = Surv(see_combine_73_14995$time, see_combine_73_14995$status == 1) #give note that we take time=days and status=1 is death
        kmfit_Q_73_14995 = survfit(Q_73_14995 ~ see_combine_73_14995$x_MIMAT0000073) # it doesn't matter since the number is same as we analyzed MIMAT0000073
        model_fit_Q_73_14995 <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = see_combine_73_14995)
        autoplot(model_fit_Q_73_14995) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-19a-3p & hsa-miR-3130-5p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        # surv_pvalue(model_fit_Q_73_14995) #Pval = 0.4

#Primary MIMAT0000073 & MIMAT0015066 [31.03.2023]   
        see_73_15066= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==2 & miRNA_Primary$x_MIMAT0015066 == 2 ,]  #Get the set fo 
        see2_73_15066_2= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==1 & miRNA_Primary$x_MIMAT0015066 == 1 ,] #
        see_combine_73_15066=rbind(see_73_15066,see2_73_15066_2) #combined base on col names
        
        Q_73_15066 = Surv(see_combine_73_15066$time, see_combine_73_15066$status == 1) #give note that we take time=days and status=1 is death
        kmfit_Q_73_15066 = survfit(Q_73_15066 ~ see_combine_73_15066$x_MIMAT0000073) # it doesn't matter since the number is same as we analyzed MIMAT0000073
        model_fit_Q_73_15066 <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = see_combine_73_15066)
        autoplot(model_fit_Q_73_15066) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-19a-3p & hsa-miR-3065-5p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        # surv_pvalue(model_fit_Q_73_15066) #Pval = 0.024

# Primary MIMAT0000099 &  MIMAT0014995 [31.03.2023] 
        see_99_14995= miRNA_Primary[miRNA_Primary$x_MIMAT0000099 ==2 & miRNA_Primary$x_MIMAT0014995 == 2 ,]  #Get the set fo 
        see2_99_14995_2= miRNA_Primary[miRNA_Primary$x_MIMAT0000099 ==1 & miRNA_Primary$x_MIMAT0014995 == 1 ,] #
        see_combine_99_14995=rbind(see_99_14995,see2_99_14995_2) #combined base on col names
        
        Q_99_14995 = Surv(see_combine_99_14995$time, see_combine_99_14995$status == 1) #give note that we take time=days and status=1 is death
        kmfit_Q_99_14995 = survfit(Q_99_14995 ~ see_combine_99_14995$x_MIMAT0000099) # it doesn't matter since the number is same as we analyzed MIMAT0000099
        model_fit_Q_99_14995 <- survfit(Surv(time, status) ~ x_MIMAT0000099, data = see_combine_99_14995)
        autoplot(model_fit_Q_99_14995) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-101-3p & hsa-miR-3130-5p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        # surv_pvalue(model_fit_Q_99_14995) #Find real p-value
        # Result : 0.57
        
# Primary MIMAT0000099 & MIMAT0015066 [31.03.2023]
        see_99_15066= miRNA_Primary[miRNA_Primary$x_MIMAT0000099 ==2 & miRNA_Primary$x_MIMAT0015066 == 2 ,]  #Get the set fo 
        see2_99_15066_2= miRNA_Primary[miRNA_Primary$x_MIMAT0000099 ==1 & miRNA_Primary$x_MIMAT0015066 == 1 ,] #
        see_combine_99_15066=rbind(see_99_15066,see2_99_15066_2) #combined base on col names
        
        Q_99_15066 = Surv(see_combine_99_15066$time, see_combine_99_15066$status == 1) #give note that we take time=days and status=1 is death
        kmfit_Q_99_15066 = survfit(Q_99_15066 ~ see_combine_99_15066$x_MIMAT0000099) # it doesn't matter since the number is same as we analyzed MIMAT0000099
        model_fit_Q_99_15066 <- survfit(Surv(time, status) ~ x_MIMAT0000099, data = see_combine_99_15066)
        surv_pvalue(model_fit_Q_99_15066) #Find real p-value
        #Result : p = 0.059
        autoplot(model_fit_Q_99_15066) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-101-3p & hsa-miR-3065-5p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        # surv_pvalue(model_fit_Q_99_15066) #Find real p-value = 0.059
        
# Primary MIMAT0014995 & MIMAT0015066 [31.03.2023] 
        see_14995_15066= miRNA_Primary[miRNA_Primary$x_MIMAT0014995 ==2 & miRNA_Primary$x_MIMAT0015066 == 2 ,]  #Get the set fo 
        see2_14995_15066_2= miRNA_Primary[miRNA_Primary$x_MIMAT0014995 ==1 & miRNA_Primary$x_MIMAT0015066 == 1 ,] #
        see_combine_14995_15066=rbind(see_14995_15066,see2_14995_15066_2) #combined base on col names
        
        Q_14995_15066 = Surv(see_combine_14995_15066$time, see_combine_14995_15066$status == 1) #give note that we take time=days and status=1 is death
        kmfit_Q_14995_15066 = survfit(Q_14995_15066 ~ see_combine_14995_15066$x_MIMAT0014995) # it doesn't matter since the number is same as we analyzed MIMAT0000099
        model_fit_Q_14995_15066 <- survfit(Surv(time, status) ~ x_MIMAT0014995, data = see_combine_14995_15066)
        surv_pvalue(model_fit_Q_14995_15066) #Find real p-value
        # Results : p = 0.43
                
#Primary MIMAT0000073 & MIMAT0000099 & MIMAT0015066
        see_3mir= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==2 & miRNA_Primary$x_MIMAT0000099 == 2 & miRNA_Primary$x_MIMAT0015066 == 2,]  #Get the set fo 
        see2_3mir= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==1 & miRNA_Primary$x_MIMAT0000099 == 1 & miRNA_Primary$x_MIMAT0015066 == 1,] #
        see_combine_3mir=rbind(see_3mir,see2_3mir) #combined base on col names
        
        Q3 = Surv(see_combine_3mir$time, see_combine_3mir$status == 1)
        kmfit_Q3 = survfit(Q3 ~ see_combine_3mir$x_MIMAT0000073) # it doesn't matter since the number is same; patietns who has characterstics High expression in 3 miRNA
        model_fit_Q3 <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = see_combine_3mir)
        autoplot(model_fit_Q3) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-19a-3p, hsa-miR-101-3p, hsa-miR-3065-5p  KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        # surv_pvalue(model_fit_Q3) #Find real p-value = 0.014

# MIMAT0000073 & MIMAT0000099 & MIMAT0014995 [31.03.2023] 
        see_73_99_14995 = miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==2 & miRNA_Primary$x_MIMAT0000099 ==2 & miRNA_Primary$x_MIMAT0014995 ==2,]  #Get the set fo 
        see2_73_99_14995_2= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==1 & miRNA_Primary$x_MIMAT0000099 ==1 & miRNA_Primary$x_MIMAT0014995 ==1,]  #
        see_combine_73_99_14995=rbind(see_73_99_14995,see2_73_99_14995_2) #combined base on col names
        
        Q_73_99_14995 = Surv(see_combine_73_99_14995$time, see_combine_73_99_14995$status == 1) #give note that we take time=days and status=1 is death
        kmfit_Q_73_99_14995 = survfit(Q_73_99_14995 ~ see_combine_73_99_14995$x_MIMAT0000073) # it doesn't matter since the number is same as we analyzed MIMAT0000099
        model_fit_Q_73_99_14995 <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = see_combine_73_99_14995)
        surv_pvalue(model_fit_Q_73_99_14995) #Find real p-value
        #Results p = 0.42 #if se set model as MIMAT0000073
        model_fit_Q_73_99_14995_2 <- survfit(Surv(time, status) ~ x_MIMAT0000099, data = see_combine_73_99_14995)
        surv_pvalue(model_fit_Q_73_99_14995_2) #Find real p-value
        # Results p-val = 0.42

# MIMAT0000073 & MIMAT0014995 & MIMAT0015066 [31.03.2023] 
        see_73_15066_14995 = miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==2 & miRNA_Primary$x_MIMAT0015066 ==2 & miRNA_Primary$x_MIMAT0014995 ==2,]  #Get the set fo 
        see2_73_15066_14995_2= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==1 & miRNA_Primary$x_MIMAT0015066 ==1 & miRNA_Primary$x_MIMAT0014995 ==1,]  #
        see_combine_73_15066_14995=rbind(see_73_15066_14995,see2_73_15066_14995_2) #combined base on col names
        
        Q_73_15066_14995 = Surv(see_combine_73_15066_14995$time, see_combine_73_15066_14995$status == 1) #give note that we take time=days and status=1 is death
        kmfit_Q_73_15066_14995 = survfit(Q_73_15066_14995 ~ see_combine_73_15066_14995$x_MIMAT0000073) # it doesn't matter since the number is same as we analyzed MIMAT0000099
        model_fit_Q_73_15066_14995 <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = see_combine_73_15066_14995)
        surv_pvalue(model_fit_Q_73_15066_14995) #Find real p-value
        #Results p = 0.42
        
#Primary MIMAT0000073 & MIMAT0000099 & MIMAT0015066 & MIMAT0014995	
        see_4mir= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==2 & miRNA_Primary$x_MIMAT0000099 == 2 
                                & miRNA_Primary$x_MIMAT0015066 == 2 & miRNA_Primary$x_MIMAT0014995 == 2,]  #Get the set fo 
        see2_4mir= miRNA_Primary[miRNA_Primary$x_MIMAT0000073 ==1 & miRNA_Primary$x_MIMAT0000099 == 1 
                                 & miRNA_Primary$x_MIMAT0015066 == 1 & miRNA_Primary$x_MIMAT0014995 == 1,] #
        see_combine_4mir=rbind(see_4mir,see2_4mir) #combined base on col names
        
        Q4 = Surv(see_combine_4mir$time, see_combine_4mir$status == 1)
        kmfit_Q4 = survfit(Q4 ~ see_combine_4mir$x_MIMAT0000073) # it doesn't matter since the number is same; patietns who has characterstics High expression in 3 miRNA
        model_fit_Q4 <- survfit(Surv(time, status) ~ x_MIMAT0000073, data = see_combine_4mir)
        autoplot(model_fit_Q4) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-19a-3p, hsa-miR-101-3p, hsa-miR-3065-5p, & hsa-miR-3130-5p  KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        surv_pvalue(model_fit_Q4) #Find real p-value
        #Results : p = 0.46
#NOTE : AS THE COMBINATION INCREASE THE SAMPLE IS LOWER, WHICH MAY CAUSING INNACURATE IN THE RESULTS


#Metastasis
#Note we dont know which combination is the best for KM-Plot
        library(gtools)
        list_miR_metas = c("MIMAT0004494", "MIMAT0000450", "MIMAT0002171", "MIMAT0000280", "MIMAT0019880", 
                        "MIMAT0002807", "MIMAT0000318", "MIMAT0003218", "MIMAT0000254", "MIMAT0000689")
        combination_miR_metas_4key=combinations(n=10,r=4,v=list_miR_metas,set=TRUE, repeats.allowed=F) # 4 combination key without replacement
        combination_miR_metas_2key=combinations(n=10,r=2,v=list_miR_metas,set=TRUE,repeats.allowed=F) # 2 combination key without replacement [ideal combination otherwise sample is too small]
        combination_miR_metas_3key=combinations(n=10,r=3,v=list_miR_metas,set=TRUE,repeats.allowed=F) # 3 combination key
        print(nrow(combination_miR_metas_4key)) # result 210 combination for fixed 4 key
        print(nrow(combination_miR_metas_2key)) # result 45 combination for fixed 2 key
        print(nrow(combination_miR_metas_3key)) # 120 combination for fixed 3 key
        
        combination_miR_metas_5key=combinations(n=10,r=5,v=list_miR_metas,set=TRUE, repeats.allowed=F) # 4 combination key without replacement
        combination_miR_metas_6key=combinations(n=10,r=6,v=list_miR_metas,set=TRUE,repeats.allowed=F) # 2 combination key without replacement [ideal combination otherwise sample is too small]
        combination_miR_metas_7key=combinations(n=10,r=7,v=list_miR_metas,set=TRUE,repeats.allowed=F) # 3 combination key
        combination_miR_metas_8key=combinations(n=10,r=8,v=list_miR_metas,set=TRUE,repeats.allowed=F)
        combination_miR_metas_9key=combinations(n=10,r=9,v=list_miR_metas,set=TRUE,repeats.allowed=F)
        combination_miR_metas_10key=combinations(n=10,r=10,v=list_miR_metas,set=TRUE,repeats.allowed=F)
        
        #library(combinat)
        #res <- combn(list_miR_metas)
        
#Metastasis	MIMAT0004494 & MIMAT0000450 
        seeM= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==2 & miRNA_Metastasis$x_MIMAT0000450 == 2 ,]  #Get the set fo 
        seeM2= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==1 & miRNA_Metastasis$x_MIMAT0000450 == 1 ,] #
        see_combineM=rbind(seeM,seeM2) #combined base on col names
        
        R = Surv(see_combineM$time, see_combineM$status == 1)
        kmfit_R = survfit(R ~ see_combineM$x_MIMAT0004494) # it doesn't matter since the number is same 
        model_fit_R <- survfit(Surv(time, status) ~ x_MIMAT0004494, data = see_combineM)
        autoplot(model_fit_R) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-21-3p & hsa-miR-149-5p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        # surv_pvalue(model_fit_R) #pval = 0.93
        
#Metastasis	MIMAT0004494 & MIMAT0000450 & MIMAT0002171
        seeM_3mir= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==2 & miRNA_Metastasis$x_MIMAT0000450 == 2 & miRNA_Metastasis$x_MIMAT0002171 == 2,]  #Get the set fo 
        seeM2_3mir= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==1 & miRNA_Metastasis$x_MIMAT0000450 == 1 & miRNA_Metastasis$x_MIMAT0002171 == 1 ,] #
        see_combineM_3mir=rbind(seeM_3mir,seeM2_3mir) #combined base on col names
        
        R2 = Surv(see_combineM_3mir$time, see_combineM_3mir$status == 1)
        kmfit_R2 = survfit(R2 ~ see_combineM_3mir$x_MIMAT0004494) # it doesn't matter since the number is same 
        model_fit_R2 <- survfit(Surv(time, status) ~ x_MIMAT0004494, data = see_combineM_3mir)
        autoplot(model_fit_R2) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-21-3p, hsa-miR-149-5p, & hsa-miR-410-3p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        surv_pvalue(model_fit_R2) #0.66

#Metastasis	MIMAT0004494 & MIMAT0000450 & MIMAT0002171 &  MIMAT0000280
        seeM_4mir= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==2 & miRNA_Metastasis$x_MIMAT0000450 == 2 
                                    & miRNA_Metastasis$x_MIMAT0002171 == 2 & miRNA_Metastasis$x_MIMAT0000280==2,]  #Get the set fo 
        seeM2_4mir= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==1 & miRNA_Metastasis$x_MIMAT0000450 == 1 
                                     & miRNA_Metastasis$x_MIMAT0002171 == 1 & miRNA_Metastasis$x_MIMAT0000280==1,] #
        see_combineM_4mir=rbind(seeM_4mir,seeM2_4mir) #combined base on col names
        
        R3 = Surv(see_combineM_4mir$time, see_combineM_4mir$status == 1)
        kmfit_R3 = survfit(R3 ~ see_combineM_4mir$x_MIMAT0004494) # it doesn't matter since the number is same 
        model_fit_R3 <- survfit(Surv(time, status) ~ x_MIMAT0004494, data = see_combineM_4mir)
        autoplot(model_fit_R3) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-21-3p, hsa-miR-149-5p, hsa-miR-410-3p, & hsa-miR-223-3p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))        
        surv_pvalue(model_fit_R3) #0.47
        
#Metastasis	MIMAT0004494 & MIMAT0000450 & MIMAT0002171 &  MIMAT0000280 & MIMAT0019880 [low sample!!! : only 17]
        seeM_5mir= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==2 & miRNA_Metastasis$x_MIMAT0000450 == 2 
                                    & miRNA_Metastasis$x_MIMAT0002171 == 2 & miRNA_Metastasis$x_MIMAT0000280==2
                                    & miRNA_Metastasis$x_MIMAT0019880 == 2,]  #Get the set fo 
        seeM2_5mir= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0004494 ==1 & miRNA_Metastasis$x_MIMAT0000450 == 1 
                                     & miRNA_Metastasis$x_MIMAT0002171 == 1 & miRNA_Metastasis$x_MIMAT0000280==1
                                     & miRNA_Metastasis$x_MIMAT0019880 == 1,] #
        see_combineM_5mir=rbind(seeM_5mir,seeM2_5mir) #combined base on col names
        
        R4 = Surv(see_combineM_5mir$time, see_combineM_5mir$status == 1)
        kmfit_R4 = survfit(R4 ~ see_combineM_5mir$x_MIMAT0004494) # it doesn't matter since the number is same 
        model_fit_R4 <- survfit(Surv(time, status) ~ x_MIMAT0004494, data = see_combineM_5mir)
        autoplot(model_fit_R4) + 
          labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
               title = "Survival Times Of \n hsa-miR-21-3p, hsa-miR-149-5p, hsa-miR-410-3p, & hsa-miR-223-3p KIRC  \n 1=low expression 2=high expression") + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
                axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
                legend.title = element_text(face="bold", size = 10))
        surv_pvalue(model_fit_R4) #0.72
        
# To see best combination need to check again combination of HR per miRNA [HR=? Pval=?]; ANSWER: using Survminer to get p-valu3
# 4 miRNA in Primary & 10miRNA metastasis were chosen due to significant in multiple cox-regression

        
        
        
##Special Code#------------------------------------Final Function combined--Loop ---------------------------------------------------###        
# Input: combination of miRNA lists
# Function-1 : Variable setting
# Function-2 : Survfit
# Output : P-value [in colwise]
# Required library : dplyr
        
New_miRNA_Metastasis=data.frame(miRNA_Metastasis[,-c(2:11)] ) #focus only value 1 or 2       
pval_2combination=data.frame(Pval_col = rep(NA, length(combination_miR_metas_2key[,1]))) #Store in here

#------[2]------ Two miRNA combination
x=0
for(x in 1:length(combination_miR_metas_2key[,1])){ 
      a=combination_miR_metas_2key[x,1] #grab first miRNA
      b=combination_miR_metas_2key[x,2] #grab second miRNA
      #grab colomn based on miRNA name combination
        a=select(New_miRNA_Metastasis,contains(a))
        b=select(New_miRNA_Metastasis,contains(b))
        up= New_miRNA_Metastasis[a == 2 & b == 2,]  #Get the set fo 
        down= New_miRNA_Metastasis[a ==1 & b == 1,] #
        comb=rbind(up,down) #combined base on col names
        
        # Main Function
        f=paste("survfit(Surv(time, status) ~", colnames(a),",data=comb)",sep="") #create character function
        model = eval(parse(text=paste(f))) #coverting text(chr) function into R-function
        k=data.frame(surv_pvalue(model)) #Find real p-value
      
      # Store
        pval_2combination[x,] = k$pval
        x=x+1
}
  
#------[3]------ Three miRNA combination
pval_3combination=data.frame(Pval_col = rep(NA, length(combination_miR_metas_3key[,1]))) #Store in here
x=0
for(x in 1:length(combination_miR_metas_3key[,1])){ 
  a=combination_miR_metas_3key[x,1] #grab first miRNA
  b=combination_miR_metas_3key[x,2] #grab second miRNA
  c=combination_miR_metas_3key[x,3] #grab second miRNA
  #grab colomn based on miRNA name combination
  a=select(New_miRNA_Metastasis,contains(a))
  b=select(New_miRNA_Metastasis,contains(b))
  c=select(New_miRNA_Metastasis,contains(c))
  up= New_miRNA_Metastasis[a == 2 & b == 2 & c==2,]  #Get the set fo 
  down= New_miRNA_Metastasis[a ==1 & b == 1 & c==1,] #
  comb=rbind(up,down) #combined base on col names
  
  # Main Function
  f=paste("survfit(Surv(time, status) ~", colnames(a),",data=comb)",sep="") #create character function
  model = eval(parse(text=paste(f))) #coverting text(chr) function into R-function
  k=data.frame(surv_pvalue(model)) #Find real p-value
  
  # Store
  pval_3combination[x,] = k$pval
  x=x+1
}

#------[4]------ Four miRNA combination
pval_4combination=data.frame(Pval_col = rep(NA, length(combination_miR_metas_4key[,1]))) #Store in here
x=0
for(x in 1:length(combination_miR_metas_4key[,1])){ 
  a=combination_miR_metas_4key[x,1] #grab first miRNA #CHECK THE INPUT
  b=combination_miR_metas_4key[x,2] #grab second miRNA
  c=combination_miR_metas_4key[x,3] #grab second miRNA
  d=combination_miR_metas_4key[x,4] #grab second miRNA
  #grab colomn based on miRNA name combination
  a=select(New_miRNA_Metastasis,contains(a))
  b=select(New_miRNA_Metastasis,contains(b))
  c=select(New_miRNA_Metastasis,contains(c))
  d=select(New_miRNA_Metastasis,contains(d))
  up= New_miRNA_Metastasis[a == 2 & b == 2 & c==2 & d==2,]  #Get the set fo 
  down= New_miRNA_Metastasis[a ==1 & b == 1 & c==1 & d==1,] #
  comb=rbind(up,down) #combined base on col names
  
  # Main Function
  f=paste("survfit(Surv(time, status) ~", colnames(a),",data=comb)",sep="") #create character function
  model = eval(parse(text=paste(f))) #coverting text(chr) function into R-function
  k=data.frame(surv_pvalue(model)) #Find real p-value
  
  # Store
  pval_4combination[x,] = k$pval
  x=x+1
}

#------[5]------ FIVE miRNA combination
pval_5combination=data.frame(Pval_col = rep(NA, length(combination_miR_metas_5key[,1]))) #Store in here
x=0
for(x in 1:length(combination_miR_metas_5key[,1])){ 
  a=combination_miR_metas_5key[x,1] #grab first miRNA #CHECK THE INPUT
  b=combination_miR_metas_5key[x,2] #grab second miRNA
  c=combination_miR_metas_5key[x,3] #grab Third miRNA
  d=combination_miR_metas_5key[x,4] #grab Fourth miRNA
  e=combination_miR_metas_5key[x,5] #grab Five miRNA
  
  #grab colomn based on miRNA name combination
  a=select(New_miRNA_Metastasis,contains(a))
  b=select(New_miRNA_Metastasis,contains(b))
  c=select(New_miRNA_Metastasis,contains(c))
  d=select(New_miRNA_Metastasis,contains(d))
  e=select(New_miRNA_Metastasis,contains(e))
  
  up= New_miRNA_Metastasis[a == 2 & b == 2 & c==2 & d==2 & e==2,]  #Get the set fo 
  down= New_miRNA_Metastasis[a ==1 & b == 1 & c==1 & d==1 & e==1,] #
  comb=rbind(up,down) #combined base on col names
  
  # Main Function
  f=paste("survfit(Surv(time, status) ~", colnames(a),",data=comb)",sep="") #create character function
  model = eval(parse(text=paste(f))) #coverting text(chr) function into R-function
  k=data.frame(surv_pvalue(model)) #Find real p-value
  
  # Store
  pval_5combination[x,] = k$pval
  x=x+1
}
#got Warning : 1: In pchisq(chi, df, lower.tail = FALSE) : NaNs produced

#------[6]------ SIX miRNA combination
pval_6combination=data.frame(Pval_col = rep(NA, length(combination_miR_metas_6key[,1]))) #Store in here
x=0
for(x in 1:length(combination_miR_metas_6key[,1])){ 
  a=combination_miR_metas_6key[x,1] #grab first miRNA #CHECK THE INPUT
  b=combination_miR_metas_6key[x,2] #grab second miRNA
  c=combination_miR_metas_6key[x,3] #grab Third miRNA
  d=combination_miR_metas_6key[x,4] #grab Fourth miRNA
  e=combination_miR_metas_6key[x,5] #grab Five miRNA
  f=combination_miR_metas_6key[x,6] #grab Five miRNA
  
  #grab colomn based on miRNA name combination
  a=select(New_miRNA_Metastasis,contains(a))
  b=select(New_miRNA_Metastasis,contains(b))
  c=select(New_miRNA_Metastasis,contains(c))
  d=select(New_miRNA_Metastasis,contains(d))
  e=select(New_miRNA_Metastasis,contains(e))
  f=select(New_miRNA_Metastasis,contains(f))
  
  up= New_miRNA_Metastasis[a == 2 & b == 2 & c==2 & d==2 & e==2 & f==2,]  #Get the set fo 
  down= New_miRNA_Metastasis[a ==1 & b == 1 & c==1 & d==1 & e==1 & f==1,] #
  comb=rbind(up,down) #combined base on col names
  
  # Main Function
  f=paste("survfit(Surv(time, status) ~", colnames(a),",data=comb)",sep="") #create character function
  model = eval(parse(text=paste(f))) #coverting text(chr) function into R-function
  k=data.frame(surv_pvalue(model)) #Find real p-value
  
  # Store
  pval_6combination[x,] = k$pval
  x=x+1
}
# note :There were 45 warnings (use warnings() to see them)
# note : number of sample after filter = 5 samples
# FURTHER IN HERE WILL CAUSING NO SIGNIFICANT RESULTS
########################### STOP HERE ####################################################################
# SAVE 
KM_surv_2mirna = data.frame(combination_miR_metas_2key,pval_2combination)
KM_surv_3mirna = data.frame(combination_miR_metas_3key,pval_3combination)
KM_surv_4mirna = data.frame(combination_miR_metas_4key,pval_4combination)
KM_surv_5mirna = data.frame(combination_miR_metas_5key,pval_5combination)
KM_surv_6mirna = data.frame(combination_miR_metas_6key,pval_6combination)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Combination 2 miRNA")
addWorksheet(wb, "Combination 3 miRNA")
addWorksheet(wb, "Combination 4 miRNA")
addWorksheet(wb, "Combination 5 miRNA")
addWorksheet(wb, "Combination 6 miRNA")

writeData(wb, 1, KM_surv_2mirna)
writeData(wb, 2, KM_surv_3mirna)
writeData(wb, 3, KM_surv_4mirna)
writeData(wb, 4, KM_surv_5mirna)
writeData(wb, 5, KM_surv_6mirna)

saveWorkbook(wb, "D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\Combination_miRNA_KManalysis.xlsx")

# ERROR : difficult to call column names inside function; rejected non variable name
#TRY-1 : use paste() function to construct
        #use parse() function to coverting "text type of funcgtion" into real function
        
                
### note : MOST OF THESE KM PLOT WERE ANALYZE BASED ON : NORMAL VS PRIMARY AND PRIMARY VS METASTASIS
### SOMETIMES IT DOESNT EFFECT THAT MUCH
### OR MEDIAN IS NOT GOOD MEASUREMENT TO CUT OFF DATA


#QC cek data types
#str(CCLE_value)

############*****# Plot for Publication Purposes ##############################

# MIMAT0000450 &	MIMAT0003218
seeM_pubx= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000450 ==2 & miRNA_Metastasis$x_MIMAT0003218 == 2 ,]  #Get the set fo 
seeM2_pubx= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000450 ==1 & miRNA_Metastasis$x_MIMAT0003218 == 1 ,] #
see_combineM_pubx=rbind(seeM_pubx,seeM2_pubx) #combined base on col names

R_pubx = Surv(see_combineM_pubx$time, see_combineM_pubx$status == 1)
kmfit_R_pubx = survfit(R_pubx ~ see_combineM_pubx$x_MIMAT0000450) # it doesn't matter since the number is same 
model_fit_R_pubx <- survfit(Surv(time, status) ~ x_MIMAT0000450, data = see_combineM_pubx)
autoplot(model_fit_R_pubx) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-149-5p & hsa-miR-92b-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))

# MIMAT0000318 &	MIMAT0003218
seeM_pub= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000318 ==2 & miRNA_Metastasis$x_MIMAT0003218 == 2 ,]  #Get the set fo 
seeM2_pub= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000318 ==1 & miRNA_Metastasis$x_MIMAT0003218 == 1 ,] #
see_combineM_pub=rbind(seeM_pub,seeM2_pub) #combined base on col names

R_pub = Surv(see_combineM_pub$time, see_combineM_pub$status == 1)
kmfit_R_pub = survfit(R_pub ~ see_combineM_pub$x_MIMAT0000318) # it doesn't matter since the number is same 
model_fit_R_pub <- survfit(Surv(time, status) ~ x_MIMAT0000318, data = see_combineM_pub)
autoplot(model_fit_R_pub) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-200b-3p & hsa-miR-92b-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))

# MIMAT0000318 &	MIMAT0000689	& MIMAT0003218
seeM_pub3= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000318 ==2 & miRNA_Metastasis$x_MIMAT0000689 == 2 & miRNA_Metastasis$x_MIMAT0003218 == 2,]  #Get the set fo 
seeM2_pub3= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000318 ==1 & miRNA_Metastasis$x_MIMAT0000689 == 1 & miRNA_Metastasis$x_MIMAT0003218 == 1 ,] #
see_combineM_3mir_pub=rbind(seeM_pub3,seeM2_pub3) #combined base on col names

R2_pub3 = Surv(see_combineM_3mir_pub$time, see_combineM_3mir_pub$status == 1)
kmfit_R2_pub3 = survfit(R2_pub3 ~ see_combineM_3mir_pub$x_MIMAT0000318) # it doesn't matter since the number is same 
model_fit_R2_pub3 <- survfit(Surv(time, status) ~ x_MIMAT0000318, data = see_combineM_3mir_pub)
autoplot(model_fit_R2_pub3) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-200b-3p, hsa-miR-99b-5p, & hsa-miR-92b-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))


# MIMAT0000318 &	MIMAT0000689 &	MIMAT0002171 &	MIMAT0003218
seeM_pub4= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000318 ==2 & miRNA_Metastasis$x_MIMAT0000689 == 2 & miRNA_Metastasis$x_MIMAT0002171 == 2 & miRNA_Metastasis$x_MIMAT0003218 == 2,]  #Get the set fo 
seeM2_pub4= miRNA_Metastasis[miRNA_Metastasis$x_MIMAT0000318 ==1 & miRNA_Metastasis$x_MIMAT0000689 == 1 & miRNA_Metastasis$x_MIMAT0002171 == 1 & miRNA_Metastasis$x_MIMAT0003218 == 1 ,] #
see_combineM_4mir_pub=rbind(seeM_pub4,seeM2_pub4) #combined base on col names

R2_pub4 = Surv(see_combineM_4mir_pub$time, see_combineM_4mir_pub$status == 1)
kmfit_R2_pub4 = survfit(R2_pub4 ~ see_combineM_4mir_pub$x_MIMAT0000318) # it doesn't matter since the number is same 
model_fit_R2_pub4 <- survfit(Surv(time, status) ~ x_MIMAT0000318, data = see_combineM_4mir_pub)
autoplot(model_fit_R2_pub4) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of \n hsa-miR-200b-3p, hsa-miR-99b-5p, hsa-miR-410-3p & hsa-miR-92b-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))


# Adjusting 
#colnames(CCLE_data)[1]=""
#rownames(CCLE_data)=CCLE_data[,1]; 
#CCLE_data <- mCCLE_data[, -c(1)] 

#1. miRNA Cox's Group in Metastasis (10)					miRNA name
#	MIMAT0004494				hsa-miR-21-3p
#	MIMAT0000450				hsa-miR-149-5p
#	MIMAT0002171				hsa-miR-410-3p
#	MIMAT0000280				hsa-miR-223-3p
#	MIMAT0019880				hsa-miR-4746-5p
#	MIMAT0002807				hsa-miR-491-5p
#	MIMAT0000318				hsa-miR-200b-3p
#	MIMAT0003218				hsa-miR-92b-3p
#	MIMAT0000254				hsa-miR-10b-5p
#	MIMAT0000689				hsa-miR-99b-5p
					
#2. miRNA Cox's Group in Primary (4)					miRNA name
#MIMAT0000073				hsa-miR-19a-3p
#MIMAT0000099				hsa-miR-101-3p
#MIMAT0015066				hsa-miR-3065-5p
#MIMAT0014995				hsa-miR-3130-5p

