#Validation by KM plot GSE/GEO Study
#


library(ggfortify) #need to install
library(haven) #For importing STATA files.
library(survival)
library(survminer)

library("readxl")
library(dplyr)
library(stringr)


# Load Information from GEO study
GSE131960=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\GSE131960 MIR SURV PMID 32534961\\GSE131960_miR_surv_data.xlsx")

# MIMAT0000073 OR hsa-miR-10b-5p
# 10-5p
Y = Surv(GSE131960$time_DAYS, GSE131960$status == 1) # 1 = dead/deceased
kmfit_Y = survfit(Y ~ GSE131960$Code_10b) #
model_fit_Y <- survfit(Surv(time_DAYS, status) ~ Code_10b, data = GSE131960)
autoplot(model_fit_Y) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of GSE131960 \n hsa-10b-5p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
  surv_pvalue(model_fit_Y) # p = 0.11

#21-5p in our data is 21-3p [function maybe different]
kmfit_Y1 = survfit(Y ~ GSE131960$Code_21) #
model_fit_Y1 <- survfit(Surv(time_DAYS, status) ~ Code_21, data = GSE131960)
autoplot(model_fit_Y1) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of GSE131960 \n hsa-21-5p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
  surv_pvalue(model_fit_Y1) #p = 0.71
  
#223-3p 
  kmfit_Y2 = survfit(Y ~ GSE131960$Code_223) #
  model_fit_Y2 <- survfit(Surv(time_DAYS, status) ~ Code_223, data = GSE131960)
  autoplot(model_fit_Y2) + 
    labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
         title = "Survival Times Of GSE131960 \n hsa-223-3p KIRC  \n 1=low expression 2=high expression") + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
          axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
          legend.title = element_text(face="bold", size = 10))
  surv_pvalue(model_fit_Y2) #p = 0.66

# 10b & 21
  see= GSE131960[GSE131960$Code_10b ==2 & GSE131960$Code_21 == 2 ,]  #Get the set fo 
  see2= GSE131960[GSE131960$Code_10b ==1 & GSE131960$Code_21 == 1 ,] #
  see_combine=rbind(see,see2) #combined base on col names
  
  Q = Surv(see_combine$time_DAYS, see_combine$status == 1) #give note that we take time=days and status=1 is death
  kmfit_Q = survfit(Q ~ see_combine$Code_10b) # it doesn't matter since the number is same as we analyzed MIMAT0000073
  model_fit_Q <- survfit(Surv(time_DAYS, status) ~ Code_10b, data = see_combine)
  autoplot(model_fit_Q) + 
    labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
         title = "Survival Times Of GSE131960 \n hsa-miR-10b-5p & hsa-miR-21-5p KIRC  \n 1=low expression 2=high expression") + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
          axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
          legend.title = element_text(face="bold", size = 10))
  surv_pvalue(model_fit_Q) #p = 0.25

# 10b & 223
  see_1= GSE131960[GSE131960$Code_10b ==2 & GSE131960$Code_223 == 2 ,]  #Get the set fo 
  see2_1= GSE131960[GSE131960$Code_10b ==1 & GSE131960$Code_223 == 1 ,] #
  see_combine_1=rbind(see_1,see2_1) #combined base on col names
  
  Q1 = Surv(see_combine_1$time_DAYS, see_combine_1$status == 1) #give note that we take time=days and status=1 is death
  kmfit_Q1 = survfit(Q1 ~ see_combine_1$Code_10b) # it doesn't matter since the number is same as we analyzed MIMAT0000073
  model_fit_Q1 <- survfit(Surv(time_DAYS, status) ~ Code_10b, data = see_combine_1)
  autoplot(model_fit_Q1) + 
    labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
         title = "Survival Times Of GSE131960 \n hsa-miR-10b-5p & hsa-miR-223-3p KIRC  \n 1=low expression 2=high expression") + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
          axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
          legend.title = element_text(face="bold", size = 10))
  surv_pvalue(model_fit_Q1) #p = 0.28
  
  
  
  
  
  
# GSE131959
# Load Information from GEO study [this study have missing data : manual impute by mean value]
library(gtools)
GSE131959=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\GSE131959 marray MIR\\GSE131959_marray_filtered.xlsx")

#Testing Primary
GSE131959_P=data.frame(GSE131959[,24:25])

Y = Surv(GSE131959$time_DAYS, GSE131959$status == 1) # 1 = dead/deceased
#miR-19a
model_fit_19a <- survfit(Surv(time_DAYS, status) ~ hsa_miR_19a_3p_exp, data = GSE131959)
p_val_19a=surv_pvalue(model_fit_19a) 

#miR-101
model_fit_101 <- survfit(Surv(time_DAYS, status) ~ hsa_miR_101_3p_exp, data = GSE131959)
p_val_101=surv_pvalue(model_fit_101) 

#Testing Metastasis
GSE131959_M=data.frame(GSE131959[,15:23])
combination_miR_metas_2key=combinations(n=9,r=2,v=names(GSE131959_M),set=TRUE,repeats.allowed=F) # n=possible item, r=2 combination key without replacement
combination_miR_metas_3key=combinations(n=9,r=3,v=names(GSE131959_M),set=TRUE, repeats.allowed=F) #set= duplicate will be removed
combination_miR_metas_4key=combinations(n=9,r=4,v=names(GSE131959_M),set=TRUE, repeats.allowed=F)

# Render 2 combination
pval_2combination=data.frame(Pval_col = rep(NA, length(combination_miR_metas_2key[,1]))) #Store in here
x=0
for(x in 1:length(combination_miR_metas_2key[,1])){ 
  a=combination_miR_metas_2key[x,1] #grab first miRNA
  b=combination_miR_metas_2key[x,2] #grab second miRNA
  #grab colomn based on miRNA name combination
  a=select(GSE131959,contains(a))
  b=select(GSE131959,contains(b))
  up= GSE131959[a == 2 & b == 2,]  #Get the set fo 
  down= GSE131959[a ==1 & b == 1,] #
  comb=rbind(up,down) #combined base on col names
  
  # Main Function
  f=paste("survfit(Surv(time_DAYS, status) ~", colnames(a),",data=comb)",sep="") #create character function
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
  a=select(GSE131959,contains(a))
  b=select(GSE131959,contains(b))
  c=select(GSE131959,contains(c))
  up= GSE131959[a == 2 & b == 2 & c==2,]  #Get the set fo 
  down= GSE131959[a ==1 & b == 1 & c==1,] #
  comb=rbind(up,down) #combined base on col names
  
  # Main Function
  f=paste("survfit(Surv(time_DAYS, status) ~", colnames(a),",data=comb)",sep="") #create character function
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
  a=select(GSE131959,contains(a))
  b=select(GSE131959,contains(b))
  c=select(GSE131959,contains(c))
  d=select(GSE131959,contains(d))
  up= GSE131959[a == 2 & b == 2 & c==2 & d==2,]  #Get the set fo 
  down= GSE131959[a ==1 & b == 1 & c==1 & d==1,] #
  comb=rbind(up,down) #combined base on col names
  
  # Main Function
  f=paste("survfit(Surv(time_DAYS, status) ~", colnames(a),",data=comb)",sep="") #create character function
  model = eval(parse(text=paste(f))) #coverting text(chr) function into R-function
  k=data.frame(surv_pvalue(model)) #Find real p-value
  
  # Store
  pval_4combination[x,] = k$pval
  x=x+1
}

## Combination of Compiled File
two_combination_results = data.frame(combination_miR_metas_2key,pval_2combination)
three_combination_results = data.frame(combination_miR_metas_3key,pval_3combination)
four_combination_results = data.frame(combination_miR_metas_4key,pval_4combination)

filter_1=four_combination_results$X1["hsa_miR_200b_3p_exp",]


## fOR rendering purposes ##
#miR-149 & 92b
seeM_149_92b= GSE131959[GSE131959$hsa_miR_149_5p_exp ==2 & GSE131959$hsa_miR_92b_3p_exp == 2 ,]  #Get the set fo 
seeM2_149_92b= GSE131959[GSE131959$hsa_miR_149_5p_exp ==1 & GSE131959$hsa_miR_92b_3p_exp == 1 ,] #
see_combineM_149_92b=rbind(seeM_149_92b,seeM2_149_92b) #combined base on col names

model_fit_149_92b <- survfit(Surv(time_DAYS, status) ~ hsa_miR_149_5p_exp, data = see_combineM_149_92b)
autoplot(model_fit_149_92b ) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of GSE131959 \n hsa-miR-149-5p & hsa-miR-92-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))


#miR-410 & 92b
seeM_410_92b= GSE131959[GSE131959$hsa_miR_410_3p_exp ==2 & GSE131959$hsa_miR_92b_3p_exp == 2,]  #Get the set fo 
seeM2_410_92b= GSE131959[GSE131959$hsa_miR_410_3p_exp ==1 & GSE131959$hsa_miR_92b_3p_exp == 1,] #
see_combineM_410_92b=rbind(seeM_410_92b,seeM2_410_92b) #combined base on col names

model_fit_410_92b <- survfit(Surv(time_DAYS, status) ~ hsa_miR_410_3p_exp, data = see_combineM_410_92b)
autoplot(model_fit_410_92b  ) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of GSE131959 \n hsa-miR-410-3p & hsa-miR-92-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))

#miR-491 & 92b
seeM_491_92b= GSE131959[GSE131959$hsa_miR_491_5p_exp ==2 & GSE131959$hsa_miR_92b_3p_exp == 2,]  #Get the set fo 
seeM2_491_92b= GSE131959[GSE131959$hsa_miR_491_5p_exp ==1 & GSE131959$hsa_miR_92b_3p_exp == 1,] #
see_combineM_491_92b=rbind(seeM_491_92b,seeM2_491_92b) #combined base on col names

model_fit_491_92b <- survfit(Surv(time_DAYS, status) ~ hsa_miR_491_5p_exp, data = see_combineM_491_92b)
autoplot(model_fit_491_92b  ) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of GSE131959 \n hsa-miR-491-5p & hsa-miR-92-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))


#miR-200b & 99b & 92b
seeM_200b_99b_92b= GSE131959[GSE131959$hsa_miR_200b_3p_exp ==2 & GSE131959$hsa_miR_92b_3p_exp == 2 & GSE131959$hsa_miR_99b_5p_exp ==2,]  #Get the set fo 
seeM2_200b_99b_92b= GSE131959[GSE131959$hsa_miR_200b_3p_exp ==1 & GSE131959$hsa_miR_92b_3p_exp == 1 & GSE131959$hsa_miR_99b_5p_exp ==1,] #
see_combineM_200b_99b_92b=rbind(seeM_200b_99b_92b,seeM2_200b_99b_92b) #combined base on col names

model_fit_200b_99b_92b <- survfit(Surv(time_DAYS, status) ~ hsa_miR_200b_3p_exp, data = see_combineM_200b_99b_92b)
autoplot(model_fit_200b_99b_92b  ) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of GSE131959 \n hsa-miR-200b-3p, hsa-miR-99b-5p & hsa-miR-92-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))
#small sample size unable to plot KM correctly

#miR-99b & 410 & 92b
seeM_410_99b_92b= GSE131959[GSE131959$hsa_miR_410_3p_exp ==2 & GSE131959$hsa_miR_92b_3p_exp == 2 & GSE131959$hsa_miR_99b_5p_exp ==2,]  #Get the set fo 
seeM2_410_99b_92b= GSE131959[GSE131959$hsa_miR_410_3p_exp ==1 & GSE131959$hsa_miR_92b_3p_exp == 1 & GSE131959$hsa_miR_99b_5p_exp ==1,] #
see_combineM_410_99b_92b=rbind(seeM_410_99b_92b,seeM2_410_99b_92b) #combined base on col names

model_fit_410_99b_92b <- survfit(Surv(time_DAYS, status) ~ hsa_miR_410_3p_exp, data = see_combineM_410_99b_92b)
autoplot(model_fit_410_99b_92b  ) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times Of GSE131959 \n hsa-miR-410-3p, hsa-miR-99b-5p & hsa-miR-92-3p KIRC  \n 1=low expression 2=high expression") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
        axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
        legend.title = element_text(face="bold", size = 10))








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
  