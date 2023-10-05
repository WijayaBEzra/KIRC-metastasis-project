# Source : https://stats.stackexchange.com/questions/506941/cox-regression-with-multiple-factors-r
# Source : https://www.r-bloggers.com/2016/12/cox-proportional-hazards-model/

library("survival")
library("survminer")
library("readxl")
library(dplyr)
library(stringr)
data("lung")
head(lung)

options(scipen = 999) #avoid unit as exponential number as 1e+100

#Load Data
Metastasis=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Metastasis.xlsx")
Primary=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Primary.xlsx")
Normal=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Normal.xlsx")

d_Metastasis=t(Metastasis) #transpose matrix
d_Primary=t(Primary) #transpose matrix
d_Normal=t(Normal) #transpose matrix

colnames(d_Metastasis)=d_Metastasis[1,] ; d_Metastasis <- d_Metastasis[-1, ] # first row into col names
colnames(d_Primary)=d_Primary[1,] ; d_Primary <- d_Primary[-1, ] # first row into col names
colnames(d_Normal)=d_Normal[1,] ; d_Normal <- d_Normal[-1, ] # first row into col names

myFile2 = read.delim(file.choose(), header = TRUE) #warning: check your tab, choose your clinical file manually in your file directory
#[!] UPgrade : API address can be use
# Instruction : https://www.r-bloggers.com/2016/12/cox-proportional-hazards-model/
time=data.frame(myFile2$days_to_death,myFile2$days_to_last_followup) #Combined time days_to_death & Days to follow up by choosing the bigger dayss
time[is.na(time)] <- 0 #NA set as zero
time=pmax(time$myFile2.days_to_death, time$myFile2.days_to_last_followup)
time=data.frame(time)
patient_id=str_replace_all(myFile2$sampleID,"-",".")
vital_stat=str_replace_all(myFile2$vital_status,c(LIVING="1", DECEASED = "2"))
vital_stat=data.frame(vital_stat)

Censored_=data.frame(patient_id,time,vital_stat) #final clinical survival time
colnames(Censored_)=c("patient_ID","time","")

### Filtered Patient_ID with miRNA_HIseq's Patient's ID ###
# Set Cancer data in Clinical Censored variable
Censored_c=Censored_[grep(".01",Censored_$patient_ID ),] #Eliminate Normal Sample by choosing name with .01 pattern or primary

#Find Metastasis Patient ID
metas_patient=(paste(row.names(d_Metastasis), collapse="|"))
Metastasis_test2=Censored_c[grep(metas_patient,Censored_c$patient_ID),] #Reverse we chech only in clinical data prior to metastasis sample

# QC double check its similarities
#test=data.frame(Metastasis_test2$patient_ID, rownames(d_Metastasis))
#library(writexl)
#write_xlsx(test,"D:\\Dissertation 4.08.2022\\test.xlsx") # note : correct Result!

#Combined into one frame
row.names(Metastasis_test2)=Metastasis_test2$patient_ID
Metastasis_test2=Metastasis_test2[,-c(1)]
d_Metastasis_x=merge(d_Metastasis,Metastasis_test2, by=0) #combined based on rownames

#QC
#see=d_Metastasis_x[,630:633] #check either it already merged correctly
#write_xlsx(d_Metastasis_x,"D:\\Dissertation 4.08.2022\\checkcombination.xlsx") # note: correct it matched

#Find Primary Patient ID
Primary_patient=(paste(row.names(d_Primary), collapse="|"))
Primary_test2=Censored_c[grep(Primary_patient,Censored_c$patient_ID),]
row.names(Primary_test2)=Primary_test2$patient_ID
Primary_test2=Primary_test2[,-c(1)]
d_Primary_x=merge(d_Primary,Primary_test2, by=0) #combined based on rownames

#Prepare data structure for analysis set as NUMERIC from CHAR Problems
row.names(d_Primary_x)=d_Primary_x[,1]; d_Primary_x=d_Primary_x[,-c(1)]
row.names(d_Metastasis_x)=d_Metastasis_x[,1]; d_Metastasis_x=d_Metastasis_x[,-c(1)]

d_Primary_xm=as.numeric(d_Primary_x) # Error: 'list' object cannot be coerced to type 'double'
d_Metastasis_xm=as.numeric(d_Metastasis_x) # Error: 'list' object cannot be coerced to type 'double'

d_Primary_xm <- sapply(d_Primary_x,as.numeric) #Error :'data' must be a data.frame, not a matrix or an array
d_Metastasis_xm <- sapply(d_Metastasis_x,as.numeric) #Error :'data' must be a data.frame, not a matrix or an array

d_Primary_xm <- data.frame(d_Primary_xm) # Return to DF worked!
d_Metastasis_xm <- data.frame(d_Metastasis_xm) # Return to DF worked!

rm(d_Metastasis_xm)
rm(d_Primary_xm)



#test cox function in PRIMARY miRNA-Patients (note: rowname not established)
# WARNING : value still in RPM 
covariates <- colnames(d_Primary_x[,2:630]) #grab only miRNA names to be analyzed
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = d_Primary_xm)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
summary(res)
HR_Result_Primary_miR=data.frame(rownames(res),res)

#HR in Metastasis each miRNA
univ_models_m <- lapply(univ_formulas, function(x){coxph(x, data = d_Metastasis_xm)})
univ_results_m <- lapply(univ_models_m,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res2<-c(beta, HR, wald.test, p.value)
                         names(res2)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res2)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

res2 <- t(as.data.frame(univ_results_m, check.names = FALSE))
as.data.frame(res2)
HR_Result_Metastasis_miR=data.frame(rownames(res2),res2) 

#Only Significant
HR_Metastasis_filtered <- filter(HR_Result_Metastasis_miR,p.value < 0.05)
HR_Primary_filtered <- filter(HR_Result_Primary_miR,p.value < 0.05)

#Save Results (no Filtered p-val)
write_xlsx(HR_Result_Metastasis_miR,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\HR_Result_Metastasis_miR.xlsx")
write_xlsx(HR_Result_Primary_miR,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\HR_Result_Primary_miR.xlsx")

write_xlsx(HR_Metastasis_filtered,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\HR_Metastasis_filtered.xlsx")
write_xlsx(HR_Primary_filtered,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\HR_Primary_filtered.xlsx")
#NOTE : above results based on Individual features of MIRNA hazard to survival

#Multiple CLINICAL INFO cox
# For application refer to link  : https://www.tutorialspoint.com/r/r_multiple_regression.htm
# Another multivariate link : https://cran.r-project.org/web/packages/survivalAnalysis/vignettes/multivariate.html
# Glimpse of Model Prediction : https://www.geeksforgeeks.org/logistic-regression-in-r-programming/

###----------------First Run ----------------------------------------------------###
#use filtered miRNA as New Covarities
Sig_miRNAp <- HR_Primary_filtered$rownames.res.
Sig_miRNAm <- HR_Metastasis_filtered$rownames.res2.

# Try use All miRNA sig Cox (Metastasis)
ALL_miRm = paste(Sig_miRNAm[1:27],"+") #EDITED in NOTEPAD++ remove the(") symbol
all_miRNA_cox <- coxph(Surv(time, status) ~ 
                         MIMAT0004819 + MIMAT0004494 + MIMAT0022287 + MIMAT0000455 + MIMAT0000450 + MIMAT0019729 + MIMAT0002171 +
                         MIMAT0000280 + MIMAT0019880 + MIMAT0003880 + MIMAT0002807 + MIMAT0000318 + MIMAT0004703 + MIMAT0000076 +
                         MIMAT0000094 + MIMAT0004571 + MIMAT0004985 + MIMAT0003218 + MIMAT0000254 + MIMAT0004556 + MIMAT0000689 +
                         MIMAT0000686 + MIMAT0006790 + MIMAT0000690 + MIMAT0009199 + MIMAT0000100 + MIMAT0004699          
                         , data =  d_Metastasis_xm) # Manually inputted miRNA
summary(all_miRNA_cox)
test=summary(all_miRNA_cox)
Multicox_Metas_results=data.frame(test$coefficients[,1:2],test$conf.int[,3:4],test$coefficients[,3:5])

# Try use All miRNA sig Cox (Primary)
ALL_miRp = paste(Sig_miRNAp[1:length(Sig_miRNAp)],"+") #EDITED in NOTEPAD++ remove the(") symbol
all_miRNA_cox_p=coxph(Surv(time, status) ~ 
       MIMAT0004766 + MIMAT0004768 + MIMAT0004808 + MIMAT0004482 + MIMAT0003284 + MIMAT0000448 + MIMAT0000449 +
       MIMAT0000445 + MIMAT0000446 + MIMAT0000447 + MIMAT0004813 + MIMAT0004495 + MIMAT0026472 + MIMAT0019208 +
       MIMAT0002173 + MIMAT0002174 + MIMAT0022833 + MIMAT0000733 + MIMAT0000280 + MIMAT0000281 + MIMAT0002820 +
       MIMAT0000462 + MIMAT0003163 + MIMAT0003885 + MIMAT0004604 + MIMAT0004601 + MIMAT0004518 + MIMAT0004515 +
       MIMAT0003393 + MIMAT0000265 + MIMAT0000263 + MIMAT0002808 + MIMAT0002809 + MIMAT0000064 + MIMAT0000063 +
       MIMAT0004501 + MIMAT0004503 + MIMAT0000279 + MIMAT0000278 + MIMAT0000070 + MIMAT0000073 + MIMAT0000072 +
       MIMAT0000074 + MIMAT0000076 + MIMAT0000099 + MIMAT0000098 + MIMAT0004570 + MIMAT0004571 + MIMAT0000646 +
       MIMAT0000427 + MIMAT0000425 + MIMAT0000421 + MIMAT0003218 + MIMAT0004563 + MIMAT0004569 + MIMAT0022693 +
       MIMAT0000254 + MIMAT0000250 + MIMAT0004776 + MIMAT0004556 + MIMAT0004552 + MIMAT0004551 + MIMAT0000688 +
       MIMAT0015066 + MIMAT0001636 + MIMAT0001638 + MIMAT0001639 + MIMAT0000693 + MIMAT0027513 + MIMAT0000232 +
       MIMAT0026480 + MIMAT0000754 + MIMAT0000758 + MIMAT0004921 + MIMAT0004926 + MIMAT0003258 + MIMAT0009198 +
       MIMAT0004927 + MIMAT0000100 + MIMAT0000103 + MIMAT0003888 + MIMAT0004680 + MIMAT0004688 + MIMAT0014995 +
       MIMAT0004674 + MIMAT0015020 + MIMAT0004589 + MIMAT0004588
       , data =  d_Primary_xm)

summary(all_miRNA_cox_p)
test2=summary(all_miRNA_cox_p)
Multicox_Primary_results=data.frame(test2$coefficients[,1:2],test2$conf.int[,3:4],test2$coefficients[,3:5])

options(scipen = 9) #hold only 9 sig figure
# NOTE: WALD TEST :used to compare models on best fit criteria in case of logistic regression

#QC use 2 miRNA as multivariate analysis
#try = coxph(Surv(time, status) ~  ALL_miR , data =  d_Primary_xm) 
#summary(try)
#Note : the p-value will changed each time we create different miRNA combination
#NOTE : exp(coef) is Hazard Ratio

###-----------------Second RUN---------------------------###
# In This Run, we filtered miRNA again
#use filtered miRNA as New Covarities
Multicox_m_filtered <- filter(Multicox_Metas_results,Multicox_Metas_results$Pr...z.. < 0.05)
Multicox_p_filtered <- filter(Multicox_Primary_results,Multicox_Primary_results$Pr...z.. < 0.05)
All_miR_m=paste(row.names(Multicox_m_filtered),"+")
All_miR_p=paste(row.names(Multicox_p_filtered),"+")

all_miRNA_coxm_2nd <- coxph(Surv(time, status) ~ 
                         MIMAT0004494 + MIMAT0022287 + MIMAT0000450 + MIMAT0002171 + MIMAT0000280 + MIMAT0019880 + MIMAT0002807 +
                         MIMAT0000318 + MIMAT0003218 + MIMAT0000254 + MIMAT0000689 + MIMAT0000686 + MIMAT0000690 + MIMAT0009199 
                       , data =  d_Metastasis_xm) # Manually inputted miRNA
summary(all_miRNA_coxm_2nd)
test3=summary(all_miRNA_coxm_2nd)
Multicox_Metas_results_2nd=data.frame(test3$coefficients[,1:2],test3$conf.int[,3:4],test3$coefficients[,3:5])

#Primary 2nd round
all_miRNA_coxp_2nd=coxph(Surv(time, status) ~ 
                        MIMAT0004766 + MIMAT0004808 + MIMAT0004482 + MIMAT0000448 + MIMAT0000449 + MIMAT0000446 + MIMAT0000447 +
                        MIMAT0004813 + MIMAT0004495 + MIMAT0026472 + MIMAT0019208 + MIMAT0002173 + MIMAT0022833 + MIMAT0000733 +
                        MIMAT0000280 + MIMAT0000281 + MIMAT0002820 + MIMAT0000462 + MIMAT0003163 + MIMAT0003885 + MIMAT0004604 +
                        MIMAT0004601 + MIMAT0004515 + MIMAT0003393 + MIMAT0000265 + MIMAT0000263 + MIMAT0002808 + MIMAT0002809 +
                        MIMAT0000064 + MIMAT0000063 + MIMAT0004501 + MIMAT0004503 + MIMAT0000279 + MIMAT0000278 + MIMAT0000070 +
                        MIMAT0000073 + MIMAT0000072 + MIMAT0000074 + MIMAT0000076 + MIMAT0000099 + MIMAT0000098 + MIMAT0004570 +
                        MIMAT0004571 + MIMAT0000646 + MIMAT0000427 + MIMAT0000425 + MIMAT0003218 + MIMAT0004563 + MIMAT0004569 +
                        MIMAT0022693 + MIMAT0000254 + MIMAT0000250 + MIMAT0004776 + MIMAT0004556 + MIMAT0004552 + MIMAT0004551 +
                        MIMAT0000688 + MIMAT0015066 + MIMAT0001636 + MIMAT0001638 + MIMAT0001639 + MIMAT0000693 + MIMAT0000232 +
                        MIMAT0026480 + MIMAT0000754 + MIMAT0004921 + MIMAT0004926 + MIMAT0003258 + MIMAT0009198 + MIMAT0004927 +
                        MIMAT0000100 + MIMAT0000103 + MIMAT0004680 + MIMAT0004688 + MIMAT0014995 + MIMAT0004674 + MIMAT0015020 +
                        MIMAT0004589 + MIMAT0004588
                      , data =  d_Primary_xm)

summary(all_miRNA_coxp_2nd)
test4=summary(all_miRNA_coxp_2nd)
Multicox_Primary_results_2nd=data.frame(test4$coefficients[,1:2],test4$conf.int[,3:4],test4$coefficients[,3:5])

###-----------------Third RUN---------------------------###
#use filtered miRNA as New Covarities
Multicox_m_filtered3 <- filter(Multicox_Metas_results_2nd,Multicox_Metas_results_2nd$Pr...z.. < 0.05)
Multicox_p_filtered3 <- filter(Multicox_Primary_results_2nd,Multicox_Primary_results_2nd$Pr...z.. < 0.05)
All_miR_m3=paste(row.names(Multicox_m_filtered3),"+")
All_miR_p3=paste(row.names(Multicox_p_filtered3),"+")

all_miRNA_coxm_3rd <- coxph(Surv(time, status) ~ 
                              MIMAT0004494 + MIMAT0000450 + MIMAT0002171 + MIMAT0000280 + MIMAT0019880 + MIMAT0002807 + MIMAT0000318 +
                              MIMAT0003218 + MIMAT0000254 + MIMAT0000689 + MIMAT0009199
                            , data =  d_Metastasis_xm) # Manually inputted miRNA
summary(all_miRNA_coxm_3rd)
test5=summary(all_miRNA_coxm_3rd)
Multicox_Metas_results_3rd=data.frame(test5$coefficients[,1:2],test5$conf.int[,3:4],test5$coefficients[,3:5])

#Primary 3rd round
all_miRNA_coxp_3rd=coxph(Surv(time, status) ~ 
                           MIMAT0004766 + MIMAT0004808 + MIMAT0004482 + MIMAT0000448 + MIMAT0000449 + MIMAT0000446 + MIMAT0000447 +
                           MIMAT0004813 + MIMAT0004495 + MIMAT0026472 + MIMAT0019208 + MIMAT0002173 + MIMAT0022833 + MIMAT0000733 +
                           MIMAT0000281 + MIMAT0002820 + MIMAT0000462 + MIMAT0003163 + MIMAT0003885 + MIMAT0004604 + MIMAT0004601 +
                           MIMAT0004515 + MIMAT0000265 + MIMAT0000263 + MIMAT0002808 + MIMAT0002809 + MIMAT0000064 + MIMAT0000063 +
                           MIMAT0004501 + MIMAT0004503 + MIMAT0000279 + MIMAT0000278 + MIMAT0000070 + MIMAT0000073 + MIMAT0000072 +
                           MIMAT0000074 + MIMAT0000076 + MIMAT0000099 + MIMAT0000098 + MIMAT0004570 + MIMAT0004571 + MIMAT0000646 +
                           MIMAT0000425 + MIMAT0004563 + MIMAT0004569 + MIMAT0022693 + MIMAT0000254 + MIMAT0000250 + MIMAT0004552 +
                           MIMAT0004551 + MIMAT0000688 + MIMAT0015066 + MIMAT0001636 + MIMAT0001638 + MIMAT0001639 + MIMAT0000232 +
                           MIMAT0026480 + MIMAT0000754 + MIMAT0004921 + MIMAT0004926 + MIMAT0003258 + MIMAT0009198 + MIMAT0004927 +
                           MIMAT0000100 + MIMAT0000103 + MIMAT0004680 + MIMAT0004688 + MIMAT0014995 + MIMAT0004674 + MIMAT0015020 +
                           MIMAT0004589 + MIMAT0004588
                         , data =  d_Primary_xm)

summary(all_miRNA_coxp_3rd)
test6=summary(all_miRNA_coxp_3rd)
Multicox_Primary_results_3rd=data.frame(test6$coefficients[,1:2],test6$conf.int[,3:4],test6$coefficients[,3:5])

###-----------------Fourth RUN---------------------------###
Multicox_m_filtered4 <- filter(Multicox_Metas_results_3rd,Multicox_Metas_results_3rd$Pr...z.. < 0.05)
Multicox_p_filtered4 <- filter(Multicox_Primary_results_3rd,Multicox_Primary_results_3rd$Pr...z.. < 0.05)
All_miR_m4=paste(row.names(Multicox_m_filtered4),"+")
All_miR_p4=paste(row.names(Multicox_p_filtered4),"+")

all_miRNA_coxm_4th <- coxph(Surv(time, status) ~ 
                              MIMAT0004494 + MIMAT0000450 + MIMAT0002171 + MIMAT0000280 + MIMAT0019880 + MIMAT0002807 + MIMAT0000318 +
                              MIMAT0003218 + MIMAT0000254 + MIMAT0000689
                            , data =  d_Metastasis_xm) # Manually inputted miRNA
summary(all_miRNA_coxm_4th)
test7=summary(all_miRNA_coxm_4th)
Multicox_Metas_results_4th=data.frame(test7$coefficients[,1:2],test7$conf.int[,3:4],test7$coefficients[,3:5])

all_miRNA_coxp_4th=coxph(Surv(time, status) ~ 
                           MIMAT0004766 + MIMAT0004808 + MIMAT0004482 + MIMAT0000448 + MIMAT0000449 + MIMAT0000446 + MIMAT0000447 +
                           MIMAT0004495 + MIMAT0026472 + MIMAT0019208 + MIMAT0002173 + MIMAT0022833 + MIMAT0000733 + MIMAT0000281 +
                           MIMAT0000462 + MIMAT0003885 + MIMAT0004604 + MIMAT0004601 + MIMAT0004515 + MIMAT0000265 + MIMAT0000263 +
                           MIMAT0002808 + MIMAT0002809 + MIMAT0000064 + MIMAT0000063 + MIMAT0004501 + MIMAT0004503 + MIMAT0000279 +
                           MIMAT0000278 + MIMAT0000070 + MIMAT0000073 + MIMAT0000072 + MIMAT0000074 + MIMAT0000076 + MIMAT0000099 +
                           MIMAT0004570 + MIMAT0004571 + MIMAT0000646 + MIMAT0004563 + MIMAT0004569 + MIMAT0022693 + MIMAT0000254 +
                           MIMAT0000250 + MIMAT0004552 + MIMAT0004551 + MIMAT0000688 + MIMAT0015066 + MIMAT0001636 + MIMAT0001638 +
                           MIMAT0001639 + MIMAT0000232 + MIMAT0026480 + MIMAT0000754 + MIMAT0004921 + MIMAT0004926 + MIMAT0003258 +
                           MIMAT0009198 + MIMAT0004927 + MIMAT0000103 + MIMAT0004680 + MIMAT0004688 + MIMAT0014995 + MIMAT0004674 +
                           MIMAT0015020 + MIMAT0004589 + MIMAT0004588
                         , data =  d_Primary_xm)

summary(all_miRNA_coxp_4th)
test8=summary(all_miRNA_coxp_4th)
Multicox_Primary_results_4th=data.frame(test8$coefficients[,1:2],test8$conf.int[,3:4],test8$coefficients[,3:5])

###-----------------Fifth RUN---------------------------###
Multicox_m_filtered5 <- filter(Multicox_Metas_results_4th,Multicox_Metas_results_4th$Pr...z.. < 0.05)
Multicox_p_filtered5 <- filter(Multicox_Primary_results_4th,Multicox_Primary_results_4th$Pr...z.. < 0.05)
All_miR_m5=paste(row.names(Multicox_m_filtered5),"+")
All_miR_p5=paste(row.names(Multicox_p_filtered5),"+")

all_miRNA_coxm_5th <- coxph(Surv(time, status) ~ 
                              MIMAT0004494 + MIMAT0000450 + MIMAT0002171 + MIMAT0000280 + MIMAT0019880 + MIMAT0002807 + MIMAT0000318 +
                              MIMAT0003218 + MIMAT0000254 + MIMAT0000689
                            , data =  d_Metastasis_xm) # Manually inputted miRNA
summary(all_miRNA_coxm_5th)
test9=summary(all_miRNA_coxm_5th)
Multicox_Metas_results_5th=data.frame(test9$coefficients[,1:2],test9$conf.int[,3:4],test9$coefficients[,3:5])
####NOTIFICATION!!!!!########## METASTASIS ALREADY STABLE ,see : filtered miRNA cox running 4 =10 and rungging 5=10 and same results

all_miRNA_coxp_5th=coxph(Surv(time, status) ~ 
                           MIMAT0004766 + MIMAT0000448 + MIMAT0000449 + MIMAT0000446 + MIMAT0000447 + MIMAT0004495 + MIMAT0002173 +
                           MIMAT0022833 + MIMAT0000733 + MIMAT0000281 + MIMAT0000462 + MIMAT0003885 + MIMAT0004604 + MIMAT0004515 +
                           MIMAT0000265 + MIMAT0002808 + MIMAT0002809 + MIMAT0000278 + MIMAT0000073 + MIMAT0000099 + MIMAT0004570 +
                           MIMAT0004571 + MIMAT0000646 + MIMAT0004563 + MIMAT0022693 + MIMAT0000688 + MIMAT0015066 + MIMAT0001636 +
                           MIMAT0001638 + MIMAT0000232 + MIMAT0004926 + MIMAT0009198 + MIMAT0004680 + MIMAT0014995 + MIMAT0004674 +
                           MIMAT0004589
                         , data =  d_Primary_xm)

summary(all_miRNA_coxp_5th)
test10=summary(all_miRNA_coxp_5th)
Multicox_Primary_results_5th=data.frame(test10$coefficients[,1:2],test10$conf.int[,3:4],test10$coefficients[,3:5])

###-----------------Sixth RUN {Only Primary} ---------------------------###
Multicox_m_filtered6 <- filter(Multicox_Metas_results_5th,Multicox_Metas_results_5th$Pr...z.. < 0.05) #same as above
Multicox_p_filtered6 <- filter(Multicox_Primary_results_5th,Multicox_Primary_results_5th$Pr...z.. < 0.05)
All_miR_p6=paste(row.names(Multicox_p_filtered6),"+")

all_miRNA_coxp_6th=coxph(Surv(time, status) ~ 
                           MIMAT0000449 + MIMAT0002173 + MIMAT0004515 + MIMAT0000265 + MIMAT0002808 + MIMAT0000073 + MIMAT0000099 +
                           MIMAT0004571 + MIMAT0000646 + MIMAT0022693 + MIMAT0015066 + MIMAT0001638 + MIMAT0004926 + MIMAT0009198 +
                           MIMAT0014995 + MIMAT0004674
                         ,data =  d_Primary_xm)

summary(all_miRNA_coxp_6th)
test11=summary(all_miRNA_coxp_6th)
Multicox_Primary_results_6th=data.frame(test11$coefficients[,1:2],test11$conf.int[,3:4],test11$coefficients[,3:5])

###-----------------Seventh RUN {Only Primary} ---------------------------###
Multicox_p_filtered7 <- filter(Multicox_Primary_results_6th,Multicox_Primary_results_6th$Pr...z.. < 0.05)
All_miR_p7=paste(row.names(Multicox_p_filtered7),"+")

all_miRNA_coxp_7th=coxph(Surv(time, status) ~ 
                           MIMAT0002173 + MIMAT0004515 + MIMAT0000073 + MIMAT0000099 + MIMAT0022693 + MIMAT0015066 + MIMAT0014995 +
                           MIMAT0004674
                         ,data =  d_Primary_xm)

summary(all_miRNA_coxp_7th)
test12=summary(all_miRNA_coxp_7th)
Multicox_Primary_results_7th=data.frame(test12$coefficients[,1:2],test12$conf.int[,3:4],test12$coefficients[,3:5])

###-----------------Eighth RUN {Only Primary} ---------------------------###
Multicox_p_filtered8 <- filter(Multicox_Primary_results_7th,Multicox_Primary_results_7th$Pr...z.. < 0.05)
All_miR_p8=paste(row.names(Multicox_p_filtered8),"+")

all_miRNA_coxp_8th=coxph(Surv(time, status) ~ 
                           MIMAT0002173 + MIMAT0000073 + MIMAT0000099 + MIMAT0022693 + MIMAT0015066 + MIMAT0014995 + MIMAT0004674
                         ,data =  d_Primary_xm)

summary(all_miRNA_coxp_8th)
test13=summary(all_miRNA_coxp_8th)
Multicox_Primary_results_8th=data.frame(test13$coefficients[,1:2],test13$conf.int[,3:4],test13$coefficients[,3:5])

###-----------------Ninth RUN {Only Primary} ---------------------------###
Multicox_p_filtered9 <- filter(Multicox_Primary_results_8th,Multicox_Primary_results_8th$Pr...z.. < 0.05)
All_miR_p9=paste(row.names(Multicox_p_filtered9),"+")

all_miRNA_coxp_9th=coxph(Surv(time, status) ~ 
                           MIMAT0000073 + MIMAT0000099 + MIMAT0022693 + MIMAT0015066 + MIMAT0014995
                         ,data =  d_Primary_xm)

summary(all_miRNA_coxp_9th)
test14=summary(all_miRNA_coxp_9th)
Multicox_Primary_results_9th=data.frame(test14$coefficients[,1:2],test14$conf.int[,3:4],test14$coefficients[,3:5])

###-----------------Tenth RUN {Only Primary} ---------------------------###
Multicox_p_filtered10 <- filter(Multicox_Primary_results_9th,Multicox_Primary_results_9th$Pr...z.. < 0.05)
All_miR_p10=paste(row.names(Multicox_p_filtered10),"+")

all_miRNA_coxp_10th=coxph(Surv(time, status) ~ 
                            MIMAT0000073 + MIMAT0000099 + MIMAT0015066 + MIMAT0014995
                          ,data =  d_Primary_xm)

summary(all_miRNA_coxp_10th)
test15=summary(all_miRNA_coxp_10th)
Multicox_Primary_results_10th=data.frame(test15$coefficients[,1:2],test15$conf.int[,3:4],test15$coefficients[,3:5])

####END OF RUN COX RESULTS######

#SAVE AND STORED INFORMATIONS!!!!#
library(writexl)
Multicox_Primary_results_10th_=data.frame(row.names(Multicox_Primary_results_10th),Multicox_Primary_results_10th)
Multicox_Metas_results_5th_=data.frame(row.names(Multicox_Metas_results_5th),Multicox_Metas_results_5th)
ggforest(all_miRNA_coxp_10th) #forest Plot HR Primary
ggforest(all_miRNA_coxm_5th) #forest Plot HR Metastasis
write_xlsx(Multicox_Primary_results_10th_,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Multicox_Primary_results_10th.xlsx")
write_xlsx(Multicox_Metas_results_5th_,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Multicox_Metas_results_5th.xlsx")
