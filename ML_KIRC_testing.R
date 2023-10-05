# Machine Learning Features miRNA
# Purpose : Find Feature importance and accuracy between ML 
# Adding feature : merging 2 or more miRNA into single variable
# In this code : 
# 1.Feature selection based on significant cox miRNA 
# 2.Try Correlation as additional feature selection (complement/additional)
# 3.Test Performace of each ML method

library("readxl")
library(caret) # sometimes error due to 'cli' is old version
library(dplyr)
library(mlbench)
#library(UBL) #for data imbalance
#library(themis) #for data imbalance
#read_excel("Path where your Excel file is stored\\File Name.xlsx")

options(scipen = 9) #avoid unit as exponential number as 1e+100

#Get miRNA name pick any data
# We can check as its condition OR use combined data
Metastasis=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Metastasis.xlsx")
Primary=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Primary.xlsx")
Normal=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\miRNA_Normal.xlsx")

cox_mirna_metastasis=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Multicox_Metas_results_5th.xlsx")
cox_mirna_primary=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Multicox_Primary_results_10th.xlsx")

##### PREPARING ML Comparison #####
# Filtered based on Cox miRNA #
id_cox_m=data.frame(cox_mirna_metastasis$row.names.Multicox_Metas_results_5th.) #get miRNA from cox metastasis
id_cox_p=data.frame(cox_mirna_primary$row.names.Multicox_Primary_results_10th.) #get miRNA from cox Primary
names(id_cox_m)="miRNA"
names(id_cox_p)="miRNA"
all_cox_mir=rbind(id_cox_m,id_cox_p)
dap=(paste(all_cox_mir$miRNA, collapse="|")) #To make it as "OR"
# To Filter, miRNA needs to be in Row positions (caution if the code run 2 times it will return back to original transpose)

#Back Up Variable (if FAIL): recovery code variable
#Metastasis2=Metastasis 
#Primary2=Primary
#Normal2=Normal
#Metastasis=Metastasis2
#Primary=Primary2
#Normal=Normal2
#d_Metastasis=data.frame(d_Metastasis)
#Metastasis_test=d_Metastasis %>% dplyr::select(matches(dap)) #take miRNA in cox list

d_Metastasis=t(Metastasis) #transpose matrix 
d_Primary=t(Primary) #transpose matrix
d_Normal=t(Normal) #transpose matrix

colnames(d_Metastasis)=d_Metastasis[1,] ; d_Metastasis <- d_Metastasis[-1, ] # first row into col names
colnames(d_Primary)=d_Primary[1,] ; d_Primary <- d_Primary[-1, ] # first row into col names
colnames(d_Normal)=d_Normal[1,] ; d_Normal <- d_Normal[-1, ] # first row into col names

d_Metastasis=data.frame(d_Metastasis) #Return to data frame so we can apply dplyr filter
d_Primary=data.frame(d_Primary) #Return to data frame so we can apply dplyr filter
d_Normal=data.frame(d_Normal) #Return to data frame so we can apply dplyr filter

#FILTER CODE
d_Metastasis=d_Metastasis %>% dplyr::select(matches(dap)) #take miRNA in cox in COLOM name
d_Primary=d_Primary%>% dplyr::select(matches(dap)) #take miRNA in cox list in COLOM name
d_Normal=d_Normal %>% dplyr::select(matches(dap)) #take miRNA in cox list in COLOM name

design_m=data.frame(row.names(d_Metastasis),code=2) ; colnames(design_m)=c("patient_ID","ML_tags")
design_p=data.frame(row.names(d_Primary),code=1) ; colnames(design_p)=c("patient_ID","ML_tags")
design_n=data.frame(row.names(d_Normal),code=0) ; colnames(design_n)=c("patient_ID","ML_tags")

d_Metastasis=data.frame(d_Metastasis,design_m$ML_tags) #combine miRNA-Patient data with Tags of disease stage
d_Primary=data.frame(d_Primary,design_p$ML_tags) #combine miRNA-Patient data with Tags of disease stage
d_Normal=data.frame(d_Normal,design_n$ML_tags) #combine miRNA-Patient data with Tags of disease stage

#standarized col names for ML_tags
colnames(d_Metastasis)[15] <- "ML_tags"
colnames(d_Primary)[15] <- "ML_tags"
colnames(d_Normal)[15] <- "ML_tags"

# Combine PvN & MvP (possible since miRNA name is the same/col names is same)
d_MvP=rbind(d_Metastasis,d_Primary)
d_PvN=rbind(d_Primary,d_Normal)

##TRY standar ML testing ##
# Source : https://machinelearningmastery.com/compare-the-performance-of-machine-learning-algorithms-in-r/
control <- trainControl(method="repeatedcv", number=10, repeats=3) #Default
# tags is : ML_tags
d_MvP$ML_tags =as.factor(d_MvP$ML_tags)
d_PvN$ML_tags =as.factor(d_PvN$ML_tags)

#TEST using smaller miRNA variable 
#Test=d_PvN[,1:10]
#Test=data.matrix(Test) #Input need data matrix
#Test=log2((Test)+1) #Need to Log2 some number to big
#Test=data.frame(Test,d_PvN$ML_tags)
#colnames(Test)[11] <- "ML_tags"

#WARNING : MAYBE NEED TO BE NUMERIC FOR MIRNA EXPRESSION
#QC testing data -1 (ALL NUMERIC)
str(d_PvN)
str(d_MvP)
d_PvN1=d_PvN #back up data
d_MvP1=d_MvP #back up data
#d_PvN1 <- data.frame(data.matrix(d_PvN1)) # Resest to data matrix then data frame [ REvised this code 25.07.2023]
d_PvN1  <- sapply(d_PvN1,as.numeric) #Successfully change character CHR into NUM; NOTE : although miRNA name still in char! maybe unstable
#d_MvP1 <- data.frame(data.matrix(d_MvP1)) # Resest to data matrix then data frame [ REvised this code 25.07.2023]
d_MvP1 <- sapply(d_MvP1,as.numeric) #Successfully change character CHR into NUM; NOTE : although miRNA name still in char! maybe unstable
#Save back up file
#library(writexl)
#d_PvN1=data.frame(row.names(d_PvN1),d_PvN1)
#d_MvP1=data.frame(row.names(d_MvP1),d_MvP1)
#write_xlsx(d_PvN1,"D:\\Dissertation 4.08.2022\\Pan-Cancer miRNA\\d_PvN1.xlsx")
#write_xlsx(d_MvP1,"D:\\Dissertation 4.08.2022\\Pan-Cancer miRNA\\d_MvP1.xlsx")

###EXPERIMENT-1# ML TESTING ## Data in RPM & ALL NUMERIC###
# Note: takes quite a minute to test ML -> move to PC computer
# CART
set.seed(7)
fit.cart <- train(ML_tags~., data=d_PvN1, method="rpart", trControl=control) #got warning due to not factor for ML_tags
# LDA
#set.seed(7)
#fit.lda <- train(ML_tags~., data=d_PvN1, method="lda", trControl=control) #ERROR SYSTEM STOP
# SVM
set.seed(7)
fit.svm <- train(ML_tags~., data=d_PvN1, method="svmRadial", trControl=control) #ERROR due to kernlab cannot be installed [say 'no' to download from source]
# kNN
set.seed(7)
fit.knn <- train(ML_tags~., data=d_PvN1, method="knn", trControl=control)
# Random Forest
set.seed(7)
fit.rf <- train(ML_tags~., data=d_PvN1, method="rf", trControl=control) #takes a long time to process
# collect resamples
#results <- resamples(list(CART=fit.cart, LDA=fit.lda, SVM=fit.svm, KNN=fit.knn, RF=fit.rf))
results <- resamples(list(CART=fit.cart, SVM=fit.svm, KNN=fit.knn, RF=fit.rf))

# Try Metastasis vs Primary with cox miRNA
# CART
set.seed(7)
fit.cart2 <- train(ML_tags~., data=d_MvP1, method="rpart", trControl=control)
# LDA
set.seed(7)
fit.lda2 <- train(ML_tags~., data=d_MvP1, method="lda", trControl=control) 
# SVM
set.seed(7)
fit.svm2 <- train(ML_tags~., data=d_MvP1, method="svmRadial", trControl=control)
# kNN
set.seed(7)
fit.knn2 <- train(ML_tags~., data=d_MvP1, method="knn", trControl=control)
# Random Forest
set.seed(7)
fit.rf2 <- train(ML_tags~., data=d_MvP1, method="rf", trControl=control)
# collect resamples
results2 <- resamples(list(CART=fit.cart2, SVM=fit.svm2, KNN=fit.knn2, RF=fit.rf2))

# dot plots of accuracy
# more about kappa statistic https://thedatascientist.com/performance-measures-cohens-kappa-statistic/
# Cohen’s kappa statistic is a very good measure that can handle very well both multi-class and imbalanced class problems.

# summarize differences between modes
summary(results)
# summarize differences between modes
summary(results2)

print(fit.cart) 
print(fit.svm)
print(fit.knn)
print(fit.rf)

# estimate variable importance
importance_fit.cart <- varImp(fit.cart,scale=FALSE) 
importance_fit.knn <- varImp(fit.knn, scale=FALSE) #Error in modifyList(data, lapply(data[, fc], as.numeric)) : 
#importance_fit.lda <- varImp(fit.lda, scale=FALSE)
importance_fit.rf <- varImp(fit.rf, scale=FALSE)
importance_fit.svm <- varImp(fit.svm, scale=FALSE) #error: Error in auc3_(actual, predicted, ranks) :

importance_fit.cart2 <- varImp(fit.cart2,scale=FALSE) 
importance_fit.knn2 <- varImp(fit.knn2, scale=FALSE) #Error in modifyList(data, lapply(data[, fc], as.numeric)) : 
#importance_fit.lda <- varImp(fit.lda2, scale=FALSE)
importance_fit.rf2 <- varImp(fit.rf2, scale=FALSE)
importance_fit.svm2 <- varImp(fit.svm2, scale=FALSE) #error: Error in auc3_(actual, predicted, ranks) :

print(importance_fit.cart)
print(importance_fit.knn) 
print(importance_fit.lda)
print(importance_fit.rf)
print(importance_fit.svm)

plot(importance_fit.cart)
plot(importance_fit.knn)
plot(importance_fit.lda)
plot(importance_fit.rf)
plot(importance_fit.svm)

print(importance_fit.cart2)
print(importance_fit.knn2) 
print(importance_fit.lda2)
print(importance_fit.rf2)
print(importance_fit.svm2)

plot(importance_fit.cart2)
plot(importance_fit.knn2)
plot(importance_fit.lda2)
plot(importance_fit.rf2)
plot(importance_fit.svm2)

###WARNING : MAYBE NEED TO LOG2(RPM) TO GET BETTER RESULT###
#QC testing data -1 (combined NUMERIC, Factor, as.table + Log2(RPM+1))
d_PvN2=d_PvN[,1:14] #back up data
d_MvP2=d_MvP[,1:14] #back up data

d_PvN2 <- data.frame(data.matrix(d_PvN2)) # Resest to data matrix then data frame
d_PvN2  <- sapply(d_PvN2,as.numeric) #Successfully change character CHR into NUM; NOTE : although miRNA name still in char! maybe unstable
d_MvP2 <- data.frame(data.matrix(d_MvP2)) # Resest to data matrix then data frame
d_MvP2 <- sapply(d_MvP2,as.numeric) #Successfully change character CHR into NUM; NOTE : although miRNA name still in char! maybe unstable

d_PvN2=log2(d_PvN2+1) #Converting to log2(RPM+1)
d_MvP2=log2(d_MvP2+1) #Converting to log2(RPM+1)

d_PvN2=data.frame(d_PvN2) #converting back to data.frame
d_MvP2=data.frame(d_MvP2) #converting back to data.frame

ML_tags =as.factor(d_PvN$ML_tags) #adding Tags as factor
ML_tags_ =as.factor(d_MvP$ML_tags) #adding tags as factor

d_PvN2=data.frame(d_PvN2,ML_tags) #combined both numeric + factor
d_MvP2=data.frame(d_MvP2,ML_tags_) #combined both numeric + factor
colnames(d_MvP2)[15] <- "ML_tags"

#EXPERIMENT-2# ML TESTING ## Data in LOG2RPM & COMBINED NUMERIC & fACTOR
# Note: takes quite a minute to test ML -> move to PC computer
# note: MOST STABLE WHEN CHANGED THE VARIABLE TYPE ACCORIDING TO INSTRUCTION FROM : https://machinelearningmastery.com/compare-the-performance-of-machine-learning-algorithms-in-r/
# CART
set.seed(7)
fit.cart3 <- train(ML_tags~., data=d_PvN2, method="rpart", trControl=control) #got warning due to not factor for ML_tags
# LDA
set.seed(7)
fit.lda3 <- train(ML_tags~., data=d_PvN2, method="lda", trControl=control) #ERROR SYSTEM STOP
# SVM
set.seed(7)
fit.svm3 <- train(ML_tags~., data=d_PvN2, method="svmRadial", trControl=control) #ERROR due to kernlab cannot be installed [say 'no' to download from source]
# kNN
set.seed(7)
fit.knn3 <- train(ML_tags~., data=d_PvN2, method="knn", trControl=control)
# Random Forest
set.seed(7)
fit.rf3 <- train(ML_tags~., data=d_PvN2, method="rf", trControl=control,) #takes a long time to process
# collect resamples
results3 <- resamples(list(CART=fit.cart3, LDA=fit.lda3, SVM=fit.svm3, KNN=fit.knn3, RF=fit.rf3))
#results3 <- resamples(list(CART=fit.cart, SVM=fit.svm, KNN=fit.knn, RF=fit.rf))

#REVIEW XGboost 23.06.2023 [Primary Vs Normal]
# Code from : https://towardsdatascience.com/a-guide-to-using-caret-in-r-71dec0bda208
# Error/ Warning you can check : https://stackoverflow.com/questions/70453557/caret-xgbtree-warning-ntree-limit-is-deprecated-use-iteration-range-instea
set.seed(7)
fit.xgboost3 <- train(ML_tags~., data=d_PvN2, method="xgbTree", trControl=control, verbose = FALSE, verbosity = 0)
results_23062023 <- resamples(list(CART=fit.cart3, LDA=fit.lda3, SVM=fit.svm3, KNN=fit.knn3, RF=fit.rf3, XGB=fit.xgboost3))
summary(results_23062023)

# Try Metastasis vs Primary with cox miRNA
# CART
set.seed(7)
fit.cart4 <- train(ML_tags~., data=d_MvP2, method="rpart", trControl=control)
# LDA
set.seed(7)
fit.lda4 <- train(ML_tags~., data=d_MvP2, method="lda", trControl=control) 
# SVM
set.seed(7)
fit.svm4 <- train(ML_tags~., data=d_MvP2, method="svmRadial", trControl=control)
# kNN
set.seed(7)
fit.knn4 <- train(ML_tags~., data=d_MvP2, method="knn", trControl=control)
# Random Forest
set.seed(7)
fit.rf4 <- train(ML_tags~., data=d_MvP2, method="rf", trControl=control)
# collect resamples
results4 <- resamples(list(CART=fit.cart4, LDA=fit.lda4, SVM=fit.svm4, KNN=fit.knn4, RF=fit.rf4))

#REVIEW XGboost 23.06.2023 [Primary Vs Normal]
# Code from : https://towardsdatascience.com/a-guide-to-using-caret-in-r-71dec0bda208
# Error/ Warning you can check : https://stackoverflow.com/questions/70453557/caret-xgbtree-warning-ntree-limit-is-deprecated-use-iteration-range-instea
set.seed(7)
fit.xgboost4 <- train(ML_tags~., data=d_MvP2, method="xgbTree", trControl=control, verbose = FALSE, verbosity = 0)
results_23062023_mvp <- resamples(list(CART=fit.cart4, LDA=fit.lda4, SVM=fit.svm4, KNN=fit.knn4, RF=fit.rf4, XGB=fit.xgboost4))
summary(results_23062023_mvp)


# dot plots of accuracy
# more about kappa statistic https://thedatascientist.com/performance-measures-cohens-kappa-statistic/
# Cohen’s kappa statistic is a very good measure that can handle very well both multi-class and imbalanced class problems.

# summarize differences between modes Primary vs Normal 
summary(results3)
# summarize differences between modes Metastasis vs Primary
summary(results4)

# show individual results for ML
#print(fit.cart3) 
#print(fit.svm)
#print(fit.knn)
#print(fit.rf)

# estimate variable importance
importance_fit.cart3 <- varImp(fit.cart3,scale=FALSE) 
importance_fit.knn3 <- varImp(fit.knn3, scale=FALSE) #Error in modifyList(data, lapply(data[, fc], as.numeric)) : 
importance_fit.lda3 <- varImp(fit.lda3, scale=FALSE)
importance_fit.rf3 <- varImp(fit.rf3, scale=FALSE)
importance_fit.svm3 <- varImp(fit.svm3, scale=FALSE) #error: Error in auc3_(actual, predicted, ranks) :
importance_fit.xgboost3 <- varImp(fit.xgboost3, scale=FALSE) #error: Error in auc3_(actual, predicted, ranks) :

importance_fit.cart4 <- varImp(fit.cart4,scale=FALSE) 
importance_fit.knn4 <- varImp(fit.knn4, scale=FALSE) #Error in modifyList(data, lapply(data[, fc], as.numeric)) : 
importance_fit.lda4 <- varImp(fit.lda4, scale=FALSE)
importance_fit.rf4 <- varImp(fit.rf4, scale=FALSE)
importance_fit.svm4 <- varImp(fit.svm4, scale=FALSE) #error: Error in auc3_(actual, predicted, ranks) :
importance_fit.xgboost4 <- varImp(fit.xgboost4, scale=FALSE) #error: Error in auc3_(actual, predicted, ranks) :

print(importance_fit.cart3)
print(importance_fit.knn3) 
print(importance_fit.lda3)
print(importance_fit.rf3)
print(importance_fit.svm3)

plot(importance_fit.cart3)
plot(importance_fit.knn3)
plot(importance_fit.lda3)
plot(importance_fit.rf3)
plot(importance_fit.svm3)
plot(importance_fit.xgboost3)

print(importance_fit.cart4)
print(importance_fit.knn4) 
print(importance_fit.lda4)
print(importance_fit.rf4)
print(importance_fit.svm4)

plot(importance_fit.cart4)
plot(importance_fit.knn4)
plot(importance_fit.lda4)
plot(importance_fit.rf4)
plot(importance_fit.svm4)
plot(importance_fit.xgboost4)

# box and whisker plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results3, scales=scales)
bwplot(results4, scales=scales)
bwplot(results_23062023, scales=scales)
bwplot(results_23062023_mvp, scales=scales) 

densityplot(results3, scales=scales, pch = "|")
densityplot(results4, scales=scales, pch = "|")
densityplot(results_23062023, scales=scales, pch = "|")
densityplot(results_23062023_mvp, scales=scales, pch = "|")

#Save FILE as TXT [the LOG ML one; note :13.04.2023]
sink("MachineLearning_PvN.txt")
print(summary(results3))
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning_MvP.txt")
print(summary(results4))
sink()  # returns output to the console in Documents [working document]

# NOTE : Kappa values of 0.4 to 0.75 are considered moderate to good and a kappa of >0.75 represents excellent agreement.
# NOTE : The kappa coefficient measures the agreement between classification and truth values

#######################################
#         13.04.2023                  #
#          Validation                 #
#######################################
# 13.04.2023 Confusion Matrix in TCGA datasets Primary Vs Normal
# Handling with NA's Prediction: https://github.com/topepo/caret/issues/809
# using LDA [More Stable for Probability]
predictions_TCGA_KIRC_LDA <- predict(fit.lda3 , d_PvN2[,-15], type= "prob") #Testing Accuracy
#predictions_TCGA_KIRC=data.frame(predictions_TCGA_KIRC); 
pred_lda3 = ifelse(predictions_TCGA_KIRC_LDA$`1` > 0.5, "1", "0")
pred_lda3_factor= as.factor(pred_lda3)
lda3_matrix=confusionMatrix(d_PvN2$ML_tags,pred_lda3_factor) # check inside if we want to know the results

library("readxl")
GSE151423=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\....GSE151423 miRNA raw\\GSE151423_miRNARAW_filtered_renamed.xlsx")
#GSE151423_miR=GSE151423[,1]  #Get First Col = miRNA name
#GSE151423=GSE151423[,-1] # Removing miRNA names
GSE151423_MLtags <- as.factor(GSE151423$ML_tags)
GSE151423 <- data.frame(data.matrix(GSE151423)) # Set as data matrix
GSE151423 = sapply(GSE151423,as.numeric)  # set as Numeric

#Try without CPM conversion function
GSE151423_log2 = log2((GSE151423)+1)  # set as Numeric
GSE151423_log2_= GSE151423_log2[,3:16] # Similar to TCGA
predictions_GSE151423_LDA <- predict(fit.lda3 , GSE151423_log2_, type= "prob") #Testing Accuracy using our LDA model
pred_lda3_GSE151423 = ifelse(predictions_GSE151423_LDA$`1` > 0.5, "1", "0")
pred_lda3_GSE151423_factor= as.factor(pred_lda3_GSE151423 )
lda3_matrix_GSE151423=confusionMatrix(GSE151423_MLtags,pred_lda3_GSE151423_factor) # check inside if we want to know the results
Check=data.frame(GSE151423_MLtags,pred_lda3_GSE151423_factor)
# Conclusion: 100% True

## 13.04.2023 Confusion Matrix in TCGA datasets Metastasis vs Primary
# LDA
predictions_TCGA_KIRC_LDA_MVP <- predict(fit.lda4 , d_MvP2[,-15], type= "prob") #Testing Accuracy
#predictions_TCGA_KIRC=data.frame(predictions_TCGA_KIRC); 
pred_lda4 = ifelse(predictions_TCGA_KIRC_LDA_MVP$`2` > 0.5, "2", "1")
pred_lda4_factor= as.factor(pred_lda4)
lda4_matrix=confusionMatrix(d_MvP2$ML_tags,pred_lda4_factor) # check inside if we want to know the results
# Conclusion : Using all miRNA cox we able to predict 80% effect between Primary and Metastasis

# RF
predictions_TCGA_KIRC_RF_MVP <- predict(fit.rf4 , d_MvP2[,-15], type= "prob") #Testing Accuracy
#predictions_TCGA_KIRC=data.frame(predictions_TCGA_KIRC); 
pred_rf4 = ifelse(predictions_TCGA_KIRC_RF_MVP$`2` > 0.5, "2", "1")
pred_rf4_factor= as.factor(pred_rf4)
rf4_matrix_coxmiR=confusionMatrix(d_MvP2$ML_tags,pred_rf4_factor) # check inside if we want to know the results
# Conclusion : Using all miRNA cox RF able to predict 100% effect between Primary and  but Still Kappa is Low


# Try to Remove certain Feature in Metastasis vs Primary case
# Link : https://topepo.github.io/caret/feature-selection-using-univariate-filters.html
# Link : https://rdrr.io/cran/caret/man/sbfControl.html
# Univariate Feature
filterCtrl <- sbfControl(functions = rfSBF, method="repeatedcv", number=10, repeats=3) #Default
set.seed(10)
a=as.numeric(d_MvP2$ML_tags)
rfWithFilter <- sbf(d_MvP2[-15], a, sbfControl = filterCtrl) #Error in { : task 1 failed - "package gam is required"
#Solution : create separate predictor, and remove tags from train datasets
# check : https://www.rdocumentation.org/packages/caret/versions/6.0-92/topics/sbf
# Results :
#   During resampling, the top 5 selected variables (out of a possible 9):
#     MIMAT0000254 (100%), MIMAT0000450 (100%), MIMAT0000280 (83.3%), MIMAT0000073 (56.7%), MIMAT0004494 (56.7%)

d_MvP2_Filtered = data.frame(d_MvP2$MIMAT0000254,d_MvP2$MIMAT0000450,d_MvP2$MIMAT0000280,d_MvP2$MIMAT0000073,d_MvP2$MIMAT0004494,d_MvP2$ML_tags)
colnames(d_MvP2_Filtered)= c( "MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "MIMAT0000073" , "MIMAT0004494", "ML_tags")   

set.seed(7)
fit.cart4_fil <- train(ML_tags~., data=d_MvP2_Filtered, method="rpart", trControl=control)
# LDA
set.seed(7)
fit.lda4_fil <- train(ML_tags~., data=d_MvP2_Filtered, method="lda", trControl=control) 
# SVM
set.seed(7)
fit.svm4_fil <- train(ML_tags~., data=d_MvP2_Filtered, method="svmRadial", trControl=control)
# kNN
set.seed(7)
fit.knn4_fil <- train(ML_tags~., data=d_MvP2_Filtered, method="knn", trControl=control)
# Random Forest
set.seed(7)
fit.rf4_fil <- train(ML_tags~., data=d_MvP2_Filtered, method="rf", trControl=control)
results5 <- resamples(list(CART=fit.cart4_fil, LDA=fit.lda4_fil, SVM=fit.svm4_fil, KNN=fit.knn4_fil, RF=fit.rf4_fil))
summary(results5) #Kappa Still Low 

#Xgboost using filtered (23.06.2023); testing for validation in LINE 547
set.seed(7)
fit.xgboost4_fil <- train(ML_tags~., data=d_MvP2_Filtered, method="xgbTree", trControl=control, verbose = FALSE, verbosity = 0)
results5_xg <- resamples(list(CART=fit.cart4_fil, LDA=fit.lda4_fil, SVM=fit.svm4_fil, KNN=fit.knn4_fil, RF=fit.rf4_fil, XGB=fit.xgboost4_fil))
summary(results5_xg) #Kappa Still Low 

# Testing Prediction using RF
predictions_TCGA_KIRC_RF_MVP_5miR <- predict(fit.rf4_fil  , d_MvP2_Filtered[,-6], type= "prob") #Testing Accuracy
#predictions_TCGA_KIRC=data.frame(predictions_TCGA_KIRC); 
pred_rf4_fil = ifelse(predictions_TCGA_KIRC_RF_MVP_5miR$`2` > 0.5, "2", "1")
pred_rf4_fil_factor= as.factor(pred_rf4_fil)
rf4_matrix_5miR=confusionMatrix( d_MvP2_Filtered$ML_tags,pred_rf4_fil_factor) # check inside if we want to know the results



# Experiment 2 for MvP Kappa using TOP 3 miRNA 100% 100% 83%
d_MvP2_Filtered2 = data.frame(d_MvP2$MIMAT0000254,d_MvP2$MIMAT0000450,d_MvP2$MIMAT0000280,d_MvP2$ML_tags)
colnames(d_MvP2_Filtered2)= c( "MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "ML_tags") 

set.seed(7)
fit.cart4_fil2 <- train(ML_tags~., data=d_MvP2_Filtered2, method="rpart", trControl=control)
# LDA
set.seed(7)
fit.lda4_fil2 <- train(ML_tags~., data=d_MvP2_Filtered2, method="lda", trControl=control) 
# SVM
set.seed(7)
fit.svm4_fil2 <- train(ML_tags~., data=d_MvP2_Filtered2, method="svmRadial", trControl=control)
# kNN
set.seed(7)
fit.knn4_fil2 <- train(ML_tags~., data=d_MvP2_Filtered2, method="knn", trControl=control)
# Random Forest
set.seed(7)
fit.rf4_fil2 <- train(ML_tags~., data=d_MvP2_Filtered2, method="rf", trControl=control)
results6 <- resamples(list(CART=fit.cart4_fil2, LDA=fit.lda4_fil2, SVM=fit.svm4_fil2, KNN=fit.knn4_fil2, RF=fit.rf4_fil2))
summary(results6) #only RF is increase

# add xgboost for 3miRNA 23.06.2023
set.seed(7)
fit.xgboost4_fil2 <- train(ML_tags~., data=d_MvP2_Filtered2, method="xgbTree", trControl=control, verbose = FALSE, verbosity = 0)
results6_xg <- resamples(list(CART=fit.cart4_fil2, LDA=fit.lda4_fil2, SVM=fit.svm4_fil2, KNN=fit.knn4_fil2, RF=fit.rf4_fil2, XGB=fit.xgboost4_fil2))
summary(results6_xg)

# Testing Prediction using RF
predictions_TCGA_KIRC_RF_MVP <- predict(fit.rf4_fil2  , d_MvP2_Filtered2[,-4], type= "prob") #Testing Accuracy
#predictions_TCGA_KIRC=data.frame(predictions_TCGA_KIRC); 
pred_rf4_fil2 = ifelse(predictions_TCGA_KIRC_RF_MVP$`2` > 0.5, "2", "1")
pred_rf4_fil2_factor= as.factor(pred_rf4_fil2)
rf4_matrix=confusionMatrix( d_MvP2_Filtered2$ML_tags,pred_rf4_fil2_factor) # check inside if we want to know the results

#Save FILE as TXT [14.04.2023]
sink("MachineLearning_MvP_5miR.txt")
print(summary(results5))
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning_MvP_3miR.txt")
print(summary(results6))
sink()  # returns output to the console in Documents [working document]

importance_fit.cart4_fil <- varImp(fit.cart4_fil,scale=FALSE) 
importance_fit.knn4_fil <- varImp(fit.knn4_fil, scale=FALSE) #Error in modifyList(data, lapply(data[, fc], as.numeric)) : 
importance_fit.lda4_fil <- varImp(fit.lda4_fil, scale=FALSE)
importance_fit.rf4_fil <- varImp(fit.rf4_fil, scale=FALSE)
importance_fit.svm4_fil <- varImp(fit.svm4_fil, scale=FALSE) #error: Error in auc3_(actual, predicted, ranks) :
plot(importance_fit.cart4_fil)
plot(importance_fit.knn4_fil)
plot(importance_fit.lda4_fil)
plot(importance_fit.rf4_fil) # score significant changed
plot(importance_fit.svm4_fil)

#######################################
# Validation Metastasis vs Primary    #
#######################################
## Validation GSE207557 [use to predict 5 year metastasis development]
#Testing prediction using RF [all cox-miRNA]
GSE207557_all=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\XXXXX GSE207557 PMID 35931808\\GSE207557_ALLcoxmiR_logCPM.xlsx")
GSE207557_all_MLtags <- as.factor(GSE207557_all$group)
GSE207557_all <- data.frame(data.matrix(GSE207557_all)) # Set as data 
GSE207557_all_=data.frame(GSE207557_all[3:16])
colnames(GSE207557_all_)=c( "MIMAT0004494",
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

#testing all miRNA cox# WARNING!!!!!!!!! FAIL to OBTAIN RESULTS
predictions_GSE207557_RF_MVP_allmiR <- predict(fit.rf4  , GSE207557_all_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_rf4_GSE207557_allmiR = ifelse(predictions_GSE207557_RF_MVP_allmiR$`2` > 0.5, "2", "1")
pred_rf4_GSE207557_allmiR_factor= as.factor(pred_rf4_GSE207557_allmiR)
rf4_matrix_GSE207557_allmiR=confusionMatrix(GSE207557_all_MLtags,pred_rf4_GSE207557_allmiR_factor) # check inside if we want to know the results

#[23.06.2023] XGboost for all miRNA ---Fail to predict
predictions_GSE207557_Xgboost_MVP_allmiR <- predict(fit.xgboost4, GSE207557_all_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_xgboost_GSE207557_allmiR = ifelse(predictions_GSE207557_Xgboost_MVP_allmiR$`2` > 0.5, "2", "1")
pred_xgboost_GSE207557_allmiR_factor= as.factor(pred_xgboost_GSE207557_allmiR)
xgboost_matrix_GSE207557_allmiR=confusionMatrix(GSE207557_all_MLtags,pred_xgboost_GSE207557_allmiR_factor) # check inside if we want to know the results

#testing 3 & 5 miRNA cox
GSE207557=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\XXXXX GSE207557 PMID 35931808\\GSE207557_5miR_logCPM.xlsx")
GSE207557_MLtags <- as.factor(GSE207557$group)
GSE207557 <- data.frame(data.matrix(GSE207557)) # Set as data matrix
GSE207557_ <- data.frame(GSE207557[,3:7])
colnames(GSE207557_)=c( "MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "MIMAT0000073" , "MIMAT0004494")

# Testing Prediction using RF GSE207557 [5-miR]
predictions_GSE207557_RF_MVP_5miR <- predict(fit.rf4_fil  , GSE207557_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_rf4_GSE207557 = ifelse(predictions_GSE207557_RF_MVP_5miR$`2` > 0.5, "2", "1")
pred_rf4_GSE207557_factor= as.factor(pred_rf4_GSE207557)
rf4_matrix_GSE207557_5miR=confusionMatrix(GSE207557_MLtags,pred_rf4_GSE207557_factor) # check inside if we want to know the results

#[23.06.2023] Testing using Xgboost GSE207557  for 5 miRNA ---Fail to predict
predictions_GSE207557_Xgboost_MVP_5miR <- predict(fit.xgboost4_fil, GSE207557_, type= "prob")
pred_xgboost_GSE207557 = ifelse(predictions_GSE207557_Xgboost_MVP_5miR$`2` > 0.5, "2", "1")
pred_xgboost_GSE207557_factor= as.factor(pred_xgboost_GSE207557)
xgboost_matrix_GSE207557_5miR=confusionMatrix(GSE207557_MLtags,pred_xgboost_GSE207557_factor) # check inside if we want to know the results

#*****testing Prediction using LDA GSE207557 [5-miR] 20.05.2023
#predictions_GSE207557_LDA_MVP_5miR <- predict(fit.lda4_fil  , GSE207557_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
#pred_lda4_GSE207557 = ifelse(predictions_GSE207557_LDA_MVP_5miR$`2` > 0.5, "2", "1") #check if the factor "level" only show 1, meaning this is not good indicator
# TO make the level keep as 2 level, we edited the first one [causing inaccuracy!]
#pred_lda4_GSE207557 [1]="2" #change the first value into to 2; so at least there are 2 level to compare in confusion matrix
#pred_lda4_GSE207557_factor= as.factor(pred_lda4_GSE207557)
#lda4_matrix_GSE207557_5miR=confusionMatrix(GSE207557_MLtags,pred_lda4_GSE207557_factor) # check inside if we want to know the results
# error : the data cannot have more levels than the reference; meaning all predicted as cathegory 1 only
    #solution, just change one of the value, so it at least have 2 level "1" and "2"; it works but thats consider as data fabrication    

#testing Prediction using KNN GSE207557 [5-miR] 20.05.2023
predictions_GSE207557_KNN_MVP_5miR <- predict(fit.knn4_fil  , GSE207557_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_knn4_GSE207557 = ifelse(predictions_GSE207557_KNN_MVP_5miR$`2` > 0.5, "2", "1") #check if the factor "level" only show 1, meaning this is not good indicator
pred_knn4_GSE207557_factor= as.factor(pred_knn4_GSE207557)
knn4_matrix_GSE207557_5miR=confusionMatrix(GSE207557_MLtags,as.factor(pred_knn4_GSE207557_factor))
          

# Testing Prediction using RF [3-miR] # # WARNING!!!!!!!!! FAIL to OBTAIN RESULTS
predictions_GSE207557_RF_MVP_3miR <- predict(fit.rf4_fil2  , GSE207557_[1:3], type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_rf4_GSE207557_3mir = ifelse(predictions_GSE207557_RF_MVP_3miR$`2` > 0.5, "2", "1")
pred_rf4_GSE207557_3mir_factor= as.factor(pred_rf4_GSE207557_3mir)
rf4_matrix_GSE207557_3miR=confusionMatrix(GSE207557_MLtags,pred_rf4_GSE207557_3mir_factor) # check inside if we want to know the results

#[23.06.2023] Testing using Xgboost GSE207557  for 3 miRNA ---Fail to predict
predictions_GSE207557_Xgboost_MVP_3miR <- predict(fit.xgboost4_fil2, GSE207557_, type= "prob")
pred_xgboost_GSE207557_3miR = ifelse(predictions_GSE207557_Xgboost_MVP_3miR$`2` > 0.5, "2", "1")
pred_xgboost_GSE207557_factor_3miR= as.factor(pred_xgboost_GSE207557_3miR)
xgboost_matrix_GSE207557_3miR=confusionMatrix(GSE207557_MLtags,pred_xgboost_GSE207557_factor_3miR) # check inside if we want to know the results


#ONLY TESTING GSE155209 [microarray : use to predict 5 year metastasis development]
GSE155209_all=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\GSE155209 miRNA PMID 32850445\\GSE155209_ARRAY_5miR.xlsx")
GSE155209_all_MLtags <- as.factor(GSE155209_all$group)
GSE155209_all <- data.frame(data.matrix(GSE155209_all)) # Set as data 
GSE155209_all_=data.frame(GSE155209_all[3:7])
colnames(GSE155209_all_)=c( "MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "MIMAT0000073" , "MIMAT0004494")

# Testing Prediction using RF GSE155209 [5-miR]
GSE155209_all_=data.frame(GSE155209_all[2:6])
predictions_GSE155209_RF_MVP_5miR <- predict(fit.rf4_fil_GSE155209, GSE155209_all_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_rf4_GSE155209 = ifelse(predictions_GSE155209_RF_MVP_5miR$`2` > 0.5, "2", "1")
pred_rf4_GSE155209_factor= as.factor(pred_rf4_GSE155209)
rf4_matrix_GSE155209_5miR=confusionMatrix(GSE155209_all_MLtags,pred_rf4_GSE155209_factor) # check inside if we want to know the results


## Validation GSE93175 [all metastasis from plasma rcc] ##unfortunately No Results
GSE93175_all=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\XXno comparatorXX GSE93175 miRNA PMID  28969022\\GSE93175_miRNA_log2_selected.xlsx")
GSE93175_all_MLtags <- as.factor(GSE93175_all$group)
GSE93175_all <- data.frame(data.matrix(GSE93175_all)) # Set as data 
GSE93175_all_=data.frame(GSE93175_all[3:7])
colnames(GSE93175_all_)=c( "MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "MIMAT0000073" , "MIMAT0004494")

# Testing Prediction using RF GSE93175 [5-miR]
predictions_GSE93175_RF_MVP_5miR <- predict(fit.rf4_fil  , GSE93175_all_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_rf4_GSE93175 = ifelse(predictions_GSE93175_RF_MVP_5miR$`2` > 0.5, "2", "1")
pred_rf4_GSE93175_factor= as.factor(pred_rf4_GSE93175)
rf4_matrix_GSE93175_5miR=confusionMatrix(GSE93175_all_MLtags,pred_rf4_GSE93175_factor) 

#Save Conf Matix [14.04.2023]
# Save
sink("MachineLearning_lda_PvN_Predict_GSE151423.txt")
print(lda3_matrix_GSE151423)
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning__lda_PvN_Predict_TCGA.txt")
print(lda3_matrix)
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning__rf_MvP_Predict_TCGA_5miR.txt")
print(rf4_matrix_5miR)
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning__rf_MvP_Predict_TCGA.txt")
print(rf4_matrix_coxmiR)
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning__rf_MvP_Predict_GSE207557_progression.txt")
print(rf4_matrix_GSE207557_5miR)
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning__rf_MvP_Predict_GSE155209_progression.txt")
print(rf4_matrix_GSE155209_5miR)
sink()  # returns output to the console in Documents [working document]

#[23.06.2023] xgboost update SAVE_New
# Save
sink("MachineLearning_PvN_23.06.2023.txt")
print(summary(results_23062023))
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning_MvP_23.06.2023.txt")
print(summary(results_23062023_mvp))
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning_MvP_5miR_23.06.2023.txt")
print(summary(results5_xg))
sink()  # returns output to the console in Documents [working document]

sink("MachineLearning_MvP_3miR_23.06.2023.txt")
print(summary(results6_xg))
sink()  # returns output to the console in Documents [working document]

###########################
## MICROARRAY VALIDATION : Metastasis vs Primary ##
###########################

# constructing ML usig GSE155209 [microarray : use to predict 5 year metastasis development] FAIL 
#GSE155209_all=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\GSE155209 miRNA PMID 32850445\\GSE155209_ARRAY_5miR.xlsx")
#GSE155209_all_MLtags <- as.factor(GSE155209_all$group)
#GSE155209_all$group= as.factor(GSE155209_all$group) 
#GSE155209_all <- data.frame(data.matrix(GSE155209_all)) # Set as data 
#colnames(GSE155209_all)=c( "ID","group","MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "MIMAT0000073" , "MIMAT0004494")
#GSE155209_all = GSE155209_all [,-1]#removing Patient's ID

#ONLY TESTING GSE155209 [microarray : use to predict 5 year metastasis development]
GSE155209_all=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\GSE155209 miRNA PMID 32850445\\GSE155209_ARRAY_5miR.xlsx")
GSE155209_all_MLtags <- as.factor(GSE155209_all$group)
GSE155209_all <- data.frame(data.matrix(GSE155209_all)) # Set as data 
GSE155209_all_=data.frame(GSE155209_all[3:7])
colnames(GSE155209_all_)=c( "MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "MIMAT0000073" , "MIMAT0004494")

# Testing Prediction using RF GSE155209 [5-miR]
GSE155209_all_=data.frame(GSE155209_all[2:6])
predictions_GSE155209_RF_MVP_5miR <- predict(fit.rf4_fil_GSE155209, GSE155209_all_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_rf4_GSE155209 = ifelse(predictions_GSE155209_RF_MVP_5miR$`2` > 0.5, "2", "1")
pred_rf4_GSE155209_factor= as.factor(pred_rf4_GSE155209)
rf4_matrix_GSE155209_5miR=confusionMatrix(GSE155209_all_MLtags,pred_rf4_GSE155209_factor) # check inside if we want to know the results


# Validation GSE131959 [microarray : use to predict 5 year metastasis development] NA replace with mean value
GSE131959_all=read_excel("D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\Validation\\GSE Study Validation\\GSE131959 marray MIR\\GSE131959_miRNA_nonNormalize_array.xlsx")
GSE131959_all_MLtags <- as.factor(GSE131959_all$group)
GSE131959_all <- data.frame(data.matrix(GSE131959_all)) # Set as data 
GSE131959_all_=data.frame(GSE131959_all[3:7])
colnames(GSE131959_all_)=c( "MIMAT0000254" , "MIMAT0000450" , "MIMAT0000280" , "MIMAT0000073" , "MIMAT0004494")

GSE131959_all_=log2(GSE131959_all_)
# Testing Prediction using RF GSE131959 [5-miR]
predictions_GSE131959_RF_MVP_5miR <- predict(fit.rf4_fil  , GSE131959_all_, type= "prob") #Testing Accuracy, note that the score is low due to real expression to train is bigger
pred_rf4_GSE131959 = ifelse(predictions_GSE131959_RF_MVP_5miR$`2` > 0.5, "2", "1")
pred_rf4_GSE131959_factor= as.factor(pred_rf4_GSE131959)
rf4_matrix_GSE131959_5miR=confusionMatrix(GSE131959_all_MLtags,pred_rf4_GSE131959_factor) # check inside if we want to know the results

# 25 07 2023
# Save SCV for datasets
write.csv(d_MvP,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\datasets_MvP.csv",row.names=TRUE)
write.csv(d_MvP1,"D:\\Dissertation 4.08.2022\\Publication - 1 Bioinformatics miRNA\\datasets_log2_MvP.csv",row.names=TRUE)





#------------ Preparing datasets for CORRELATION miRNA  ---------------- #
#Set Metastasis
Metastasis_=t(Metastasis)  #transpose data frame
colnames(Metastasis_) <- Metastasis_[1,] #set first row into col names
Metastasis_ <- Metastasis_[-1, ] #remove first row (old ones)
row.names(Metastasis_) <- NULL

#Set Primary
Primary_=t(Primary)  #transpose data frame
colnames(Primary_) <- Primary_[1,] #set first row into col names
Primary_ <- Primary_[-1, ] #remove first row (old ones)
row.names(Primary_) <- NULL

Metastasis_ <- data.frame(data.matrix(Metastasis_)) # Resest to data matrix then data frame
Metastasis_  <- sapply(Metastasis_,as.numeric) #Successfully change character CHR into NUM; NOTE : although miRNA name still in char! maybe unstable

Primary_ <- data.frame(data.matrix(Primary_)) # Resest to data matrix then data frame
Primary_ <- sapply(Primary_,as.numeric) #Successfully change character CHR into NUM; NOTE : although miRNA name still in char! maybe unstable

# JUNK FAIL CODE : converting character to numeric
#Metastasis_  <- as.numeric(Metastasis_) 
#Metastasis_  <- lapply(Metastasis_, as.numeric) #from char to Numeric

#QC check data strtucture type
#str(Metastasis_)

#+++++++++# [ADDITIIONAL] : Feature Selection using Correlation #+++++++++# 
#correlation calculation intepretaion : http://www.statstutor.ac.uk/resources/uploaded/spearmans.pdf
correlationMatrixm <- cor(Metastasis_)
highlyCorrelatedm <- findCorrelation(correlationMatrixm, cutoff = .8) #very strong #verbos function to print in console high correlated miRNA
correlationMatrixp <- cor(Primary_)
highlyCorrelatedp <- findCorrelation(correlationMatrixp, cutoff = .8)
#highlyCorrelated <- findCorrelation(correlationMatrix, cutoff = .6, verbose = TRUE) #strong correlation
#highcorr=data.frame(highlyCorrelated)
print(highlyCorrelatedm)
hcm = sort(highlyCorrelatedm)
hc_mirna_m = Metastasis_[,c(hcm)] #its left you with 284/630 miRNA highly correlated keep
print(highlyCorrelatedp)
hcp = sort(highlyCorrelatedp)
hc_mirna_p = Metastasis_[,c(hcp)] #its left you with 128/630 miRNA highly correlated keep
#show in pheatmap to see which miRNA is highly correlated
library(pheatmap)
pheatmap(correlationMatrixm,cutree_rows = 3,
         cutree_cols = 3) #show cluster of highcorrelation matrix wait few min
pheatmap(correlationMatrixp,cutree_rows = 3,
         cutree_cols = 3) #show cluster of highcorrelation matrix wait few min
#install.packages("corrplot")
#library(corrplot)
#corrplot(correlationMatrix, method="circle") #WARNING!!! takes a lot of time




## QC code to ensure tag in place ##
#d_Primary[,630:631] #Check last col 
#d_Metastasis[,630:631] #Check last col 
#d_Normal[,630:631] #Check last col 
#d_PvN[,630:631]

# Combine PvN & MvP (possible since miRNA name is the same/col names is same)
d_MvP=rbind(d_Metastasis,d_Primary)
d_PvN=rbind(d_Primary,d_Normal)


#Error: no applicable method for 'select' applied to an object of class "c('matrix', 'array', 'character')"
#Note : it will remove ML tags!


## ML testing Optimization FIND importance feature from each ML method 

## Clustering in patients

## Molecular Pathway analysis